# Nonparametric Spatio-Temporal Hawkes (ETAS-type) fitting
# --------------------------------------------------------
# How to run:
#   - R >= 4.2 recommended
#   - Install required packages listed in `required_pkgs` (use renv or install.packages)
#   - Set `data_path_events*` and `output_path` below
#   - Source this file and run from R or Rscript



# ---- Configuration ------------------------------------------------------------

options(stringsAsFactors = FALSE)
set.seed(1)

required_pkgs <- c(
  "doParallel", "foreach", "dplyr", "lubridate", "sf",
  "spatstat", "fields", "RANN", "np", "sparr", "ks"
)

# fail fast if pkgs missing (better for reproducibility on GitHub)
check_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Please install the following packages before running:\n  ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}
check_packages(required_pkgs)

# ---- Helper utilities ---------------------------------------------------------

# Midpoints of bin edges
midpoints <- function(breaks) breaks[-1] - diff(breaks) / 2

# ---- Step-function estimators -------------------------------------------------
# All step estimators return densities (intensity per unit of bin size) when
# `intensity = TRUE`; otherwise they normalize by total weights to return a PDF.

#' Univariate step-function estimator
#' @param data_vec vector of data (e.g., distances or times)
#' @param weights  vector of nonnegative weights (same length as data_vec)
#' @param breaks   strictly increasing vector of bin boundaries
#' @param intensity logical; TRUE returns intensity per unit width
mystepfun_uni <- function(data_vec, weights, breaks, intensity = TRUE) {
  keep_idx <- which(data_vec >= breaks[1] & data_vec <= breaks[length(breaks)])
  if (!length(keep_idx)) {
    return(list(funvals = numeric(length(breaks) - 1L), breaks = breaks))
  }
  x <- data_vec[keep_idx]
  w <- weights[keep_idx]
  
  dat_idx <- findInterval(x, breaks, rightmost.closed = TRUE)
  vals <- numeric(length(breaks) - 1L)
  
  for (i in seq_len(length(breaks) - 1L)) {
    bin_w <- sum(w[dat_idx == i])
    vals[i] <- bin_w / diff(breaks[i:(i + 1)])
  }
  if (!intensity) vals <- vals / sum(w)
  list(funvals = vals, breaks = breaks)
}

#' Bivariate step-function estimator on a rectangle grid
#' @param data_mat  N x 2 matrix (first margin breaks1, second margin breaks2)
#' @param weights   length-N weights
#' @param breaks1   x-axis bin boundaries (increasing)
#' @param breaks2   y-axis bin boundaries (increasing)
#' @param intensity TRUE for intensity; FALSE for PDF over area
mystepfun_bi <- function(data_mat, weights, breaks1, breaks2, intensity = TRUE) {
  keep_idx <- which(
    data_mat[, 1] >= breaks1[1] & data_mat[, 1] <= breaks1[length(breaks1)] &
      data_mat[, 2] >= breaks2[1] & data_mat[, 2] <= breaks2[length(breaks2)]
  )
  if (!length(keep_idx)) {
    funval_mat <- matrix(0, nrow = length(breaks2) - 1L, ncol = length(breaks1) - 1L)
    return(list(funvals = funval_mat, breaks1 = breaks1, breaks2 = breaks2))
  }
  
  x <- data_mat[keep_idx, , drop = FALSE]
  w <- weights[keep_idx]
  
  idx1 <- findInterval(x[, 1], breaks1, rightmost.closed = TRUE)
  idx2 <- findInterval(x[, 2], breaks2, rightmost.closed = TRUE)
  
  nrow_mat <- length(breaks2) - 1L
  ncol_mat <- length(breaks1) - 1L
  
  lin_idx <- (idx1 - 1L) * nrow_mat + idx2
  w_sum <- tapply(w, lin_idx, sum, default = 0)
  
  funval_mat <- matrix(0, nrow = nrow_mat, ncol = ncol_mat)
  funval_mat[as.integer(names(w_sum))] <- w_sum
  
  # normalize by bin area
  area_mat <- outer(diff(breaks2), diff(breaks1), "*")
  funval_mat <- funval_mat / area_mat
  if (!intensity) funval_mat <- funval_mat / sum(w)
  
  list(funvals = funval_mat, breaks1 = breaks1, breaks2 = breaks2)
}

#' Trivariate step-function estimator on a 3D grid (y,x,t ordering)
#' @param loc_mat   N x 2 matrix of (x,y)
#' @param time_vec  length-N numeric or POSIXct time
#' @param weights   length-N weights
#' @param breaks1   x bin edges
#' @param breaks2   y bin edges
#' @param breaks3   t bin edges
#' @param intensity TRUE for intensity; FALSE for PDF over volume
mystepfun_tri <- function(loc_mat, time_vec, weights,
                          breaks1, breaks2, breaks3, intensity = TRUE) {
  keep_idx <- which(
    loc_mat[, 1] >= breaks1[1] & loc_mat[, 1] <= breaks1[length(breaks1)] &
      loc_mat[, 2] >= breaks2[1] & loc_mat[, 2] <= breaks2[length(breaks2)] &
      time_vec      >= breaks3[1] & time_vec      <= breaks3[length(breaks3)]
  )
  nrow_arr <- length(breaks2) - 1L  # rows  (2nd margin)
  ncol_arr <- length(breaks1) - 1L  # cols  (1st margin)
  nlay_arr <- length(breaks3) - 1L  # layers(3rd margin)
  
  if (!length(keep_idx)) {
    funval_arr <- array(0, dim = c(nrow_arr, ncol_arr, nlay_arr))
    return(list(funvals = funval_arr, breaks1 = breaks1, breaks2 = breaks2, breaks3 = breaks3))
  }
  
  x   <- loc_mat[keep_idx, , drop = FALSE]
  tsv <- time_vec[keep_idx]
  w   <- weights[keep_idx]
  
  idx1 <- findInterval(x[, 1], breaks1, rightmost.closed = TRUE)
  idx2 <- findInterval(x[, 2], breaks2, rightmost.closed = TRUE)
  idx3 <- findInterval(tsv,    breaks3, rightmost.closed = TRUE)
  
  lin_idx <- idx2 + (idx1 - 1L) * nrow_arr + (idx3 - 1L) * nrow_arr * ncol_arr
  w_sum   <- tapply(w, lin_idx, sum, default = 0)
  
  funval_arr <- array(0, dim = c(nrow_arr, ncol_arr, nlay_arr))
  funval_arr[as.integer(names(w_sum))] <- w_sum
  
  # normalize by bin volume: Δy × Δx × Δt
  base_area <- outer(diff(breaks2), diff(breaks1), "*")
  vol_arr   <- array(base_area, dim = c(nrow_arr, ncol_arr, nlay_arr))
  vol_arr   <- sweep(vol_arr, 3, diff(breaks3), "*")
  
  funval_arr <- funval_arr / vol_arr
  if (!intensity) funval_arr <- funval_arr / sum(w)
  
  list(funvals = funval_arr, breaks1 = breaks1, breaks2 = breaks2, breaks3 = breaks3)
}

# ---- Data loading -------------------------------------------------------------

data_path_events   <- "~/path/to/lightning/data.RDS"
output_path        <- "~/path/to/output/folder"

lightning2022         <- readRDS(data_path_events)

# Count events per window in the metadata table
size_vec <- numeric(nrow(lightning2022$info))
for (i in seq_len(nrow(lightning2022$info))) {
  dt_tmp <- dplyr::filter(
    lightning2022$data,
    일시 >= lightning2022$info$start[i],
    일시 <= lightning2022$info$end[i]
  )
  size_vec[i] <- nrow(dt_tmp)
}

# ---- Model settings -----------------------------------------------------------

bgr_est_type  <- "nonconst"  # "const" or "nonconst"
sden_est_type <- "1d"        # "1d" (radial) or "2d" (x,y)
bw_mu         <- 1           # (placeholder; not used in current step-function impl.)

# Spatial domain (convex hull polygon) and its area
w0           <- lightning2022$chull
domain_sarea <- sf::st_area(w0) |> as.numeric()

lon_lim <- c(116, 138)
lat_lim <- c(28, 46)

# Background grid
resol_bgr    <- 0.5
bgr_breaks1  <- seq(lon_lim[1], lon_lim[2], by = resol_bgr)
bgr_breaks2  <- seq(lat_lim[1], lat_lim[2], by = resol_bgr)
plane_grid   <- cbind(rep(bgr_breaks1, each = length(bgr_breaks2)),
                      rep(bgr_breaks2,       length(bgr_breaks1)))

# Mark grid points falling inside polygon w0 (to correct for outer bins)
within_vec <- numeric(nrow(plane_grid))
for (i in seq_len(nrow(plane_grid))) {
  is_within <- sf::st_within(sf::st_point(plane_grid[i, ]), w0)
  if (length(is_within[[1]]) > 0) within_vec[i] <- 1
}

# Distance / time binning
blk          <- 0.005 * 2^(0:13)
ds_breaks    <- c(0, blk)
ds_bin_mid   <- midpoints(ds_breaks)

ll_max       <- 29
dxy_breaks1  <- c(rev(-blk), 0, blk)
dxy_breaks1  <- dxy_breaks1[dxy_breaks1 >= -ll_max & dxy_breaks1 <= ll_max]
dxy_breaks1  <- c(-ll_max, dxy_breaks1, ll_max)
dxy_breaks2  <- dxy_breaks1

dt_breaks    <- c(0, 1e-5 * 4^(4:14), 300000)

# ---- Parallel setup -----------------------------------------------------------

core_num <- 2
cl <- parallel::makeCluster(core_num)
doParallel::registerDoParallel(cl)
on.exit(parallel::stopCluster(cl), add = TRUE)

# ---- Main loop over windows ---------------------------------------------------

result_total <- foreach::foreach(
  data_idx = 1:34,
  .packages = c("spatstat", "lubridate", "dplyr", "sf")
) %dopar% {
  
  tryCatch({
    time_start <- Sys.time()
    
    # Subset events for the current window
    dt_tmp <- dplyr::filter(
      lightning2022$data,
      시분초 >= lightning2022$info$start[data_idx],
      시분초 <= lightning2022$info$end[data_idx]
    )
    
    # Construct training catalog (rename Korean columns to English)
    catalog <- dplyr::tibble(
      time = dt_tmp$시분초,
      lon  = dt_tmp$경도,
      lat  = dt_tmp$위도
    ) |> dplyr::arrange(time)
    
    period_begin   <- as.POSIXct(lightning2022$info$start[data_idx], format = "%Y-%m-%d")
    period_end     <- as.POSIXct(lightning2022$info$end[data_idx],   format = "%Y-%m-%d")
    domain_tlength <- as.numeric(difftime(period_end, period_begin, units = "secs"))
    
    # Background time bins (hourly)
    bgr_breaks3 <- seq(period_begin, period_end, by = 3600)
    
    # ---- Initialize EM --------------------------------------------------------
    
    n <- nrow(catalog)
    if (n < 2L) {
      return(list(error = TRUE, message = "Not enough events in this window."))
    }
    
    iter_num <- 200
    
    # Initial triggering probability matrix (lower triangle only matters)
    p_mat <- matrix(0, nrow = n, ncol = n)
    for (i in seq_len(n)) p_mat[i, 1:i] <- 1 / i
    
    # Pairwise lags (lower triangle vectorization)
    time_diff_mat <- matrix(0, n, n)
    for (i in 2:n) time_diff_mat[i, 1:(i - 1)] <- as.numeric(difftime(catalog$time[i], catalog$time[1:(i - 1)], "secs"))
    time_diff_vec <- time_diff_mat[lower.tri(time_diff_mat)]
    
    lon_diff_mat <- matrix(0, n, n)
    for (i in 2:n) lon_diff_mat[i, 1:(i - 1)] <- catalog$lon[i] - catalog$lon[1:(i - 1)]
    lon_diff_vec <- lon_diff_mat[lower.tri(lon_diff_mat)]
    
    lat_diff_mat <- matrix(0, n, n)
    for (i in 2:n) lat_diff_mat[i, 1:(i - 1)] <- catalog$lat[i] - catalog$lat[1:(i - 1)]
    lat_diff_vec <- lat_diff_mat[lower.tri(lat_diff_mat)]
    
    lonlat_diff <- cbind(lon_diff_vec, lat_diff_vec)
    
    diff_vec <- numeric(iter_num)
    t0 <- Sys.time()
    
    for (iter in seq_len(iter_num)) {
      # ---- M-step: background rate ------------------------------------------
      if (bgr_est_type == "const") {
        # Homogeneous background: total expected mains / (space x time)
        bgr_rate_vec <- rep(sum(diag(p_mat)) / (domain_sarea * domain_tlength), n)
        bgr_rate <- NULL
      } else {
        # Inhomogeneous background over an s-t grid
        bgr_rate <- mystepfun_tri(
          loc_mat  = cbind(catalog$lon, catalog$lat),
          time_vec = catalog$time,
          weights  = diag(p_mat),
          breaks1  = bgr_breaks1,
          breaks2  = bgr_breaks2,
          breaks3  = bgr_breaks3,
          intensity = TRUE
        )
        # Correct for grid points falling outside spatial polygon
        bgr_rate$funvals <- bgr_rate$funvals * length(within_vec) / sum(within_vec == 1)
        
        # Look up background rate at each event
        xidx <- findInterval(catalog$lon,  bgr_breaks1, rightmost.closed = TRUE)
        yidx <- findInterval(catalog$lat,  bgr_breaks2, rightmost.closed = TRUE)
        tidx <- findInterval(catalog$time, bgr_breaks3, rightmost.closed = TRUE)
        
        bgr_rate_vec <- numeric(n)
        for (j in seq_len(n)) bgr_rate_vec[j] <- bgr_rate$funvals[yidx[j], xidx[j], tidx[j]]
      }
      
      # ---- M-step: productivity (scalar) ------------------------------------
      p_offdiag <- p_mat; diag(p_offdiag) <- 0
      alpha_hat <- sum(p_offdiag) / (n - 1)  # average expected offspring per parent
      alpha_vec <- rep(alpha_hat, n)
      
      # ---- M-step: triggering kernel density --------------------------------
      pij <- p_mat[lower.tri(p_mat)]
      pij[pij < 0] <- 0
      
      if (sden_est_type == "2d") {
        xyden <- mystepfun_bi(lonlat_diff, weights = pij,
                              breaks1 = dxy_breaks1, breaks2 = dxy_breaks2,
                              intensity = FALSE)
        dxidx <- findInterval(lonlat_diff[, 1], dxy_breaks1, rightmost.closed = TRUE)
        dyidx <- findInterval(lonlat_diff[, 2], dxy_breaks2, rightmost.closed = TRUE)
      } else {
        dist_vec <- sqrt(lonlat_diff[, 1]^2 + lonlat_diff[, 2]^2)
        xyden <- mystepfun_uni(dist_vec, weights = pij, breaks = ds_breaks, intensity = TRUE)
        dsidx <- findInterval(dist_vec, ds_breaks, rightmost.closed = TRUE)
      }
      
      tden   <- mystepfun_uni(time_diff_vec, weights = pij, breaks = dt_breaks, intensity = TRUE)
      dtidx  <- findInterval(time_diff_vec, dt_breaks, rightmost.closed = TRUE)
      
      dens_vec <- numeric(nrow(lonlat_diff))
      for (j in seq_along(dens_vec)) {
        if (sden_est_type == "2d") {
          dens_vec[j] <- xyden$funvals[dyidx[j], dxidx[j]] * tden$funvals[dtidx[j]]
        } else {
          # radial × temporal density; divide by 2πr to convert radial shell to areal density
          dens_vec[j] <- xyden$funvals[dsidx[j]] * tden$funvals[dtidx[j]] / (2 * pi * ds_bin_mid[dsidx[j]])
        }
      }
      trig_st <- matrix(0, n, n); trig_st[lower.tri(trig_st)] <- dens_vec
      
      # ---- E-step: update p_mat ---------------------------------------------
      p_diff <- 0
      for (i in 2:n) {
        trig_vec <- alpha_vec[1:(i - 1)] * trig_st[i, 1:(i - 1)]
        denom    <- bgr_rate_vec[i] + sum(trig_vec)
        
        main_contrib <- bgr_rate_vec[i] / denom
        p_diff <- max(p_diff, abs(p_mat[i, i] - main_contrib))
        p_mat[i, i] <- main_contrib
        
        after_contrib <- trig_vec / denom
        p_diff <- max(p_diff, max(abs(p_mat[i, 1:(i - 1)] - after_contrib)))
        p_mat[i, 1:(i - 1)] <- after_contrib
      }
      
      diff_vec[iter] <- p_diff
      cat(sprintf("Iteration %d: maxΔ = %.4f | elapsed = %s\n",
                  iter, p_diff, format(Sys.time() - t0)))
      if (is.nan(p_diff) || p_diff < 0.01) break
    }
    
    time_end <- Sys.time()
    
    list(
      data              = catalog,
      lon_lim           = lon_lim,
      lat_lim           = lat_lim,
      p_mat             = p_mat,
      mainshock_num     = sum(diag(p_mat)),
      evt_prod          = colSums(p_mat),
      bgr_rate_vec      = bgr_rate_vec,
      bgr_rate          = if (exists("bgr_rate")) bgr_rate else NULL,
      offspring_productivity = alpha_vec,
      p_mat_diff        = diff_vec[diff_vec > 0],
      xyden             = xyden,
      tden              = tden,
      iter              = iter,
      bgr_est_type      = bgr_est_type,
      calcTime          = time_end - time_start,
      domain_sarea      = domain_sarea,
      domain_tlength    = domain_tlength,
      w0                = w0,
      period_begin      = period_begin,
      period_end        = period_end
    )
    
  }, error = function(e) {
    list(error = TRUE, message = e$message)
  })
}

# ---- Save results -------------------------------------------------------------

dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(
  output_path,
  sprintf("nonpar_results2022_st_inhomo_bgr_sden_%s.rds", sden_est_type)
)
saveRDS(result_total, outfile)
cat("Saved results to: ", outfile, "\n")
