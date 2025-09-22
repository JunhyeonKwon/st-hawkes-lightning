# Super-thinning residuals for nonparametric ST-Hawkes fits
# ---------------------------------------------------------
# How to run:
#   - Install packages listed in `required_pkgs`
#   - Set data/output paths in the CONFIG section
#   - Source this file and run from R / Rscript
#
# Notes:
#   - k defaults to average rate N / |S×T| (or set a target residual count).

# ---- Dependencies -------------------------------------------------------------

required_pkgs <- c("stpp", "sf", "spatstat.geom", "spatstat.random", "dplyr", "fields")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Please install: ", paste(missing, collapse = ", "), call. = FALSE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ---- CONFIG ------------------------------------------------------------------

# path to fitted objects (list of length 34; each element is your `result_list`)
fits_path <- "~/path/to/nonpar_results2022_st_inhomo_bgr_sden_1d.rds"

# where to save outputs
out_dir    <- "~/path/to/output/folder"
thin_out   <- file.path(out_dir, "thinning_results.rds")
pcf_out    <- file.path(out_dir, "pcf_results.rds")

# super-thinning rate choice:
#   (a) leave k = NULL to use N / |S×T|
#   (b) or set k_target_count to an integer M; we'll use k = M / |S×T|
k_default        <- NULL
k_target_count   <- NULL   # e.g., 2000  (set to NULL to disable)

# PCF evaluation lags: we will compute from each fit's step bins (midpoints)
# and use first 8 temporal lags as in your script.

# ---- Small helpers ------------------------------------------------------------

midpoints <- function(breaks) breaks[-1] - diff(breaks) / 2

# Convert an sf POLYGON/MULTIPOLYGON to spatstat's owin
sf_to_owin <- function(sfg) {
  # sfg: sf geometry column (e.g., fit$w0$geometry)
  # works for single polygon; multipolygons are combined
  poly <- sf::st_geometry(sf::st_as_sf(sfg))
  poly <- sf::st_cast(poly, "POLYGON")
  coords_list <- lapply(poly, function(pg) sf::st_coordinates(pg)[, 1:2, drop = FALSE])
  spatstat.geom::owin(poly = lapply(coords_list, function(x) list(x = x[,1], y = x[,2])))
}

# Build 3D stpp object from vectors; time must be numeric (seconds)
as_stpp_3dpoints <- function(x, y, t) {
  pts <- cbind(x = x, y = y, t = t)
  pts <- pts[order(pts[, 3]), , drop = FALSE]
  class(pts) <- "stpp"
  pts
}

# Conditional intensity at either event locations or supplied newdata
compute_ci <- function(fit, newdata = NULL) {
  data_df        <- fit$data
  bgr_rate       <- fit$bgr_rate
  alpha          <- fit$offspring_productivity[1]
  tden           <- fit$tden
  xyden          <- fit$xyden     # radial step density when sden_est_type == "1d"
  domain_tlength <- fit$domain_tlength
  
  # 1) evaluation coordinates
  if (is.null(newdata)) {
    x <- data_df$lon; y <- data_df$lat; tt <- data_df$time
  } else {
    x <- newdata$lon; y <- newdata$lat; tt <- newdata$time
  }
  n <- length(tt)
  
  # 2) pairwise lags (vectorized over lower triangle)
  low <- lower.tri(matrix(0, n, n), diag = FALSE)
  dt  <- outer(tt, tt, "-")[low]                    # positive time lags
  dx  <- outer(x , x , "-")[low]
  dy  <- outer(y , y , "-")[low]
  ds  <- sqrt(dx^2 + dy^2)
  
  dti <- findInterval(dt, tden$breaks, rightmost.closed = TRUE)
  dsi <- findInterval(ds, xyden$breaks, rightmost.closed = TRUE)
  
  # radial × temporal; divide by 2πr to convert shell to areal density
  kern_vec <- tden$funvals[dti] * xyden$funvals[dsi] / (2 * pi * pmax(ds, .Machine$double.eps))
  
  trig_mat <- matrix(0, n, n)
  trig_mat[low] <- kern_vec
  
  # background at each point (take care of indices at boundaries)
  loni <- pmax(1L, pmin(findInterval(x,  bgr_rate$breaks1, rightmost.closed = TRUE),
                        length(bgr_rate$breaks1) - 1L))
  lati <- pmax(1L, pmin(findInterval(y,  bgr_rate$breaks2, rightmost.closed = TRUE),
                        length(bgr_rate$breaks2) - 1L))
  tti  <- pmax(1L, pmin(findInterval(tt, bgr_rate$breaks3, rightmost.closed = TRUE),
                        length(bgr_rate$breaks3) - 1L))
  
  bgr_rate_vec <- bgr_rate$funvals[cbind(lati, loni, tti)] / domain_tlength
  
  # 3) conditional intensity
  lambda <- bgr_rate_vec + alpha * rowSums(trig_mat)
  lambda
}

# Super-thinning for a single fitted window
superthin <- function(fit, k = NULL, k_target_count = NULL) {
  data_df        <- fit$data
  domain_sarea   <- fit$domain_sarea
  domain_tlength <- fit$domain_tlength
  
  # choose k
  if (!is.null(k_target_count)) k <- k_target_count / (domain_sarea * domain_tlength)
  if (is.null(k))               k <- nrow(data_df) / (domain_sarea * domain_tlength)
  if (k <= 0) stop("k must be > 0")
  
  # intensity at the observed events
  lambda_obs <- compute_ci(fit)
  
  # (A) simulate a homogeneous Poisson field with rate k over S×T
  expected_resid <- k * domain_sarea * domain_tlength
  n_add          <- rpois(1, expected_resid)
  
  residual_add <- data.frame()
  if (n_add > 0) {
    # spatial: homogeneous k * T over S
    win_owin <- sf_to_owin(fit$w0$geometry)
    new_pts  <- spatstat.random::rpoispp(lambda = k * domain_tlength, win = win_owin)
    
    # time: uniform over [T0, T1]
    t0 <- as.numeric(fit$period_begin)
    t1 <- as.numeric(fit$period_end)
    new_t <- as.POSIXct(runif(new_pts$n, t0, t1), origin = "1970-01-01", tz = attr(fit$period_begin, "tzone"))
    
    residual_add <- dplyr::tibble(lon = new_pts$x, lat = new_pts$y, time = new_t)
    
    # keep added points with probability (k - lambda(new))/k where lambda(new) computed at added points
    lambda_new <- compute_ci(fit, newdata = residual_add)
    keep_add   <- runif(length(lambda_new)) <= pmax(0, (k - lambda_new) / k)
    residual_add <- residual_add[keep_add, , drop = FALSE]
  }
  
  # (B) thin the observed points where lambda > k
  keep_obs1 <- data_df[lambda_obs <= k, c("경도","위도","시분초")]
  if (nrow(keep_obs1)) names(keep_obs1) <- c("lon","lat","time")
  
  thin_df   <- data_df[lambda_obs >  k, c("경도","위도","시분초")]
  thin_lmb  <- lambda_obs[lambda_obs > k]
  names(thin_df) <- c("lon","lat","time")
  
  keep_prob <- k / thin_lmb
  keep_idx  <- runif(length(keep_prob)) <= keep_prob
  keep_obs2 <- thin_df[keep_idx, , drop = FALSE]
  deleted   <- thin_df[!keep_idx, , drop = FALSE]
  
  residuals_df <- rbind(keep_obs1, keep_obs2, residual_add)
  residuals_df <- residuals_df[order(residuals_df$time), , drop = FALSE]
  
  list(
    k          = k,
    X          = data_df,
    residuals  = residuals_df,
    super      = residual_add,
    kept_obs   = rbind(keep_obs1, keep_obs2),
    deleted    = deleted
  )
}

# ---- Load fits ---------------------------------------------------------------

fits <- readRDS(fits_path)

# quick sanity check
stopifnot(is.list(fits), length(fits) >= 1)

# ---- Run super-thinning over all windows -------------------------------------

thinning_results <- vector("list", length(fits))
for (i in seq_along(fits)) {
  t1 <- Sys.time()
  cat(sprintf("Super-thinning for dataset %d started.\n", i))
  
  thinning_results[[i]] <- superthin(
    fits[[i]],
    k = k_default,
    k_target_count = k_target_count
  )
  
  cat(sprintf("Super-thinning for dataset %d ended. Elapsed: %s\n\n",
              i, format(Sys.time() - t1)))
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(thinning_results, thin_out)
cat("Saved thinning results to: ", thin_out, "\n")

# ---- PCF estimation on residuals ---------------------------------------------

pcf_results <- vector("list", length(fits))

for (i in seq_along(fits)) {
  t1 <- Sys.time(); cat(sprintf("PCFhat for dataset %d started.\n", i))
  
  fit_i   <- fits[[i]]
  thin_i  <- thinning_results[[i]]$residuals
  
  # build stpp object (time as numeric seconds)
  stpp_i <- as_stpp_3dpoints(
    x = thin_i$lon,
    y = thin_i$lat,
    t = as.numeric(thin_i$time)
  )
  
  # spatial/temporal lag grids from fit bins (midpoints)
  # spatial uses the radial step bins stored in fit$xyden$breaks
  dist_grid  <- midpoints(fit_i$xyden$breaks)
  time_grid  <- midpoints(fit_i$tden$breaks)[1:8]
  
  # spatial region boundary as matrix of coordinates
  s_region <- sf::st_coordinates(fit_i$w0$geometry)
  
  pcf_i <- stpp::PCFhat(
    xyt       = stpp_i,
    s.region  = s_region,
    t.region  = as.numeric(c(fit_i$period_begin, fit_i$period_end)),
    correction = "border",
    dist       = dist_grid,
    times      = time_grid
  )
  pcf_results[[i]] <- pcf_i
  
  # optional quick plot
  # stpp::plotPCF(pcf_i, which = "border", main = i)
  
  cat(sprintf("PCFhat for dataset %d ended. Elapsed: %s\n\n",
              i, format(Sys.time() - t1)))
}

saveRDS(pcf_results, pcf_out)
cat("Saved PCF results to: ", pcf_out, "\n")
