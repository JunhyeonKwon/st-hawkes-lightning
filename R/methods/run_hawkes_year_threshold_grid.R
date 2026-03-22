rm(list = ls())
gc()


# Load all required packages up front so runtime failures happen early.
load_required_packages <- function(packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0L) {
    stop(
      "Install the required packages before running this script: ",
      paste(missing_packages, collapse = ", ")
    )
  }
  
  invisible(lapply(packages, library, character.only = TRUE))
}


# Estimate a 1D step density on pre-defined bins.
mystepfun_uni <- function(data_vec, weights, breaks) {
  keep_idx <- which(data_vec >= breaks[1] & data_vec <= breaks[length(breaks)])
  data_vec <- data_vec[keep_idx]
  w_tmp <- weights[keep_idx]
  
  dat_idx <- findInterval(data_vec, breaks, rightmost.closed = TRUE)
  funval_vec <- numeric(length(breaks) - 1L)
  
  weight_sum <- tapply(w_tmp, dat_idx, sum, default = 0)
  if (length(weight_sum) > 0L) {
    funval_vec[as.integer(names(weight_sum))] <- as.numeric(weight_sum)
  }
  
  funval_vec <- funval_vec / diff(breaks) / sum(w_tmp)
  
  list(funvals = funval_vec, breaks = breaks)
}

mystepfun_bi <- function(data_mat, time_vec, weights, breaks1, breaks2, intensity = TRUE) {
  # Estimate a 2D step density using the same row/column convention as the original script.
  keep_idx <- which(
    data_mat[, 1] >= breaks1[1] & data_mat[, 1] <= breaks1[length(breaks1)] &
      data_mat[, 2] >= breaks2[1] & data_mat[, 2] <= breaks2[length(breaks2)]
  )
  
  data_mat <- data_mat[keep_idx, , drop = FALSE]
  w_tmp <- weights[keep_idx]
  
  dat_idx1 <- findInterval(data_mat[, 1], breaks1, rightmost.closed = TRUE)
  dat_idx2 <- findInterval(data_mat[, 2], breaks2, rightmost.closed = TRUE)
  
  nrow_mat <- length(breaks2) - 1L
  ncol_mat <- length(breaks1) - 1L
  linear_indices <- (dat_idx1 - 1L) * nrow_mat + dat_idx2
  
  weight_sum <- tapply(w_tmp, linear_indices, sum, default = 0)
  funval_mat <- matrix(0, nrow = nrow_mat, ncol = ncol_mat)
  
  if (length(weight_sum) > 0L) {
    funval_mat[as.integer(names(weight_sum))] <- as.numeric(weight_sum)
  }
  
  funval_mat <- funval_mat / outer(diff(breaks2), diff(breaks1), "*")
  if (!intensity) {
    funval_mat <- funval_mat / sum(w_tmp)
  }
  
  list(funvals = funval_mat, breaks1 = breaks1, breaks2 = breaks2)
}

mystepfun_tri <- function(loc_mat, time_vec, weights, breaks1, breaks2, breaks3, intensity = TRUE) {
  # Estimate a 3D step intensity over longitude, latitude, and time bins.
  keep_idx <- which(
    loc_mat[, 1] >= breaks1[1] & loc_mat[, 1] <= breaks1[length(breaks1)] &
      loc_mat[, 2] >= breaks2[1] & loc_mat[, 2] <= breaks2[length(breaks2)] &
      time_vec >= breaks3[1] & time_vec <= breaks3[length(breaks3)]
  )
  
  if (length(keep_idx) == 0L) {
    funval_arr <- array(
      0,
      dim = c(length(breaks2) - 1L, length(breaks1) - 1L, length(breaks3) - 1L)
    )
    
    return(list(funvals = funval_arr, breaks1 = breaks1, breaks2 = breaks2, breaks3 = breaks3))
  }
  
  loc_mat <- loc_mat[keep_idx, , drop = FALSE]
  time_vec <- time_vec[keep_idx]
  w_tmp <- weights[keep_idx]
  
  dat_idx1 <- findInterval(loc_mat[, 1], breaks1, rightmost.closed = TRUE)
  dat_idx2 <- findInterval(loc_mat[, 2], breaks2, rightmost.closed = TRUE)
  dat_idx3 <- findInterval(time_vec, breaks3, rightmost.closed = TRUE)
  
  nrow_arr <- length(breaks2) - 1L
  ncol_arr <- length(breaks1) - 1L
  nlay_arr <- length(breaks3) - 1L
  
  linear_indices <- dat_idx2 +
    (dat_idx1 - 1L) * nrow_arr +
    (dat_idx3 - 1L) * nrow_arr * ncol_arr
  
  weight_sum <- tapply(w_tmp, linear_indices, sum, default = 0)
  funval_arr <- array(0, dim = c(nrow_arr, ncol_arr, nlay_arr))
  
  if (length(weight_sum) > 0L) {
    funval_arr[as.integer(names(weight_sum))] <- as.numeric(weight_sum)
  }
  
  base_area <- outer(diff(breaks2), diff(breaks1), "*")
  vol_arr <- array(base_area, dim = c(nrow_arr, ncol_arr, nlay_arr))
  vol_arr <- sweep(vol_arr, 3, diff(breaks3), "*")
  
  funval_arr <- funval_arr / vol_arr
  if (!intensity) {
    funval_arr <- funval_arr / sum(w_tmp)
  }
  
  list(funvals = funval_arr, breaks1 = breaks1, breaks2 = breaks2, breaks3 = breaks3)
}

pairwise_lag_vector <- function(values) {
  # Store pairwise lags directly in lower-triangular order to avoid full lag matrices.
  n_values <- length(values)
  
  if (n_values < 2L) {
    return(numeric(0))
  }
  
  lag_vec <- numeric(n_values * (n_values - 1L) / 2L)
  start_idx <- 1L
  
  for (j in seq_len(n_values - 1L)) {
    end_idx <- start_idx + (n_values - j) - 1L
    lag_vec[start_idx:end_idx] <- values[(j + 1L):n_values] - values[j]
    start_idx <- end_idx + 1L
  }
  
  lag_vec
}

initialize_probability_matrix <- function(n_events) {
  # Keep the original initialization: row i gets probability 1 / i over its admissible parents.
  prob_mat <- matrix(0, nrow = n_events, ncol = n_events)
  for (i in seq_len(n_events)) {
    prob_mat[i, seq_len(i)] <- 1 / i
  }
  
  prob_mat
}

calculate_event_sizes <- function(lightning_year) {
  # Count the number of strikes in each analysis window.
  info <- lightning_year$info
  event_time <- lightning_year$data$일시
  
  vapply(
    seq_len(nrow(info)),
    function(i) {
      sum(event_time >= info$start[i] & event_time <= info$end[i])
    },
    numeric(1)
  )
}

create_background_setup <- function(w0, lon_lim, lat_lim, resol_bgr) {
  # Build the spatial grid once per year and precompute the in-domain correction factor.
  bgr_breaks1 <- seq(lon_lim[1], lon_lim[2], by = resol_bgr)
  bgr_breaks2 <- seq(lat_lim[1], lat_lim[2], by = resol_bgr)
  
  plane_grid <- cbind(
    rep(bgr_breaks1, each = length(bgr_breaks2)),
    rep(bgr_breaks2, length(bgr_breaks1))
  )
  
  grid_points <- sf::st_as_sf(
    data.frame(lon = plane_grid[, 1], lat = plane_grid[, 2]),
    coords = c("lon", "lat"),
    crs = sf::st_crs(w0)
  )
  
  within_vec <- as.integer(lengths(sf::st_within(grid_points, w0)) > 0L)
  
  list(
    bgr_breaks1 = bgr_breaks1,
    bgr_breaks2 = bgr_breaks2,
    within_vec = within_vec,
    within_ratio = length(within_vec) / sum(within_vec == 1L)
  )
}

make_batch_plan <- function(size_vec, size_breaks, core_by_group) {
  # Allocate cores by event-count thresholds to keep memory use predictable across years.
  ordered_indices <- order(size_vec)
  batch_plan <- list()
  batch_id <- 1L
  
  if (length(ordered_indices) == 0L) {
    return(batch_plan)
  }
  
  for (group_id in seq_along(core_by_group)) {
    lower_bound <- size_breaks[group_id]
    upper_bound <- size_breaks[group_id + 1L]
    
    if (is.infinite(upper_bound)) {
      group_indices <- ordered_indices[size_vec[ordered_indices] > lower_bound]
      label <- sprintf("gt_%d", lower_bound)
    } else if (lower_bound == 0L) {
      group_indices <- ordered_indices[size_vec[ordered_indices] <= upper_bound]
      label <- sprintf("le_%d", upper_bound)
    } else {
      group_indices <- ordered_indices[
        size_vec[ordered_indices] > lower_bound & size_vec[ordered_indices] <= upper_bound
      ]
      label <- sprintf("%d_%d", lower_bound, upper_bound)
    }
    
    if (length(group_indices) > 0L) {
      batch_plan[[batch_id]] <- list(
        label = label,
        indices = group_indices,
        cores = core_by_group[group_id]
      )
      batch_id <- batch_id + 1L
    }
  }
  
  batch_plan
}


# Run the EM-like estimation for one event window while preserving the original calculations.
estimate_event <- function(data_idx, lightning_year, year, w0, domain_sarea, background_setup, config) {
  tryCatch(
    {
      time_start <- Sys.time()
      
      interval_start <- lightning_year$info$start[data_idx]
      interval_end <- lightning_year$info$end[data_idx]
      
      keep_idx <- lightning_year$data$시분초 >= interval_start &
        lightning_year$data$시분초 <= interval_end
      data_tmp <- lightning_year$data[keep_idx, , drop = FALSE]
      data_tmp$위도 <- round(data_tmp$위도, 3)
      
      period_begin <- as.POSIXct(lightning_year$info$start[data_idx], format = "%Y-%m-%d")
      period_end <- as.POSIXct(lightning_year$info$end[data_idx], format = "%Y-%m-%d")
      bgr_breaks3 <- seq(period_begin, period_end, by = 3600)
      
      catalog_training <- data.frame(
        time = data_tmp$시분초,
        lon = data_tmp$경도,
        lat = data_tmp$위도
      )
      catalog_training <- catalog_training[order(catalog_training$time), , drop = FALSE]
      
      domain_tlength <- as.numeric(difftime(period_end, period_begin, units = "secs"))
      n_catalog <- nrow(catalog_training)
      
      p_mat_old <- initialize_probability_matrix(n_catalog)
      
      time_num <- as.numeric(catalog_training$time)
      timediff_vec <- pairwise_lag_vector(time_num)
      londiff_vec <- pairwise_lag_vector(catalog_training$lon)
      latdiff_vec <- pairwise_lag_vector(catalog_training$lat)
      lonlatdiff_mat <- cbind(londiff_vec, latdiff_vec)
      
      diff_vec <- numeric(config$iter_num)
      
      for (iter in seq_len(config$iter_num)) {
        if (config$bgr_est_type == "const") {
          bgr_rate_vec <- rep(sum(diag(p_mat_old)) / (domain_sarea * domain_tlength), n_catalog)
          bgr_rate <- NULL
        } else {
          # Background rate is evaluated on fixed space-time bins and then mapped back to events.
          bgr_rate <- mystepfun_tri(
            loc_mat = cbind(catalog_training$lon, catalog_training$lat),
            time_vec = catalog_training$time,
            weights = diag(p_mat_old),
            breaks1 = background_setup$bgr_breaks1,
            breaks2 = background_setup$bgr_breaks2,
            breaks3 = bgr_breaks3,
            intensity = TRUE
          )
          
          bgr_rate$funvals <- bgr_rate$funvals * background_setup$within_ratio
          
          xindex <- findInterval(catalog_training$lon, background_setup$bgr_breaks1, rightmost.closed = TRUE)
          yindex <- findInterval(catalog_training$lat, background_setup$bgr_breaks2, rightmost.closed = TRUE)
          tindex <- findInterval(catalog_training$time, bgr_breaks3, rightmost.closed = TRUE)
          
          bgr_rate_vec <- bgr_rate$funvals[cbind(yindex, xindex, tindex)]
        }
        
        p_mat_old_offdiag <- p_mat_old
        diag(p_mat_old_offdiag) <- 0
        
        alpha_tmp <- sum(p_mat_old_offdiag) / (n_catalog - 1)
        alpha_vec <- rep(alpha_tmp, n_catalog)
        
        pij <- p_mat_old[lower.tri(p_mat_old)]
        pij <- ifelse(pij < 0, 0, pij)
        
        if (config$sden_est_type == "2d") {
          xyden <- mystepfun_bi(
            lonlatdiff_mat,
            time_vec = NULL,
            weights = pij,
            breaks1 = config$dxy_breaks1,
            breaks2 = config$dxy_breaks2,
            intensity = FALSE
          )
          
          dxindex <- findInterval(lonlatdiff_mat[, 1], config$dxy_breaks1, rightmost.closed = TRUE)
          dyindex <- findInterval(lonlatdiff_mat[, 2], config$dxy_breaks2, rightmost.closed = TRUE)
        } else {
          distdiff_vec <- sqrt(lonlatdiff_mat[, 1]^2 + lonlatdiff_mat[, 2]^2)
          xyden <- mystepfun_uni(distdiff_vec, weights = pij, breaks = config$ds_breaks)
          dsindex <- findInterval(distdiff_vec, config$ds_breaks, rightmost.closed = TRUE)
        }
        
        tden <- mystepfun_uni(timediff_vec, weights = pij, breaks = config$dt_breaks)
        dtindex <- findInterval(timediff_vec, config$dt_breaks, rightmost.closed = TRUE)
        
        if (config$sden_est_type == "2d") {
          dens_tmp <- xyden$funvals[cbind(dyindex, dxindex)] * tden$funvals[dtindex]
        } else {
          dens_tmp <- xyden$funvals[dsindex] * tden$funvals[dtindex] /
            (2 * pi * config$ds_bin_mid[dsindex])
        }
        
        trig_st <- matrix(0, nrow = n_catalog, ncol = n_catalog)
        trig_st[lower.tri(trig_st)] <- dens_tmp
        
        p_mat_diff <- 0
        for (i in 2:n_catalog) {
          # Update the diagonal background term and the triggering terms exactly as in v2.
          trig_vec <- alpha_vec[1:(i - 1)] * trig_st[i, 1:(i - 1)]
          
          main_contrib <- bgr_rate_vec[i] / (bgr_rate_vec[i] + sum(trig_vec))
          p_mat_diff <- max(c(p_mat_diff, abs(p_mat_old[i, i] - main_contrib)))
          p_mat_old[i, i] <- main_contrib
          
          if (is.nan(p_mat_diff)) {
            stop("main")
          }
          
          after_contrib <- trig_vec / (bgr_rate_vec[i] + sum(trig_vec))
          p_mat_diff <- max(c(p_mat_diff, abs(p_mat_old[i, 1:(i - 1)] - after_contrib)))
          p_mat_old[i, 1:(i - 1)] <- after_contrib
          
          if (is.nan(p_mat_diff)) {
            stop("after")
          }
        }
        
        diff_vec[iter] <- p_mat_diff
        
        if (p_mat_diff < 0.01) {
          break
        }
      }
      
      time_end <- Sys.time()
      
      list(
        data = catalog_training,
        lon_lim = config$lon_lim,
        lat_lim = config$lat_lim,
        p_mat = p_mat_old,
        mainshock_num = sum(diag(p_mat_old)),
        evt_prod = colSums(p_mat_old),
        bgr_rate_vec = bgr_rate_vec,
        bgr_rate = bgr_rate,
        offspring_productivity = alpha_vec,
        p_mat_diff = diff_vec[seq_len(iter)],
        xyden = xyden,
        tden = tden,
        iter = iter,
        bgr_est_type = config$bgr_est_type,
        calcTime = time_end - time_start,
        domain_sarea = domain_sarea,
        domain_tlength = domain_tlength,
        w0 = w0,
        period_begin = period_begin,
        period_end = period_end,
        year = year,
        data_idx = data_idx
      )
    },
    error = function(e) {
      list(error = TRUE, message = e$message, year = year, data_idx = data_idx)
    }
  )
}


run_parallel_batch <- function(lightning_year, year, batch, w0, domain_sarea, background_setup, config) {
  # Each batch uses its own cluster so large windows can run with a smaller memory footprint.
  available_cores <- parallel::detectCores(logical = TRUE)
  core_num <- batch$cores
  
  if (!is.na(available_cores)) {
    core_num <- min(core_num, available_cores)
  }
  
  cl <- parallel::makeCluster(core_num)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  worker_exports <- c(
    "mystepfun_uni",
    "mystepfun_bi",
    "mystepfun_tri",
    "pairwise_lag_vector",
    "initialize_probability_matrix",
    "estimate_event"
  )
  
  foreach::foreach(
    data_idx = batch$indices,
    .packages = c("spatstat", "sf"),
    .export = worker_exports
  ) %dopar% {
    estimate_event(
      data_idx = data_idx,
      lightning_year = lightning_year,
      year = year,
      w0 = w0,
      domain_sarea = domain_sarea,
      background_setup = background_setup,
      config = config
    )
  }
}


build_case_label <- function(year, threshold) {
  sprintf("%d | thsd=%d", year, threshold)
}

build_input_path <- function(data_dir, year, threshold) {
  file.path(data_dir, sprintf("Lightningevent_%d_%d.RDS", year, threshold))
}

build_output_path <- function(output_dir, year, threshold, sden_est_type) {
  file.path(
    output_dir,
    sprintf(
      "nonpar_results%d_thsd%d_st_inhomo_bgr_sden_%s_v3.rds",
      year,
      threshold,
      sden_est_type
    )
  )
}


run_year_analysis <- function(year, threshold, data_dir, output_dir, config) {
  # Run all windows for one year-threshold combination and save a distinct RDS file.
  case_label <- build_case_label(year, threshold)
  input_path <- build_input_path(data_dir, year, threshold)
  output_path <- build_output_path(output_dir, year, threshold, config$sden_est_type)
  
  message(sprintf("[%s] Loading data from %s", case_label, input_path))
  lightning_year <- readRDS(input_path)
  
  size_vec <- calculate_event_sizes(lightning_year)
  batch_plan <- make_batch_plan(
    size_vec = size_vec,
    size_breaks = config$size_breaks,
    core_by_group = config$core_by_group
  )
  
  w0 <- lightning_year$chull
  domain_sarea <- area(w0)
  background_setup <- create_background_setup(
    w0 = w0,
    lon_lim = config$lon_lim,
    lat_lim = config$lat_lim,
    resol_bgr = config$resol_bgr
  )
  
  result_total <- vector("list", nrow(lightning_year$info))
  
  for (batch_id in seq_along(batch_plan)) {
    batch <- batch_plan[[batch_id]]
    message(
      sprintf(
        "[%s] Batch %d/%d (%s): %d windows on %d cores.",
        case_label,
        batch_id,
        length(batch_plan),
        batch$label,
        length(batch$indices),
        batch$cores
      )
    )
    
    batch_result <- run_parallel_batch(
      lightning_year = lightning_year,
      year = year,
      batch = batch,
      w0 = w0,
      domain_sarea = domain_sarea,
      background_setup = background_setup,
      config = config
    )
    
    result_total[batch$indices] <- batch_result
    gc()
  }
  
  names(result_total) <- rownames(lightning_year$info)
  attr(result_total, "year") <- year
  attr(result_total, "threshold") <- threshold
  attr(result_total, "case_label") <- case_label
  attr(result_total, "size_vec") <- size_vec
  attr(result_total, "processing_order") <- order(size_vec)
  attr(result_total, "batch_plan") <- data.frame(
    batch = seq_along(batch_plan),
    label = vapply(batch_plan, `[[`, character(1), "label"),
    cores = vapply(batch_plan, `[[`, numeric(1), "cores"),
    n_windows = vapply(batch_plan, function(x) length(x$indices), numeric(1))
  )
  attr(result_total, "input_path") <- input_path
  attr(result_total, "output_path") <- output_path
  
  saveRDS(result_total, output_path)
  message(sprintf("[%s] Saved results to %s", case_label, output_path))
  
  invisible(output_path)
}


required_packages <- c("doParallel", "foreach", "spatstat", "sf")
load_required_packages(required_packages)

data_dir <- path.expand("~/path/to/data/folder")
output_dir <- file.path(data_dir, "more", "path", "to", "output", "subfolder")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Core analysis settings shared across all years.
config <- list(
  bgr_est_type = "nonconst",
  sden_est_type = "1d",
  lon_lim = c(116, 138),
  lat_lim = c(28, 46),
  resol_bgr = 0.5,
  ll_max = 29,
  iter_num = 200L,
  size_breaks = c(0L, 5000L, 15000L, Inf),
  core_by_group = c(12L, 2L, 1L)
)

blk <- 0.005 * 2^(0:13)
config$ds_breaks <- c(0, blk)
config$ds_bin_mid <- config$ds_breaks[-1] - diff(config$ds_breaks) / 2
config$dxy_breaks1 <- c(rev(-blk), 0, blk)
config$dxy_breaks1 <- config$dxy_breaks1[
  config$dxy_breaks1 >= -config$ll_max & config$dxy_breaks1 <= config$ll_max
]
config$dxy_breaks1 <- c(-config$ll_max, config$dxy_breaks1, config$ll_max)
config$dxy_breaks2 <- config$dxy_breaks1
config$dt_breaks <- c(0, 1e-5 * 4^(4:14), 300000)

years_to_analyze <- c(2022L, 2023L)
thresholds_to_analyze <- c(40L, 45L, 55L, 60L)
analysis_grid <- expand.grid(
  year = years_to_analyze,
  threshold = thresholds_to_analyze,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
analysis_grid <- analysis_grid[order(analysis_grid$year, analysis_grid$threshold), , drop = FALSE]

output_paths <- lapply(
  seq_len(nrow(analysis_grid)),
  function(i) {
    run_year_analysis(
      year = analysis_grid$year[i],
      threshold = analysis_grid$threshold[i],
      data_dir = data_dir,
      output_dir = output_dir,
      config = config
    )
  }
)

names(output_paths) <- sprintf(
  "%d_thsd%d",
  analysis_grid$year,
  analysis_grid$threshold
)

invisible(output_paths)
