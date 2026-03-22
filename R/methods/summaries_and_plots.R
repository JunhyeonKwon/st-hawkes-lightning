# Summaries & plots for nonparametric ST-Hawkes fits
# --------------------------------------------------
# What this script does
#   1) Loads the list of fitted objects (length = 34)
#   2) Computes quick summaries (N, mains fraction, offspring productivity)
#   3) Plots background rate slices (quilt plots)
#   4) Plots distributions of offspring productivity across datasets
#   5) Plots spatial & temporal triggering densities (boxplots over bins)


# ---- CONFIG -------------------------------------------------------------------

fits_path <- "~/path/to/nonpar_results2022_st_inhomo_bgr_sden_1d.rds"

out_dir   <- "~/path/to/output/folder"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Background-rate figure (quilt plots)
bgr_png   <- file.path(out_dir, "group_3_bgr.png")

# Offspring productivity boxplot (PDF)
alpha_pdf <- file.path(out_dir, "hawkes_fit_alpha.pdf")

# Spatial triggering density boxplot (PDF)
h_pdf     <- file.path(out_dir, "hawkes_fit_h.pdf")

# Temporal triggering density boxplot (PDF)
g_pdf     <- file.path(out_dir, "hawkes_fit_g1.pdf")

# Dataset index to showcase background rate tiles
example_idx <- 3L

# Background-rate time slice indices to display (will be clipped safely)
bgr_time_slices <- (1:4) * 4 - 2   # same intent as your original: 2, 6, 10, 14

# How many leading spatial/temporal bins to show in boxplots
n_spatial_bins  <- 7L
n_temporal_bins <- 7L

# ---- Dependencies -------------------------------------------------------------

pkgs <- c(
  "dplyr", "tidyr", "ggplot2",
  "fields", "maps", "viridis",
  "sf"
)
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Please install: ", paste(missing, collapse = ", "), call. = FALSE)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Helpers ------------------------------------------------------------------

midpoints <- function(breaks) breaks[-1] - diff(breaks) / 2

safe_idx <- function(idx, max_len) {
  idx <- as.integer(idx)
  idx[idx < 1L] <- 1L
  idx[idx > max_len] <- max_len
  unique(idx)
}

# ---- Load fits ---------------------------------------------------------------

fits <- readRDS(fits_path)
stopifnot(is.list(fits), length(fits) >= 1)

# ---- Quick summaries ----------------------------------------------------------

# number of events per dataset
n_events <- vapply(fits, function(x) nrow(x$data), integer(1))

# fraction of expected mains (diag p_mat) among events
mains_frac <- vapply(
  seq_along(fits),
  function(i) fits[[i]]$mainshock_num / max(1L, nrow(fits[[i]]$data)),
  numeric(1)
)

# offspring productivity (scalar) per dataset
alpha_vec <- vapply(fits, function(x) x$offspring_productivity[1], numeric(1))

print(n_events)
print(mains_frac)
print(round(alpha_vec, 2))

# ---- Background rate slices (quilt plots) ------------------------------------

fit_ex <- fits[[example_idx]]
bgr_ex <- fit_ex$bgr_rate
stopifnot(!is.null(bgr_ex), is.array(bgr_ex$funvals))

lon_mid <- midpoints(bgr_ex$breaks1)
lat_mid <- midpoints(bgr_ex$breaks2)

# pick four time slices, clipped to available range
slice_idx <- safe_idx(bgr_time_slices, dim(bgr_ex$funvals)[3])

ppi <- 300
png(bgr_png, width = 16 * ppi, height = 4 * ppi, res = ppi)
par(mfrow = c(1, length(slice_idx)), oma = c(0, 0, 0, 1))
zlim_all <- range(log1p(bgr_ex$funvals[, , slice_idx]), finite = TRUE)

for (t in slice_idx) {
  fields::quilt.plot(
    x  = rep(lon_mid, each = length(lat_mid)),
    y  = rep(lat_mid,       length(lon_mid)),
    z  = log1p(bgr_ex$funvals[, , t]),
    nx = length(lon_mid),
    ny = length(lat_mid),
    col = viridis::turbo(128),
    main = as.character(bgr_ex$breaks3[t]),
    zlim = zlim_all,
    xlab = "Longitude", ylab = "Latitude"
  )
  maps::map("world", col = "white", add = TRUE)
}
par(mfrow = c(1, 1))
dev.off()
cat("Saved background-rate tiles to:", bgr_png, "\n")

# ---- Offspring productivity: scatter & boxplot -------------------------------

# quick scatter (optional)
plot(n_events, alpha_vec,
     xlab = "Number of events", ylab = "Offspring productivity (alpha)",
     main = "Alpha vs. catalog size")

# ggplot boxplot (publication)
alpha_df <- data.frame(value = alpha_vec)

p_alpha <- ggplot(alpha_df, aes(x = "", y = value)) +
  geom_boxplot(width = 0.4, outlier.size = 0.8) +
  labs(x = NULL, y = "Offspring productivity") +
  theme_bw(base_size = 11) +
  coord_flip()

ggsave(alpha_pdf, plot = p_alpha, width = 8, height = 3)
cat("Saved alpha boxplot to:", alpha_pdf, "\n")

# ---- Spatial triggering density (radial) -------------------------------------

# build bin table (labels use [left-right])
s_breaks <- fits[[1]]$xyden$breaks
bins_s <- tibble::tibble(
  bin_id = seq_len(length(s_breaks) - 1L),
  left   = head(s_breaks, -1L),
  right  = tail(s_breaks, -1L)
) |>
  dplyr::mutate(
    mid   = (left + right) / 2,
    width = right - left,
    label = paste0(format(left, digits = 3, trim = TRUE), "-",
                   format(right, digits = 3, trim = TRUE))
  )

# collect funvals across datasets (rows: bins, cols: dataset)
all_h <- do.call(cbind, lapply(fits, function(res) res$xyden$funvals))
stopifnot(nrow(all_h) == nrow(bins_s))

df_h <- as.data.frame(all_h) |>
  dplyr::mutate(bin_id = bins_s$bin_id) |>
  tidyr::pivot_longer(-bin_id, names_to = "dataset", values_to = "funval") |>
  dplyr::left_join(bins_s, by = "bin_id") |>
  dplyr::mutate(funval_log = log1p(funval))

# first n_spatial_bins bins only
keep_s <- safe_idx(seq_len(n_spatial_bins), nrow(bins_s))
df_h_sub <- dplyr::filter(df_h, bin_id %in% keep_s)
bins_s_sub <- dplyr::slice(bins_s, keep_s)

p_h <- ggplot(df_h_sub, aes(x = factor(bin_id), y = funval_log)) +
  geom_boxplot(outlier.size = 0.7) +
  scale_x_discrete(labels = bins_s_sub$label) +
  labs(x = "Spatial lags", y = "log(Estimated value + 1)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(h_pdf, plot = p_h, width = 6, height = 4)
cat("Saved spatial density boxplot to:", h_pdf, "\n")

# ---- Temporal triggering density ---------------------------------------------

t_breaks <- fits[[1]]$tden$breaks * 1000  # ms for labeling
bins_t <- tibble::tibble(
  bin_id = seq_len(length(t_breaks) - 1L),
  left   = head(t_breaks, -1L),
  right  = tail(t_breaks, -1L)
) |>
  dplyr::mutate(
    label = paste0(
      formatC(left,  format = "f", digits = 2, drop0trailing = TRUE), "-",
      formatC(right, format = "f", digits = 2, drop0trailing = TRUE)
    )
  )

all_g <- do.call(cbind, lapply(fits, function(res) res$tden$funvals))
stopifnot(nrow(all_g) == nrow(bins_t))

df_g <- as.data.frame(all_g) |>
  dplyr::mutate(bin_id = bins_t$bin_id) |>
  tidyr::pivot_longer(-bin_id, names_to = "dataset", values_to = "funval") |>
  dplyr::left_join(bins_t, by = "bin_id") |>
  dplyr::mutate(funval_log = log1p(funval))

keep_t <- safe_idx(seq_len(n_temporal_bins), nrow(bins_t))
df_g_sub  <- dplyr::filter(df_g, bin_id %in% keep_t)
bins_t_sub <- dplyr::slice(bins_t, keep_t)

p_g <- ggplot(df_g_sub, aes(x = factor(bin_id), y = funval_log)) +
  geom_boxplot(outlier.size = 0.7) +
  scale_x_discrete(labels = bins_t_sub$label) +
  labs(x = "Temporal lags (ms)", y = "log(Estimated value + 1)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(g_pdf, plot = p_g, width = 6, height = 4)
cat("Saved temporal density boxplot to:", g_pdf, "\n")
