#!/usr/bin/env Rscript

# ============================================================
# TASSEL LD CSV -> LD decay plot + LD decay table
#
# Input expected from TASSEL LD export. Required columns:
#   Locus1, Locus2, Dist_bp, R^2
#
# Optional columns used for filtering/modeling if present:
#   pDiseq, N
#
# Implemented LD decay methods:
#   rolling_mean  : binned mean r2 + moving average
#   loess         : LOESS fitted curve through binned mean r2
#   raw_bin_mean  : binned mean r2 without smoothing
#   hill_weir     : Hill & Weir / Remington et al. expected r2 model
#
# Background:
#   Pairwise LD values are calculated by TASSEL. LD decay distance is then
#   estimated here as the physical distance where the selected fitted/smoothed
#   LD curve crosses a chosen r2 threshold.
#
# References for the Hill-Weir / Remington model:
#   Hill WG, Weir BS. 1988. Variances and covariances of squared linkage
#   disequilibria in finite populations. Theoretical Population Biology.
#   Remington DL et al. 2001. Structure of linkage disequilibrium and
#   phenotypic associations in the maize genome. PNAS 98:11479-11484.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(zoo)
  library(readr)
  library(tidyr)
})

# -----------------------------
# USER SETTINGS
# -----------------------------

input_csv <- "LD_maize_TASSEL.csv"
out_prefix <- "LD_decay"

# Distance/binning settings
max_dist_kb <- 2000               # plot and calculation maximum distance, kb
bin_size_kb <- 1                  # 1 = detailed but noisy; 5/10/20 = smoother/faster
rolling_window <- 7               # used only for rolling_mean
loess_span <- 0.25                # used only for loess; try 0.15-0.40
loess_degree <- 2                 # 1 or 2
thresholds <- c(0.2, 0.1)

# Choose one method, or run all methods below.
# Options: "rolling_mean", "loess", "raw_bin_mean", "hill_weir"
decay_method <- "hill_weir"
compare_all_methods <- FALSE

# Optional filtering. TASSEL LD export often includes pDiseq and N.
# Set use_pDiseq_filter <- TRUE if you want to keep only significant LD pairs.
use_pDiseq_filter <- FALSE
pDiseq_max <- 0.001

# Set min_N to a positive value if you want to keep only SNP pairs with enough samples.
# Example: min_N <- 150. Use NA or 0 to disable.
min_N <- NA

# Remove bins with very few marker pairs. This can reduce unstable spikes.
# Use 1 to disable this effectively.
min_pairs_per_bin <- 1

# If TRUE and the curve is already below a threshold at the first evaluable point,
# the decay distance is returned as NA with status "below_first_point".
return_na_if_already_below_threshold <- TRUE

# TASSEL sometimes uses chromosomes as A, B, C...
# For maize-like 10 chromosome data: A -> Chr1, B -> Chr2, ..., J -> Chr10.
convert_chr_names <- TRUE

# Hill-Weir settings
hill_weir_max_points_per_chr <- 100000   # use NA to fit all points; sampling speeds up nls
hill_weir_seed <- 123
hill_weir_start_values <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.1)
hill_weir_prediction_step_kb <- 1

# Mean row is only calculated if at least this many chromosomes have values.
min_valid_chr_for_mean <- 5

# Plot settings
show_threshold_lines <- TRUE
plot_width <- 7.5
plot_height <- 7
plot_dpi <- 300

# -----------------------------
# FUNCTIONS
# -----------------------------

chr_label <- function(x) {
  x <- as.character(x)

  if (convert_chr_names) {
    letter_map <- setNames(paste0("Chr", 1:26), LETTERS)
    x <- ifelse(x %in% names(letter_map), letter_map[x], x)
  }

  ifelse(grepl("^Chr", x, ignore.case = TRUE), x, paste0("Chr", x))
}

natural_chr_levels <- function(chr_vec) {
  tibble(chr = unique(as.character(chr_vec))) %>%
    mutate(chr_num = suppressWarnings(as.numeric(gsub("[^0-9]", "", chr)))) %>%
    arrange(is.na(chr_num), chr_num, chr) %>%
    pull(chr)
}

mean_or_na <- function(x, min_valid = 5) {
  x <- x[!is.na(x)]
  if (length(x) < min_valid) {
    return(NA_real_)
  }
  mean(x)
}

threshold_label <- function(x) {
  paste0("LD_decay_kb_r2_", gsub("\\.", ".", as.character(x), fixed = FALSE))
}

safe_parse_number <- function(x) {
  suppressWarnings(readr::parse_number(as.character(x)))
}

rollmean_safe <- function(x, k) {
  if (length(x) < k) {
    return(rep(NA_real_, length(x)))
  }
  zoo::rollmean(x, k = k, fill = NA, align = "center")
}

get_decay_from_curve <- function(curve_df, threshold) {
  curve_df <- curve_df %>%
    arrange(dist_kb) %>%
    filter(!is.na(dist_kb), !is.na(r2_curve))

  if (nrow(curve_df) < 2) {
    return(list(distance = NA_real_, status = "too_few_points"))
  }

  first_r2 <- curve_df$r2_curve[1]

  if (!is.na(first_r2) && first_r2 <= threshold && return_na_if_already_below_threshold) {
    return(list(distance = NA_real_, status = "below_first_point"))
  }

  crossing_idx <- which(
    curve_df$r2_curve[-nrow(curve_df)] > threshold &
      curve_df$r2_curve[-1] <= threshold
  )

  if (length(crossing_idx) == 0) {
    if (all(curve_df$r2_curve > threshold, na.rm = TRUE)) {
      return(list(distance = NA_real_, status = "not_reached_within_max_dist"))
    }
    if (all(curve_df$r2_curve <= threshold, na.rm = TRUE)) {
      return(list(distance = NA_real_, status = "always_below_threshold"))
    }
    return(list(distance = NA_real_, status = "no_downward_crossing"))
  }

  i <- crossing_idx[1]

  x1 <- curve_df$dist_kb[i]
  x2 <- curve_df$dist_kb[i + 1]
  y1 <- curve_df$r2_curve[i]
  y2 <- curve_df$r2_curve[i + 1]

  if (is.na(y1) || is.na(y2) || y1 == y2) {
    return(list(distance = x2, status = "ok_no_interpolation"))
  }

  x_cross <- x1 + (threshold - y1) * (x2 - x1) / (y2 - y1)
  list(distance = x_cross, status = "ok")
}

hill_weir_expected_r2 <- function(dist_bp, N, C) {
  ((10 + C * dist_bp) / ((2 + C * dist_bp) * (11 + C * dist_bp))) *
    (1 + ((3 + C * dist_bp) *
            (12 + 12 * C * dist_bp + (C * dist_bp)^2)) /
       (2 * N * (2 + C * dist_bp) * (11 + C * dist_bp)))
}

fit_hill_weir_one_chr <- function(df, max_dist_kb) {
  df <- df %>%
    filter(
      !is.na(r2),
      !is.na(dist_bp),
      !is.na(N),
      r2 >= 0,
      r2 <= 1,
      dist_bp > 0,
      dist_bp <= max_dist_kb * 1000,
      N > 0
    ) %>%
    select(dist_bp, r2, N)

  if (nrow(df) < 100) {
    return(list(
      curve = tibble(dist_bp = numeric(), dist_kb = numeric(), N = numeric(), r2_curve = numeric()),
      C = NA_real_,
      status = "too_few_points"
    ))
  }

  if (!is.na(hill_weir_max_points_per_chr) && nrow(df) > hill_weir_max_points_per_chr) {
    set.seed(hill_weir_seed)
    df <- df %>% slice_sample(n = hill_weir_max_points_per_chr)
  }

  best_fit <- NULL
  best_rss <- Inf
  best_start <- NA_real_

  for (start_C in hill_weir_start_values) {
    fit <- tryCatch(
      suppressWarnings(
        nls(
          r2 ~ hill_weir_expected_r2(dist_bp, N, C),
          data = df,
          start = list(C = start_C),
          algorithm = "port",
          lower = c(C = 0),
          control = nls.control(maxiter = 200, warnOnly = TRUE)
        )
      ),
      error = function(e) NULL
    )

    if (!is.null(fit)) {
      rss <- sum(residuals(fit)^2, na.rm = TRUE)
      if (is.finite(rss) && rss < best_rss) {
        best_fit <- fit
        best_rss <- rss
        best_start <- start_C
      }
    }
  }

  if (is.null(best_fit)) {
    return(list(
      curve = tibble(dist_bp = numeric(), dist_kb = numeric(), N = numeric(), r2_curve = numeric()),
      C = NA_real_,
      status = "nls_failed"
    ))
  }

  C_est <- as.numeric(coef(best_fit)[["C"]])

  pred_grid <- tibble(
    dist_kb = seq(hill_weir_prediction_step_kb, max_dist_kb, by = hill_weir_prediction_step_kb),
    dist_bp = dist_kb * 1000,
    N = mean(df$N, na.rm = TRUE)
  ) %>%
    mutate(r2_curve = hill_weir_expected_r2(dist_bp = dist_bp, N = N, C = C_est))

  list(
    curve = pred_grid,
    C = C_est,
    status = paste0("ok_start_", best_start)
  )
}

make_binned_ld <- function(ld) {
  ld %>%
    mutate(bin_kb = floor(dist_kb / bin_size_kb) * bin_size_kb + bin_size_kb / 2) %>%
    group_by(chr, bin_kb) %>%
    summarise(
      mean_r2 = mean(r2, na.rm = TRUE),
      median_r2 = median(r2, na.rm = TRUE),
      n_pairs = n(),
      .groups = "drop"
    ) %>%
    filter(n_pairs >= min_pairs_per_bin) %>%
    arrange(chr, bin_kb) %>%
    rename(dist_kb = bin_kb)
}

make_curve_for_method <- function(method, ld, ld_binned) {
  if (method == "raw_bin_mean") {
    return(ld_binned %>%
             transmute(chr, dist_kb, r2_curve = mean_r2, method = method))
  }

  if (method == "rolling_mean") {
    return(ld_binned %>%
             arrange(chr, dist_kb) %>%
             group_by(chr) %>%
             mutate(r2_curve = rollmean_safe(mean_r2, rolling_window)) %>%
             ungroup() %>%
             transmute(chr, dist_kb, r2_curve, method = method))
  }

  if (method == "loess") {
    return(ld_binned %>%
             group_by(chr) %>%
             group_modify(~{
               df <- .x %>% filter(!is.na(mean_r2), !is.na(dist_kb))
               if (nrow(df) < 20) {
                 return(tibble(dist_kb = numeric(), r2_curve = numeric(), method = character()))
               }

               fit <- tryCatch(
                 loess(
                   mean_r2 ~ dist_kb,
                   data = df,
                   weights = pmax(df$n_pairs, 1),
                   span = loess_span,
                   degree = loess_degree,
                   control = loess.control(surface = "direct")
                 ),
                 error = function(e) NULL
               )

               if (is.null(fit)) {
                 return(tibble(dist_kb = numeric(), r2_curve = numeric(), method = character()))
               }

               pred_grid <- tibble(
                 dist_kb = seq(min(df$dist_kb, na.rm = TRUE), max(df$dist_kb, na.rm = TRUE), by = 1)
               )

               pred_grid$r2_curve <- as.numeric(predict(fit, newdata = pred_grid))
               pred_grid$r2_curve <- pmin(pmax(pred_grid$r2_curve, 0), 1)
               pred_grid$method <- method
               pred_grid
             }) %>%
             ungroup())
  }

  if (method == "hill_weir") {
    return(ld %>%
             group_by(chr) %>%
             group_modify(~{
               hw <- fit_hill_weir_one_chr(.x, max_dist_kb = max_dist_kb)
               hw$curve %>%
                 mutate(method = method, C = hw$C, fit_status = hw$status)
             }) %>%
             ungroup())
  }

  stop("Unknown decay_method: ", method)
}

make_decay_table_and_diagnostics <- function(curve_df, method) {
  chr_levels <- natural_chr_levels(curve_df$chr)

  diagnostics <- curve_df %>%
    mutate(chr = factor(chr, levels = chr_levels)) %>%
    group_by(chr) %>%
    group_modify(~{
      rows <- lapply(thresholds, function(thr) {
        ans <- get_decay_from_curve(.x, thr)
        tibble(
          threshold = thr,
          decay_distance_kb = ans$distance,
          status = ans$status
        )
      })
      bind_rows(rows)
    }) %>%
    ungroup() %>%
    mutate(method = method)

  chr_table <- diagnostics %>%
    mutate(col = paste0("LD_decay_kb_r2_", threshold)) %>%
    select(chr, col, decay_distance_kb) %>%
    tidyr::pivot_wider(names_from = col, values_from = decay_distance_kb) %>%
    mutate(chr = as.character(chr))

  # Ensure expected threshold columns exist even if all values are missing
  for (thr in thresholds) {
    nm <- paste0("LD_decay_kb_r2_", thr)
    if (!nm %in% names(chr_table)) {
      chr_table[[nm]] <- NA_real_
    }
  }

  # Add Hill-Weir model parameters if available
  if (method == "hill_weir" && "C" %in% names(curve_df)) {
    hw_params <- curve_df %>%
      group_by(chr) %>%
      summarise(
        C = unique(na.omit(C))[1],
        fit_status = unique(na.omit(fit_status))[1],
        .groups = "drop"
      ) %>%
      mutate(
        C = as.numeric(C),
        fit_status = as.character(fit_status)
      )

    chr_table <- chr_table %>% left_join(hw_params, by = "chr")

    half_diag <- curve_df %>%
      group_by(chr) %>%
      group_modify(~{
        maxld <- max(.x$r2_curve, na.rm = TRUE)
        if (!is.finite(maxld)) {
          return(tibble(LD_half_decay_kb = NA_real_, half_decay_threshold = NA_real_))
        }
        half_thr <- maxld * 0.5
        ans <- get_decay_from_curve(.x, half_thr)
        tibble(LD_half_decay_kb = ans$distance, half_decay_threshold = half_thr)
      }) %>%
      ungroup()

    chr_table <- chr_table %>% left_join(half_diag, by = "chr")
  }

  # Natural chromosome order
  chr_table <- chr_table %>%
    mutate(chr_num = suppressWarnings(as.numeric(gsub("[^0-9]", "", chr)))) %>%
    arrange(is.na(chr_num), chr_num, chr) %>%
    select(-chr_num)

  # Mean row
  mean_data <- list(chr = "Mean")
  for (thr in thresholds) {
    nm <- paste0("LD_decay_kb_r2_", thr)
    mean_data[[nm]] <- mean_or_na(chr_table[[nm]], min_valid = min_valid_chr_for_mean)
  }
  if ("LD_half_decay_kb" %in% names(chr_table)) {
    mean_data[["LD_half_decay_kb"]] <- mean_or_na(chr_table$LD_half_decay_kb, min_valid = min_valid_chr_for_mean)
  }
  if ("half_decay_threshold" %in% names(chr_table)) {
    mean_data[["half_decay_threshold"]] <- mean_or_na(chr_table$half_decay_threshold, min_valid = min_valid_chr_for_mean)
  }
  if ("C" %in% names(chr_table)) {
    mean_data[["C"]] <- mean(chr_table$C, na.rm = TRUE)
  }
  if ("fit_status" %in% names(chr_table)) {
    mean_data[["fit_status"]] <- NA_character_
  }

  mean_row <- as_tibble(mean_data)

  out_table <- bind_rows(chr_table, mean_row) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))

  list(table = out_table, diagnostics = diagnostics)
}

plot_curve <- function(curve_df, method, output_prefix) {
  chr_levels <- natural_chr_levels(curve_df$chr)
  curve_df <- curve_df %>%
    mutate(chr = factor(chr, levels = chr_levels))

  p <- ggplot(curve_df, aes(x = dist_kb, y = r2_curve, color = chr)) +
    geom_line(linewidth = 0.55, na.rm = TRUE) +
    scale_x_continuous(
      limits = c(0, max_dist_kb),
      breaks = seq(0, max_dist_kb, by = 500)
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      title = paste("LD decay -", method),
      x = "Distance (Kb)",
      y = expression(r^2),
      color = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 11)
    )

  if (show_threshold_lines) {
    for (thr in thresholds) {
      p <- p + geom_hline(yintercept = thr, linetype = "dashed", linewidth = 0.3)
    }
  }

  ggsave(
    filename = paste0(output_prefix, ".png"),
    plot = p,
    width = plot_width,
    height = plot_height,
    dpi = plot_dpi
  )

  # Requires svglite in most R installations. Install with: install.packages("svglite")
  ggsave(
    filename = paste0(output_prefix, ".svg"),
    plot = p,
    width = plot_width,
    height = plot_height,
    device = "svg"
  )

  invisible(p)
}

# -----------------------------
# READ TASSEL LD FILE
# -----------------------------

message("Reading TASSEL LD file...")

header <- names(fread(input_csv, nrows = 0, check.names = FALSE))

required_cols <- c("Locus1", "Locus2", "Dist_bp", "R^2")
missing_required <- setdiff(required_cols, header)

if (length(missing_required) > 0) {
  stop("Missing required columns in TASSEL LD file: ", paste(missing_required, collapse = ", "))
}

optional_cols <- intersect(c("pDiseq", "N"), header)
read_cols <- c(required_cols, optional_cols)

ld <- fread(
  input_csv,
  select = read_cols,
  check.names = FALSE,
  showProgress = TRUE
)

new_names <- c("chr1", "chr2", "dist_bp", "r2", optional_cols)
names(ld) <- new_names

if (!"pDiseq" %in% names(ld)) {
  ld$pDiseq <- NA_real_
  message("Column pDiseq was not found. pDiseq filtering will be ignored.")
}

if (!"N" %in% names(ld)) {
  ld$N <- NA_real_
  message("Column N was not found. Hill-Weir method and N filtering require N.")
}

message("Initial column structure:")
str(ld)

# -----------------------------
# CLEAN AND FILTER DATA
# -----------------------------

message("Cleaning and filtering LD data...")

ld <- ld %>%
  mutate(
    chr1 = as.character(chr1),
    chr2 = as.character(chr2),
    dist_bp = safe_parse_number(dist_bp),
    r2 = safe_parse_number(r2),
    pDiseq = safe_parse_number(pDiseq),
    N = safe_parse_number(N)
  ) %>%
  filter(chr1 == chr2) %>%
  mutate(
    chr = chr_label(chr1),
    dist_kb = dist_bp / 1000
  ) %>%
  filter(
    !is.na(r2),
    !is.na(dist_bp),
    !is.na(dist_kb),
    r2 >= 0,
    r2 <= 1,
    dist_bp >= 0,
    dist_kb <= max_dist_kb
  )

if (use_pDiseq_filter) {
  if (all(is.na(ld$pDiseq))) {
    warning("use_pDiseq_filter is TRUE, but pDiseq is missing or all NA. Skipping pDiseq filter.")
  } else {
    ld <- ld %>% filter(!is.na(pDiseq), pDiseq <= pDiseq_max)
    message("Applied pDiseq filter: pDiseq <= ", pDiseq_max)
  }
}

if (!is.na(min_N) && min_N > 0) {
  if (all(is.na(ld$N))) {
    warning("min_N is set, but N is missing or all NA. Skipping N filter.")
  } else {
    ld <- ld %>% filter(!is.na(N), N >= min_N)
    message("Applied N filter: N >= ", min_N)
  }
}

message("Rows after filtering: ", nrow(ld))

message("Distance summary, bp:")
print(summary(ld$dist_bp))

message("r2 summary:")
print(summary(ld$r2))

if (!all(is.na(ld$N))) {
  message("N summary:")
  print(summary(ld$N))
}

if (nrow(ld) == 0) {
  stop(
    "No LD records remained after filtering. ",
    "Check chromosome names, Dist_bp/R^2 columns, pDiseq/N filters, or max_dist_kb."
  )
}

# -----------------------------
# BINNING
# -----------------------------

message("Binning LD values...")

ld_binned <- make_binned_ld(ld)
chr_levels <- natural_chr_levels(ld_binned$chr)

ld_binned <- ld_binned %>%
  mutate(chr = factor(chr, levels = chr_levels))

write_csv(ld_binned, paste0(out_prefix, "_binned_values.csv"))

# -----------------------------
# RUN SELECTED LD DECAY METHOD(S)
# -----------------------------

all_methods <- c("rolling_mean", "loess", "raw_bin_mean", "hill_weir")
methods_to_run <- if (compare_all_methods) all_methods else decay_method

bad_methods <- setdiff(methods_to_run, all_methods)
if (length(bad_methods) > 0) {
  stop("Unknown decay method(s): ", paste(bad_methods, collapse = ", "))
}

created_files <- c(paste0(out_prefix, "_binned_values.csv"))

for (method in methods_to_run) {
  message("------------------------------------------------------------")
  message("Running LD decay method: ", method)

  if (method == "hill_weir" && all(is.na(ld$N))) {
    warning("Skipping hill_weir because N column is missing or all NA.")
    next
  }

  output_prefix <- paste0(out_prefix, "_", method)

  curve_df <- make_curve_for_method(method, ld = ld, ld_binned = ld_binned)

  if (nrow(curve_df) == 0) {
    warning("No curve values were generated for method: ", method)
    next
  }

  curve_df <- curve_df %>%
    mutate(chr = factor(chr, levels = chr_levels))

  table_diag <- make_decay_table_and_diagnostics(curve_df, method)

  write_csv(curve_df, paste0(output_prefix, "_curve_values.csv"))
  write_csv(table_diag$table, paste0(output_prefix, "_table.csv"))
  write_csv(table_diag$diagnostics, paste0(output_prefix, "_decay_diagnostics.csv"))

  message("LD decay table for method: ", method)
  print(table_diag$table)

  message("Making LD decay plot for method: ", method)
  plot_curve(curve_df, method = method, output_prefix = output_prefix)

  created_files <- c(
    created_files,
    paste0(output_prefix, ".png"),
    paste0(output_prefix, ".svg"),
    paste0(output_prefix, "_curve_values.csv"),
    paste0(output_prefix, "_table.csv"),
    paste0(output_prefix, "_decay_diagnostics.csv")
  )
}

message("------------------------------------------------------------")
message("Done.")
message("Created files:")
for (f in created_files) {
  message("  ", f)
}
