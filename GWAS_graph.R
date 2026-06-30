#!/usr/bin/env Rscript

# ============================================================
# TASSEL GWAS output -> Manhattan + QQ plots
# Robust version v3
#
# Updates in v3:
#   - QQ plot now includes a null-distribution confidence envelope.
#   - Envelope is based on beta distribution of ordered p-values.
#   - Manhattan plot uses colour = chromosome, not shape.
#
# Input:
#   TASSEL pipeline CSV-like output files
#
# Uses:
#   Only chromosome-assigned markers:
#     Marker != None
#     Chr not missing
#     Pos not missing
#     p valid: 0 < p <= 1
#
# Output:
#   plots/
#     CTI2018_manhattan.svg
#     CTI2018_manhattan.png
#     CTI2018_qq.svg
#     CTI2018_qq.png
#     CTI2022_manhattan.svg
#     CTI2022_manhattan.png
#     CTI2022_qq.svg
#     CTI2022_qq.png
#
#   tables/
#     CTI2018_cleaned_chromosome_markers.csv
#     CTI2018_top_hits.csv
#     CTI2018_chromosome_marker_summary.csv
#     CTI2022_cleaned_chromosome_markers.csv
#     CTI2022_top_hits.csv
#     CTI2022_chromosome_marker_summary.csv
#     CTI_summary_statistics.csv
# ============================================================


# ----------------------------
# 1. Package handling
# ----------------------------

required_packages <- c(
  "tidyverse",
  "readr",
  "ggplot2",
  "svglite",
  "ragg",
  "scales"
)

install_if_missing <- function(pkgs) {
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

install_if_missing(required_packages)

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(svglite)
  library(ragg)
  library(scales)
})


# ----------------------------
# 2. User settings
# ----------------------------

input_files <- c(
  SPAD2018 = "SPAD_2018.csv",
  SPAD2022 = "SPAD_2022.csv"
)

output_plot_dir <- "plots"
output_table_dir <- "tables"

dir.create(output_plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_table_dir, showWarnings = FALSE, recursive = TRUE)

rename_chr_to_numbers <- FALSE

chr_order_letters <- LETTERS[1:10]

chr_label_map <- tibble(
  Chr = chr_order_letters,
  Chr_plot = if (rename_chr_to_numbers) as.character(1:10) else chr_order_letters
)

chromosome_palette_letters <- c(
  "A" = "#1b9e77",
  "B" = "#d95f02",
  "C" = "#7570b3",
  "D" = "#e7298a",
  "E" = "#66a61e",
  "F" = "#e6ab02",
  "G" = "#a6761d",
  "H" = "#666666",
  "I" = "#1f78b4",
  "J" = "#b2df8a"
)

chromosome_palette_numbers <- c(
  "1"  = "#1b9e77",
  "2"  = "#d95f02",
  "3"  = "#7570b3",
  "4"  = "#e7298a",
  "5"  = "#66a61e",
  "6"  = "#e6ab02",
  "7"  = "#a6761d",
  "8"  = "#666666",
  "9"  = "#1f78b4",
  "10" = "#b2df8a"
)

chromosome_palette <- if (rename_chr_to_numbers) {
  chromosome_palette_numbers
} else {
  chromosome_palette_letters
}

# QQ-plot konfidencia-boríték
qq_confidence_level <- 0.95

# Ábraméretek
manhattan_width <- 13
manhattan_height <- 5.5

qq_width <- 6
qq_height <- 6

png_dpi <- 300

# Top találatok száma
n_top_hits <- 50


# ----------------------------
# 3. Helper functions
# ----------------------------

check_required_columns <- function(df, file_label) {
  required_cols <- c("Trait", "Marker", "Chr", "Pos", "p")
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", file_label, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
}


read_tassel_csv <- function(file_path, file_label) {
  if (!file.exists(file_path)) {
    stop("Input file not found: ", file_path)
  }
  
  message("Reading: ", file_path)
  
  df <- read_csv(
    file = file_path,
    show_col_types = FALSE,
    na = c("", "NA", "NaN", "nan", "None")
  )
  
  check_required_columns(df, file_label)
  
  df
}


clean_tassel_gwas <- function(df, file_label) {
  
  cleaned <- df %>%
    mutate(
      Trait = as.character(Trait),
      Marker = as.character(Marker),
      Chr = as.character(Chr),
      Pos = suppressWarnings(as.numeric(Pos)),
      p = suppressWarnings(as.numeric(p)),
      F = if ("F" %in% names(.)) suppressWarnings(as.numeric(F)) else NA_real_,
      MarkerR2 = if ("MarkerR2" %in% names(.)) suppressWarnings(as.numeric(MarkerR2)) else NA_real_
    ) %>%
    filter(
      !is.na(Marker),
      Marker != "",
      Marker != "None",
      !is.na(Chr),
      Chr != "",
      !is.na(Pos),
      is.finite(Pos),
      Pos >= 0,
      !is.na(p),
      is.finite(p),
      p > 0,
      p <= 1
    ) %>%
    mutate(
      Chr = toupper(Chr)
    ) %>%
    filter(
      Chr %in% chr_order_letters
    ) %>%
    distinct(Chr, Pos, Marker, .keep_all = TRUE)
  
  if (nrow(cleaned) == 0) {
    stop("No valid chromosome-assigned markers remained after filtering for ", file_label)
  }
  
  cleaned
}


prepare_manhattan_data <- function(df) {
  
  df2 <- df %>%
    mutate(
      Chr = factor(Chr, levels = chr_order_letters),
      log10p = -log10(p)
    ) %>%
    arrange(Chr, Pos)
  
  chr_info <- df2 %>%
    group_by(Chr) %>%
    summarise(
      chr_len = max(Pos, na.rm = TRUE),
      n_markers_chr = n(),
      .groups = "drop"
    ) %>%
    arrange(Chr) %>%
    mutate(
      offset = lag(cumsum(chr_len), default = 0)
    )
  
  plot_df <- df2 %>%
    left_join(chr_info, by = "Chr") %>%
    mutate(
      pos_cum = Pos + offset
    ) %>%
    left_join(chr_label_map, by = "Chr") %>%
    mutate(
      Chr_plot = factor(
        Chr_plot,
        levels = if (rename_chr_to_numbers) as.character(1:10) else chr_order_letters
      )
    )
  
  axis_df <- plot_df %>%
    group_by(Chr, Chr_plot) %>%
    summarise(
      center = mean(range(pos_cum, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(Chr)
  
  list(
    plot_df = plot_df,
    axis_df = axis_df,
    chr_info = chr_info
  )
}


calculate_thresholds <- function(df) {
  n_markers <- nrow(df)
  
  tibble(
    n_markers = n_markers,
    bonferroni_p = 0.05 / n_markers,
    suggestive_p = 1 / n_markers,
    bonferroni_log10p = -log10(0.05 / n_markers),
    suggestive_log10p = -log10(1 / n_markers)
  )
}


calculate_lambda_gc <- function(p_values) {
  p_values <- p_values[is.finite(p_values) & p_values > 0 & p_values <= 1]
  
  if (length(p_values) < 10) {
    return(NA_real_)
  }
  
  chisq_values <- qchisq(1 - p_values, df = 1)
  lambda_gc <- median(chisq_values, na.rm = TRUE) / qchisq(0.5, df = 1)
  
  lambda_gc
}


make_manhattan_plot <- function(plot_df, axis_df, thresholds, plot_title) {
  
  ggplot(plot_df, aes(x = pos_cum, y = log10p)) +
    geom_point(
      aes(colour = Chr_plot),
      size = 0.85,
      alpha = 0.85,
      show.legend = TRUE
    ) +
    scale_colour_manual(
      name = "Chromosome",
      values = chromosome_palette,
      drop = FALSE
    ) +
    geom_hline(
      yintercept = thresholds$suggestive_log10p,
      linewidth = 0.4,
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = thresholds$bonferroni_log10p,
      linewidth = 0.4,
      linetype = "solid"
    ) +
    annotate(
      "text",
      x = max(plot_df$pos_cum, na.rm = TRUE),
      y = thresholds$bonferroni_log10p,
      label = "Bonferroni",
      hjust = 1.02,
      vjust = -0.4,
      size = 3
    ) +
    annotate(
      "text",
      x = max(plot_df$pos_cum, na.rm = TRUE),
      y = thresholds$suggestive_log10p,
      label = "Suggestive",
      hjust = 1.02,
      vjust = -0.4,
      size = 3
    ) +
    scale_x_continuous(
      breaks = axis_df$center,
      labels = axis_df$Chr_plot,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.10))
    ) +
    labs(
      title = plot_title,
      subtitle = paste0(
        "Chromosome-assigned markers only; N = ",
        comma(thresholds$n_markers),
        " | Bonferroni p = ",
        scientific(thresholds$bonferroni_p, digits = 3),
        " | Suggestive p = ",
        scientific(thresholds$suggestive_p, digits = 3)
      ),
      x = "Chromosome",
      y = expression(-log[10](p))
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(size = 3, alpha = 1),
        nrow = 1
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9),
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 9)
    )
}


# ----------------------------
# 4. QQ plot with confidence envelope
# ----------------------------

prepare_qq_data <- function(df, confidence_level = 0.95) {
  
  p_values <- df$p
  p_values <- p_values[is.finite(p_values) & p_values > 0 & p_values <= 1]
  p_values <- sort(p_values)
  
  n <- length(p_values)
  
  if (n < 10) {
    stop("Too few valid p-values for QQ plot.")
  }
  
  alpha <- 1 - confidence_level
  rank_i <- seq_len(n)
  
  expected_p <- ppoints(n)
  
  # Under the null hypothesis, the i-th ordered p-value follows:
  # p_(i) ~ Beta(i, n - i + 1)
  #
  # Since y = -log10(p), lower p-values correspond to larger y-values.
  ci_p_lower <- qbeta(alpha / 2, rank_i, n - rank_i + 1)
  ci_p_upper <- qbeta(1 - alpha / 2, rank_i, n - rank_i + 1)
  
  qq_df <- tibble(
    rank = rank_i,
    expected = -log10(expected_p),
    observed = -log10(p_values),
    ci_lower = -log10(ci_p_upper),
    ci_upper = -log10(ci_p_lower)
  ) %>%
    arrange(expected)
  
  qq_df
}


make_qq_plot <- function(qq_df, lambda_gc, plot_title, n_markers, confidence_level = 0.95) {
  
  max_val <- max(
    c(
      qq_df$expected,
      qq_df$observed,
      qq_df$ci_lower,
      qq_df$ci_upper
    ),
    na.rm = TRUE
  )
  
  envelope_label <- paste0(round(confidence_level * 100), "% confidence envelope")
  
  ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_ribbon(
      aes(
        x = expected,
        ymin = ci_lower,
        ymax = ci_upper
      ),
      alpha = 0.25
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_point(
      size = 1.1,
      alpha = 0.75
    ) +
    coord_equal(
      xlim = c(0, max_val * 1.03),
      ylim = c(0, max_val * 1.03)
    ) +
    labs(
      title = plot_title,
      subtitle = paste0(
        "Chromosome-assigned markers only; N = ",
        comma(n_markers),
        " | lambda GC = ",
        round(lambda_gc, 3),
        " | ",
        envelope_label
      ),
      x = expression(Expected ~ -log[10](p)),
      y = expression(Observed ~ -log[10](p))
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9)
    )
}


save_plot_svg_png <- function(plot_object, base_file_name, width, height) {
  
  svg_path <- file.path(output_plot_dir, paste0(base_file_name, ".svg"))
  png_path <- file.path(output_plot_dir, paste0(base_file_name, ".png"))
  
  ggsave(
    filename = svg_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = svglite
  )
  
  ggsave(
    filename = png_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = png_dpi,
    device = ragg::agg_png
  )
  
  message("Saved: ", svg_path)
  message("Saved: ", png_path)
}


write_top_hits <- function(df, file_label) {
  
  top_hits <- df %>%
    mutate(
      log10p = -log10(p)
    ) %>%
    arrange(p) %>%
    select(
      any_of(c(
        "Trait",
        "Marker",
        "Chr",
        "Pos",
        "p",
        "log10p",
        "F",
        "MarkerR2",
        "df",
        "errordf",
        "Genetic Var",
        "Residual Var",
        "-2LnLikelihood"
      ))
    ) %>%
    slice_head(n = n_top_hits)
  
  out_path <- file.path(output_table_dir, paste0(file_label, "_top_hits.csv"))
  write_csv(top_hits, out_path)
  message("Saved: ", out_path)
  
  top_hits
}


write_cleaned_table <- function(df, file_label) {
  
  out_path <- file.path(output_table_dir, paste0(file_label, "_cleaned_chromosome_markers.csv"))
  write_csv(df, out_path)
  message("Saved: ", out_path)
}


write_chr_marker_summary <- function(df, file_label) {
  
  chr_summary <- df %>%
    group_by(Chr) %>%
    summarise(
      n_markers = n(),
      min_pos = min(Pos, na.rm = TRUE),
      max_pos = max(Pos, na.rm = TRUE),
      min_p = min(p, na.rm = TRUE),
      max_log10p = max(-log10(p), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(factor(Chr, levels = chr_order_letters))
  
  out_path <- file.path(output_table_dir, paste0(file_label, "_chromosome_marker_summary.csv"))
  write_csv(chr_summary, out_path)
  message("Saved: ", out_path)
  
  chr_summary
}


process_gwas_file <- function(file_path, file_label) {
  
  raw_df <- read_tassel_csv(file_path, file_label)
  clean_df <- clean_tassel_gwas(raw_df, file_label)
  
  message(file_label, ": valid chromosome-assigned markers = ", nrow(clean_df))
  
  thresholds <- calculate_thresholds(clean_df)
  lambda_gc <- calculate_lambda_gc(clean_df$p)
  
  manhattan_data <- prepare_manhattan_data(clean_df)
  
  manhattan_plot <- make_manhattan_plot(
    plot_df = manhattan_data$plot_df,
    axis_df = manhattan_data$axis_df,
    thresholds = thresholds,
    plot_title = paste0(file_label, " Manhattan plot")
  )
  
  qq_df <- prepare_qq_data(
    clean_df,
    confidence_level = qq_confidence_level
  )
  
  qq_plot <- make_qq_plot(
    qq_df = qq_df,
    lambda_gc = lambda_gc,
    plot_title = paste0(file_label, " QQ plot"),
    n_markers = nrow(clean_df),
    confidence_level = qq_confidence_level
  )
  
  save_plot_svg_png(
    plot_object = manhattan_plot,
    base_file_name = paste0(file_label, "_manhattan"),
    width = manhattan_width,
    height = manhattan_height
  )
  
  save_plot_svg_png(
    plot_object = qq_plot,
    base_file_name = paste0(file_label, "_qq"),
    width = qq_width,
    height = qq_height
  )
  
  write_cleaned_table(clean_df, file_label)
  top_hits <- write_top_hits(clean_df, file_label)
  chr_summary <- write_chr_marker_summary(clean_df, file_label)
  
  summary_row <- thresholds %>%
    mutate(
      dataset = file_label,
      input_file = file_path,
      lambda_gc = lambda_gc,
      min_p = min(clean_df$p, na.rm = TRUE),
      max_log10p = max(-log10(clean_df$p), na.rm = TRUE),
      n_chr = n_distinct(clean_df$Chr),
      top_marker = top_hits$Marker[1],
      top_chr = top_hits$Chr[1],
      top_pos = top_hits$Pos[1],
      top_p = top_hits$p[1],
      qq_confidence_level = qq_confidence_level
    ) %>%
    select(
      dataset,
      input_file,
      n_markers,
      n_chr,
      min_p,
      max_log10p,
      lambda_gc,
      bonferroni_p,
      suggestive_p,
      bonferroni_log10p,
      suggestive_log10p,
      qq_confidence_level,
      top_marker,
      top_chr,
      top_pos,
      top_p
    )
  
  summary_row
}


# ----------------------------
# 5. Run analysis
# ----------------------------

summary_list <- imap(
  input_files,
  ~ process_gwas_file(file_path = .x, file_label = .y)
)

summary_table <- bind_rows(summary_list)

summary_path <- file.path(output_table_dir, "CTI_summary_statistics.csv")
write_csv(summary_table, summary_path)

message("Saved: ", summary_path)
message("Analysis completed successfully.")
print(summary_table)
