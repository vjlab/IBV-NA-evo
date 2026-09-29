library(tidyverse)
library(Biostrings)
library(ggplot2)
library(mgcv)

setwd("/Users/biol0375/Library/CloudStorage/Dropbox-vjlab/Ruopeng Xie/flub-evo/analysis/")

REF_STRAIN <- "B/Lee/1940|1940"
LINEAGE_COLORS <- c("A" = "#636363", "Y" = "#e41a1c", "V" = "#377eb8")
FIG_DIR <- "figures/Lee"

aa_diff <- function(seq1, seq2) {
  s1 <- strsplit(as.character(seq1), "", fixed = TRUE)[[1]]
  s2 <- strsplit(as.character(seq2), "", fixed = TRUE)[[1]]
  if (length(s1) != length(s2)) {
    stop("Sequence length mismatch")
  }
  valid <- !(s1 %in% c("X", "-") | s2 %in% c("X", "-"))
  sum(s1[valid] != s2[valid])
}

parse_year <- function(strain) {
  as.numeric(sub("-.*", "", sub(".*\\|", "", strain)))
}

compute_aa_diffs <- function(seqs, ref_seq) {
  ids <- names(seqs)
  setNames(
    vapply(ids, function(id) aa_diff(seqs[[id]], ref_seq), integer(1)),
    ids
  )
}

linear_subtitle <- function(df, y_var) {
  fit <- lm(as.formula(paste(y_var, "~ year")), data = df)
  slope <- coef(fit)[2]
  r2 <- summary(fit)$r.squared
  paste0("slope = ", round(slope, 4), " AA/year, R² = ", round(r2, 3))
}

linear_by_lineage_subtitle <- function(df, y_var) {
  fits <- df %>%
    group_by(lineage) %>%
    group_modify(~ {
      fit <- lm(as.formula(paste(y_var, "~ year")), data = .x)
      tibble(
        subtitle = paste0(
          lineage[1],
          ": slope = ",
          round(coef(fit)[2], 4),
          ", R² = ",
          round(summary(fit)$r.squared, 3)
        )
      )
    })
  paste(fits$subtitle, collapse = " | ")
}

gam_subtitle <- function(df, y_var) {
  fit <- gam(as.formula(paste(y_var, "~ s(year, bs = 'cs')")), data = df)
  edf <- summary(fit)$s.table[1, "edf"]
  paste0("pooled GAM edf = ", round(edf, 2))
}

gam_by_lineage_subtitle <- function(df, y_var) {
  fits <- df %>%
    group_by(lineage) %>%
    group_modify(~ {
      if (nrow(.x) < 4) {
        return(tibble(subtitle = paste0(lineage[1], ": insufficient data")))
      }
      fit <- gam(as.formula(paste(y_var, "~ s(year, bs = 'cs')")), data = .x)
      edf <- summary(fit)$s.table[1, "edf"]
      tibble(subtitle = paste0(lineage[1], ": edf = ", round(edf, 2)))
    })
  paste(fits$subtitle, collapse = " | ")
}

plot_lee40_drift <- function(df, y_var, y_label, fit_type, out_file) {
  y_sym <- rlang::sym(y_var)

  subtitle <- switch(
    fit_type,
    linear_pooled = linear_subtitle(df, y_var),
    linear_by_lineage = linear_by_lineage_subtitle(df, y_var),
    gam_pooled = gam_subtitle(df, y_var),
    gam_by_lineage = gam_by_lineage_subtitle(df, y_var)
  )

  p <- ggplot(df, aes(x = year, y = !!y_sym, color = lineage)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = LINEAGE_COLORS) +
    scale_x_continuous(breaks = seq(1940, 2020, 10)) +
    xlab("Year") +
    ylab(y_label) +
    labs(subtitle = subtitle) +
    theme_classic() +
    theme(legend.position = c(0.12, 0.82))

  if (fit_type == "linear_pooled") {
    p <- p + geom_smooth(
      aes(x = year, y = !!y_sym),
      method = "lm",
      se = TRUE,
      color = "black",
      inherit.aes = FALSE
    )
  } else if (fit_type == "linear_by_lineage") {
    p <- p + geom_smooth(method = "lm", se = FALSE)
  } else if (fit_type == "gam_pooled") {
    p <- p + geom_smooth(
      aes(x = year, y = !!y_sym),
      method = "gam",
      formula = y ~ s(x, bs = "cs"),
      se = TRUE,
      color = "black",
      inherit.aes = FALSE
    )
  } else if (fit_type == "gam_by_lineage") {
    p <- p + geom_smooth(
      method = "gam",
      formula = y ~ s(x, bs = "cs"),
      se = FALSE
    )
  }

  ggsave(out_file, p, width = 20, height = 9, units = "cm")
  invisible(p)
}

coords <- read.csv("table/NA_tips_coordinates.csv")
lineage <- read_tsv("Nextstrain/NA/data/IBV_NA_meta.tsv") %>%
  dplyr::select(strain, lineage)

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (!REF_STRAIN %in% coords$strain) {
  stop("Reference strain not found in NA_tips_coordinates.csv: ", REF_STRAIN)
}

lee_x <- coords$x_coordinate[coords$strain == REF_STRAIN]
lee_y <- coords$y_coordinate[coords$strain == REF_STRAIN]

ag_dist <- coords %>%
  mutate(
    ag_dist_lee40 = sqrt((x_coordinate - lee_x)^2 + (y_coordinate - lee_y)^2),
    ag_dist_lee40_x = abs(x_coordinate - lee_x)
  ) %>%
  dplyr::select(strain, ag_dist_lee40, ag_dist_lee40_x)

na_seqs <- readAAStringSet("regression/IBV_NA_only_NI.translated.fasta")

if (!REF_STRAIN %in% names(na_seqs)) {
  stop("Reference strain not found in translated FASTA file: ", REF_STRAIN)
}

na_diffs <- compute_aa_diffs(na_seqs, na_seqs[[REF_STRAIN]])

aa_df <- tibble(
  strain = names(na_diffs),
  na_aa_diff = as.integer(na_diffs)
)

df <- coords %>%
  left_join(ag_dist, by = "strain") %>%
  left_join(aa_df, by = "strain") %>%
  left_join(lineage, by = "strain") %>%
  mutate(
    year = parse_year(strain),
    lineage = factor(lineage, levels = c("A", "Y", "V"))
  )

missing_ag <- df$strain[is.na(df$ag_dist_lee40)]
missing_aa <- df$strain[is.na(df$na_aa_diff)]
missing_lineage <- df$strain[is.na(df$lineage)]

if (length(missing_ag) > 0) {
  warning("Strains missing antigenic distance to Lee/40: ", paste(missing_ag, collapse = ", "))
}
if (length(missing_aa) > 0) {
  warning("Strains missing AA differences: ", paste(missing_aa, collapse = ", "))
}
if (length(missing_lineage) > 0) {
  warning("Strains missing lineage annotation: ", paste(missing_lineage, collapse = ", "))
}

df_plot <- df %>%
  filter(
    !is.na(ag_dist_lee40),
    !is.na(ag_dist_lee40_x),
    !is.na(na_aa_diff),
    !is.na(lineage),
    !is.na(year)
  )

write.csv(df_plot, "table/lee40_drift_NA_only.csv", row.names = FALSE)

df_fig <- df_plot %>% filter(strain != REF_STRAIN)

message("Plotting ", nrow(df_fig), " of ", nrow(coords), " NI-only strains (excluding ", REF_STRAIN, ").")

plot_lee40_drift(
  df_fig,
  "na_aa_diff",
  "Full NA amino acid substitutions vs B/Lee/1940",
  "linear_pooled",
  file.path(FIG_DIR, "only_NI_NA_AA_vs_year_lee40_linear_pooled.pdf")
)
plot_lee40_drift(
  df_fig,
  "na_aa_diff",
  "Full NA amino acid substitutions vs B/Lee/1940",
  "linear_by_lineage",
  file.path(FIG_DIR, "only_NI_NA_AA_vs_year_lee40_linear_by_lineage.pdf")
)
plot_lee40_drift(
  df_fig,
  "ag_dist_lee40",
  "NA antigenic distance from B/Lee/1940 (2D)",
  "gam_pooled",
  file.path(FIG_DIR, "only_NI_antigenic_dist_lee40_vs_year_gam_pooled.pdf")
)
plot_lee40_drift(
  df_fig,
  "ag_dist_lee40",
  "NA antigenic distance from B/Lee/1940 (2D)",
  "gam_by_lineage",
  file.path(FIG_DIR, "only_NI_antigenic_dist_lee40_vs_year_gam_by_lineage.pdf")
)
plot_lee40_drift(
  df_fig,
  "ag_dist_lee40",
  "NA antigenic distance from B/Lee/1940 (2D)",
  "linear_by_lineage",
  file.path(FIG_DIR, "only_NI_antigenic_dist_lee40_vs_year_linear_by_lineage.pdf")
)
plot_lee40_drift(
  df_fig,
  "ag_dist_lee40_x",
  "NA antigenic distance from B/Lee/1940 (antigenic dim 1 only)",
  "gam_pooled",
  file.path(FIG_DIR, "only_NI_antigenic_dist_lee40_x_vs_year_gam_pooled.pdf")
)
plot_lee40_drift(
  df_fig,
  "ag_dist_lee40_x",
  "NA antigenic distance from B/Lee/1940 (antigenic dim 1 only)",
  "gam_by_lineage",
  file.path(FIG_DIR, "only_NI_antigenic_dist_lee40_x_vs_year_gam_by_lineage.pdf")
)
plot_lee40_drift(
  df_fig,
  "ag_dist_lee40_x",
  "NA antigenic distance from B/Lee/1940 (antigenic dim 1 only)",
  "linear_by_lineage",
  file.path(FIG_DIR, "only_NI_antigenic_dist_lee40_x_vs_year_linear_by_lineage.pdf")
)
