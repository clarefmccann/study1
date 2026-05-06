## 03_gamms.R
## GAMMs for PDS composite trajectories, per sex × reporter.
## Model: pds_comp ~ s(age, k=6) + s(id_fac, bs="re")
## Derives individual-level timing (onset age) and tempo (growth velocity
## at onset) for use as LCA inputs in 04_lca.R.

library(mgcv)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gratia)

set.seed(90025)

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
root_path <- Sys.getenv("HOME_DIR")
if (!nzchar(root_path)) {
  root_path <- Sys.getenv("HOME")
}

data_dir <- file.path(
  root_path,
  "projects/abcd-projs/abcd-data-release-6.0/cfm/physical-health/puberty"
)
if (!dir.exists(data_dir)) {
  data_dir <- file.path(
    root_path,
    "University of Oregon Dropbox/Clare McCann/mine/projects",
    "abcd-projs/abcd-data-release-6.0/cfm/physical-health/puberty"
  )
}
if (!dir.exists(data_dir)) {
  stop("Cannot locate puberty data directory: ", data_dir)
}

out_dir <- file.path(data_dir, "gamm_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------
datasets <- list(
  female_parent = read.csv(file.path(data_dir, "female_parent_long.csv")),
  female_youth = read.csv(file.path(data_dir, "female_youth_long.csv")),
  male_parent = read.csv(file.path(data_dir, "male_parent_long.csv")),
  male_youth = file.path(data_dir, "male_youth_long.csv") |> read.csv()
)

# Retain only rows with observed pds_comp and age; factor id for RE
datasets <- lapply(datasets, function(df) {
  df %>%
    filter(!is.na(pds_comp), !is.na(age)) %>%
    mutate(id_fac = factor(id))
})

cat("Dataset sizes (rows with complete pds_comp + age):\n")
print(sapply(datasets, nrow))

# ---------------------------------------------------------------------------
# FIT GAMMs
# bam() + fREML + discrete = TRUE for speed at ABCD sample sizes.
# k = 6 for s(age): adequate for ≤7 annual waves; avoid over-smoothing.
# s(id_fac, bs="re"): random participant intercept.
# ---------------------------------------------------------------------------
fit_gamm <- function(df, label) {
  cat(
    "\nFitting GAMM:",
    label,
    " (n =",
    nrow(df),
    "rows,",
    n_distinct(df$id),
    "participants)\n"
  )
  mgcv::bam(
    pds_comp ~ s(age, k = 6) + s(id_fac, bs = "re"),
    data = df,
    method = "fREML",
    discrete = TRUE
  )
}

gamm_fits <- mapply(fit_gamm, datasets, names(datasets), SIMPLIFY = FALSE)

saveRDS(gamm_fits, file.path(out_dir, "gamm_fits.rds"))
cat("\nModel objects saved.\n")

# Print summaries
for (nm in names(gamm_fits)) {
  cat("\n\n===", nm, "===\n")
  print(summary(gamm_fits[[nm]]))
}

# ---------------------------------------------------------------------------
# POPULATION SMOOTH + DERIVATIVE  (gratia)
# gratia >= 0.7: smooth_estimates() → age, .estimate, .se
#                derivatives()      → age, .derivative, .se
# ---------------------------------------------------------------------------
AGE_GRID_N <- 300

safe_col <- function(df, primary, fallback) {
  if (primary %in% names(df)) primary else fallback
}

smooth_list <- lapply(names(gamm_fits), function(nm) {
  sm <- gratia::smooth_estimates(
    gamm_fits[[nm]],
    smooth = "s(age)",
    n = AGE_GRID_N
  )
  sm$dataset <- nm
  sm
})
smooth_df <- bind_rows(smooth_list) %>%
  tidyr::separate(
    dataset,
    into = c("sex", "reporter"),
    sep = "_",
    extra = "merge"
  ) %>%
  mutate(
    sex = factor(sex, levels = c("female", "male")),
    reporter = factor(reporter, levels = c("parent", "youth"))
  )

deriv_list <- lapply(names(gamm_fits), function(nm) {
  d1 <- gratia::derivatives(
    gamm_fits[[nm]],
    term = "s(age)",
    n = AGE_GRID_N,
    type = "central",
    interval = "simultaneous"
  )
  d1$dataset <- nm
  d1
})
deriv_df <- bind_rows(deriv_list) %>%
  tidyr::separate(
    dataset,
    into = c("sex", "reporter"),
    sep = "_",
    extra = "merge"
  ) %>%
  mutate(
    sex = factor(sex, levels = c("female", "male")),
    reporter = factor(reporter, levels = c("parent", "youth"))
  )

# Detect column names (gratia version-robust)
sm_age <- safe_col(smooth_df, "age", ".smooth_covar")
sm_est <- safe_col(smooth_df, ".estimate", "est")
sm_se <- safe_col(smooth_df, ".se", "se")
d_age <- safe_col(deriv_df, "age", ".smooth_covar")
d_deriv <- safe_col(deriv_df, ".derivative", "derivative")
d_se <- safe_col(deriv_df, ".se", "se")

# Population trajectory plot
p_smooth <- ggplot(smooth_df, aes(x = .data[[sm_age]], y = .data[[sm_est]])) +
  geom_ribbon(
    aes(
      ymin = .data[[sm_est]] - 2 * .data[[sm_se]],
      ymax = .data[[sm_est]] + 2 * .data[[sm_se]]
    ),
    alpha = 0.20,
    fill = "steelblue"
  ) +
  geom_line(linewidth = 1, colour = "steelblue4") +
  facet_grid(reporter ~ sex) +
  labs(
    title = "Population PDS composite trajectory (GAMM smooth ± 2 SE)",
    x = "Age (years)",
    y = "PDS composite (population smooth)"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(out_dir, "population_trajectories.png"),
  p_smooth,
  width = 10,
  height = 6,
  dpi = 150
)

# Growth velocity (first derivative) plot
p_deriv <- ggplot(deriv_df, aes(x = .data[[d_age]], y = .data[[d_deriv]])) +
  geom_ribbon(
    aes(
      ymin = .data[[d_deriv]] - 2 * .data[[d_se]],
      ymax = .data[[d_deriv]] + 2 * .data[[d_se]]
    ),
    alpha = 0.20,
    fill = "darkorange"
  ) +
  geom_line(linewidth = 1, colour = "darkorange3") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_grid(reporter ~ sex) +
  labs(
    title = "Pubertal velocity (first derivative of GAMM smooth)",
    x = "Age (years)",
    y = "d(PDS)/d(age) per year"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(out_dir, "growth_velocity.png"),
  p_deriv,
  width = 10,
  height = 6,
  dpi = 150
)

write.csv(
  smooth_df,
  file.path(out_dir, "population_smooth.csv"),
  row.names = FALSE
)
write.csv(
  deriv_df,
  file.path(out_dir, "population_derivative.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# INDIVIDUAL TIMING + TEMPO
#
# In this random-intercept GAMM, each participant's trajectory is a vertical
# shift of the population smooth: fitted_i(age) = smooth(age) + RE_i.
#
# Timing: age (on a fine grid) at which fitted_i first reaches
#   PDS_ONSET_THRESH = 2.0 (midpoint between prepubertal 1 and adult 4).
#   Kids with a large positive RE onset earlier; those with negative RE onset
#   later — so timing varies across individuals.
#
# Tempo: the derivative of the population smooth evaluated at each person's
#   timing age. Because the smooth is not perfectly linear, the growth rate
#   at the threshold crossing differs across timing ages, giving each
#   participant a (partially) distinct tempo value.
#
# Note: tempo here is population-smooth-derived (not fully individual). For
#   truly person-specific tempo variation, extend to a factor-smooth:
#   s(age, id_fac, bs="fs", m=1) — computationally intensive at ABCD N.
# ---------------------------------------------------------------------------
PDS_ONSET_THRESH <- 2.0
AGE_FINE_MIN <- 7
AGE_FINE_MAX <- 20
AGE_FINE_N <- 1000

derive_timing_tempo <- function(gam_obj, df, ds_label) {
  cat("\nDeriving timing/tempo for:", ds_label, "\n")

  age_fine <- seq(AGE_FINE_MIN, AGE_FINE_MAX, length.out = AGE_FINE_N)
  id_levels <- levels(df$id_fac)

  # Population smooth at fine age grid (RE excluded)
  pop_newdata <- data.frame(
    age = age_fine,
    id_fac = factor(id_levels[1], levels = id_levels)
  )
  pop_fit <- predict(
    gam_obj,
    newdata = pop_newdata,
    exclude = "s(id_fac)",
    type = "response"
  )

  # Population smooth derivative via central finite differences
  n <- length(age_fine)
  pop_deriv <- rep(NA_real_, n)
  if (n >= 3) {
    pop_deriv[2:(n - 1)] <-
      (pop_fit[3:n] - pop_fit[1:(n - 2)]) /
      (age_fine[3:n] - age_fine[1:(n - 2)])
  }

  # Extract random intercepts — mgcv stores s(id_fac, bs="re") coefs as
  # s(id_fac).1 … s(id_fac).K, in the same order as levels(id_fac)
  all_coefs <- coef(gam_obj)
  re_idx <- grep("^s\\(id_fac\\)", names(all_coefs))
  re_coefs <- all_coefs[re_idx]

  if (length(re_coefs) != length(id_levels)) {
    warning(
      ds_label,
      ": RE coef count (",
      length(re_coefs),
      ") != id levels (",
      length(id_levels),
      "). Timing/tempo will use population fit for all participants."
    )
    re_map <- setNames(rep(0, length(id_levels)), id_levels)
  } else {
    re_map <- setNames(re_coefs, id_levels)
  }

  # Participants present in data
  participants <- df %>%
    distinct(id, id_str = as.character(id_fac))

  result <- participants %>%
    mutate(
      re = re_map[id_str],
      timing = mapply(
        function(re_i) {
          ind_fit <- pop_fit + re_i
          idx <- which(ind_fit >= PDS_ONSET_THRESH)
          if (length(idx) > 0) age_fine[idx[1]] else NA_real_
        },
        re
      ),
      tempo = mapply(
        function(t_age) {
          if (is.na(t_age)) {
            return(NA_real_)
          }
          pop_deriv[which.min(abs(age_fine - t_age))]
        },
        timing
      ),
      dataset = ds_label
    ) %>%
    select(id, re, timing, tempo, dataset)

  n_na <- sum(is.na(result$timing))
  cat(
    "  N with timing:",
    sum(!is.na(result$timing)),
    "  N without (never crossed threshold):",
    n_na,
    "\n"
  )
  result
}

tt_list <- mapply(
  derive_timing_tempo,
  gamm_fits,
  datasets,
  names(datasets),
  SIMPLIFY = FALSE
)

timing_tempo <- bind_rows(tt_list) %>%
  tidyr::separate(
    dataset,
    into = c("sex", "reporter"),
    sep = "_",
    extra = "merge"
  ) %>%
  mutate(
    sex = factor(sex, levels = c("female", "male")),
    reporter = factor(reporter, levels = c("parent", "youth"))
  )

write.csv(
  timing_tempo,
  file.path(out_dir, "individual_timing_tempo.csv"),
  row.names = FALSE
)
cat("\nSummary of individual timing/tempo:\n")
print(summary(timing_tempo %>% select(sex, reporter, re, timing, tempo)))

# ---------------------------------------------------------------------------
# DISTRIBUTION PLOTS: timing + tempo
# ---------------------------------------------------------------------------
p_timing <- ggplot(
  timing_tempo %>% filter(!is.na(timing)),
  aes(x = timing)
) +
  geom_histogram(binwidth = 0.25, fill = "steelblue", colour = "white") +
  facet_grid(reporter ~ sex) +
  labs(
    title = "Individual pubertal onset age (timing)",
    subtitle = paste0("Threshold: PDS composite = ", PDS_ONSET_THRESH),
    x = "Age at onset (years)",
    y = "Count"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(out_dir, "timing_distribution.png"),
  p_timing,
  width = 10,
  height = 6,
  dpi = 150
)

p_tempo <- ggplot(
  timing_tempo %>% filter(!is.na(tempo)),
  aes(x = tempo)
) +
  geom_histogram(bins = 40, fill = "darkorange", colour = "white") +
  facet_grid(reporter ~ sex) +
  labs(
    title = "Individual pubertal growth rate (tempo)",
    subtitle = "d(PDS)/d(age) of population smooth evaluated at onset age",
    x = "Growth rate at onset (PDS units/year)",
    y = "Count"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(out_dir, "tempo_distribution.png"),
  p_tempo,
  width = 10,
  height = 6,
  dpi = 150
)

p_tt <- ggplot(
  timing_tempo %>% filter(!is.na(timing), !is.na(tempo)),
  aes(x = timing, y = tempo)
) +
  geom_point(alpha = 0.10, size = 0.6) +
  geom_smooth(method = "loess", se = TRUE, colour = "firebrick") +
  facet_grid(reporter ~ sex) +
  labs(
    title = "Timing vs. tempo of pubertal development",
    x = "Onset age (timing, years)",
    y = "Growth rate at onset (tempo)"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(out_dir, "timing_vs_tempo.png"),
  p_tt,
  width = 10,
  height = 6,
  dpi = 150
)

# ---------------------------------------------------------------------------
# TIMING + TEMPO CROSS-REPORTER AGREEMENT  (within sex)
# ---------------------------------------------------------------------------
for (sx in c("female", "male")) {
  wide <- timing_tempo %>%
    filter(sex == sx, !is.na(timing)) %>%
    select(id, reporter, timing, tempo) %>%
    pivot_wider(
      names_from = reporter,
      values_from = c(timing, tempo)
    )

  if (all(c("timing_parent", "timing_youth") %in% names(wide))) {
    cor_t <- cor(
      wide$timing_parent,
      wide$timing_youth,
      use = "pairwise.complete.obs"
    )
    cor_p <- cor(
      wide$tempo_parent,
      wide$tempo_youth,
      use = "pairwise.complete.obs"
    )
    cat(
      "\n",
      sx,
      " — timing parent–youth r =",
      round(cor_t, 3),
      "  tempo parent–youth r =",
      round(cor_p, 3),
      "\n"
    )

    p_agree <- ggplot(wide, aes(x = timing_parent, y = timing_youth)) +
      geom_point(alpha = 0.15, size = 0.6) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        colour = "grey40"
      ) +
      geom_smooth(method = "lm", se = TRUE) +
      labs(
        title = paste0(sx, ": parent vs. youth onset age"),
        subtitle = paste0("r = ", round(cor_t, 3)),
        x = "Parent-report onset age",
        y = "Youth-report onset age"
      ) +
      theme_minimal(base_size = 13)

    ggsave(
      file.path(out_dir, paste0(sx, "_timing_agreement.png")),
      p_agree,
      width = 6,
      height = 5,
      dpi = 150
    )
  }
}

cat("\n\nAll GAMM outputs written to:", out_dir, "\n")
