## 03_gamms_hpc.R
## HPC version: one dataset per job, factor-smooth GAMM for truly individual
## timing and tempo.
##
## Usage:
##   Rscript 03_gamms_hpc.R <dataset_name> [nthreads]
##   dataset_name: female_parent | female_youth | male_parent | male_youth
##   nthreads: integer, default 4
##
## Model:
##   pds_comp ~ s(age, k=6) + s(age, id_fac, bs="fs", k=4, m=2)
##   s(age): population mean trajectory
##   s(age, id_fac, bs="fs"): individual deviations (random smooth slopes)
##   m=2 penalizes second derivative of individual deviation → allows each
##   person's curve to deviate linearly from population (random intercepts
##   AND slopes). This is critical: m=1 would penalize first derivative,
##   collapsing individuals to random intercepts only and making tempo a
##   mathematical function of timing.

## to do: add model for each item (other than pete for females, since its binary. that would need to be a different model?)

suppressPackageStartupMessages({
  library(mgcv)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gratia)
})

set.seed(90025)

# ---------------------------------------------------------------------------
# ARGUMENTS
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
ds_name <- args[1]
n_threads <- if (length(args) >= 2) as.integer(args[2]) else 4L

valid_ds <- c("female_parent", "female_youth", "male_parent", "male_youth")
if (!ds_name %in% valid_ds) {
  stop("dataset_name must be one of: ", paste(valid_ds, collapse = ", "))
}
cat("Dataset:", ds_name, "  nthreads:", n_threads, "\n")

# ---------------------------------------------------------------------------
# PATHS
# Checks: DATA_DIR / OUT_DIR env vars (set in .sh) → Box local fallback
# ---------------------------------------------------------------------------
root_path <- Sys.getenv("HOME_DIR")
if (!nzchar(root_path)) {
  root_path <- Sys.getenv("HOME")
}

data_dir <- Sys.getenv("DATA_DIR")
if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  data_dir <- file.path(
    root_path,
    "projects/abcd-projs/abcd-data-release-6.0/cfm/physical-health/puberty"
  )
}
if (!dir.exists(data_dir)) {
  stop("Cannot locate data directory: ", data_dir)
}

out_base <- Sys.getenv("OUT_DIR")
if (!nzchar(out_base)) {
  out_base <- file.path(
    root_path,
    "Library/CloudStorage/Box-Box/everything/projects/abcd-projs",
    "dissertation/study1/outputs"
  )
}
out_dir <- file.path(out_base, "gamm")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------
csv_path <- file.path(data_dir, paste0(ds_name, "_long.csv"))
if (!file.exists(csv_path)) {
  stop("File not found: ", csv_path)
}

df <- read.csv(csv_path) %>%
  filter(!is.na(pds_comp), !is.na(age)) %>%
  mutate(id_fac = factor(id))

cat(
  "Loaded:",
  nrow(df),
  "rows,",
  n_distinct(df$id),
  "participants,",
  n_distinct(df$wave),
  "waves\n"
)

# ---------------------------------------------------------------------------
# FIT FACTOR-SMOOTH GAMM
#
# s(age, k=6):                        population mean trajectory
# s(age, id_fac, bs="fs", k=6, m=1): person-specific deviations from mean
#   bs="fs" = factor-smooth interaction
#   m=1     = first-derivative penalty on individual deviations
#             (encourages individual curves to be smooth perturbations of mean)
#
# discrete=TRUE + nthreads speeds up bam() substantially at ABCD N.
# Expect ~30–90 min per dataset on Hoffman2 depending on N.
# ---------------------------------------------------------------------------
rds_path <- file.path(out_dir, paste0(ds_name, "_gamm_fs.rds"))

t0 <- proc.time()

if (file.exists(rds_path)) {
  cat("\nLoading existing model from:", rds_path, "\n")
  m <- readRDS(rds_path)
} else {
  cat("\nFitting factor-smooth GAMM (this takes a while)...\n")
  m <- mgcv::bam(
    pds_comp ~ s(age, k = 6) +
      s(age, id_fac, bs = "fs", k = 4, m = 2),
    data = df,
    method = "fREML",
    discrete = TRUE,
    nthreads = n_threads
  )
  elapsed <- (proc.time() - t0)[["elapsed"]]
  cat("Fit complete in", round(elapsed / 60, 1), "min\n")
  saveRDS(m, rds_path)
  cat("Model saved:", rds_path, "\n")
}

# ---------------------------------------------------------------------------
# INDIVIDUAL TIMING + TEMPO
#
# Saved immediately after the model so this CSV is written even if the
# gratia or per-item sections fail (they are more memory-intensive).
# ---------------------------------------------------------------------------
PDS_ONSET_THRESH <- 3.0
AGE_FINE_MIN <- 7
AGE_FINE_MAX <- 20
AGE_FINE_N <- 300

cat(
  "\nBuilding prediction grid (",
  n_distinct(df$id),
  "×",
  AGE_FINE_N,
  "rows)...\n"
)

age_fine <- seq(AGE_FINE_MIN, AGE_FINE_MAX, length.out = AGE_FINE_N)
id_levels <- levels(df$id_fac)
n_ids <- length(id_levels)

pred_grid <- expand.grid(
  age = age_fine,
  id_fac = id_levels,
  stringsAsFactors = FALSE
) %>%
  mutate(id_fac = factor(id_fac, levels = id_levels))

cat("Predicting", nrow(pred_grid), "values...\n")
t1 <- proc.time()

pred_grid$pred <- predict(m, newdata = pred_grid, type = "response")

cat(
  "Prediction complete in",
  round((proc.time() - t1)[["elapsed"]] / 60, 1),
  "min\n"
)

pred_mat <- matrix(pred_grid$pred, nrow = AGE_FINE_N, ncol = n_ids)

da <- diff(age_fine)[1]
deriv_mat <- matrix(NA_real_, nrow = AGE_FINE_N, ncol = n_ids)
deriv_mat[2:(AGE_FINE_N - 1), ] <-
  (pred_mat[3:AGE_FINE_N, ] - pred_mat[1:(AGE_FINE_N - 2), ]) / (2 * da)

cat("Deriving timing and tempo...\n")

# Completion = when predicted PDS reaches 95% of the individual's predicted max.
# tempo = (PDS_at_completion - PDS_at_onset) / (age_at_completion - age_at_onset)
#         i.e. average PDS gain per year during the active pubertal window.
# With m=2 in the fs smooth, individual deviations can be linear in age
# (random intercepts + slopes), so peak_velocity and tempo genuinely vary
# across individuals rather than being fixed by the population curve shape.
PDS_COMPLETION_PCT <- 0.95

timing_vec        <- rep(NA_real_, n_ids)
peak_velocity_vec <- rep(NA_real_, n_ids)
duration_vec      <- rep(NA_real_, n_ids)
tempo_vec         <- rep(NA_real_, n_ids)

for (j in seq_len(n_ids)) {
  traj <- pred_mat[, j]
  d    <- deriv_mat[, j]
  cross <- which(traj >= PDS_ONSET_THRESH)
  if (length(cross) > 0) {
    onset_idx  <- cross[1]
    timing_vec[j] <- age_fine[onset_idx]

    # completion: first point at/after onset where PDS >= 95% of individual max
    ind_max   <- max(traj, na.rm = TRUE)
    compl     <- which(traj[onset_idx:length(traj)] >= PDS_COMPLETION_PCT * ind_max)
    if (length(compl) > 0) {
      compl_idx <- onset_idx + compl[1] - 1L
      dur <- age_fine[compl_idx] - age_fine[onset_idx]
      duration_vec[j] <- dur
      if (dur > 0) {
        tempo_vec[j] <- (traj[compl_idx] - traj[onset_idx]) / dur
      }
    }
  }
  peak_velocity_vec[j] <- max(d, na.rm = TRUE)
}

# Empirical tempo: OLS slope of observed pds_comp ~ age per person.
# This is independent of the GAMM penalty and serves as a direct check.
obs_slopes <- df %>%
  group_by(id) %>%
  filter(n() >= 2L) %>%
  summarise(
    obs_tempo = tryCatch(
      coef(lm(pds_comp ~ age))[["age"]],
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

timing_tempo <- data.frame(
  id            = id_levels,
  timing        = timing_vec,
  peak_velocity = peak_velocity_vec,
  duration      = duration_vec,
  tempo         = tempo_vec,
  dataset       = ds_name,
  stringsAsFactors = FALSE
) %>%
  left_join(obs_slopes, by = "id")

cat(
  "N with timing:",
  sum(!is.na(timing_tempo$timing)),
  "  N without:",
  sum(is.na(timing_tempo$timing)),
  "\n"
)

tt_path <- file.path(out_dir, paste0(ds_name, "_timing_tempo.csv"))
write.csv(timing_tempo, tt_path, row.names = FALSE)
cat("Timing/tempo saved:", tt_path, "\n")

# ---------------------------------------------------------------------------
# INDIVIDUAL TRAJECTORY PLOT
# Sample up to 300 participants and plot their fitted curves as thin lines,
# with the population mean overlaid. Colour by timing (early / on-time / late
# tertiles) so the reader can see how the spread relates to onset age.
# ---------------------------------------------------------------------------
set.seed(90025)
n_sample <- min(300L, n_ids)
sample_idx <- sample(seq_len(n_ids), n_sample)

# long-form data for sampled trajectories
traj_df <- data.frame(
  age    = rep(age_fine, times = n_sample),
  pred   = as.vector(pred_mat[, sample_idx]),
  id     = rep(id_levels[sample_idx], each = AGE_FINE_N)
)

# attach timing for colour coding
timing_lookup <- timing_tempo[, c("id", "timing")]
traj_df <- merge(traj_df, timing_lookup, by = "id", all.x = TRUE)

# tertile-based colour label
tert <- quantile(timing_tempo$timing, probs = c(1/3, 2/3), na.rm = TRUE)
traj_df$timing_group <- cut(
  traj_df$timing,
  breaks = c(-Inf, tert[1], tert[2], Inf),
  labels = c("Early", "On-time", "Late"),
  right  = TRUE
)
traj_df$timing_group[is.na(traj_df$timing_group)] <- "No onset"

# population mean from pred_mat (mean across ALL participants, not just sample)
pop_mean_df <- data.frame(
  age  = age_fine,
  pred = rowMeans(pred_mat)
)

p_traj <- ggplot(traj_df, aes(x = age, y = pred, group = id,
                               colour = timing_group)) +
  geom_line(alpha = 0.18, linewidth = 0.3) +
  geom_line(data = pop_mean_df, aes(x = age, y = pred, group = NULL),
            colour = "black", linewidth = 1.2, inherit.aes = FALSE) +
  geom_hline(yintercept = PDS_ONSET_THRESH,
             linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  scale_colour_manual(
    values = c(Early = "#d73027", `On-time` = "#4575b4",
               Late  = "#1a9850", `No onset` = "grey60"),
    na.value = "grey60"
  ) +
  labs(
    title    = paste("Individual pubertal trajectories:", ds_name),
    subtitle = paste0("n = ", n_sample, " sampled; black = population mean; ",
                      "dashed = onset threshold (PDS = ", PDS_ONSET_THRESH, ")"),
    x        = "Age (years)",
    y        = "Fitted PDS composite",
    colour   = "Timing"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(
  file.path(out_dir, paste0(ds_name, "_individual_trajectories.png")),
  p_traj,
  width  = 8,
  height = 5,
  dpi    = 150
)
cat("Individual trajectory plot saved.\n")

p_t <- ggplot(timing_tempo %>% filter(!is.na(timing)), aes(x = timing)) +
  geom_histogram(binwidth = 0.25, fill = "steelblue", colour = "white") +
  labs(
    title = paste("Onset age distribution:", ds_name),
    subtitle = paste0("Threshold: PDS = ", PDS_ONSET_THRESH),
    x = "Age at onset (years)",
    y = "Count"
  ) +
  theme_minimal(base_size = 13)
ggsave(
  file.path(out_dir, paste0(ds_name, "_timing_dist.png")),
  p_t,
  width = 6,
  height = 4,
  dpi = 150
)

p_v <- ggplot(
  timing_tempo %>% filter(!is.na(peak_velocity)),
  aes(x = peak_velocity)
) +
  geom_histogram(bins = 40, fill = "darkorange", colour = "white") +
  labs(
    title = paste("Peak velocity distribution:", ds_name),
    x = "Max d(PDS)/d(age) per year",
    y = "Count"
  ) +
  theme_minimal(base_size = 13)
ggsave(
  file.path(out_dir, paste0(ds_name, "_peak_velocity_dist.png")),
  p_v,
  width = 6,
  height = 4,
  dpi = 150
)

p_tv <- ggplot(
  timing_tempo %>% filter(!is.na(timing), !is.na(peak_velocity)),
  aes(x = timing, y = peak_velocity)
) +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_smooth(method = "loess", se = TRUE, colour = "firebrick") +
  labs(
    title = paste("Timing vs. peak velocity:", ds_name),
    x = "Onset age (years)",
    y = "Peak d(PDS)/d(age)"
  ) +
  theme_minimal(base_size = 13)
ggsave(
  file.path(out_dir, paste0(ds_name, "_timing_vs_velocity.png")),
  p_tv,
  width = 6,
  height = 4,
  dpi = 150
)

# ---------------------------------------------------------------------------
# POPULATION SMOOTH + DERIVATIVE  (gratia)
# Exclude the fs term to get the mean trajectory only.
# Wrapped in tryCatch so an OOM or gratia error does not abort per-item models.
# ---------------------------------------------------------------------------
tryCatch({

AGE_GRID_N <- 300

# Try column names in order; error with diagnostics if none found
safe_col <- function(df, ...) {
  candidates <- c(...)
  found <- intersect(candidates, names(df))
  if (length(found) == 0) {
    stop(
      "None of (",
      paste(candidates, collapse = ", "),
      ") found. Columns present: ",
      paste(names(df), collapse = ", ")
    )
  }
  found[1]
}

sm <- gratia::smooth_estimates(m, smooth = "s(age)", n = AGE_GRID_N)
cat("smooth_estimates columns:", paste(names(sm), collapse = ", "), "\n")

# gratia 0.8.x: covariate in "age" or "data"; estimate/se without dots
# gratia 0.9.x+: covariate in ".smooth_covar"; estimate/se with dots
sm_age <- safe_col(sm, "age", "data", ".smooth_covar")
sm_est <- safe_col(sm, "est", ".estimate")
sm_se <- safe_col(sm, "se", ".se")

d1 <- gratia::derivatives(
  m,
  term = "s(age)",
  n = AGE_GRID_N,
  type = "central",
  interval = "simultaneous"
)
cat("derivatives columns:", paste(names(d1), collapse = ", "), "\n")

# gratia 0.8.x: covariate in "data" (with "var" holding name); 0.9.x+: actual name
d_age <- safe_col(d1, "data", "age", ".smooth_covar")
d_deriv <- safe_col(d1, "derivative", ".derivative")
d_se <- safe_col(d1, "se", ".se")

write.csv(
  sm %>% mutate(dataset = ds_name),
  file.path(out_dir, paste0(ds_name, "_population_smooth.csv")),
  row.names = FALSE
)
write.csv(
  d1 %>% mutate(dataset = ds_name),
  file.path(out_dir, paste0(ds_name, "_population_derivative.csv")),
  row.names = FALSE
)

# Population smooth plot
p_smooth <- ggplot(sm, aes(x = .data[[sm_age]], y = .data[[sm_est]])) +
  geom_ribbon(
    aes(
      ymin = .data[[sm_est]] - 2 * .data[[sm_se]],
      ymax = .data[[sm_est]] + 2 * .data[[sm_se]]
    ),
    alpha = 0.2,
    fill = "steelblue"
  ) +
  geom_line(linewidth = 1, colour = "steelblue4") +
  labs(
    title = paste("Population PDS trajectory:", ds_name),
    x = "Age (years)",
    y = "PDS composite (population smooth)"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(out_dir, paste0(ds_name, "_population_smooth.png")),
  p_smooth,
  width = 7,
  height = 5,
  dpi = 150
)

# Velocity plot
p_deriv <- ggplot(d1, aes(x = .data[[d_age]], y = .data[[d_deriv]])) +
  geom_ribbon(
    aes(
      ymin = .data[[d_deriv]] - 2 * .data[[d_se]],
      ymax = .data[[d_deriv]] + 2 * .data[[d_se]]
    ),
    alpha = 0.2,
    fill = "darkorange"
  ) +
  geom_line(linewidth = 1, colour = "darkorange3") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(
    title = paste("Pubertal velocity:", ds_name),
    x = "Age (years)",
    y = "d(PDS)/d(age) per year"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(out_dir, paste0(ds_name, "_growth_velocity.png")),
  p_deriv,
  width = 7,
  height = 5,
  dpi = 150
)

}, error = function(e) {
  cat("\n[WARNING] Gratia section failed for", ds_name, "—", conditionMessage(e), "\n")
  cat("Skipping population smooth/derivative plots. Continuing to per-item models.\n\n")
})

# ---------------------------------------------------------------------------
# PER-ITEM MODELS
#
# Edit the item lists below to match your CSV column names.
# Any item in `binary_items` gets family = binomial(); all others gaussian().
# Timing/tempo is skipped for binary items (log-odds scale not meaningful).
# ---------------------------------------------------------------------------

female_items <- c("peta", "petb", "petc", "petd", "fpete")
male_items <- c("peta", "petb", "petc", "petd", "mpete")
binary_items <- c("fpete")

sex <- if (startsWith(ds_name, "female")) "female" else "male"
items <- if (sex == "female") female_items else male_items

for (item in items) {
  if (!item %in% names(df)) {
    cat("\nSkipping", item, "— column not found in data\n")
    next
  }

  df_item <- df %>% filter(!is.na(.data[[item]]))
  is_binary <- item %in% binary_items
  fam <- if (is_binary) binomial() else gaussian()

  item_rds <- file.path(out_dir, paste0(ds_name, "_", item, "_gamm_fs.rds"))

  cat("\n---", item, "| n =", nrow(df_item), "| family:", fam$family, "---\n")

  if (file.exists(item_rds)) {
    cat("Loading existing model from:", item_rds, "\n")
    m_item <- readRDS(item_rds)
  } else {
    t_item <- proc.time()
    f_item <- as.formula(paste0(
      item,
      " ~ s(age, k = 6) + s(age, id_fac, bs = 'fs', k = 4, m = 1)"
    ))
    m_item <- mgcv::bam(
      f_item,
      data = df_item,
      family = fam,
      method = "fREML",
      discrete = TRUE,
      nthreads = n_threads
    )
    cat(
      item,
      "fit in",
      round((proc.time() - t_item)[["elapsed"]] / 60, 1),
      "min\n"
    )
    saveRDS(m_item, item_rds)
  }

  # Population smooth
  sm_i <- gratia::smooth_estimates(m_item, smooth = "s(age)", n = AGE_GRID_N)
  write.csv(
    sm_i %>% mutate(dataset = ds_name, item = item),
    file.path(out_dir, paste0(ds_name, "_", item, "_smooth.csv")),
    row.names = FALSE
  )

  sm_i_age <- safe_col(sm_i, "age", "data", ".smooth_covar")
  sm_i_est <- safe_col(sm_i, "est", ".estimate")
  sm_i_se <- safe_col(sm_i, "se", ".se")

  y_lab <- if (is_binary) {
    "Log-odds (population smooth)"
  } else {
    "Score (population smooth)"
  }

  p_i <- ggplot(sm_i, aes(x = .data[[sm_i_age]], y = .data[[sm_i_est]])) +
    geom_ribbon(
      aes(
        ymin = .data[[sm_i_est]] - 2 * .data[[sm_i_se]],
        ymax = .data[[sm_i_est]] + 2 * .data[[sm_i_se]]
      ),
      alpha = 0.2,
      fill = "steelblue"
    ) +
    geom_line(linewidth = 1, colour = "steelblue4") +
    labs(
      title = paste("PDS item:", item, "|", ds_name),
      x = "Age (years)",
      y = y_lab
    ) +
    theme_minimal(base_size = 13)

  ggsave(
    file.path(out_dir, paste0(ds_name, "_", item, "_smooth.png")),
    p_i,
    width = 7,
    height = 5,
    dpi = 150
  )

  # Derivatives + velocity plot (continuous items only)
  if (!is_binary) {
    d_i <- gratia::derivatives(
      m_item,
      term = "s(age)",
      n = AGE_GRID_N,
      type = "central",
      interval = "simultaneous"
    )
    write.csv(
      d_i %>% mutate(dataset = ds_name, item = item),
      file.path(out_dir, paste0(ds_name, "_", item, "_derivative.csv")),
      row.names = FALSE
    )

    d_i_age <- safe_col(d_i, "data", "age", ".smooth_covar")
    d_i_deriv <- safe_col(d_i, "derivative", ".derivative")
    d_i_se <- safe_col(d_i, "se", ".se")

    p_d_i <- ggplot(d_i, aes(x = .data[[d_i_age]], y = .data[[d_i_deriv]])) +
      geom_ribbon(
        aes(
          ymin = .data[[d_i_deriv]] - 2 * .data[[d_i_se]],
          ymax = .data[[d_i_deriv]] + 2 * .data[[d_i_se]]
        ),
        alpha = 0.2,
        fill = "darkorange"
      ) +
      geom_line(linewidth = 1, colour = "darkorange3") +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      labs(
        title = paste("Velocity:", item, "|", ds_name),
        x = "Age (years)",
        y = "d(item)/d(age) per year"
      ) +
      theme_minimal(base_size = 13)

    ggsave(
      file.path(out_dir, paste0(ds_name, "_", item, "_velocity.png")),
      p_d_i,
      width = 7,
      height = 5,
      dpi = 150
    )
  }
}

cat("\nPer-item models complete.\n")

cat(
  "\nDone. All outputs written to:",
  out_dir,
  "\nElapsed total:",
  round((proc.time() - t0)[["elapsed"]] / 60, 1),
  "min\n"
)
