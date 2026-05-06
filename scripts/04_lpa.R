## 04_lpa.R
## Latent Profile Analysis on GAMM-derived pubertal timing and tempo.
## Uses tidyLPA (wrapping mclust) for model comparison and profile extraction.
##
## Inputs:  {ds_name}_timing_tempo.csv  (from 03_gamms_hpc.R)
## Outputs: model comparison tables, profile assignments, plots
##
## Runs 6 LPA solutions:
##   female_parent, female_youth, male_parent, male_youth  (sex × reporter)
##   female_averaged, male_averaged                        (reporter-averaged)
##
## Model selection: compare covariance structures × K profiles; BIC-optimal
## solution is saved; comparison plots allow manual override.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tidyLPA)
})

set.seed(90025)

LPA_VARS <- c("timing", "peak_velocity", "tempo_at_onset")
N_PROFILES <- 1:6
# tidyLPA model numbers (mclust parameterisation):
#   1 = EEI  equal variance, zero covariance
#   2 = VVI  varying variance, zero covariance
#   3 = EEE  equal variance + covariance
#   6 = VVV  varying variance + covariance  (most flexible)
LPA_MODELS <- c(1, 2, 3, 6)

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
    "projects/abcd-projs",
    "dissertation/study1/outputs"
  )
}
gamm_dir <- file.path(out_base, "gamm")
out_dir <- file.path(out_base, "lpa")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# LOAD timing/tempo from all four datasets
# ---------------------------------------------------------------------------
ds_names <- c("female_parent", "female_youth", "male_parent", "male_youth")

tt_list <- lapply(ds_names, function(ds) {
  path <- file.path(gamm_dir, paste0(ds, "_timing_tempo.csv"))
  if (!file.exists(path)) {
    warning("Not found: ", path)
    return(NULL)
  }
  read.csv(path)
})
names(tt_list) <- ds_names

# ---------------------------------------------------------------------------
# BUILD AVERAGED DATASETS (parent + youth averaged per person per sex)
# ---------------------------------------------------------------------------
build_averaged <- function(sex) {
  p <- tt_list[[paste0(sex, "_parent")]]
  y <- tt_list[[paste0(sex, "_youth")]]

  if (is.null(p) && is.null(y)) {
    return(NULL)
  }
  if (is.null(p)) {
    return(y %>% mutate(reporter = "youth_only"))
  }
  if (is.null(y)) {
    return(p %>% mutate(reporter = "parent_only"))
  }

  full_join(
    p %>%
      select(id, all_of(LPA_VARS)) %>%
      rename_with(~ paste0(., "_p"), all_of(LPA_VARS)),
    y %>%
      select(id, all_of(LPA_VARS)) %>%
      rename_with(~ paste0(., "_y"), all_of(LPA_VARS)),
    by = "id"
  ) %>%
    mutate(
      timing = rowMeans(cbind(timing_p, timing_y), na.rm = TRUE),
      peak_velocity = rowMeans(
        cbind(peak_velocity_p, peak_velocity_y),
        na.rm = TRUE
      ),
      tempo_at_onset = rowMeans(
        cbind(tempo_at_onset_p, tempo_at_onset_y),
        na.rm = TRUE
      ),
      dataset = paste0(sex, "_averaged"),
      reporter = "averaged"
    ) %>%
    select(id, dataset, reporter, all_of(LPA_VARS))
}

tt_list[["female_averaged"]] <- build_averaged("female")
tt_list[["male_averaged"]] <- build_averaged("male")

all_labels <- c(ds_names, "female_averaged", "male_averaged")

# ---------------------------------------------------------------------------
# LPA HELPER
# ---------------------------------------------------------------------------
run_lpa <- function(df, label) {
  df_cc <- df %>% filter(complete.cases(.[LPA_VARS]))
  cat("\n", strrep("=", 60), "\n")
  cat("LPA:", label, "| complete cases:", nrow(df_cc), "\n")
  cat(strrep("=", 60), "\n")

  if (nrow(df_cc) < 50) {
    warning("Too few complete cases for ", label, " — skipping.")
    return(invisible(NULL))
  }

  # Standardise within this sample (LPA on z-scores; raw means saved separately)
  df_z <- df_cc %>%
    mutate(across(all_of(LPA_VARS), ~ as.numeric(scale(.))))

  # --- Model comparison ---------------------------------------------------
  # Estimate each model/k individually so failed mclust fits don't crash the
  # whole comparison (compare_solutions() chokes on NULL fits with empty names).
  cat("Comparing models (this may take a few minutes)...\n")
  successful_fits <- list()
  for (.m in LPA_MODELS) {
    for (.k in N_PROFILES) {
      .key <- paste0("m", .m, "_k", .k)
      .fit <- tryCatch(
        suppressWarnings(
          df_z %>%
            select(all_of(LPA_VARS)) %>%
            estimate_profiles(n_profiles = .k, models = .m)
        ),
        error = function(e) NULL
      )
      if (!is.null(.fit)) {
        .stats <- tryCatch(get_fit(.fit), error = function(e) NULL)
        if (!is.null(.stats) && nrow(.stats) > 0) {
          # Column names vary by tidyLPA/tibble version; match by pattern.
          .bic_col <- grep("^BIC", names(.stats), value = TRUE)[1]
          .ent_col <- grep(
            "(?i)^entropy",
            names(.stats),
            value = TRUE,
            perl = TRUE
          )[1]
          .bic_val <- if (!is.na(.bic_col)) {
            as.numeric(.stats[[.bic_col]][1])
          } else {
            NA_real_
          }
          .ent_val <- if (!is.na(.ent_col)) {
            as.numeric(.stats[[.ent_col]][1])
          } else {
            NA_real_
          }
          successful_fits[[.key]] <- list(
            fit = .fit,
            model = .m,
            k = .k,
            BIC = .bic_val,
            Entropy = .ent_val
          )
        }
      }
    }
  }

  if (length(successful_fits) == 0) {
    warning("No models converged for ", label, " — skipping.")
    return(invisible(NULL))
  }

  fits_df <- bind_rows(lapply(successful_fits, function(x) {
    data.frame(
      Model = x$model,
      Classes = x$k,
      BIC = x$BIC,
      Entropy = x$Entropy,
      stringsAsFactors = FALSE
    )
  }))

  cat("\nFit statistics:\n")
  print(fits_df)

  write.csv(
    fits_df,
    file.path(out_dir, paste0(label, "_model_comparison.csv")),
    row.names = FALSE
  )

  # BIC plot (standard BIC: lower = better)
  p_bic <- ggplot(
    fits_df,
    aes(x = Classes, y = BIC, colour = factor(Model), group = factor(Model))
  ) +
    geom_line() +
    geom_point(size = 2) +
    scale_x_continuous(breaks = N_PROFILES) +
    labs(
      title = paste("BIC by profiles and model —", label),
      x = "Number of profiles",
      y = "BIC (lower = better)",
      colour = "Model"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  # Entropy plot
  p_ent <- ggplot(
    fits_df %>% filter(Classes > 1),
    aes(x = Classes, y = Entropy, colour = factor(Model), group = factor(Model))
  ) +
    geom_line() +
    geom_point(size = 2) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(breaks = N_PROFILES) +
    labs(
      title = paste("Entropy by profiles and model —", label),
      x = "Number of profiles",
      y = "Entropy (>0.80 = good separation)",
      colour = "Model"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(
    file.path(out_dir, paste0(label, "_bic_plot.png")),
    p_bic,
    width = 8,
    height = 5,
    dpi = 150
  )
  ggsave(
    file.path(out_dir, paste0(label, "_entropy_plot.png")),
    p_ent,
    width = 8,
    height = 5,
    dpi = 150
  )

  # --- Best solution (BIC-optimal) ----------------------------------------
  # get_fit() returns standard BIC (lower = better); exclude K=1 since a
  # single-component solution is not a meaningful LPA result.
  valid_fits <- Filter(function(x) !is.na(x$BIC) && x$k > 1, successful_fits)
  if (length(valid_fits) == 0) {
    warning("No valid multi-profile BIC values for ", label, " — skipping.")
    return(invisible(NULL))
  }
  best_key <- names(which.min(sapply(valid_fits, `[[`, "BIC")))
  best_k <- valid_fits[[best_key]]$k
  best_model <- valid_fits[[best_key]]$model
  fit_best <- valid_fits[[best_key]]$fit
  cat("\nBIC-optimal solution: K =", best_k, "| model =", best_model, "\n")

  assigns <- get_data(fit_best) %>%
    select(Class, starts_with("CPROB")) %>%
    bind_cols(df_cc %>% select(id, all_of(LPA_VARS)))

  write.csv(
    assigns,
    file.path(out_dir, paste0(label, "_k", best_k, "_assignments.csv")),
    row.names = FALSE
  )

  # --- Profile means plot (raw scale) -------------------------------------
  profile_means <- assigns %>%
    group_by(Class) %>%
    summarise(
      n = n(),
      across(
        all_of(LPA_VARS),
        list(mean = mean, se = ~ sd(.) / sqrt(n())),
        .names = "{.col}__{.fn}"
      ),
      .groups = "drop"
    )

  write.csv(
    profile_means,
    file.path(out_dir, paste0(label, "_k", best_k, "_profile_means.csv")),
    row.names = FALSE
  )

  # Pivot for plotting
  pm_long <- profile_means %>%
    pivot_longer(
      cols = ends_with("__mean") | ends_with("__se"),
      names_to = c("variable", ".value"),
      names_sep = "__"
    ) %>%
    mutate(
      variable = factor(
        variable,
        levels = LPA_VARS,
        labels = c("Timing (age at onset)", "Peak velocity", "Tempo at onset")
      ),
      Class = factor(Class)
    )

  p_prof <- ggplot(
    pm_long,
    aes(x = variable, y = mean, colour = Class, group = Class)
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = mean - 2 * se, ymax = mean + 2 * se),
      width = 0.15,
      linewidth = 0.6
    ) +
    labs(
      title = paste0(
        "Pubertal profiles — ",
        label,
        " (K = ",
        best_k,
        ", model ",
        best_model,
        ")"
      ),
      x = NULL,
      y = "Mean (raw scale)",
      colour = "Profile",
      caption = paste0(
        "n per profile: ",
        paste(profile_means$n, collapse = " / ")
      )
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      legend.position = "bottom"
    )

  ggsave(
    file.path(out_dir, paste0(label, "_k", best_k, "_profile_plot.png")),
    p_prof,
    width = 8,
    height = 5,
    dpi = 150
  )

  # --- Also save K = 2,3,4 of best model type for inspection --------------
  for (k in setdiff(2:4, best_k)) {
    .key_k <- paste0("m", best_model, "_k", k)
    fit_k <- if (!is.null(successful_fits[[.key_k]])) {
      successful_fits[[.key_k]]$fit
    } else {
      tryCatch(
        suppressWarnings(
          df_z %>%
            select(all_of(LPA_VARS)) %>%
            estimate_profiles(n_profiles = k, models = best_model)
        ),
        error = function(e) NULL
      )
    }
    if (is.null(fit_k)) {
      next
    }

    asgn_k <- get_data(fit_k) %>%
      select(Class, starts_with("CPROB")) %>%
      bind_cols(df_cc %>% select(id, all_of(LPA_VARS)))

    write.csv(
      asgn_k,
      file.path(out_dir, paste0(label, "_k", k, "_assignments.csv")),
      row.names = FALSE
    )
  }

  cat("Outputs written for", label, "\n")
  invisible(list(
    comparison = fits_df,
    best_fit = fit_best,
    best_k = best_k,
    best_model = best_model
  ))
}

# ---------------------------------------------------------------------------
# RUN ALL SIX LPAs
# ---------------------------------------------------------------------------
results <- list()
for (label in all_labels) {
  df <- tt_list[[label]]
  if (is.null(df)) {
    next
  }
  results[[label]] <- run_lpa(df, label)
}

# ---------------------------------------------------------------------------
# SUMMARY TABLE across all solutions
# ---------------------------------------------------------------------------
summary_rows <- lapply(names(results), function(label) {
  r <- results[[label]]
  if (is.null(r)) {
    return(NULL)
  }
  data.frame(
    label = label,
    best_K = r$best_k,
    best_model = r$best_model,
    stringsAsFactors = FALSE
  )
})
summary_tbl <- bind_rows(Filter(Negate(is.null), summary_rows))
cat("\n=== BIC-optimal solutions across all runs ===\n")
print(summary_tbl)
write.csv(
  summary_tbl,
  file.path(out_dir, "lpa_best_solutions_summary.csv"),
  row.names = FALSE
)

cat("\nDone. All outputs written to:", out_dir, "\n")
