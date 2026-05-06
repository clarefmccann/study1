## 02_psychometrics.R
## Psychometric evaluation of PDS items across waves, reporters, and sex.
## Answers two questions before downstream modelling:
##   (1) Which items are informative? (GRM discrimination + item info)
##   (2) Do we need both parent and youth items, or are they redundant?
##       (cross-reporter bifactor CFA)
## Also runs longitudinal measurement invariance and invariance across
## age tertiles, race/ethnicity, and BMI tertiles.
## Requires outputs from 00_data_foundation.R.

# export DATA_DIR="/u/project/silvers/data/ABCD/ABCD-release-6.0/cfm/physical-health/puberty"
# export OUT_DIR="/u/home/c/clarefmc/projects/abcd-projs/dissertation/study1/outputs"
# Rscript 02_psychometrics.R

pacman::p_load(
  dplyr,
  tidyr,
  tibble,
  ggplot2,
  psych,
  mirt,
  lavaan,
  semTools,
  install = TRUE
)

set.seed(90025)

# safe print: ensures tibble dispatch so n = Inf works
print_all <- function(x) print(as_tibble(x), n = Inf)

# safe fit-measures row: avoids t() matrix issues with vctrs bind_rows
fit_to_row <- function(fit, fi_names) {
  if (is.null(fit)) {
    return(NULL)
  }
  fi <- tryCatch(lavaan::fitMeasures(fit, fi_names), error = function(e) NULL)
  if (is.null(fi) || length(fi) == 0) {
    return(NULL)
  }
  do.call(data.frame, c(as.list(fi), list(stringsAsFactors = FALSE)))
}

# ---------------------------------------------------------------------------
# PATHS
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
pub_root <- data_dir

out_base <- Sys.getenv("OUT_DIR")
if (!nzchar(out_base)) {
  out_base <- file.path(
    root_path,
    "projects/abcd-projs",
    "dissertation/study1/outputs"
  )
}
out_dir <- file.path(out_base, "psychometrics")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# LOAD
# ---------------------------------------------------------------------------
female_parent <- read.csv(file.path(pub_root, "female_parent_long.csv"))
female_youth <- read.csv(file.path(pub_root, "female_youth_long.csv"))
male_parent <- read.csv(file.path(pub_root, "male_parent_long.csv"))
male_youth <- read.csv(file.path(pub_root, "male_youth_long.csv"))

wave_order <- c("bl", "fu1", "fu2", "fu3", "fu4", "fu5", "fu6")

# fpete (menarche, binary 1/2) is female-only; mpete (genital dev, ordinal 1-4) is male-only.
# Invariance sections and bifactor use the shared 4-item ordinal set only.
ordinal_items <- c("peta", "petb", "petc", "petd")
female_items <- c("peta", "petb", "petc", "petd", "fpete")
female_binary_items <- "fpete"
male_items <- c("peta", "petb", "petc", "petd", "mpete")
male_binary_items <- character(0)

for (nm in c("female_parent", "female_youth", "male_parent", "male_youth")) {
  d <- get(nm)
  d$wave <- factor(d$wave, levels = wave_order)
  assign(nm, d)
}

# ---------------------------------------------------------------------------
# SECTION 1: PER-WAVE PSYCHOMETRICS
# For each sex × reporter × wave:
#   polychoric correlations, alpha, omega, EFA (1-factor), CFA (WLSMV), GRM
# ---------------------------------------------------------------------------

datasets <- list(
  female_parent = list(
    df = female_parent,
    sex = "female",
    reporter = "parent",
    items = female_items,
    binary = female_binary_items
  ),
  female_youth = list(
    df = female_youth,
    sex = "female",
    reporter = "youth",
    items = female_items,
    binary = female_binary_items
  ),
  male_parent = list(
    df = male_parent,
    sex = "male",
    reporter = "parent",
    items = male_items,
    binary = male_binary_items
  ),
  male_youth = list(
    df = male_youth,
    sex = "male",
    reporter = "youth",
    items = male_items,
    binary = male_binary_items
  )
)

wave_results <- list()
disc_rows <- list() # GRM discrimination params, collected for summary table
iif_plots <- list() # item information function plots

for (ds_name in names(datasets)) {
  ds <- datasets[[ds_name]]
  df <- ds$df
  sex_lab <- ds$sex
  rep_lab <- ds$reporter
  items_bin <- ds$binary
  items_ord <- setdiff(ds$items, items_bin)

  for (wv in wave_order) {
    key <- paste(ds_name, wv, sep = "_")

    sub_raw <- df %>%
      filter(wave == wv) %>%
      select(any_of(ds$items)) %>%
      mutate(across(everything(), as.integer))

    items_ord_avail <- intersect(items_ord, names(sub_raw))
    items_bin_avail <- intersect(items_bin, names(sub_raw))
    items_avail <- c(items_ord_avail, items_bin_avail)

    sub <- sub_raw %>%
      filter(if_all(all_of(items_ord_avail), ~ !is.na(.) & . %in% 1:4))
    if (length(items_bin_avail) > 0) {
      sub <- sub %>%
        filter(if_all(all_of(items_bin_avail), ~ !is.na(.) & . %in% 1:2))
    }
    sub <- select(sub, all_of(items_avail))

    if (nrow(sub) < 100) {
      message("Skipping ", key, ": n = ", nrow(sub))
      next
    }

    result <- list(n = nrow(sub))

    # --- polychoric + reliability ---
    pc <- tryCatch(
      psych::polychoric(sub, na.rm = TRUE, correct = 0),
      error = function(e) NULL
    )
    if (!is.null(pc)) {
      result$poly_rho <- pc$rho
      result$alpha <- tryCatch(
        psych::alpha(sub, na.rm = TRUE, check.keys = TRUE)$total,
        error = function(e) NULL
      )
      # omega_h requires >= 2 factors; skip for per-wave (used in bifactor section)
    }

    # --- EFA 1-factor (polychoric) ---
    result$efa <- tryCatch(
      psych::fa(
        sub,
        nfactors = 1,
        fm = "minres",
        rotate = "none",
        cor = "poly"
      ),
      error = function(e) NULL
    )

    # --- CFA 1-factor WLSMV ---
    # WLSMV treats all ordered vars correctly regardless of number of categories.
    cfa_syntax <- paste("pub =~", paste(items_avail, collapse = " + "))
    result$cfa <- tryCatch(
      lavaan::cfa(
        cfa_syntax,
        data = sub,
        ordered = items_avail,
        estimator = "WLSMV",
        missing = "pairwise"
      ),
      error = function(e) NULL
    )

    # --- GRM (ordinal items only; binary items such as fpete are excluded) ---
    grm_fit <- tryCatch(
      {
        mat <- as.matrix(sub[items_ord_avail]) - 1L # mirt expects 0-based
        mirt::mirt(mat, 1, itemtype = "graded", method = "EM", verbose = FALSE)
      },
      error = function(e) {
        message("GRM failed [", key, "]: ", e$message)
        NULL
      }
    )

    if (!is.null(grm_fit)) {
      result$grm <- grm_fit
      coefs <- tryCatch(
        mirt::coef(grm_fit, IRTpars = TRUE, simplify = TRUE)$items,
        error = function(e) NULL
      )
      if (!is.null(coefs)) {
        result$grm_coefs <- coefs
        coef_df <- as.data.frame(coefs)
        # mirt names discrimination 'a' with IRTpars=TRUE (1D); fallback to 'a1'
        a_col <- intersect(c("a", "a1"), names(coef_df))[1]
        if (is.na(a_col)) {
          a_col <- names(coef_df)[1]
        }
        for (col in c("b1", "b2", "b3")) {
          if (!col %in% names(coef_df)) coef_df[[col]] <- NA_real_
        }
        disc_rows[[key]] <- data.frame(
          item = rownames(coef_df),
          a = coef_df[[a_col]],
          b1 = coef_df[["b1"]],
          b2 = coef_df[["b2"]],
          b3 = coef_df[["b3"]],
          wave = wv,
          sex = sex_lab,
          reporter = rep_lab,
          dataset = ds_name,
          stringsAsFactors = FALSE
        )
      }

      # item information function plot
      iif_plots[[key]] <- tryCatch(
        {
          p <- plot(
            grm_fit,
            type = "infotrace",
            main = paste("IIF:", ds_name, wv)
          )
          p
        },
        error = function(e) NULL
      )
    }

    wave_results[[key]] <- result
  }
}

# --- GRM discrimination summary table ---
disc_table <- bind_rows(disc_rows)
cat("\n=== GRM item discrimination (a) by wave × reporter × sex ===\n")
print_all(disc_table)
write.csv(
  disc_table,
  file.path(out_dir, "grm_discrimination.csv"),
  row.names = FALSE
)

# flag items with consistently low discrimination (a < 0.7 at majority of waves)
if (nrow(disc_table) > 0) {
  low_disc <- disc_table %>%
    group_by(item, sex, reporter) %>%
    summarise(
      mean_a = mean(a, na.rm = TRUE),
      pct_low = mean(a < 0.7, na.rm = TRUE),
      n_waves = n(),
      .groups = "drop"
    ) %>%
    arrange(sex, reporter, item)
} else {
  warning(
    "No GRM results — mirt may not be installed correctly. Skipping discrimination summary."
  )
  low_disc <- data.frame()
}

cat(
  "\n=== Item discrimination summary (flag if a < 0.7 at majority of waves) ===\n"
)
print_all(low_disc)
write.csv(
  low_disc,
  file.path(out_dir, "item_discrimination_summary.csv"),
  row.names = FALSE
)

# --- CFA fit summary table ---
fit_indices <- c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr")
cfa_fits <- lapply(names(wave_results), function(key) {
  fit <- wave_results[[key]]$cfa
  if (is.null(fit)) {
    return(NULL)
  }
  fi <- tryCatch(lavaan::fitMeasures(fit, fit_indices), error = function(e) {
    NULL
  })
  if (is.null(fi) || length(fi) < length(fit_indices)) {
    return(NULL)
  }
  parts <- strsplit(key, "_")[[1]]
  tibble::tibble(
    dataset = paste(parts[1], parts[2], sep = "_"),
    wave = tail(parts, 1),
    cfi.scaled = unname(fi["cfi.scaled"]),
    tli.scaled = unname(fi["tli.scaled"]),
    rmsea.scaled = unname(fi["rmsea.scaled"]),
    srmr = unname(fi["srmr"])
  )
}) %>%
  dplyr::bind_rows()

cat("\n=== CFA fit indices per wave × reporter × sex ===\n")
print_all(cfa_fits)
write.csv(
  cfa_fits,
  file.path(out_dir, "cfa_fit_per_wave.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# SECTION 2: CROSS-REPORTER BIFACTOR CFA
# Pooled across waves. Per sex.
# Tests whether a general puberty factor + reporter method factors
# fits better than two separate reporter factors.
# ---------------------------------------------------------------------------

build_cr_data <- function(parent_df, youth_df, include_pete = FALSE) {
  base_items <- c("peta", "petb", "petc", "petd")
  # fpete (female) and mpete (male) are sex-specific; any_of() picks whichever exists
  all_items <- if (include_pete) c(base_items, "fpete", "mpete") else base_items

  parent_sel <- parent_df %>%
    select(id, wave, any_of(all_items)) %>%
    rename_with(~ paste0(., "_p"), any_of(all_items))

  youth_sel <- youth_df %>%
    select(id, wave, any_of(all_items)) %>%
    rename_with(~ paste0(., "_y"), any_of(all_items))

  inner_join(parent_sel, youth_sel, by = c("id", "wave")) %>%
    select(-id, -wave) %>%
    mutate(across(everything(), as.integer)) %>%
    # binary pete (1/2) passes %in% 1:4, so no special case needed here
    filter(if_all(everything(), ~ !is.na(.) & . %in% 1:4))
}

run_bifactor_section <- function(
  parent_df,
  youth_df,
  sex_label,
  include_pete = FALSE
) {
  dat <- build_cr_data(parent_df, youth_df, include_pete = include_pete)
  cat(
    "\n\n=== Cross-reporter bifactor | sex:",
    sex_label,
    "| n =",
    nrow(dat),
    "===\n"
  )

  if (nrow(dat) < 200) {
    message("Insufficient data for bifactor analysis (", sex_label, ")")
    return(NULL)
  }

  pc <- tryCatch(
    psych::polychoric(dat, na.rm = TRUE, correct = 0)$rho,
    error = function(e) {
      message("polychoric failed")
      NULL
    }
  )
  if (is.null(pc)) {
    return(NULL)
  }

  # Schmid-Leiman / omega to get omega_h
  sl <- tryCatch(
    psych::omega(
      pc,
      nfactors = 2,
      fm = "minres",
      plot = FALSE,
      sl = TRUE,
      key = NULL
    ),
    error = function(e) {
      message("omega/SL failed: ", e$message)
      NULL
    }
  )
  if (!is.null(sl)) {
    cat("Omega_h (general factor):", round(sl$omega_h, 3), "\n")
    cat("Omega_t (total)         :", round(sl$omega.tot, 3), "\n")
    cat("ECV (explained common variance by general):", round(sl$ECV, 3), "\n")
  }

  # Derive item lists from the joined data (automatically includes pete if present)
  cr_items <- names(dat)
  cr_items_p <- grep("_p$", cr_items, value = TRUE)
  cr_items_y <- grep("_y$", cr_items, value = TRUE)

  # Model 1: no puberty trait, freely correlated method factors
  # Items cluster only by reporter; no shared latent construct.
  mod_m1 <- paste0(
    "method_p =~ ",
    paste(cr_items_p, collapse = " + "),
    "\n",
    "method_y =~ ",
    paste(cr_items_y, collapse = " + ")
  )
  fit_m1 <- tryCatch(
    lavaan::cfa(
      mod_m1,
      data = dat,
      ordered = cr_items,
      estimator = "WLSMV",
      missing = "pairwise"
    ),
    error = function(e) {
      message("Model 1 (no trait) failed: ", e$message)
      NULL
    }
  )

  # Model 2: puberty trait + freely correlated method factors
  # General factor ⊥ each method factor (bifactor constraint); method factors
  # are free to correlate with each other — parses shared vs. reporter-unique variance.
  mod_m2 <- paste0(
    "general  =~ ",
    paste(cr_items, collapse = " + "),
    "\n",
    "method_p =~ ",
    paste(cr_items_p, collapse = " + "),
    "\n",
    "method_y =~ ",
    paste(cr_items_y, collapse = " + "),
    "\n",
    "general ~~ 0*method_p\n",
    "general ~~ 0*method_y"
  )
  fit_m2 <- tryCatch(
    lavaan::cfa(
      mod_m2,
      data = dat,
      ordered = cr_items,
      estimator = "WLSMV",
      missing = "pairwise"
    ),
    error = function(e) {
      message("Model 2 (trait + free methods) failed: ", e$message)
      NULL
    }
  )

  # Model 3: puberty trait, perfectly correlated methods (single factor)
  # Method factors contribute no unique variance beyond the trait.
  mod_m3 <- paste("pub =~", paste(cr_items, collapse = " + "))
  fit_m3 <- tryCatch(
    lavaan::cfa(
      mod_m3,
      data = dat,
      ordered = cr_items,
      estimator = "WLSMV",
      missing = "pairwise"
    ),
    error = function(e) {
      message("Model 3 (trait only) failed: ", e$message)
      NULL
    }
  )

  # fit comparison table
  fi_names <- c(
    "cfi.scaled",
    "tli.scaled",
    "rmsea.scaled",
    "srmr",
    "chisq.scaled",
    "df.scaled"
  )
  fit_table <- lapply(
    list(
      `m1_no_trait` = fit_m1,
      `m2_trait_free_meth` = fit_m2,
      `m3_trait_only` = fit_m3
    ),
    function(fit) fit_to_row(fit, fi_names)
  ) %>%
    bind_rows(.id = "model") %>%
    mutate(sex = sex_label)

  cat("\n--- Model fit comparison ---\n")
  print(fit_table)

  list(
    fit_table = fit_table,
    omega = sl,
    fit_m1 = fit_m1,
    fit_m2 = fit_m2,
    fit_m3 = fit_m3
  )
}

bf_female <- run_bifactor_section(
  female_parent,
  female_youth,
  "female",
  include_pete = TRUE
)
bf_male <- run_bifactor_section(
  male_parent,
  male_youth,
  "male",
  include_pete = TRUE
)

bifactor_fits <- bind_rows(
  bf_female$fit_table,
  bf_male$fit_table
)
write.csv(
  bifactor_fits,
  file.path(out_dir, "bifactor_model_comparison.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# SECTION 2b: CROSS-REPORTER LONGITUDINAL INVARIANCE
# Answers: does the multimethod model hold at each wave, and are the loadings
# stable enough across waves to support a consistent longitudinal model?
# Per-wave fits: M1/M2/M3 at each wave separately.
# Multigroup invariance: M1 (reliable) and M2 (attempted) with wave as group.
# ---------------------------------------------------------------------------

run_cr_longitudinal_invariance <- function(
  parent_df,
  youth_df,
  sex_label,
  include_pete = FALSE
) {
  cat("\n\n=== Cross-reporter longitudinal invariance |", sex_label, "===\n")

  base_items <- c("peta", "petb", "petc", "petd")
  all_items <- if (include_pete) c(base_items, "fpete", "mpete") else base_items

  parent_sel <- parent_df %>%
    select(id, wave, any_of(all_items)) %>%
    rename_with(~ paste0(., "_p"), any_of(all_items))

  youth_sel <- youth_df %>%
    select(id, wave, any_of(all_items)) %>%
    rename_with(~ paste0(., "_y"), any_of(all_items))

  dat_long <- inner_join(parent_sel, youth_sel, by = c("id", "wave")) %>%
    mutate(
      across(-c(id, wave), as.integer),
      wave = factor(wave, levels = wave_order)
    ) %>%
    filter(if_all(-c(id, wave), ~ !is.na(.) & . %in% 1:4))

  cr_items <- setdiff(names(dat_long), c("id", "wave"))
  cr_items_p <- grep("_p$", cr_items, value = TRUE)
  cr_items_y <- grep("_y$", cr_items, value = TRUE)

  mod_m1 <- paste0(
    "method_p =~ ",
    paste(cr_items_p, collapse = " + "),
    "\n",
    "method_y =~ ",
    paste(cr_items_y, collapse = " + ")
  )
  mod_m2 <- paste0(
    "general  =~ ",
    paste(cr_items, collapse = " + "),
    "\n",
    "method_p =~ ",
    paste(cr_items_p, collapse = " + "),
    "\n",
    "method_y =~ ",
    paste(cr_items_y, collapse = " + "),
    "\n",
    "general ~~ 0*method_p\n",
    "general ~~ 0*method_y"
  )
  mod_m3 <- paste("pub =~", paste(cr_items, collapse = " + "))

  fi_names <- c(
    "cfi.scaled",
    "tli.scaled",
    "rmsea.scaled",
    "srmr",
    "chisq.scaled",
    "df.scaled"
  )

  # --- Per-wave fits --------------------------------------------------------
  cat("\n--- Per-wave fits (M1 / M2 / M3) ---\n")
  wave_fits <- lapply(wave_order, function(wv) {
    sub <- dat_long %>% filter(wave == wv) %>% select(all_of(cr_items))
    if (nrow(sub) < 100) {
      return(NULL)
    }
    lapply(
      list(
        m1_no_trait = mod_m1,
        m2_trait_free_meth = mod_m2,
        m3_trait_only = mod_m3
      ),
      function(mod) {
        fit <- tryCatch(
          lavaan::cfa(
            mod,
            data = sub,
            ordered = cr_items,
            estimator = "WLSMV",
            missing = "pairwise"
          ),
          error = function(e) NULL
        )
        row <- fit_to_row(fit, fi_names)
        if (!is.null(row)) mutate(row, wave = wv) else NULL
      }
    ) %>%
      bind_rows(.id = "model") %>%
      mutate(sex = sex_label)
  }) %>%
    bind_rows()

  print_all(wave_fits)

  # --- Multigroup invariance across waves -----------------------------------
  dat_grp <- dat_long %>%
    select(-id) %>%
    filter(!is.na(wave)) %>%
    mutate(across(all_of(cr_items), ~ as.ordered(.)))

  fi_inv <- c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr")
  inv_levels <- list(
    configural = character(0),
    metric = "loadings",
    scalar = c("loadings", "thresholds")
  )

  run_mg_inv <- function(mod, label) {
    fits <- lapply(inv_levels, function(ge) {
      tryCatch(
        lavaan::cfa(
          mod,
          data = dat_grp,
          group = "wave",
          ordered = cr_items,
          estimator = "WLSMV",
          parameterization = "theta",
          group.equal = ge
        ),
        error = function(e) {
          message(
            sex_label,
            " ",
            label,
            " [",
            paste(ge, collapse = "+"),
            "]: ",
            e$message
          )
          NULL
        }
      )
    })
    tbl <- lapply(fits, function(fit) fit_to_row(fit, fi_inv)) %>%
      bind_rows(.id = "model") %>%
      mutate(sex = sex_label, multimethod_model = label)
    cfg_cfi <- tbl$cfi.scaled[tbl$model == "configural"]
    tbl$delta_cfi <- if (length(cfg_cfi) == 1 && !is.na(cfg_cfi)) {
      tbl$cfi.scaled - cfg_cfi
    } else {
      NA_real_
    }
    tbl
  }

  cat("\n--- Multigroup invariance across waves: M1 (no-trait baseline) ---\n")
  inv_m1 <- run_mg_inv(mod_m1, "M1_no_trait")
  print(inv_m1)

  cat(
    "\n--- Multigroup invariance across waves: M2 (trait + free methods) ---\n"
  )
  inv_m2 <- run_mg_inv(mod_m2, "M2_trait_free_meth")
  print(inv_m2)

  list(wave_fits = wave_fits, inv_m1 = inv_m1, inv_m2 = inv_m2)
}

# Run with pete items — shows full picture including per-wave failures at
# baseline due to fpete floor effects; useful for diagnosing item behaviour.
cr_long_female <- run_cr_longitudinal_invariance(
  female_parent,
  female_youth,
  "female",
  include_pete = TRUE
)
cr_long_male <- run_cr_longitudinal_invariance(
  male_parent,
  male_youth,
  "male",
  include_pete = TRUE
)

cr_wave_fits <- bind_rows(cr_long_female$wave_fits, cr_long_male$wave_fits)
write.csv(
  cr_wave_fits,
  file.path(out_dir, "cr_per_wave_fits_with_pete.csv"),
  row.names = FALSE
)

# Run without pete items — ordinal-only (peta–petd) for stable multigroup
# invariance estimates unaffected by fpete floor effects at early waves.
cr_long_female_ord <- run_cr_longitudinal_invariance(
  female_parent,
  female_youth,
  "female",
  include_pete = FALSE
)
cr_long_male_ord <- run_cr_longitudinal_invariance(
  male_parent,
  male_youth,
  "male",
  include_pete = FALSE
)

cr_wave_fits_ord <- bind_rows(
  cr_long_female_ord$wave_fits,
  cr_long_male_ord$wave_fits
)
write.csv(
  cr_wave_fits_ord,
  file.path(out_dir, "cr_per_wave_fits_ordinal.csv"),
  row.names = FALSE
)

cr_inv_table <- bind_rows(
  cr_long_female_ord$inv_m1,
  cr_long_female_ord$inv_m2,
  cr_long_male_ord$inv_m1,
  cr_long_male_ord$inv_m2
)
write.csv(
  cr_inv_table,
  file.path(out_dir, "cr_longitudinal_invariance.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# SECTION 3: LONGITUDINAL MEASUREMENT INVARIANCE
# Group = wave (7 annual timepoints). Per sex × reporter.
# Configural → metric → scalar → strict (WLSMV, theta parameterization)
# ---------------------------------------------------------------------------

run_longitudinal_invariance <- function(df, label, items = ordinal_items) {
  cat("\n\n=== Longitudinal invariance:", label, "===\n")

  dat <- df %>%
    select(wave, any_of(items)) %>%
    filter(wave %in% wave_order) %>%
    mutate(
      wave = factor(wave, levels = wave_order),
      across(any_of(items), ~ as.ordered(as.integer(.)))
    )
  items_present <- intersect(items, names(dat))
  dat <- dat %>% filter(if_all(all_of(items_present), ~ !is.na(.)))

  if (nrow(dat) < 500) {
    message("Skipping ", label, ": insufficient data")
    return(NULL)
  }

  syntax <- paste("pub =~", paste(items_present, collapse = " + "))

  models <- list(
    configural = list(group.equal = character(0)),
    metric = list(group.equal = "loadings"),
    scalar = list(group.equal = c("loadings", "thresholds")),
    strict = list(group.equal = c("loadings", "thresholds", "residuals"))
  )

  fits <- lapply(models, function(args) {
    tryCatch(
      lavaan::cfa(
        syntax,
        data = dat,
        group = "wave",
        ordered = items_present,
        estimator = "WLSMV",
        parameterization = "theta",
        group.equal = args$group.equal
      ),
      error = function(e) {
        message(label, " | ", e$message)
        NULL
      }
    )
  })

  fi_names <- c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr")
  fit_table <- lapply(fits, function(fit) fit_to_row(fit, fi_names)) %>%
    bind_rows(.id = "model") %>%
    mutate(group = label)

  if (nrow(fit_table) == 0 || !"cfi.scaled" %in% names(fit_table)) {
    return(NULL)
  }

  # delta CFI from configural (threshold: < -.010)
  fit_table <- fit_table %>%
    mutate(delta_cfi = cfi.scaled - cfi.scaled[model == "configural"])

  cat("\n--- Fit indices ---\n")
  print(fit_table)

  # LRT (DIFFTEST for WLSMV)
  non_null <- Filter(Negate(is.null), fits)
  if (length(non_null) >= 2) {
    lrt <- tryCatch(
      do.call(
        lavaan::lavTestLRT,
        c(non_null, list(model.names = names(non_null)))
      ),
      error = function(e) NULL
    )
    if (!is.null(lrt)) {
      cat("\n--- LRT ---\n")
      print(lrt)
    }
  }

  list(fits = fits, fit_table = fit_table)
}

long_inv <- list(
  female_parent = run_longitudinal_invariance(
    female_parent,
    "female_parent",
    items = female_items
  ),
  female_youth = run_longitudinal_invariance(
    female_youth,
    "female_youth",
    items = female_items
  ),
  male_parent = run_longitudinal_invariance(
    male_parent,
    "male_parent",
    items = male_items
  ),
  male_youth = run_longitudinal_invariance(
    male_youth,
    "male_youth",
    items = male_items
  )
)

long_inv_table <- bind_rows(lapply(long_inv, `[[`, "fit_table"))
write.csv(
  long_inv_table,
  file.path(out_dir, "longitudinal_invariance.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# SECTION 4: INVARIANCE ACROSS BACKGROUND COVARIATES
# Groups: age tertile, race/ethnicity, BMI tertile
# Per sex × reporter. Configural vs. scalar only (delta CFI criterion).
# ---------------------------------------------------------------------------

make_tertile <- function(x) {
  cuts <- quantile(x, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
  case_when(
    x <= cuts[1] ~ "low",
    x <= cuts[2] ~ "mid",
    TRUE ~ "high"
  )
}

run_covariate_invariance <- function(
  df,
  covariate,
  group_label,
  ds_label,
  items = ordinal_items
) {
  dat <- df %>%
    filter(!is.na(.data[[covariate]])) %>%
    mutate(
      group_var = .data[[covariate]],
      across(any_of(items), ~ as.ordered(as.integer(.)))
    )
  items_present <- intersect(items, names(dat))
  dat <- dat %>%
    filter(if_all(all_of(items_present), ~ !is.na(.))) %>%
    select(group_var, all_of(items_present))

  if (length(unique(dat$group_var)) < 2) {
    message("Skipping ", ds_label, " | ", group_label, ": < 2 groups")
    return(NULL)
  }

  syntax <- paste("pub =~", paste(items_present, collapse = " + "))

  fit_config <- tryCatch(
    lavaan::cfa(
      syntax,
      data = dat,
      group = "group_var",
      ordered = items_present,
      estimator = "WLSMV",
      parameterization = "theta"
    ),
    error = function(e) NULL
  )
  fit_scalar <- tryCatch(
    lavaan::cfa(
      syntax,
      data = dat,
      group = "group_var",
      ordered = items_present,
      estimator = "WLSMV",
      parameterization = "theta",
      group.equal = c("loadings", "thresholds")
    ),
    error = function(e) NULL
  )

  fi_names <- c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr")
  row_config <- fit_to_row(fit_config, fi_names)
  row_scalar <- fit_to_row(fit_scalar, fi_names)

  result <- bind_rows(
    if (!is.null(row_config)) {
      mutate(row_config, model = "configural")
    } else {
      NULL
    },
    if (!is.null(row_scalar)) mutate(row_scalar, model = "scalar") else NULL
  )

  if (nrow(result) == 0 || !"cfi.scaled" %in% names(result)) {
    return(NULL)
  }

  result %>%
    mutate(
      covariate = group_label,
      dataset = ds_label,
      delta_cfi = cfi.scaled - cfi.scaled[model == "configural"]
    )
}

# add tertile columns before testing
add_group_cols <- function(df) {
  df <- df %>%
    mutate(
      age_tertile = make_tertile(age),
      bmi_tertile = make_tertile(bmi),
      race_group = as.character(race)
    )
  if ("site" %in% names(df)) {
    df <- df %>% mutate(site_group = as.character(site))
  } else {
    message(
      "'site' column missing — re-run 00_data_foundation.R to enable site invariance"
    )
    df$site_group <- NA_character_
  }
  df
}

female_parent <- add_group_cols(female_parent)
female_youth <- add_group_cols(female_youth)
male_parent <- add_group_cols(male_parent)
male_youth <- add_group_cols(male_youth)

cov_inv_results <- list()
for (ds_name in names(datasets)) {
  df_aug <- get(ds_name) # already has group cols from add_group_cols above
  ds_items <- datasets[[ds_name]]$items
  for (cov in c("age_tertile", "bmi_tertile", "race_group", "site_group")) {
    res <- run_covariate_invariance(df_aug, cov, cov, ds_name, items = ds_items)
    if (!is.null(res)) cov_inv_results[[paste(ds_name, cov, sep = "_")]] <- res
  }
}

cov_inv_table <- bind_rows(cov_inv_results)
cat("\n\n=== Invariance across background covariates ===\n")
print_all(cov_inv_table)
write.csv(
  cov_inv_table,
  file.path(out_dir, "covariate_invariance.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# SUMMARY PLOT: item discrimination heatmap across waves (slide-ready)
# ---------------------------------------------------------------------------
if (nrow(disc_table) > 0) {
  disc_plot_df <- disc_table %>%
    mutate(
      wave = factor(wave, levels = wave_order),
      item = factor(
        item,
        levels = rev(c("peta", "petb", "petc", "petd", "mpete"))
      ),
      reporter_label = dplyr::recode(
        reporter,
        parent = "Parent report",
        youth = "Youth report"
      ),
      sex_label = dplyr::recode(sex, female = "Female", male = "Male"),
      panel = paste0(sex_label, "\n", reporter_label)
    ) %>%
    mutate(
      panel = factor(
        panel,
        levels = c(
          "Female\nParent report",
          "Female\nYouth report",
          "Male\nParent report",
          "Male\nYouth report"
        )
      )
    )

  a_max <- ceiling(max(disc_plot_df$a, na.rm = TRUE))

  p_disc <- ggplot(disc_plot_df, aes(x = wave, y = item, fill = a)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.1f", a)), size = 3.8, colour = "grey15") +
    scale_fill_gradient2(
      low = "#fee090",
      mid = "#1a9850",
      high = "#005a32",
      midpoint = 2,
      limits = c(0, a_max),
      name = "a",
      breaks = c(0, 1, 2, 3, 4),
      labels = c("0", "1", "2", "3", "4")
    ) +
    facet_wrap(~panel, nrow = 1) +
    labs(
      title = "GRM item discrimination across waves",
      subtitle = "a: discrimination parameter  │  higher = better item–trait separation",
      x = "Wave",
      y = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.key.width = unit(1.8, "cm"),
      plot.subtitle = element_text(size = 11, colour = "grey45"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 12, face = "italic")
    )

  ggsave(
    file.path(out_dir, "discrimination_heatmap.png"),
    p_disc,
    width = 13,
    height = 5,
    dpi = 180
  )
}

# ---------------------------------------------------------------------------
# SUMMARY PLOT: covariate invariance heatmap (slide-ready)
# ---------------------------------------------------------------------------
if (nrow(cov_inv_table) > 0 && "delta_cfi" %in% names(cov_inv_table)) {
  inv_plot_df <- cov_inv_table %>%
    filter(model == "scalar", !is.na(delta_cfi)) %>%
    mutate(
      fails = delta_cfi < -0.010,
      covariate_label = dplyr::recode(
        covariate,
        age_tertile = "Age tertile",
        bmi_tertile = "BMI tertile",
        race_group = "Race / ethnicity",
        site_group = "Site"
      ),
      dataset_label = dplyr::recode(
        dataset,
        female_parent = "Female\nparent",
        female_youth = "Female\nyouth",
        male_parent = "Male\nparent",
        male_youth = "Male\nyouth"
      ),
      dataset_label = factor(
        dataset_label,
        levels = c(
          "Female\nparent",
          "Female\nyouth",
          "Male\nparent",
          "Male\nyouth"
        )
      )
    )

  p_inv <- ggplot(
    inv_plot_df,
    aes(x = dataset_label, y = covariate_label, fill = delta_cfi)
  ) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(
      aes(
        label = sprintf("%+.3f", delta_cfi),
        colour = fails
      ),
      size = 5,
      fontface = "bold"
    ) +
    scale_fill_gradient2(
      low = "#d73027",
      mid = "#fee090",
      high = "#1a9850",
      midpoint = -0.010,
      limits = c(min(inv_plot_df$delta_cfi, na.rm = TRUE) - 0.005, 0),
      name = "ΔCFI",
      breaks = c(-0.04, -0.02, -0.01, 0),
      labels = c("−.04", "−.02", "−.01", "0")
    ) +
    scale_colour_manual(
      values = c("FALSE" = "grey25", "TRUE" = "#a50026"),
      guide = "none"
    ) +
    labs(
      title = "Measurement invariance across background covariates",
      subtitle = "ΔCFI (scalar − configural)  │  red bold = |CFI| > .010 (non-invariant)",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(1.8, "cm"),
      plot.subtitle = element_text(size = 11, colour = "grey45"),
      axis.text = element_text(size = 13)
    )

  ggsave(
    file.path(out_dir, "covariate_invariance_heatmap.png"),
    p_inv,
    width = 9,
    height = 5,
    dpi = 180
  )
  cat("\nCovariate invariance heatmap written.\n")
}

# ---------------------------------------------------------------------------
# SUMMARY PLOT: longitudinal invariance heatmap (slide-ready)
# ---------------------------------------------------------------------------
if (length(long_inv) > 0) {
  long_inv_plot_df <- long_inv_table %>%
    filter(model != "configural", !is.na(delta_cfi)) %>%
    mutate(
      model = factor(model, levels = c("metric", "scalar", "strict")),
      group_label = dplyr::recode(
        group,
        female_parent = "Female\nparent",
        female_youth = "Female\nyouth",
        male_parent = "Male\nparent",
        male_youth = "Male\nyouth"
      ),
      group_label = factor(
        group_label,
        levels = c(
          "Female\nparent",
          "Female\nyouth",
          "Male\nparent",
          "Male\nyouth"
        )
      ),
      fails = delta_cfi < -0.010
    )

  p_long_inv <- ggplot(
    long_inv_plot_df,
    aes(x = group_label, y = model, fill = delta_cfi)
  ) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(
      aes(
        label = sprintf("%+.3f", delta_cfi),
        colour = fails
      ),
      size = 5,
      fontface = "bold"
    ) +
    scale_fill_gradient2(
      low = "#d73027",
      mid = "#fee090",
      high = "#1a9850",
      midpoint = -0.010,
      limits = c(min(long_inv_plot_df$delta_cfi, na.rm = TRUE) - 0.01, 0),
      name = "ΔCFI vs. configural",
      breaks = c(-0.15, -0.10, -0.05, -0.01, 0),
      labels = c("−.15", "−.10", "−.05", "−.01", "0")
    ) +
    scale_colour_manual(
      values = c("FALSE" = "grey25", "TRUE" = "#a50026"),
      guide = "none"
    ) +
    scale_y_discrete(limits = rev) +
    labs(
      title = "Longitudinal measurement invariance across 7 waves",
      subtitle = "ΔCFI (vs. configural baseline)  │  red bold = |ΔCFI| > .010 (non-invariant)",
      x = NULL,
      y = "Invariance level"
    ) +
    theme_minimal(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(1.8, "cm"),
      plot.subtitle = element_text(size = 11, colour = "grey45"),
      axis.text = element_text(size = 13)
    )

  ggsave(
    file.path(out_dir, "longitudinal_invariance_heatmap.png"),
    p_long_inv,
    width = 9,
    height = 5,
    dpi = 180
  )
  cat("\nLongitudinal invariance heatmap written.\n")
}

cat("\n\nAll psychometrics output written to:", out_dir, "\n")
