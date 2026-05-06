pacman::p_load(
  "rstan",
  "posterior",
  "dplyr",
  "tidyr",
  "purrr",
  "readr",
  "lavaan",
  "semTools",
  "ggplot2",
  "NBDCtools"
)

set.seed(90025)

# quick-run toggle
test_mode <- TRUE
test_n_ids <- 250

root_path <- Sys.getenv("HOME_DIR")
if (!nzchar(root_path)) {
  root_path <- Sys.getenv("HOME")
}
proj_path <- here::here()
data_root <- file.path(
  root_path,
  "projects/abcd-projs/abcd-data-release-6.0/nbdc-tools-data"
)
if (!dir.exists(data_root)) {
  alt_data_root <- file.path(
    root_path,
    "Library/CloudStorage/Box-Box/everything/projects/abcd-projs/abcd-data-release-6.0/nbdc-tools-data"
  )
  if (dir.exists(alt_data_root)) {
    data_root <- alt_data_root
  }
}
if (!dir.exists(data_root)) {
  stop(
    sprintf("Could not find nbdc-tools-data directory at '%s'.", data_root),
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 2. DATA LOADING AND PREPARATION
# -----------------------------------------------------------------------------

# use nbdctools to load in data (caregiver info = )

vars <- c(
  "ab_g_dyn__visit_age",
  "ab_g_stc__cohort_sex",
  "ph_y_anthr__waist_001",
  "ph_y_anthr__height_mean",
  "ab_g_dyn__visit__day1_inform",
  "ab_g_stc__cohort_ethnrace__mhisp",
  "ph_p_pds_001",
  "ph_p_pds_002",
  "ph_p_pds_003",
  "ph_p_pds__f_001",
  "ph_p_pds__f_002",
  "ph_p_pds__m_001",
  "ph_p_pds__m_002",
  "ph_y_pds_001",
  "ph_y_pds_002",
  "ph_y_pds_003",
  "ph_y_pds__f_001",
  "ph_y_pds__f_002",
  "ph_y_pds__m_001",
  "ph_y_pds__m_002"
)


data <- create_dataset(
  dir_data = data_root,
  study = "abcd",
  vars = vars,
  value_to_na = TRUE,
  bind_shadow = TRUE
)

prepare_observed_puberty_lmnlfa_long <- function(data) {
  sex_is_male <- function(x) {
    x_chr <- tolower(trimws(as.character(x)))
    x_chr %in% c("1", "m", "male", "boy")
  }
  first_non_missing <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA_real_)
    }
    x[1]
  }

  base <- data %>%
    transmute(
      id = participant_id,
      wave = session_id,
      age = as.numeric(ab_g_dyn__visit_age),
      sex_male = sex_is_male(ab_g_stc__cohort_sex),
      race = as.numeric(ab_g_stc__cohort_ethnrace__mhisp),
      bmi = as.numeric(ph_y_anthr__waist_001) /
        as.numeric(ph_y_anthr__height_mean),
      ph_p_pds_001 = as.numeric(ph_p_pds_001),
      ph_p_pds_002 = as.numeric(ph_p_pds_002),
      ph_p_pds_003 = as.numeric(ph_p_pds_003),
      ph_p_pds__f_001 = as.numeric(ph_p_pds__f_001),
      ph_p_pds__f_002 = as.numeric(ph_p_pds__f_002),
      ph_p_pds__m_001 = as.numeric(ph_p_pds__m_001),
      ph_p_pds__m_002 = as.numeric(ph_p_pds__m_002),
      ph_y_pds_001 = as.numeric(ph_y_pds_001),
      ph_y_pds_002 = as.numeric(ph_y_pds_002),
      ph_y_pds_003 = as.numeric(ph_y_pds_003),
      ph_y_pds__f_001 = as.numeric(ph_y_pds__f_001),
      ph_y_pds__f_002 = as.numeric(ph_y_pds__f_002),
      ph_y_pds__m_001 = as.numeric(ph_y_pds__m_001),
      ph_y_pds__m_002 = as.numeric(ph_y_pds__m_002)
    ) %>%
    filter(!is.na(id), !is.na(age))

  occ <- base %>%
    mutate(
      p4 = if_else(sex_male, ph_p_pds__m_001, ph_p_pds__f_001),
      p5 = if_else(sex_male, ph_p_pds__m_002, ph_p_pds__f_002),
      y4 = if_else(sex_male, ph_y_pds__m_001, ph_y_pds__f_001),
      y5 = if_else(sex_male, ph_y_pds__m_002, ph_y_pds__f_002)
    ) %>%
    select(
      id,
      wave,
      age,
      race,
      bmi,
      ph_p_pds_001,
      ph_p_pds_002,
      ph_p_pds_003,
      p4,
      p5,
      ph_y_pds_001,
      ph_y_pds_002,
      ph_y_pds_003,
      y4,
      y5
    ) %>%
    distinct(id, wave, .keep_all = TRUE) %>%
    group_by(id) %>%
    arrange(age, .by_group = TRUE) %>%
    mutate(
      time = dplyr::dense_rank(age),
      base_age = first(age),
      race_id = first_non_missing(race),
      bmi_id = first_non_missing(bmi),
      age_c = age - 12,
      age2_c = age_c^2
    ) %>%
    ungroup()

  occ %>%
    pivot_longer(
      cols = c(
        ph_p_pds_001,
        ph_p_pds_002,
        ph_p_pds_003,
        p4,
        p5,
        ph_y_pds_001,
        ph_y_pds_002,
        ph_y_pds_003,
        y4,
        y5
      ),
      names_to = "item_name",
      values_to = "resp_raw"
    ) %>%
    mutate(
      item = recode(
        item_name,
        ph_p_pds_001 = 1L,
        ph_p_pds_002 = 2L,
        ph_p_pds_003 = 3L,
        p4 = 4L,
        p5 = 5L,
        ph_y_pds_001 = 6L,
        ph_y_pds_002 = 7L,
        ph_y_pds_003 = 8L,
        y4 = 9L,
        y5 = 10L
      ),
      resp = case_when(
        item %in%
          c(1L, 2L, 3L, 4L, 6L, 7L, 8L, 9L) &
          resp_raw %in% 1:4 ~ as.integer(resp_raw),
        item %in% c(5L, 10L) & resp_raw %in% c(0, 1) ~ as.integer(resp_raw),
        item %in% c(5L, 10L) & !is.na(resp_raw) ~ as.integer(resp_raw > 1),
        TRUE ~ NA_integer_
      )
    ) %>%
    filter(!is.na(resp)) %>%
    select(
      id, time, age, age_c, age2_c, item, resp,
      race = race_id, bmi = bmi_id, base_age
    )
}


mc_cores <- parallel::detectCores()
if (is.na(mc_cores)) {
  mc_cores <- 1L
}
options(mc.cores = mc_cores)
rstan_options(auto_write = TRUE)

# --- helpers ---
logistic_15 <- function(age, floor = 1, ceiling = 5, alpha, t0) {
  floor + (ceiling - floor) / (1 + exp(-alpha * (age - t0)))
}

get_script_dir <- function() {
  normalize_local_path <- function(path) {
    if (is.null(path) || !nzchar(path)) {
      return(NULL)
    }

    if (grepl("^file://", path, ignore.case = TRUE)) {
      path <- utils::URLdecode(sub("^file://", "", path, ignore.case = TRUE))
    }

    path <- path.expand(path)

    if (!file.exists(path)) {
      return(NULL)
    }

    normalizePath(path, mustWork = TRUE)
  }

  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    file_path <- normalize_local_path(sub("^--file=", "", file_arg[1]))
    if (!is.null(file_path)) {
      return(dirname(file_path))
    }
  }

  for (frame_idx in rev(seq_len(sys.nframe()))) {
    ofile <- sys.frame(frame_idx)$ofile
    file_path <- normalize_local_path(ofile)
    if (!is.null(file_path)) {
      return(dirname(file_path))
    }
  }

  NULL
}

default_script_dir <- get_script_dir()

resolve_input_path <- function(path) {
  if (file.exists(path)) {
    return(normalizePath(path, mustWork = TRUE))
  }

  script_dir <- default_script_dir
  if (is.null(script_dir)) {
    script_dir <- get_script_dir()
  }
  if (!is.null(script_dir)) {
    script_path <- file.path(script_dir, path)
    if (file.exists(script_path)) {
      return(normalizePath(script_path, mustWork = TRUE))
    }
  }

  stop(
    sprintf(
      paste(
        "Could not find input file '%s'.",
        "Checked relative to getwd() = '%s'%s."
      ),
      path,
      normalizePath(getwd(), mustWork = TRUE),
      if (is.null(script_dir)) {
        ""
      } else {
        sprintf(" and script_dir = '%s'", script_dir)
      }
    ),
    call. = FALSE
  )
}

simulate_puberty_lmnfa_long <- function(
  n,
  max_waves = 7,
  age_min = 9.5,
  age_max = 15.5,
  floor = 1,
  ceiling = 5,
  t0_range = c(11, 13),
  alpha_range = c(0.7, 0.95),
  p_race = 0.30,
  bmi_mu = 25,
  bmi_sd = 5,
  base_age_mu = 11.5,
  base_age_sd = 1.0,
  p_items = 10,
  ordinal_items = c(1:4, 6:9),
  binary_items = c(5, 10),
  loading_mu = 1.0,
  loading_sd = 0.15,
  resid_sd = 1.0,
  mean_waves_obs = 4
) {
  stopifnot(p_items == 10)

  id <- 1:n

  mods <- tibble(
    id = id,
    race = rbinom(n, 1, p_race),
    bmi = rnorm(n, bmi_mu, bmi_sd),
    base_age = rnorm(n, base_age_mu, base_age_sd)
  )

  pars <- tibble(
    id = id,
    t0 = runif(n, t0_range[1], t0_range[2]),
    alpha = runif(n, alpha_range[1], alpha_range[2])
  )

  w_obs <- pmin(max_waves, pmax(1, rpois(n, lambda = mean_waves_obs)))

  occ <- map2_dfr(id, w_obs, function(i, w) {
    tibble(
      id = i,
      time = 1:w,
      age = sort(runif(w, age_min, age_max))
    )
  }) %>%
    left_join(mods, by = "id") %>%
    left_join(pars, by = "id") %>%
    mutate(
      eta = logistic_15(
        age,
        floor = floor,
        ceiling = ceiling,
        alpha = alpha,
        t0 = t0
      ),
      age_c = age - 12,
      age2_c = age_c^2
    )

  lambdas <- rnorm(p_items, loading_mu, loading_sd)
  thr_14 <- c(-0.8, 0.0, 0.8)

  occ %>%
    crossing(item = 1:p_items) %>%
    mutate(
      lambda = lambdas[item],
      y_star = lambda * eta + rnorm(n(), 0, resid_sd),
      resp = case_when(
        item %in% ordinal_items ~ as.integer(cut(
          y_star,
          breaks = c(-Inf, thr_14, Inf),
          labels = FALSE
        )),
        item %in% binary_items ~ as.integer(y_star > 0),
        TRUE ~ NA_integer_
      )
    ) %>%
    select(
      id,
      time,
      age,
      age_c,
      age2_c,
      item,
      resp,
      race,
      bmi,
      base_age,
      t0,
      alpha
    )
}

make_ldf_simple <- function(
  P = 10,
  dif_time_items = c(1, 2, 3),
  dif_inv_items = c(6, 7, 8)
) {
  Ldf <- matrix(0, nrow = P, ncol = 2)
  Ldf[dif_time_items, 1] <- 1
  Ldf[dif_inv_items, 2] <- 1
  Ldf
}

build_fa_data <- function(dat_long, Ldf) {
  first_non_missing <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA_real_)
    }
    x[1]
  }

  dat_long <- dat_long %>%
    mutate(
      person = as.integer(factor(id)),
      itm = as.integer(factor(item, levels = sort(unique(item)))),
      time = as.integer(factor(time, levels = sort(unique(time))))
    )

  person_tbl <- dat_long %>%
    group_by(id, person) %>%
    summarise(
      race = first_non_missing(race),
      bmi = first_non_missing(bmi),
      base_age = first_non_missing(base_age),
      .groups = "drop"
    ) %>%
    filter(!is.na(race), !is.na(bmi), !is.na(base_age)) %>%
    arrange(person) %>%
    mutate(
      race_c = race - mean(race, na.rm = TRUE),
      bmi_z = as.numeric(scale(bmi)),
      base_age_z = as.numeric(scale(base_age))
    )

  dat_long <- dat_long %>%
    semi_join(person_tbl %>% select(id, person), by = c("id", "person")) %>%
    select(-race, -bmi, -base_age) %>%
    left_join(
      person_tbl %>%
        select(id, person, race, bmi, base_age, race_c, bmi_z, base_age_z),
      by = c("id", "person")
    )

  # Reindex after filtering so Stan indices are contiguous and bounded.
  dat_long <- dat_long %>%
    mutate(
      person = as.integer(factor(id)),
      itm = as.integer(factor(itm, levels = sort(unique(itm)))),
      time = as.integer(factor(time, levels = sort(unique(time))))
    )

  nobs <- nrow(dat_long)
  p <- length(unique(dat_long$itm))
  ni <- length(unique(dat_long$person))
  d <- length(unique(dat_long$time))

  y <- as.integer(dat_long$resp)

  is_binary <- rep(0L, p)
  is_binary[c(5, 10)] <- 1L

  k_item <- rep(4L, p)
  k_item[c(5, 10)] <- 2L
  k_max <- max(k_item)

  xf_person <- as.matrix(person_tbl[, c("race_c", "bmi_z", "base_age_z")])
  xf <- as.matrix(dat_long[, c("race_c", "bmi_z", "base_age_z")])
  xtv <- as.matrix(dat_long[, c("age_c", "age2_c")])

  nfpreds <- ncol(xf_person)
  ntvpreds <- ncol(xtv)

  # IMPORTANT: for your Stan file, mtv/mf are item counts, not multiplied by predictor count
  mtv <- sum(Ldf[, 1] == 1)
  mf <- sum(Ldf[, 2] == 1)

  list(
    nobs = nobs,
    p = p,
    ni = ni,
    d = d,
    person = dat_long$person,
    itm = dat_long$itm,
    time = dat_long$time,
    age_c = dat_long$age_c,
    age2_c = dat_long$age2_c,
    y = y,
    is_binary = is_binary,
    k_item = k_item,
    k_max = k_max,
    nfpreds = nfpreds,
    ntvpreds = ntvpreds,
    xf_person = xf_person,
    xf = xf,
    xtv = xtv,
    ldf = Ldf,
    mtv = mtv,
    mf = mf,
    sigma_l = 2,
    sigma_nu = 3,
    sigma_cor = 2,
    sigma_f = 1,
    sigma_di = 2 # ,
    # total_var = 0.76
  )
}

compile_model <- function(stan_file) {
  stan_path <- resolve_input_path(stan_file)
  invisible(stan_model(file = stan_path))
}

fit_once <- function(
  stan_m,
  fa_data,
  chains = 2,
  iter = 1000,
  warmup = 500,
  seed = 1,
  refresh = 100,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
) {
  sampling(
    object = stan_m,
    data = fa_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    seed = seed,
    refresh = refresh,
    control = control
  )
}

get_draws_df <- function(fit) {
  posterior::as_draws_df(fit)
}

posterior_prob_gt0 <- function(draws, param) mean(draws[[param]] > 0)

cri_excludes_0 <- function(draws, param, level = 0.95) {
  q <- quantile(
    draws[[param]],
    probs = c((1 - level) / 2, 1 - (1 - level) / 2),
    na.rm = TRUE
  )
  (q[1] > 0) || (q[2] < 0)
}

cri_width <- function(draws, param, level = 0.95) {
  q <- quantile(
    draws[[param]],
    probs = c((1 - level) / 2, 1 - (1 - level) / 2),
    na.rm = TRUE
  )
  unname(q[2] - q[1])
}

run_one_rep <- function(
  n,
  stan_m,
  Ldf,
  seed = 1,
  fit_iter = 800,
  fit_warmup = 400,
  success_param = "mu_quad",
  success_type = c("prob_gt0", "cri_excl0", "precision"),
  prob_cut = 0.95,
  width_cut = 0.15,
  refresh = 100
) {
  success_type <- match.arg(success_type)

  dat <- simulate_puberty_lmnfa_long(n = n)

  # range checks
  stopifnot(all(dat$resp[dat$item %in% c(1:4, 6:9)] %in% 1:4))
  stopifnot(all(dat$resp[dat$item %in% c(5, 10)] %in% 0:1))

  fa_data <- build_fa_data(dat, Ldf = Ldf)

  fit <- try(
    fit_once(
      stan_m,
      fa_data,
      chains = 1,
      iter = fit_iter,
      warmup = fit_warmup,
      seed = seed,
      refresh = refresh
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    return(list(ok = FALSE, reason = "fit_error"))
  }

  draws <- get_draws_df(fit)
  if (!success_param %in% names(draws)) {
    return(list(ok = FALSE, reason = "param_not_found"))
  }

  ok <- switch(
    success_type,
    prob_gt0 = posterior_prob_gt0(draws, success_param) > prob_cut,
    cri_excl0 = cri_excludes_0(draws, success_param),
    precision = cri_width(draws, success_param) < width_cut
  )

  list(ok = ok, reason = "ok")
}

estimate_power_curve <- function(
  n_grid,
  reps = 100,
  stan_file = "stan/lmnlfa-quad.stan",
  dif_time_items = c(1, 2, 3),
  dif_inv_items = c(6, 7, 8),
  success_param = "mu_quad",
  success_type = "cri_excl0",
  seed = 1,
  fit_iter = 800,
  fit_warmup = 400,
  refresh = 100,
  progress = interactive()
) {
  set.seed(seed)

  Ldf <- make_ldf_simple(
    P = 10,
    dif_time_items = dif_time_items,
    dif_inv_items = dif_inv_items
  )
  stan_m <- compile_model(stan_file)

  map_dfr(n_grid, function(n) {
    if (progress) {
      message(sprintf(
        "[%s] Starting n = %s with %s replications",
        Sys.time(),
        n,
        reps
      ))
    }

    n_start <- Sys.time()
    out <- vector("list", reps)

    for (rep_idx in seq_len(reps)) {
      rep_start <- Sys.time()
      out[[rep_idx]] <- run_one_rep(
        n = n,
        stan_m = stan_m,
        Ldf = Ldf,
        seed = sample.int(1e9, 1),
        fit_iter = fit_iter,
        fit_warmup = fit_warmup,
        success_param = success_param,
        success_type = success_type,
        refresh = refresh
      )

      if (progress) {
        elapsed_rep <- round(
          as.numeric(difftime(Sys.time(), rep_start, units = "mins")),
          2
        )
        message(
          sprintf(
            "[%s] Finished n = %s replicate %s/%s in %s min (%s)",
            Sys.time(),
            n,
            rep_idx,
            reps,
            elapsed_rep,
            out[[rep_idx]]$reason
          )
        )
      }
    }

    if (progress) {
      elapsed_n <- round(
        as.numeric(difftime(Sys.time(), n_start, units = "mins")),
        2
      )
      message(sprintf(
        "[%s] Completed n = %s in %s min",
        Sys.time(),
        n,
        elapsed_n
      ))
    }

    tibble(
      n = n,
      reps = reps,
      conv_rate = mean(map_lgl(out, ~ .x$reason == "ok")),
      power = mean(map_lgl(out, ~ isTRUE(.x$ok)), na.rm = TRUE)
    )
  })
}

if (sys.nframe() == 0) {
  # --- observed data fit ---
  Ldf <- make_ldf_simple(P = 10)
  dat_long <- prepare_observed_puberty_lmnlfa_long(data)
  if (isTRUE(test_mode)) {
    set.seed(90025)
    keep_ids <- sample(unique(dat_long$id), min(test_n_ids, dplyr::n_distinct(dat_long$id)))
    dat_long <- dplyr::filter(dat_long, id %in% keep_ids)
    message(sprintf("[%s] Test mode ON: fitting %s IDs", Sys.time(), dplyr::n_distinct(dat_long$id)))
  }
  fa_data <- build_fa_data(dat_long, Ldf)
  message(sprintf("[%s] Compiling Stan model", Sys.time()))
  stan_m <- compile_model("stan/lmnlfa-quad.stan")
  message(sprintf("[%s] Finished compiling Stan model", Sys.time()))
  message(sprintf("[%s] Starting observed-data fit", Sys.time()))
  fit <- fit_once(
    stan_m,
    fa_data,
    chains = if (isTRUE(test_mode)) 1 else 2,
    iter = if (isTRUE(test_mode)) 200 else 1000,
    warmup = if (isTRUE(test_mode)) 100 else 500,
    seed = 123,
    refresh = if (isTRUE(test_mode)) 20 else 50
  )
  message(sprintf("[%s] Finished observed-data fit", Sys.time()))
  print(fit, pars = c("mu_lin", "mu_quad", "sigma_eta"))
}
