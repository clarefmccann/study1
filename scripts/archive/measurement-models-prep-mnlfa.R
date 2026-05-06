############################################################
# pubertal MNLFA in OpenMx using your cleaned data structure
############################################################

library(readr)
library(dplyr)
library(tidyr)
library(OpenMx)

set.seed(90025)

# Define file paths (please update if necessary)
pub_data_root <- file.path(
  "/Users/clarefmccann/University of Oregon Dropbox",
  "Clare McCann/mine/projects/abcd-projs/abcd-data-release-6.0/physical-health/puberty"
)

# ------------------------------------------------------------------
# 1. your loader (unchanged, just keep it)
# ------------------------------------------------------------------
load_and_prep <- function(file_path, sex_code, reporter) {
  
  if (reporter == "parent") {
    name_map <- c(
      pete_p = "fpete",
      peta_p = "peta",
      petb_p = "petb",
      petc_p = "petc",
      petd_p = "petd",
      pdss_p = "PDSS"
    )
    if (sex_code == 0) name_map["pete_p"] <- "mpete"  # males use mpete
  } else {
    name_map <- c(
      pete_y = "fpete",
      peta_y = "peta",
      petb_y = "petb",
      petc_y = "petc",
      petd_y = "petd",
      pdss_y = "PDSS"
    )
    if (sex_code == 0) name_map["pete_y"] <- "mpete"
  }
  
  read_csv(file_path, show_col_types = FALSE) %>%
    rename(any_of(name_map)) %>%
    select(id, wave, age, starts_with("pet"), starts_with("pdss")) %>%
    mutate(sex = sex_code,
           reporter = ifelse(reporter == "parent", 0, 1))  # 0 = parent, 1 = youth
}

# ------------------------------------------------------------------
# 2. load all four datasets in long form
# ------------------------------------------------------------------
pub_f_p <- load_and_prep(file.path(pub_data_root, "filtered_parent_tannerstages_f.csv"), 1, "parent")
pub_m_p <- load_and_prep(file.path(pub_data_root, "filtered_parent_tannerstages_m.csv"), 0, "parent")
pub_f_y <- load_and_prep(file.path(pub_data_root, "filtered_youth_tannerstages_f.csv"), 1, "youth")
pub_m_y <- load_and_prep(file.path(pub_data_root, "filtered_youth_tannerstages_m.csv"), 0, "youth")

# bind all long
all_long <- bind_rows(pub_f_p, pub_m_p, pub_f_y, pub_m_y) %>%
  mutate(
    peta = coalesce(peta_p, peta_y),
    petb = coalesce(petb_p, petb_y),
    petc = coalesce(petc_p, petc_y),
    petd = coalesce(petd_p, petd_y),
    pete = coalesce(pete_p, pete_y),
    pdss = coalesce(pdss_p, pdss_y)
  ) 

# if PDSS is not part of the ordinal item set for mnlfa, we can drop it later

# ------------------------------------------------------------------
# 4. recode wave (you had this, keeping it)
# ------------------------------------------------------------------
map_wave <- c("ses-00A" = "bl", "ses-01A" = "fu1", "ses-02A" = "fu2", 
              "ses-03A" = "fu3", "ses-04A" = "fu4", "ses-05A" = "fu5", "ses-06A" = "fu6")

all_long <- all_long %>%
  mutate(wave = recode(wave, !!!map_wave))

# ------------------------------------------------------------------
# 5. z-score age within sex
# ------------------------------------------------------------------
all_long <- all_long %>%
  group_by(sex) %>%
  mutate(age_z = as.numeric(scale(age))) %>%
  ungroup() %>% 
  select(id, wave, age, sex, reporter, age_z,
         peta, petb, petc, petd, pete, pdss)

# ------------------------------------------------------------------
# 6. pick item sets by sex
# females have: peta, petb, petc, petd, pete
# males have:   peta, petb, petc, petd, pete
# (because we coalesced pete)
# ------------------------------------------------------------------
female_dat <- all_long %>%
  filter(sex == 1) %>%
  select(id, wave, reporter, age, age_z,
         peta, petb, petc, petd, pete)

male_dat <- all_long %>%
  filter(sex == 0) %>%
  select(id, wave, reporter, age, age_z,
         peta, petb, petc, petd, pete)

# convert to ordered factors for OpenMx
ordify <- function(df, items) {
  for (v in items) {
    df[[v]] <- as.ordered(df[[v]])
  }
  df
}

items <- c("peta", "petb", "petc", "petd", "pete")

female_dat <- ordify(female_dat, items)
male_dat   <- ordify(male_dat, items)


# females
female_dat <- all_long %>%
  dplyr::filter(sex == 1) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(items),
      ~ OpenMx::mxFactor(.x, levels = sort(unique(.x)), ordered = TRUE)
    )
  )

# males
male_dat <- all_long %>%
  dplyr::filter(sex == 0) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(items),
      ~ OpenMx::mxFactor(.x, levels = sort(unique(.x)), ordered = TRUE)
    )
  )

# ------------------------------------------------------------------
# 7. function to build mnlfa model in OpenMx
#    using definition variables: data.reporter, data.age_z
# ------------------------------------------------------------------

build_mnlfa_model <- function(data,
                              item_names,
                              sex_label = "female",
                              model_name = "mnlfa_config",
                              free_item_mods = TRUE,
                              max_thresh = 3) {
  
  # 1) make sure it's a plain data.frame
  data <- as.data.frame(data)
  
  # 2) make sure reporter/age_z are numeric
  data$reporter <- as.numeric(data$reporter)
  data$age_z    <- as.numeric(data$age_z)
  
  # 3) enforce mxFactor on every ordinal item right here
  for (itm in item_names) {
    if (!(is.factor(data[[itm]]) && isTRUE(attr(data[[itm]], "mxFactor")))) {
      levs <- sort(unique(data[[itm]]))
      data[[itm]] <- mxFactor(data[[itm]], levels = levs, ordered = TRUE)
    }
  }
  
  manifests <- item_names
  n_items   <- length(item_names)
  # count thresholds
  thr_counts <- sapply(item_names, function(v) length(levels(data[[v]])) - 1)
  
  mx_dat <- mxData(observed = data, type = "raw")
  
  # loadings ----------------------------------------------------------
  lambda_base <- mxMatrix("Full", n_items, 1,
                          free = TRUE, values = 1,
                          labels = paste0("lam_", item_names),
                          name = "lambda_base")
  
  lambda_rep <- mxMatrix("Full", n_items, 1,
                         free = free_item_mods, values = 0,
                         labels = if (free_item_mods) paste0("lamrep_", item_names) else NA,
                         name = "lambda_rep")
  
  lambda_age <- mxMatrix("Full", n_items, 1,
                         free = free_item_mods, values = 0,
                         labels = if (free_item_mods) paste0("lamage_", item_names) else NA,
                         name = "lambda_age")
  
  lambda_alg <- mxAlgebra(
    lambda_base +
      lambda_rep * data.reporter +
      lambda_age * data.age_z,
    name = "lambda"
  )
  
  # residuals fixed to 1 ----------------------------------------------
  resid <- mxMatrix("Diag", n_items, n_items,
                    free = FALSE, values = 1,
                    name = "resid")
  
  # latent mean/var moderated -----------------------------------------
  mean_eta_beta <- mxMatrix("Full", 1, 3,
                            free = TRUE, values = c(0,0,0),
                            labels = c("eta_int","eta_rep","eta_age"),
                            name = "mean_eta_beta")
  
  mean_eta_alg <- mxAlgebra(
    mean_eta_beta[1,1] +
      mean_eta_beta[1,2] * data.reporter +
      mean_eta_beta[1,3] * data.age_z,
    name = "Mean_eta"
  )
  
  var_eta_beta <- mxMatrix("Full", 1, 3,
                           free = TRUE, values = c(0,0,0),
                           labels = c("logvar_int","logvar_rep","logvar_age"),
                           name = "var_eta_beta")
  
  var_eta_alg <- mxAlgebra(
    exp(var_eta_beta[1,1] +
          var_eta_beta[1,2] * data.reporter +
          var_eta_beta[1,3] * data.age_z),
    name = "Var_eta"
  )
  
  exp_cov <- mxAlgebra(
    lambda %*% Var_eta %*% t(lambda) + resid,
    name = "expCov"
  )
  
  # must be 1 x p
  exp_mean <- mxAlgebra(
    t(lambda %*% Mean_eta),
    name = "expMean"
  )
  
  # thresholds (padded) ------------------------------------------------
  thr_mats <- list()
  for (i in seq_along(item_names)) {
    itm <- item_names[i]
    k   <- thr_counts[i]
    
    base_thr <- mxMatrix("Full", max_thresh, 1,
                         free   = c(rep(TRUE, k), rep(FALSE, max_thresh - k)),
                         values = c(seq(-1,1,length.out = k), rep(9, max_thresh - k)),
                         labels = c(paste0("tb_", itm, "_", 1:k),
                                    rep(NA, max_thresh - k)),
                         name   = paste0("thr_base_", itm))
    
    rep_thr <- mxMatrix("Full", max_thresh, 1,
                        free   = if (free_item_mods) c(rep(TRUE, k), rep(FALSE, max_thresh - k)) else FALSE,
                        values = 0,
                        labels = if (free_item_mods)
                          c(paste0("tbrep_", itm, "_", 1:k), rep(NA, max_thresh - k))
                        else NA,
                        name   = paste0("thr_rep_", itm))
    
    age_thr <- mxMatrix("Full", max_thresh, 1,
                        free   = if (free_item_mods) c(rep(TRUE, k), rep(FALSE, max_thresh - k)) else FALSE,
                        values = 0,
                        labels = if (free_item_mods)
                          c(paste0("tbage_", itm, "_", 1:k), rep(NA, max_thresh - k))
                        else NA,
                        name   = paste0("thr_age_", itm))
    
    thr_mats <- c(thr_mats, list(base_thr, rep_thr, age_thr))
  }
  
  # build algebra for all thresholds
  thr_rows <- vector("list", max_thresh)
  for (r in 1:max_thresh) {
    row_terms <- character(n_items)
    for (i in seq_along(item_names)) {
      itm <- item_names[i]
      k   <- thr_counts[i]
      if (r <= k) {
        row_terms[i] <- paste0(
          "thr_base_", itm, "[", r, ",1] + ",
          "thr_rep_",  itm, "[", r, ",1] * data.reporter + ",
          "thr_age_",  itm, "[", r, ",1] * data.age_z"
        )
      } else {
        row_terms[i] <- "9"
      }
    }
    thr_rows[[r]] <- paste0("cbind(", paste(row_terms, collapse = ","), ")")
  }
  thr_expr <- paste0("rbind(", paste(thr_rows, collapse = ","), ")")
  
  thr_alg <- mxAlgebraFromString(
    thr_expr,
    name = "thr_all",
    dimnames = list(paste0("thr", 1:max_thresh), item_names)
  )
  
  # expectation & fit -------------------------------------------------
  exp_fun <- mxExpectationNormal(
    covariance = "expCov",
    means      = "expMean",
    thresholds = "thr_all",
    dimnames   = manifests
  )
  
  fit_fun <- mxFitFunctionML()
  
  mxModel(model_name,
          mx_dat,
          lambda_base, lambda_rep, lambda_age, lambda_alg,
          resid,
          mean_eta_beta, mean_eta_alg,
          var_eta_beta, var_eta_alg,
          exp_cov, exp_mean,
          thr_mats,
          thr_alg,
          exp_fun, fit_fun,
          name = paste0(model_name, "_", sex_label))
}

# ------------------------------------------------------------------
# 8. run configural and scalar models
# ------------------------------------------------------------------

# females
mod_f_config <- build_mnlfa_model(
  data = female_dat,
  item_names = items,
  sex_label = "female",
  model_name = "mnlfa_config",
  free_item_mods = TRUE,
  max_thresh = 3
)
fit_f_config <- mxRun(mod_f_config)

mod_f_scalar <- build_mnlfa_model(female_dat,
                                  item_names = female_items,
                                  sex_label = "female",
                                  model_name = "mnlfa_scalar",
                                  free_item_mods = FALSE)
fit_f_scalar <- mxRun(mod_f_scalar)

summary(fit_f_config)
summary(fit_f_scalar)
mxCompare(fit_f_config, fit_f_scalar)

# males
mod_m_config <- build_mnlfa_model(male_dat,
                                  item_names = male_items,
                                  sex_label = "male",
                                  model_name = "mnlfa_configural",
                                  free_item_mods = TRUE)
fit_m_config <- mxRun(mod_m_config)

mod_m_scalar <- build_mnlfa_model(male_dat,
                                  item_names = male_items,
                                  sex_label = "male",
                                  model_name = "mnlfa_scalar",
                                  free_item_mods = FALSE)
fit_m_scalar <- mxRun(mod_m_scalar)

summary(fit_m_config)
summary(fit_m_scalar)
mxCompare(fit_m_config, fit_m_scalar)

# ------------------------------------------------------------------
# 9. partial invariance helper
# ------------------------------------------------------------------
build_partial_invariance <- function(fit_config,
                                     item_names,
                                     keep_free_threshold = character(0),
                                     keep_free_loading = character(0),
                                     sex_label = "female") {
  m <- fit_config
  
  for (itm in item_names) {
    # thresholds
    if (!(itm %in% keep_free_threshold)) {
      # grab number of thresholds from fitted model
      k <- length(levels(fit_config$data[[itm]]))
      # fix reporter and age slopes to 0
      rep_labels <- paste0("tbrep_", itm, "_", 1:(k - 1))
      age_labels <- paste0("tbage_", itm, "_", 1:(k - 1))
      m <- omxSetParameters(m, labels = rep_labels, free = FALSE, values = 0)
      m <- omxSetParameters(m, labels = age_labels, free = FALSE, values = 0)
    }
    
    # loadings
    if (!(itm %in% keep_free_loading)) {
      m <- omxSetParameters(m,
                            labels = c(paste0("lamrep_", itm), paste0("lamage_", itm)),
                            free = FALSE, values = 0)
    }
  }
  
  m@name <- paste0("mnlfa_partial_", sex_label)
  mxRun(m)
}

# example:
# fit_f_partial <- build_partial_invariance(
#   fit_f_config,
#   item_names = female_items,
#   keep_free_threshold = c("pete"),   # suppose pete had DIF
#   keep_free_loading   = c("pete"),
#   sex_label = "female"
# )
# mxCompare(fit_f_config, fit_f_partial)
# mxCompare(fit_f_partial, fit_f_scalar)
