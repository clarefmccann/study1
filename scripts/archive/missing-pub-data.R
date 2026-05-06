# -----------------------------------------------------------------------------
# 1. SETUP
# -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(readr)
library(lavaan)
library(semTools)
library(NBDCtools)
library(NBDCtoolsData)
library(tableone)

set.seed(90025)

root_path <- Sys.getenv("HOME_DIR")
proj_path <- here()
data_root <- file.path(paste0(root_path, "/projects/abcd-projs/abcd-data-release-6.0/nbdc-tools-data/"))

# -----------------------------------------------------------------------------
# 2. DATA LOADING AND PREPARATION
# -----------------------------------------------------------------------------

# use nbdctools to load in data (caregiver info = )

vars <- c("ab_g_dyn__visit_age",
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
          "ph_y_pds__m_002")


desc_data <- create_dataset(
  dir_data = data_root,
  study = "abcd",
  vars = vars,
  value_to_na = TRUE,
  bind_shadow = TRUE,
  value_to_label = TRUE
)

data <- create_dataset(
  dir_data = data_root,
  study = "abcd",
  vars = vars,
  value_to_na = TRUE,
  bind_shadow = TRUE
)


## quick descriptives 
table <- CreateTableOne(vars = vars, strata = c("ab_g_stc__cohort_sex", "session_id"), data = desc_data)
print(table, quote = TRUE)

# test frequency of caregiver switching 
table(desc_data$ab_g_dyn__visit__day1_inform)

data <- data[order(data$participant_id, data$session_id), ]
caregiver_changes <- data %>%
  group_by(participant_id) %>%
  summarise(
    n_caregivers = n_distinct(ab_g_dyn__visit__day1_inform, na.rm = TRUE)
  )
table(caregiver_changes$n_caregivers > 1)
sum(caregiver_changes$n_caregivers > 1)
caregiver_changes %>%
  filter(n_caregivers > 1)


# characterize missing puberty data patterns


# test correlations between items and across reporters 

# test for longitudinal measurement invariance 

# -----------------------------------------------------------------------------
# 3. LONGITUDINAL MEASUREMENT INVARIANCE (MULTI-GROUP OVER TIME)
#    (works with group.equal + WLSMV for ordered indicators)
# -----------------------------------------------------------------------------

# map ABCD session ids to readable time labels
session_map <- c(
  "ses-00A" = "BL",
  "ses-01A" = "FU1",
  "ses-02A" = "FU2",
  "ses-03A" = "FU3",
  "ses-04A" = "FU4",
  "ses-05A" = "FU5",
  "ses-06A" = "FU6"
)

items <- c(
  "ph_p_pds_001",
  "ph_p_pds_002",
  "ph_p_pds_003",
  "ph_p_pds__f_001",
  "ph_p_pds__f_002"
)

dat_long <- data %>%
  mutate(
    time = recode(session_id, !!!session_map),
    sex  = ab_g_stc__cohort_sex
  ) %>%
  filter(!is.na(time)) %>%
  select(participant_id, time, sex, all_of(items)) %>%
  # make sure items are treated as ordered categorical
  mutate(across(all_of(items), ~ as.ordered(.)))

dat_long_f <- dat_long %>% 
  filter(sex == 2)

# single-factor model (same indicators at each wave; invariance via group.equal)
model_pds_1f <- '
  pds =~ ph_p_pds_001 + ph_p_pds_002 + ph_p_pds_003 + ph_p_pds__f_001 + ph_p_pds__f_002
'

# CONFIGURAL
fit_config <- cfa(
  model_pds_1f,
  data = dat_long_f,
  group = "time",
  ordered = items,
  estimator = "WLSMV",
  parameterization = "theta"
)

# METRIC / WEAK (equal loadings)
fit_metric <- cfa(
  model_pds_1f,
  data = dat_long_f,
  group = "time",
  ordered = items,
  estimator = "WLSMV",
  parameterization = "theta",
  group.equal = c("loadings")
)

# SCALAR / STRONG (equal loadings + thresholds)
fit_scalar <- cfa(
  model_pds_1f,
  data = dat_long_f,
  group = "time",
  ordered = items,
  estimator = "WLSMV",
  parameterization = "theta",
  group.equal = c("loadings", "thresholds")
)

# STRICT (equal loadings + thresholds + residual variances)
fit_strict <- cfa(
  model_pds_1f,
  data = dat_long_f,
  group = "time",
  ordered = items,
  estimator = "WLSMV",
  parameterization = "theta",
  group.equal = c("loadings", "thresholds", "residuals")
)

# WLSMV: use DIFFTEST-style comparisons
lavTestLRT(fit_config, fit_metric, fit_scalar, fit_strict)

# handy fit summary table
invariance_fits <- rbind(
  config = fitMeasures(fit_config, c("cfi","tli","rmsea","srmr")),
  metric = fitMeasures(fit_metric, c("cfi","tli","rmsea","srmr")),
  scalar = fitMeasures(fit_scalar, c("cfi","tli","rmsea","srmr")),
  strict = fitMeasures(fit_strict, c("cfi","tli","rmsea","srmr"))
)
invariance_fits

