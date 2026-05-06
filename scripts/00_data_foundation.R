## 00_data_foundation.R
## Loads raw ABCD PDS data (annual waves only), attaches covariates
## (age, sex, BMI, race/ethnicity, caregiver), and writes clean
## long-format datasets for each sex × reporter combination.

# export DATA_DIR="/u/project/silvers/data/ABCD/ABCD-release-6.0/cfm/physical-health/puberty"
# export OUT_DIR="/u/home/c/clarefmc/projects/abcd-projs/dissertation/study1/outputs"
# Rscript 01_data_foundation.R

library(dplyr)
library(tidyr)
library(NBDCtools)

set.seed(90025)

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
root_path <- Sys.getenv("HOME_DIR")
if (!nzchar(root_path)) {
  root_path <- Sys.getenv("HOME")
}

data_root <- file.path(
  root_path,
  "projects/abcd-projs/abcd-data-release-6.0/nbdc-tools-data"
)
if (!dir.exists(data_root)) {
  stop("Cannot locate nbdc-tools-data: ", data_root)
}

out_root <- file.path(
  root_path,
  "projects/abcd-projs/abcd-data-release-6.0/cfm/physical-health/puberty"
)

# ---------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------
vars <- c(
  "ab_g_dyn__visit_age",
  "ab_g_stc__cohort_sex",
  "ab_g_stc__cohort_ethnrace__mhisp",
  "ab_g_dyn__visit__day1_inform",
  "ab_g_dyn__design_site", # site — verify name if create_dataset errors
  "ph_y_anthr__height_mean", # inches
  "ph_y_anthr__weight_mean", # lbs  — if unavailable, BMI will be NA
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

raw <- create_dataset(
  dir_data = data_root,
  study = "abcd",
  vars = vars,
  value_to_na = TRUE,
  bind_shadow = FALSE
)

# ---------------------------------------------------------------------------
# ANNUAL WAVES ONLY + WAVE LABEL
# ---------------------------------------------------------------------------
wave_labels <- c(
  "ses-00A" = "bl",
  "ses-01A" = "fu1",
  "ses-02A" = "fu2",
  "ses-03A" = "fu3",
  "ses-04A" = "fu4",
  "ses-05A" = "fu5",
  "ses-06A" = "fu6"
)

data <- raw %>%
  filter(session_id %in% names(wave_labels)) %>%
  mutate(
    wave = recode(session_id, !!!wave_labels),
    wave = factor(wave, levels = wave_labels)
  ) %>%
  rename(
    id = participant_id,
    age = ab_g_dyn__visit_age,
    sex = ab_g_stc__cohort_sex,
    race = ab_g_stc__cohort_ethnrace__mhisp,
    caregiver = ab_g_dyn__visit__day1_inform,
    site = ab_g_dyn__design_site,
    height_in = ph_y_anthr__height_mean,
    weight_lb = ph_y_anthr__weight_mean,
    peta_p = ph_p_pds_001,
    petb_p = ph_p_pds_002,
    petc_p = ph_p_pds_003,
    petdf_p = ph_p_pds__f_001,
    fpete_p = ph_p_pds__f_002,
    petdm_p = ph_p_pds__m_001,
    mpete_p = ph_p_pds__m_002,
    peta_y = ph_y_pds_001,
    petb_y = ph_y_pds_002,
    petc_y = ph_y_pds_003,
    petdf_y = ph_y_pds__f_001,
    fpete_y = ph_y_pds__f_002,
    petdm_y = ph_y_pds__m_001,
    mpete_y = ph_y_pds__m_002
  )

# ---------------------------------------------------------------------------
# DERIVED VARIABLES
# ---------------------------------------------------------------------------

# BMI
data <- data %>%
  mutate(bmi = 703 * weight_lb / height_in^2)

# age-standardized BMI z-score within sex (for use as covariate)
data <- data %>%
  group_by(sex, wave) %>%
  mutate(bmi_z = as.numeric(scale(bmi))) %>%
  ungroup()

# caregiver switch flags (relative to each participant's prior annual wave)
data <- data %>%
  arrange(id, wave) %>%
  group_by(id) %>%
  mutate(
    caregiver_switched = caregiver != lag(caregiver) &
      !is.na(caregiver) &
      !is.na(lag(caregiver)),
    n_cg_switches = sum(caregiver_switched, na.rm = TRUE),
    ever_switched_cg = n_cg_switches > 0
  ) %>%
  ungroup()

# ---------------------------------------------------------------------------
# RECODE INVALID VALUES TO NA
# (999 = don't know, 777 = refused in ABCD)
# ---------------------------------------------------------------------------
pds_items <- c(
  "peta_p",
  "petb_p",
  "petc_p",
  "petdf_p",
  "fpete_p",
  "petdm_p",
  "mpete_p",
  "peta_y",
  "petb_y",
  "petc_y",
  "petdf_y",
  "fpete_y",
  "petdm_y",
  "mpete_y"
)

data <- data %>%
  mutate(across(
    all_of(pds_items),
    ~ if_else(. %in% c(777, 999), NA_real_, as.numeric(.))
  ))

# ---------------------------------------------------------------------------
# PDS COMPOSITES
# 4-item mean (shared items peta–petd, excluding binary menarche/facial hair)
# Females: mean(peta, petb, petc, petdf)
# Males:   mean(peta, petb, petc, petdm)
# Scored separately per reporter; NA if any item missing
# ---------------------------------------------------------------------------
data <- data %>%
  mutate(
    pds_comp_f_p = rowMeans(
      cbind(peta_p, petb_p, petc_p, petdf_p),
      na.rm = FALSE
    ),
    pds_comp_f_y = rowMeans(
      cbind(peta_y, petb_y, petc_y, petdf_y),
      na.rm = FALSE
    ),
    pds_comp_m_p = rowMeans(
      cbind(peta_p, petb_p, petc_p, petdm_p),
      na.rm = FALSE
    ),
    pds_comp_m_y = rowMeans(
      cbind(peta_y, petb_y, petc_y, petdm_y),
      na.rm = FALSE
    )
  )

# ---------------------------------------------------------------------------
# SEX × REPORTER LONG DATASETS
# ---------------------------------------------------------------------------
covariate_cols <- c(
  "id",
  "wave",
  "age",
  "sex",
  "race",
  "site",
  "bmi",
  "bmi_z",
  "caregiver",
  "caregiver_switched",
  "n_cg_switches",
  "ever_switched_cg"
)

female_parent <- data %>%
  filter(sex == 2) %>%
  select(
    all_of(covariate_cols),
    peta = peta_p,
    petb = petb_p,
    petc = petc_p,
    petd = petdf_p,
    fpete = fpete_p,
    pds_comp = pds_comp_f_p
  ) %>%
  mutate(reporter = "parent")

female_youth <- data %>%
  filter(sex == 2) %>%
  select(
    all_of(covariate_cols),
    peta = peta_y,
    petb = petb_y,
    petc = petc_y,
    petd = petdf_y,
    fpete = fpete_y,
    pds_comp = pds_comp_f_y
  ) %>%
  mutate(reporter = "youth")

male_parent <- data %>%
  filter(sex == 1) %>%
  select(
    all_of(covariate_cols),
    peta = peta_p,
    petb = petb_p,
    petc = petc_p,
    petd = petdm_p,
    mpete = mpete_p,
    pds_comp = pds_comp_m_p
  ) %>%
  mutate(reporter = "parent")

male_youth <- data %>%
  filter(sex == 1) %>%
  select(
    all_of(covariate_cols),
    peta = peta_y,
    petb = petb_y,
    petc = petc_y,
    petd = petdm_y,
    mpete = mpete_y,
    pds_comp = pds_comp_m_y
  ) %>%
  mutate(reporter = "youth")

# combined long (all four) — useful for cross-reporter models
all_long <- bind_rows(female_parent, female_youth, male_parent, male_youth)

# ---------------------------------------------------------------------------
# WRITE OUTPUTS
# ---------------------------------------------------------------------------
write.csv(
  female_parent,
  file.path(out_root, "female_parent_long.csv"),
  row.names = FALSE
)
write.csv(
  female_youth,
  file.path(out_root, "female_youth_long.csv"),
  row.names = FALSE
)
write.csv(
  male_parent,
  file.path(out_root, "male_parent_long.csv"),
  row.names = FALSE
)
write.csv(
  male_youth,
  file.path(out_root, "male_youth_long.csv"),
  row.names = FALSE
)
write.csv(all_long, file.path(out_root, "all_long.csv"), row.names = FALSE)

cat("Written to:", out_root, "\n")
cat(
  "Rows per dataset — female parent:",
  nrow(female_parent),
  "| female youth:",
  nrow(female_youth),
  "| male parent:",
  nrow(male_parent),
  "| male youth:",
  nrow(male_youth),
  "\n"
)
