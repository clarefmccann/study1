## CLEANING PUBERTY DATA FOR TRANSLATION TO TS

library(dplyr)

data_root = "/Users/clarefmccann/University\ of\ Oregon\ Dropbox/Clare\ McCann/mine/projects/abcd-projs/abcd-data-release-6.0/cfm/"

age_sex <- read.csv(paste0(data_root,"abcd-general/ab_g_dyn.csv")) %>%
  select(participant_id, session_id, site, ab_g_dyn__visit_age, ab_g_stc__cohort_sex) %>%
  rename("id" = "participant_id",
         "wave" = "session_id",
         "age" = "ab_g_dyn__visit_age", # unit is years
         "sex" = "ab_g_stc__cohort_sex") # 1 = male, 2 = female

raw_PDS_p <- read.csv(paste0(data_root, "physical-health/ph_puberty.csv")) %>%
  rename("id" = "participant_id",
         "wave" = "session_id") %>%
  mutate(wave_index = case_when(
    wave == "ses-00A" ~ 0,
    wave == "ses-00M" ~ 0.5,
    wave == "ses-01A" ~ 1,
    wave == "ses-01M" ~ 1.5,
    wave == "ses-02A" ~ 2,
    wave == "ses-02M" ~ 2.5,
    wave == "ses-03A" ~ 3,
    wave == "ses-03M" ~ 3.5,
    wave == "ses-04A" ~ 4,
    wave == "ses-04M" ~ 4.5,
    wave == "ses-05A" ~ 5,
    wave == "ses-05M" ~ 5.5,
    wave == "ses-06A" ~ 6,
    TRUE ~ NA_real_ 
  )) %>% 
  select(-contains(c("language", "002__"))) %>%
  rename("peta_p" = "ph_p_pds_001", # height
         "petb_p" = "ph_p_pds_002", # body hair
         "petc_p" = "ph_p_pds_003", # skin
         "petdf_p" = "ph_p_pds__f_001", # breast growth
         "fpete_p" = "ph_p_pds__f_002", # menarche (y/n)
         "petdm_p" = "ph_p_pds__m_001", # voice
         "mpete_p" = "ph_p_pds__m_002",
         "pdss_cat_f_p" = "ph_p_pds__f_categ",
         "pdss_cont_f_p" = "ph_p_pds__f_mean",
         "pdss_cat_m_p" = "ph_p_pds__m_categ",
         "pdss_cont_m_p" = "ph_p_pds__m_mean") %>% # facial hair
  ungroup() 

# summary_PDS_p <- raw_PDS_p %>% 
#   select(id, wave, ph_p_pds__f_categ, ph_p_pds__m_categ)

raw_PDS_y <- read.csv(paste0(data_root, "physical-health/ph_puberty.csv")) %>%
  rename("id" = "participant_id",
         "wave" = "session_id") %>%
  mutate(wave_index = case_when(
    wave == "ses-00A" ~ 0,
    wave == "ses-00M" ~ 0.5,
    wave == "ses-01A" ~ 1,
    wave == "ses-01M" ~ 1.5,
    wave == "ses-02A" ~ 2,
    wave == "ses-02M" ~ 2.5,
    wave == "ses-03A" ~ 3,
    wave == "ses-03M" ~ 3.5,
    wave == "ses-04A" ~ 4,
    wave == "ses-04M" ~ 4.5,
    wave == "ses-05A" ~ 5,
    wave == "ses-05M" ~ 5.5,
    wave == "ses-06A" ~ 6,
    TRUE ~ NA_real_ 
  )) %>% 
  select(-contains(c("language", "002__")) )%>%
  rename("peta_y" = "ph_y_pds_001", # height
         "petb_y" = "ph_y_pds_002", # body hair
         "petc_y" = "ph_y_pds_003", # skin
         "petdf_y" = "ph_y_pds__f_001", # breast growth
         "fpete_y" = "ph_y_pds__f_002", # menarche (y/n)
         "petdm_y" = "ph_y_pds__m_001", # voice
         "mpete_y" = "ph_y_pds__m_002",  # facial hair
         "pdss_cat_f_y" = "ph_y_pds__f_categ",
         "pdss_cont_f_y" = "ph_y_pds__f_mean",
         "pdss_cat_m_y" = "ph_y_pds__m_categ",
         "pdss_cont_m_y" = "ph_y_pds__m_mean"
         ) %>%
  ungroup() 

# summary_PDS_y <- raw_PDS_y %>% 
#   select(id, wave, ph_y_pds__f_categ, ph_y_pds__m_categ)
# 
# summary_PDS <- full_join(summary_PDS_p, full_join(summary_PDS_y, age_sex, by = c("id", "wave"))) %>% 
#   mutate(parent_summary = coalesce(ph_p_pds__f_categ, ph_p_pds__m_categ),
#          youth_summary = coalesce(ph_y_pds__f_categ, ph_y_pds__m_categ)) %>%
#   select(-starts_with("ph_p_pds__m_category"), 
#          -starts_with("ph_p_pds__f_category"),
#          -starts_with("ph_y_pds__m_category"), 
#          -starts_with("ph_y_pds__f_category")) 
# 
# rm(summary_PDS_p, summary_PDS_y)

pet_vars_y <- c("peta_y", "petb_y", "petc_y", "petdf_y", "fpete_y", "petdm_y", "mpete_y", "pdss_cat_f_y", "pdss_cont_f_y", "pdss_cat_m_y", "pdss_cont_m_y")
pet_vars_p <- c("peta_p", "petb_p", "petc_p", "petdf_p", "fpete_p", "petdm_p", "mpete_p", "pdss_cat_f_p", "pdss_cont_f_p", "pdss_cat_m_p", "pdss_cont_m_p")

raw_PDS <- raw_PDS_y %>%
  left_join(raw_PDS_p %>%
              select(id, wave_index, all_of(pet_vars_p)),
            by = c("id", "wave_index")) %>%
  group_by(id) %>%
  mutate(across(all_of(pet_vars_y), 
                ~ if_else(wave_index <= 2, get(sub("_y", "_p", cur_column())), .), 
                .names = "{col}")) %>%
  select(-ends_with("_p")) %>%
  ungroup()

raw_PDS <- left_join(raw_PDS, age_sex, by = c("id", "wave"))
raw_PDS_y <- left_join(raw_PDS_y, age_sex, by = c("id", "wave"))
raw_PDS_p <- left_join(raw_PDS_p, age_sex, by = c("id", "wave"))


## back to sex specific raw dfs

raw_PDS_f <- raw_PDS %>%
  filter(sex == 2) %>%
  select(id, wave, peta_y, petb_y, petc_y, petdf_y, fpete_y, pdss_cat_f_y, pdss_cont_f_y, age) %>%
  rename("peta" = "peta_y", 
         "petb" = "petb_y", 
         "petc" = "petc_y", 
         "petd" = "petdf_y",
         "fpete" ="fpete_y",
         "pdss_cont" = "pdss_cont_f_y",
         "pdss_cat" = "pdss_cat_f_y") %>% 
  filter(
    if_all(
      starts_with(c("pet", "pdss", "fpete")),
      ~ !(. %in% c(999, 777)) & !is.na(.)
    ))

raw_PDS_m <- raw_PDS %>%
  filter(sex == 1) %>%
  select(id, wave, peta_y, petb_y, petc_y, petdm_y, mpete_y, pdss_cat_m_y, pdss_cont_m_y, age) %>%
  rename("peta" = "peta_y", 
         "petb" = "petb_y", 
         "petc" = "petc_y", 
         "petd" = "petdm_y",
         "mpete" ="mpete_y",
         "pdss_cont" = "pdss_cont_m_y",
         "pdss_cat" = "pdss_cat_m_y") %>% 
  filter(
    if_all(
      starts_with(c("pet", "pdss", "mpete")),
      ~ !(. %in% c(999, 777)) & !is.na(.)
      ))

write.csv(raw_PDS_f, file = paste0(data_root, "physical-health/puberty/parentyouth_raw_PDS_f.csv"))
write.csv(raw_PDS_m, file = paste0(data_root, "physical-health/puberty/parentyouth_raw_PDS_m.csv"))
#write.csv(summary_PDS, file = paste0(data_root, "physical-health/puberty/parentyouth_summary_PDS.csv"))


raw_PDS_f_y <- raw_PDS_y %>%
  filter(sex == 2) %>%
  select(id, sex, wave, peta_y, petb_y, petc_y, petdf_y, fpete_y, pdss_cont_f_y, pdss_cat_f_y, age) %>%
  rename("peta" = "peta_y", 
         "petb" = "petb_y", 
         "petc" = "petc_y", 
         "petd" = "petdf_y",
         "fpete" ="fpete_y",
         "pdss_cont" = "pdss_cont_f_y",
         "pdss_cat" = "pdss_cat_f_y") %>% 
  filter(
    if_all(
      starts_with(c("pet", "pdss", "mpete")),
      ~ !(. %in% c(999, 777)) & !is.na(.)
    ))

raw_PDS_m_y <- raw_PDS_y %>%
  filter(sex == 1) %>%
  select(id, sex, wave, peta_y, petb_y, petc_y, petdm_y, mpete_y, pdss_cont_m_y, pdss_cat_m_y, age) %>%
  rename("peta" = "peta_y", 
         "petb" = "petb_y", 
         "petc" = "petc_y", 
         "petd" = "petdm_y",
         "mpete" ="mpete_y",
         "pdss_cont" = "pdss_cont_m_y",
         "pdss_cat" = "pdss_cat_m_y") %>% 
  filter(
    if_all(
      starts_with(c("pet", "pdss", "mpete")),
      ~ !(. %in% c(999, 777)) & !is.na(.)
    ))

write.csv(raw_PDS_f_y, file = paste0(data_root, "physical-health/puberty/youth_raw_PDS_f.csv"))
write.csv(raw_PDS_m_y, file = paste0(data_root, "physical-health/puberty/youth_raw_PDS_m.csv"))

raw_PDS_f_p <- raw_PDS_p %>%
  filter(sex == 2) %>%
  select(id, wave, sex, peta_p, petb_p, petc_p, petdf_p, fpete_p, pdss_cont_f_p, pdss_cat_f_p, age) %>%
  rename("peta" = "peta_p", 
         "petb" = "petb_p", 
         "petc" = "petc_p", 
         "petd" = "petdf_p",
         "fpete" ="fpete_p",
         "pdss_cont" = "pdss_cont_f_p",
         "pdss_cat" = "pdss_cat_f_p") %>% 
  filter(
    if_all(
      starts_with(c("pet", "pdss", "mpete")),
      ~ !(. %in% c(999, 777)) & !is.na(.)
    ))

raw_PDS_m_p <- raw_PDS_p %>%
  filter(sex == 1) %>%
  select(id, wave, sex, peta_p, petb_p, petc_p, petdm_p, mpete_p, pdss_cont_m_p, pdss_cat_m_p, age) %>%
  rename("peta" = "peta_p", 
         "petb" = "petb_p", 
         "petc" = "petc_p", 
         "petd" = "petdm_p",
         "mpete" ="mpete_p",
         "pdss_cont" = "pdss_cont_m_p",
         "pdss_cat" = "pdss_cat_m_p") %>% 
  filter(
    if_all(
      starts_with(c("pet", "pdss", "mpete")),
      ~ !(. %in% c(999, 777)) & !is.na(.)
    ))

write.csv(raw_PDS_f_p, file = paste0(data_root, "physical-health/puberty/parent_raw_PDS_f.csv"))
write.csv(raw_PDS_m_p, file = paste0(data_root, "physical-health/puberty/parent_raw_PDS_m.csv"))

rm(raw_PDS_p, raw_PDS_y, age_sex, pet_vars_p, pet_vars_y)
