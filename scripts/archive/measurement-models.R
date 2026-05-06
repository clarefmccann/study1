# -----------------------------------------------------------------------------
# 1. SETUP
# -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(readr)
library(lavaan)
library(semTools)

set.seed(90025)

# define file paths (please update if necessary)
pub_data_root <- file.path(
  "/users/clarefmccann/university of oregon dropbox",
  "clare mccann/mine/projects/abcd-projs/abcd-data-release-6.0/physical-health/puberty"
)

# -----------------------------------------------------------------------------
# 2. DATA LOADING AND PREPARATION
# -----------------------------------------------------------------------------

# helper function to read, rename, and select columns
load_and_prep <- function(file_path, sex_code, reporter) {
  
  if (reporter == "parent") {
    name_map <- c(pete_p = "fpete", peta_p = "peta", petb_p = "petb", petc_p = "petc", petd_p = "petd", pdss_p = "pdss")
    if (sex_code == 0) name_map["pete_p"] <- "mpete" # use mpete for males
  } else {
    name_map <- c(pete_y = "fpete", peta_y = "peta", petb_y = "petb", petc_y = "petc", petd_y = "petd", pdss_y = "pdss")
    if (sex_code == 0) name_map["pete_y"] <- "mpete"
  }
  
  read_csv(file_path, show_col_types = FALSE) %>%
    rename(any_of(name_map)) %>%
    select(id, wave, age, starts_with("pet"), starts_with("pdss")) %>%
    mutate(sex = sex_code)
}

# load all four datasets
pub_f_p <- load_and_prep(file.path(pub_data_root, "filtered_parent_tannerstages_f.csv"), 1, "parent")
pub_m_p <- load_and_prep(file.path(pub_data_root, "filtered_parent_tannerstages_m.csv"), 0, "parent")
pub_f_y <- load_and_prep(file.path(pub_data_root, "filtered_youth_tannerstages_f.csv"), 1, "youth")
pub_m_y <- load_and_prep(file.path(pub_data_root, "filtered_youth_tannerstages_m.csv"), 0, "youth")

# define wave mapping
map_wave <- c("ses-00A" = "bl", "ses-01A" = "fu1", "ses-02A" = "fu2", 
              "ses-03A" = "fu3", "ses-04A" = "fu4", "ses-05A" = "fu5", "ses-06A" = "fu6")
waves <- c("bl", "fu1", "fu2", "fu3", "fu4", "fu5", "fu6")

# combine, join, and pivot to wide format in one pipeline
wide_data <- bind_rows(pub_f_p, pub_m_p) %>%
  full_join(bind_rows(pub_f_y, pub_m_y), by = c("id", "wave", "sex", "age"), suffix = c("_p", "_y")) %>%
  mutate(wave = recode(wave, !!!map_wave)) %>%
  pivot_wider(
    id_cols = c(id, sex),
    names_from = wave,
    values_from = c(starts_with(c("age", "pet")), starts_with("pdss")),
    names_glue = "{.value}_{wave}"
  ) %>% 
  select(-contains(c("m.x","f.x", "m.y", "f.y")))
# we remove the select() line that tried to drop .x and .y columns, 
# as suffix() fixes this.

# create data subsets for males and females, converting all puberty items to numeric
females_wide <- wide_data %>%
  filter(sex == 1) %>%
  mutate(across(starts_with("peta") | starts_with("petb") | starts_with("petc") | 
                  starts_with("petd") | starts_with("pete"), as.ordered))

males_wide <- wide_data %>%
  filter(sex == 0) %>%
  mutate(across(starts_with("peta") | starts_with("petb") | starts_with("petc") | 
                  starts_with("petd") | starts_with("pete"), as.ordered))


# -----------------------------------------------------------------------------
# 3. GENERALIZED INVARIANCE TESTING FUNCTION (WITH PARTIAL SCALAR LOOP)
# -----------------------------------------------------------------------------
# this function automates the entire process for any pair of adjacent waves,
# for any specified reporter (_p or _y).
# it now includes an iterative loop for partial scalar invariance.

test_invariance_pair <- function(data, wave1_suffix, wave2_suffix, reporter_suffix, sex_code, all_ordered_vars, max_free_intercepts = 2) {
  
  cat(paste("\n\n--- testing invariance between", wave1_suffix, "and", wave2_suffix, "---\n"))
  
  
  # dynamically create variable lists based on the reporter
  # use sex_code to determine the correct sex-specific item name
  base_vars <- paste0(c("peta", "petb", "petc", "petd", "pete"), reporter_suffix)
  wave1_vars <- paste0(base_vars, "_", wave1_suffix)
  wave2_vars <- paste0(base_vars, "_", wave2_suffix)
  all_wave_vars <- c(wave1_vars, wave2_vars) # <-- add this line to get all item names
  current_ordered_vars <- c(wave1_vars, wave2_vars)
  
  missing_cols <- setdiff(all_wave_vars, names(data))
  if (length(missing_cols) > 0) {
    cat(paste("error: the following variables are missing from the dataset for this group:\n", 
              paste(missing_cols, collapse = ", "), "\n"))
    cat("skipping this pair.\n")
    return()
  }
  
  # [EDIT] check for zero-variance items using a method safe for ordered factors
  has_variance <- function(x) { length(unique(na.omit(x))) > 1 }
  variances <- sapply(data[all_wave_vars], has_variance)
  zero_var_items <- names(variances[variances == FALSE])
  
  if (length(zero_var_items) > 0) {
    cat(paste("warning: the following items have zero variance (only one unique value) for this group and wave-pair:\n",
              paste(zero_var_items, collapse = ", "), "\n"))
    cat("model estimation is impossible. skipping this pair.\n")
    return() # skip to the next pair
  }
  # --- [END] pre-flight check ---
  
  # dynamically create model syntax
  model_config <- paste(
    paste("puberty_", wave1_suffix, " =~ ", paste(wave1_vars, collapse = " + ")),
    paste("puberty_", wave2_suffix, " =~ ", paste(wave2_vars, collapse = " + ")),
    sep = "\n"
  )
  
  print(model_config)
  
  model_metric <- paste(
    paste0("puberty_", wave1_suffix, " =~ 1*", wave1_vars[1], " + l2*", wave1_vars[2], " + l3*", wave1_vars[3], " + l4*", wave1_vars[4], " + l5*", wave1_vars[5]),
    paste0("puberty_", wave2_suffix, " =~ 1*", wave2_vars[1], " + l2*", wave2_vars[2], " + l3*", wave2_vars[3], " + l4*", wave2_vars[4], " + l5*", wave2_vars[5]),
    sep = "\n"
  )
  
  print(model_metric)
  
  fit_config <- cfa(model_config, data = data, estimator = "wlsmv", missing = "pairwise", ordered = all_ordered_vars)
  fit_metric <- cfa(model_metric, data = data, estimator = "wlsmv", missing = "pairwise", ordered = all_ordered_vars)
  
  cfi_config <- fitmeasures(fit_config, "cfi.scaled")
  cfi_metric <- fitmeasures(fit_metric, "cfi.scaled")
  
  # --- establish baseline metric model ---
  successful_metric_fit <- NULL
  successful_metric_syntax <- NULL
  metric_level <- "none"
  
  if (cfi_config - cfi_metric <= 0.01) {
    cat("\nfull metric invariance achieved.\n")
    successful_metric_fit <- fit_metric
    successful_metric_syntax <- model_metric
    metric_level <- "full metric"
  } else {
    cat("\nfull metric invariance failed. testing for partial metric invariance...\n")
    mod_indices <- modindices(fit_metric) %>% filter(op == "=~") %>% arrange(desc(mi))
    problem_item_base_name <- gsub("_[^_]+$", "", mod_indices$rhs[1])
    problem_item_index <- which(base_vars == problem_item_base_name)
    problem_item_label <- paste0("l", problem_item_index, "_free")
    cat(paste("most non-invariant item appears to be:", problem_item_base_name, "\n"))
    
    metric_lines <- strsplit(model_metric, "\n")[[1]]
    metric_lines[2] <- gsub(paste0("l", problem_item_index, "*"), paste0(problem_item_label, "*"), metric_lines[2], fixed = TRUE)
    model_partial_metric <- paste(metric_lines, collapse = "\n")
    print(model_partial_metric)
    
    fit_partial_metric <- cfa(model_partial_metric, data = data, estimator = "wlsmv", missing = "pairwise", ordered = all_ordered_vars)
    cfi_partial_metric <- fitmeasures(fit_partial_metric, "cfi.scaled")
    
    if (cfi_config - cfi_partial_metric <= 0.01) {
      cat("partial metric invariance achieved.\n")
      successful_metric_fit <- fit_partial_metric
      successful_metric_syntax <- model_partial_metric
      metric_level <- "partial metric"
    } else {
      cat("partial metric invariance failed. stopping analysis for this pair.\n")
      return()
    }
  }
  
  # --- [NEW] iterative test for scalar invariance ---
  cat("\n--- testing for scalar invariance ---\n")
  
  # create the full set of intercept constraints
  intercept_constraints <- sapply(base_vars, function(item) {
    item_index <- which(base_vars == item)
    paste0(item, "_", wave1_suffix, " ~ i", item_index, "*1; ", item, "_", wave2_suffix, " ~ i", item_index, "*1")
  })
  
  cfi_metric_baseline <- fitmeasures(successful_metric_fit, "cfi.scaled")
  
  # initialize loop variables
  current_intercept_constraints <- intercept_constraints
  items_freed <- c()
  scalar_level <- "none" # default to failure
  
  # loop from 0 (full scalar) up to max_free_intercepts
  for (i in 0:max_free_intercepts) {
    
    if (i == 0) {
      cat("testing full scalar model...\n")
    } else {
      cat(paste("testing partial scalar model, item(s) freed:", paste(items_freed, collapse = ", "), "\n"))
    }
    
    # build the scalar model syntax for this iteration
    current_scalar_syntax <- paste(
      successful_metric_syntax, 
      paste(current_intercept_constraints, collapse = "\n"), 
      sep = "\n"
    )
    
    # fit the current scalar model
    fit_current_scalar <- cfa(current_scalar_syntax, data = data, estimator = "wlsmv", missing = "pairwise", ordered = all_ordered_vars)
    
    # check for convergence
    if (!lavInspect(fit_current_scalar, "converged")) {
      cat("model failed to converge. stopping scalar tests for this pair.\n")
      scalar_level <- "convergence error"
      break
    }
    
    cfi_current_scalar <- fitmeasures(fit_current_scalar, "cfi.scaled")
    delta_cfi <- cfi_metric_baseline - cfi_current_scalar
    
    # check if invariance is met
    if (delta_cfi <= 0.01) {
      scalar_level <- ifelse(i == 0, "full scalar", "partial scalar")
      cat(paste(scalar_level, "invariance achieved (delta_cfi =", round(delta_cfi, 4), ").\n"))
      break # success! exit the loop
    }
  
    if (i < max_free_intercepts) {
      cat(paste("model failed (delta_cfi =", round(delta_cfi, 4), "). finding worst item to free...\n"))
      
      # --- [NEW] use lavTestScore() to test intercept constraints ---
      score_tests_obj <- lavTestScore(fit_current_scalar)
      pt <- parTable(fit_current_scalar)
      
      # --- [NEW] check if '$uni' or 'lhs' column exists ---
      if (is.null(score_tests_obj$uni) || !("lhs" %in% names(score_tests_obj$uni))) {
        cat("could not find '$uni' or 'lhs' column in lavTestScore() output for scalar test. stopping.\n")
        
        # add diagnostic summary on failure
        cat("\n--- diagnostic summary for failed scalar model ---\n")
        cat("this often happens due to non-positive definite matrices or other estimation issues.\n")
        cat("review the warnings and fit measures below:\n")
        print(summary(fit_current_scalar, fit.measures = TRUE))
        cat("\n--- end diagnostic summary ---\n")
        break
      }
      score_tests_raw <- score_tests_obj$uni
      # --- [END] check ---
      
      # --- [NEW] create a lookup table to map internal plabels (e.g., .p11.) to our labels (e.g., i1)
      intercept_labels <- paste0("i", 1:5)
      lookup <- pt %>% 
        filter(label %in% intercept_labels) %>% 
        select(label, plabel) %>% 
        distinct()
      
      # join the score tests with our labels
      score_tests <- left_join(score_tests_raw, lookup, by = c("lhs" = "plabel"))
      
      # filter for our intercept labels ("i1", "i2", "i3", "i4", "i5")
      mod_indices <- score_tests %>%
        filter(label %in% intercept_labels) %>%
        arrange(desc(X2)) # X2 is the chi-square (score test)
      
      if(nrow(mod_indices) == 0) {
        cat("could not find valid score tests for intercepts. stopping.\n")
        
        # add diagnostic summary on failure
        cat("\n--- diagnostic summary for failed scalar model ---\n")
        cat("this often happens due to non-positive definite matrices or zero-variance items.\n")
        cat("review the warnings and fit measures below:\n")
        print(summary(fit_current_scalar, fit.measures = TRUE))
        cat("\n--- end diagnostic summary ---\n")
        break
      }
      
      # get the label of the worst-offending intercept
      problem_item_label <- mod_indices$label[1] # e.g., "i3"
      
      # find the base_vars index from the label (e.g., "i3" -> 3)
      problem_item_index <- as.numeric(sub("i", "", problem_item_label))
      problem_item_base_name <- base_vars[problem_item_index]
      
      # prevent an item from being freed twice
      if (problem_item_base_name %in% items_freed) {
        cat("worst item already freed. stopping.\n")
        break
      }
      
      items_freed <- c(items_freed, problem_item_base_name)
      cat(paste("freeing intercept for:", problem_item_base_name, "(mi =", round(mod_indices$X2[1], 2), ")\n"))
      
      # remove the constraint line for this item
      line_to_remove <- grep(paste0("^", problem_item_base_name), current_intercept_constraints, value = TRUE)
      current_intercept_constraints <- setdiff(current_intercept_constraints, line_to_remove)
      
    } else {
      # this was the last iteration, and it failed
      cat(paste("partial scalar invariance failed after freeing", i, "item(s).\n"))
      scalar_level <- "metric"
    }
  } # end scalar loop
  
  
  # --- final conclusion for the pair ---
  cat(paste("\n--- final conclusion for", wave1_suffix, "vs.", wave2_suffix, "---\n"))
  cat(paste("highest level of metric invariance achieved:", metric_level, "\n"))
  cat(paste("highest level of scalar invariance achieved:", scalar_level, "\n"))
  
  if (scalar_level %in% c("full scalar", "partial scalar")) {
    cat("implication: latent mean comparisons are valid for this transition.\n")
  } else {
    cat("implication: latent mean comparisons are not valid for this transition.\n")
  }
}

# -----------------------------------------------------------------------------
# 4. MAIN ANALYSIS: LOOP THROUGH GROUPS AND WAVE PAIRS
# -----------------------------------------------------------------------------

# define the groups to analyze
analysis_groups <- list(
  list(name = "females - parent report", data = females_wide, reporter_suffix = "_p", sex_code = 1),
  list(name = "males - parent report",   data = males_wide,   reporter_suffix = "_p", sex_code = 0),
  list(name = "females - youth report",  data = females_wide, reporter_suffix = "_y", sex_code = 1),
  list(name = "males - youth report",    data = males_wide,   reporter_suffix = "_y", sex_code = 0)
)

# define the wave pairs to test
wave_pairs <- list(
  c("bl", "fu1"),
  c("fu1", "fu2"),
  c("fu2", "fu3"),
  c("fu3", "fu4"),
  c("fu4", "fu5"),
  c("fu5", "fu6")
)

# main loop
for (group in analysis_groups) {
  cat(paste("\n\n\n========================================================\n"))
  cat(paste("     starting analysis for:", group$name, "\n"))
  cat(paste("========================================================\n"))
  
  for (pair in wave_pairs) {
    
    # --- [NEW] dynamically define the list of ordered vars for this pair ---
    base_vars <- paste0(c("peta", "petb", "petc", "petd", "pete"), group$reporter_suffix)
    wave1_vars <- paste0(base_vars, "_", pair[1])
    wave2_vars <- paste0(base_vars, "_", pair[2])
    current_ordered_vars <- c(wave1_vars, wave2_vars)
    
    # use a try-catch block to prevent one failed model from stopping the whole loop
    tryCatch({
      test_invariance_pair(
        data = group$data, 
        wave1_suffix = pair[1], 
        wave2_suffix = pair[2],
        reporter_suffix = group$reporter_suffix,
        sex_code = group$sex_code, # pass the sex_code to the function
        all_ordered_vars = current_ordered_vars, 
        max_free_intercepts = 2 # you can adjust this limit
      )
    }, error = function(e) {
      cat(paste("\nerror during analysis for", group$name, pair[1], "vs.", pair[2], ":", e$message, "\n"))
    })
  }
}


# get all column names related to puberty items (peta, petb, etc.)
puberty_cols <- names(wide_data)[grepl("^(peta|petb|petc|petd|pete)_(p|y)_(bl|fu[1-6])$", names(wide_data))]

long_data <- wide_data %>%
  select(id, sex, all_of(puberty_cols)) %>%
  pivot_longer(
    cols = all_of(puberty_cols),
    names_to = "item_wave",
    values_to = "score"
  ) %>%
  # separate item name, reporter, and wave
  extract(item_wave, into = c("item", "reporter", "wave"), regex = "^(.*)_(p|y)_(.*)$") %>%
  # ensure wave is ordered correctly for summarization
  mutate(wave = factor(wave, levels = waves, ordered = TRUE),
         score = as.numeric(as.character(score))) # convert score to numeric for stats

cat("--- data pivoted to long format. ---\n")

# -----------------------------------------------------------------------------
# 3. calculate descriptive statistics
# -----------------------------------------------------------------------------

cat("--- calculating descriptive statistics... ---\n")

descriptives <- long_data %>%
  group_by(sex, reporter, wave, item) %>%
  summarise(
    n = sum(!is.na(score)),
    n_unique = n_distinct(score, na.rm = TRUE),
    min = if(n > 0) min(score, na.rm = TRUE) else NA,
    max = if(n > 0) max(score, na.rm = TRUE) else NA,
    mean = mean(score, na.rm = TRUE),
    sd = sd(score, na.rm = TRUE),
    # calculate counts for each possible score value (e.g., 1 to 5)
    # assuming values are integers
    count_0 = sum(score == 0, na.rm = TRUE),
    count_1 = sum(score == 1, na.rm = TRUE),
    count_2 = sum(score == 2, na.rm = TRUE),
    count_3 = sum(score == 3, na.rm = TRUE),
    count_4 = sum(score == 4, na.rm = TRUE),
    .groups = 'drop' # drop grouping for printing
  ) %>%
  # make sex and reporter more readable
  mutate(sex = factor(sex, levels = c(0, 1), labels = c("male", "female")),
         reporter = factor(reporter, levels = c("p", "y"), labels = c("parent", "youth"))) %>%
  # arrange for easier reading
  arrange(sex, reporter, item, wave)

# -----------------------------------------------------------------------------
# 4. print results
# -----------------------------------------------------------------------------

cat("\n\n\n========================================================\n")
cat("  descriptive statistics for puberty items \n")
cat("========================================================\n")
cat(paste("generated on:", format(today(), "%y-%m-%d"), "\n\n"))

# print the full table
print(as.data.frame(descriptives), row.names = FALSE)

# highlight potential issues
zero_variance_issues <- descriptives %>% filter(n_unique <= 1 & n > 0)
if (nrow(zero_variance_issues) > 0) {
  cat("\n\n--- potential zero variance issues (n_unique <= 1) ---\n")
  print(as.data.frame(zero_variance_issues), row.names = FALSE)
} else {
  cat("\n\n--- no apparent zero variance issues found (n_unique > 1 for all items with data). ---\n")
}

cat("\n--- diagnostic script complete. ---\n")

# --- [fix 7] model syntax for females - parent report cfa only ---
model_females_parent_cfa <- paste0('
  # -------------------------------------------------
  # part 1a: parent-report measurement model (females)
  # -------------------------------------------------
  pub_p_bl  =~ 1*peta_p_bl + l2_p*petb_p_bl + l3_p*petc_p_bl + l4_p*petd_p_bl # + l5_p*pete_p_bl
  pub_p_fu1 =~ 1*peta_p_fu1 + l2_p*petb_p_fu1 + l3_p*petc_p_fu1 + l4_p*petd_p_fu1 # + l5_p*pete_p_fu1
  pub_p_fu2 =~ 1*peta_p_fu2 + l2_p*petb_p_fu2 + l3_p*petc_p_fu2 + l4_p*petd_p_fu2 # + l5_p*pete_p_fu2
  pub_p_fu3 =~ 1*peta_p_fu3 + l2_p*petb_p_fu3 + l3_p*petc_p_fu3 + l4_p*petd_p_fu3 # + l5_p*pete_p_fu3
  pub_p_fu4 =~ 1*peta_p_fu4 + l2_p*petb_p_fu4 + l3_p*petc_p_fu4 + l4_p*petd_p_fu4 # + l5_p*pete_p_fu4
  pub_p_fu5 =~ 1*peta_p_fu5 + l2_p*petb_p_fu5 + l3_p*petc_p_fu5 + l4_p*petd_p_fu5 # + l5_p*pete_p_fu5
  pub_p_fu6 =~ 1*peta_p_fu6 + l2_p*petb_p_fu6 + l3_p*petc_p_fu6 + l4_p*petd_p_fu6 # + l5_p*pete_p_fu6

  # allow factors to correlate
  pub_p_bl ~~ pub_p_fu1 + pub_p_fu2 + pub_p_fu3 + pub_p_fu4 + pub_p_fu5 + pub_p_fu6
  pub_p_fu1 ~~ pub_p_fu2 + pub_p_fu3 + pub_p_fu4 + pub_p_fu5 + pub_p_fu6
  pub_p_fu2 ~~ pub_p_fu3 + pub_p_fu4 + pub_p_fu5 + pub_p_fu6
  pub_p_fu3 ~~ pub_p_fu4 + pub_p_fu5 + pub_p_fu6
  pub_p_fu4 ~~ pub_p_fu5 + pub_p_fu6
  pub_p_fu5 ~~ pub_p_fu6
')

# --- [fix 7] model syntax for females - youth report cfa only ---
model_females_youth_cfa <- paste0('
  # -------------------------------------------------
  # part 1b: youth-report measurement model (females)
  # -------------------------------------------------
  pub_y_bl  =~ 1*peta_y_bl + l2_y*petb_y_bl + l3_y*petc_y_bl + l4_y*petd_y_bl + l5_y*pete_y_bl
  pub_y_fu1 =~ 1*peta_y_fu1 + l2_y*petb_y_fu1 + l3_y*petc_y_fu1 + l4_y*petd_y_fu1 + l5_y*pete_y_fu1
  pub_y_fu2 =~ 1*peta_y_fu2 + l2_y*petb_y_fu2 + l3_y*petc_y_fu2 + l4_y*petd_y_fu2 + l5_y*pete_y_fu2
  pub_y_fu3 =~ 1*peta_y_fu3 + l2_y*petb_y_fu3 + l3_y*petc_y_fu3 + l4_y*petd_y_fu3 + l5_y*pete_y_fu3
  pub_y_fu4 =~ 1*peta_y_fu4 + l2_y*petb_y_fu4 + l3_y*petc_y_fu4 + l4_y*petd_y_fu4 + l5_y*pete_y_fu4
  pub_y_fu5 =~ 1*peta_y_fu5 + l2_y*petb_y_fu5 + l3_y*petc_y_fu5 + l4_y*petd_y_fu5 + l5_y*pete_y_fu5
  pub_y_fu6 =~ 1*peta_y_fu6 + l2_y*petb_y_fu6 + l3_y*petc_y_fu6 + l4_y*petd_y_fu6 + l5_y*pete_y_fu6

  # allow factors to correlate
  pub_y_bl ~~ pub_y_fu1 + pub_y_fu2 + pub_y_fu3 + pub_y_fu4 + pub_y_fu5 + pub_y_fu6
  pub_y_fu1 ~~ pub_y_fu2 + pub_y_fu3 + pub_y_fu4 + pub_y_fu5 + pub_y_fu6
  pub_y_fu2 ~~ pub_y_fu3 + pub_y_fu4 + pub_y_fu5 + pub_y_fu6
  pub_y_fu3 ~~ pub_y_fu4 + pub_y_fu5 + pub_y_fu6
  pub_y_fu4 ~~ pub_y_fu5 + pub_y_fu6
  pub_y_fu5 ~~ pub_y_fu6
')

# --- [fix 7] model syntax for males - parent report cfa only ---
model_males_parent_cfa <- paste0('
  # -------------------------------------------------
  # part 1a: parent-report measurement model (males)
  # -------------------------------------------------
  pub_p_bl  =~ 1*peta_p_bl + l2_p*petb_p_bl + l3_p*petc_p_bl + l4_p*petd_p_bl + l5_p*pete_p_bl
  pub_p_fu1 =~ 1*peta_p_fu1 + l2_p*petb_p_fu1 + l3_p*petc_p_fu1 + l4_p*petd_p_fu1 + l5_p*pete_p_fu1
  pub_p_fu2 =~ 1*peta_p_fu2 + l2_p*petb_p_fu2 + l3_p*petc_p_fu2 + l4_p*petd_p_fu2 + l5_p*pete_p_fu2
  pub_p_fu3 =~ 1*peta_p_fu3 + l2_p*petb_p_fu3 + l3_p*petc_p_fu3 + l4_p*petd_p_fu3 + l5_p*pete_p_fu3
  pub_p_fu4 =~ 1*peta_p_fu4 + l2_p*petb_p_fu4 + l3_p*petc_p_fu4 + l4_p*petd_p_fu4 + l5_p*pete_p_fu4
  pub_p_fu5 =~ 1*peta_p_fu5 + l2_p*petb_p_fu5 + l3_p*petc_p_fu5 + l4_p*petd_p_fu5 + l5_p*pete_p_fu5
  pub_p_fu6 =~ 1*peta_p_fu6 + l2_p*petb_p_fu6 + l3_p*petc_p_fu6 + l4_p*petd_p_fu6 + l5_p*pete_p_fu6

  # allow factors to correlate
  pub_p_bl ~~ pub_p_fu1 + pub_p_fu2 + pub_p_fu3 + pub_p_fu4 + pub_p_fu5 + pub_p_fu6
  pub_p_fu1 ~~ pub_p_fu2 + pub_p_fu3 + pub_p_fu4 + pub_p_fu5 + pub_p_fu6
  pub_p_fu2 ~~ pub_p_fu3 + pub_p_fu4 + pub_p_fu5 + pub_p_fu6
  pub_p_fu3 ~~ pub_p_fu4 + pub_p_fu5 + pub_p_fu6
  pub_p_fu4 ~~ pub_p_fu5 + pub_p_fu6
  pub_p_fu5 ~~ pub_p_fu6
')

# --- [fix 7] model syntax for males - youth report cfa only ---
model_males_youth_cfa <- paste0('
  # -------------------------------------------------
  # part 1b: youth-report measurement model (males)
  # -------------------------------------------------
  pub_y_bl  =~ 1*peta_y_bl + l2_y*petb_y_bl + l3_y*petc_y_bl + l4_y*petd_y_bl + l5_y*pete_y_bl
  pub_y_fu1 =~ 1*peta_y_fu1 + l2_y*petb_y_fu1 + l3_y*petc_y_fu1 + l4_y*petd_y_fu1 + l5_y*pete_y_fu1
  pub_y_fu2 =~ 1*peta_y_fu2 + l2_y*petb_y_fu2 + l3_y*petc_y_fu2 + l4_y*petd_y_fu2 + l5_y*pete_y_fu2
  pub_y_fu3 =~ 1*peta_y_fu3 + l2_y*petb_y_fu3 + l3_y*petc_y_fu3 + l4_y*petd_y_fu3 + l5_y*pete_y_fu3
  pub_y_fu4 =~ 1*peta_y_fu4 + l2_y*petb_y_fu4 + l3_y*petc_y_fu4 + l4_y*petd_y_fu4 + l5_y*pete_y_fu4
  pub_y_fu5 =~ 1*peta_y_fu5 + l2_y*petb_y_fu5 + l3_y*petc_y_fu5 + l4_y*petd_y_fu5 + l5_y*pete_y_fu5
  pub_y_fu6 =~ 1*peta_y_fu6 + l2_y*petb_y_fu6 + l3_y*petc_y_fu6 + l4_y*petd_y_fu6 + l5_y*pete_y_fu6

  # allow factors to correlate
  pub_y_bl ~~ pub_y_fu1 + pub_y_fu2 + pub_y_fu3 + pub_y_fu4 + pub_y_fu5 + pub_y_fu6
  pub_y_fu1 ~~ pub_y_fu2 + pub_y_fu3 + pub_y_fu4 + pub_y_fu5 + pub_y_fu6
  pub_y_fu2 ~~ pub_y_fu3 + pub_y_fu4 + pub_y_fu5 + pub_y_fu6
  pub_y_fu3 ~~ pub_y_fu4 + pub_y_fu5 + pub_y_fu6
  pub_y_fu4 ~~ pub_y_fu5 + pub_y_fu6
  pub_y_fu5 ~~ pub_y_fu6
')


# -----------------------------------------------------------------------------
# step 3: fit the separate longitudinal cfa models (diagnostic)
# -----------------------------------------------------------------------------

# --- fit female models ---
cat("\n--- fitting longitudinal cfa for females - parent report ---\n")
fit_cfa_females_parent <- cfa(model_females_parent_cfa, # [fix 7] use cfa()
                              data = females_wide,
                              estimator = "wlsmv",
                              missing = "pairwise",
                              ordered = ordered_vars_p) # use parent-only list

cat("\n\n\n========================================================\n")
cat("  longitudinal cfa summary (females - parent) \n")
cat("========================================================\n")
print(summary(fit_cfa_females_parent, fit.measures = TRUE, standardized = TRUE))

cat("\n--- fitting longitudinal cfa for females - youth report ---\n")
fit_cfa_females_youth <- cfa(model_females_youth_cfa, # [fix 7] use cfa()
                             data = females_wide,
                             estimator = "wlsmv",
                             missing = "pairwise",
                             ordered = ordered_vars_y) # use youth-only list

cat("\n\n\n========================================================\n")
cat("  longitudinal cfa summary (females - youth) \n")
cat("========================================================\n")
print(summary(fit_cfa_females_youth, fit.measures = TRUE, standardized = TRUE))


# --- fit male models ---
cat("\n--- fitting longitudinal cfa for males - parent report ---\n")
fit_cfa_males_parent <- cfa(model_males_parent_cfa, # [fix 7] use cfa()
                            data = males_wide,
                            estimator = "wlsmv",
                            missing = "pairwise",
                            ordered = ordered_vars_p) # use parent-only list

cat("\n\n\n========================================================\n")
cat("  longitudinal cfa summary (males - parent) \n")
cat("========================================================\n")
print(summary(fit_cfa_males_parent, fit.measures = TRUE, standardized = TRUE))


cat("\n--- fitting longitudinal cfa for males - youth report ---\n")
fit_cfa_males_youth <- cfa(model_males_youth_cfa, # [fix 7] use cfa()
                           data = males_wide,
                           estimator = "wlsmv",
                           missing = "pairwise",
                           ordered = ordered_vars_y) # use youth-only list

cat("\n\n\n========================================================\n")
cat("  longitudinal cfa summary (males - youth) \n")
cat("========================================================\n")
print(summary(fit_cfa_males_youth, fit.measures = TRUE, standardized = TRUE))


cat("\n--- diagnostic script complete. ---\n")

# -----------------------------------------------------------------------------
# 5. OVERALL CONCLUSION
# -----------------------------------------------------------------------------
# review the complete output above for the specific results of each group and
# wave-pair comparison. this automated analysis provides a comprehensive test
# of both metric (weak) and scalar (strong) invariance.

# final interpretation guide:
# for each group and wave transition, the output will now tell you:
# 1. if full or partial metric invariance was achieved.
# 2. if full, partial (and with which items freed), or no scalar invariance
#    was achieved.

# -----------------------------------------------------------------------------
# (your single-wave cfa functions would go here, unchanged)
# -----------------------------------------------------------------------------

c# -----------------------------------------------------------------------------
# step 2: one-step nonlinear latent growth curve model (gold standard)
#
# this script implements the "one-step" approach detailed in
# 'analysis_plan_v2.md'. this is possible because the invariance
# tests (using wlsmv) confirmed full metric and scalar invariance.
#
# this model estimates the measurement (cfa) and structural (lgcm)
# components simultaneously, using the raw ordinal items.
#
# [fix] this version fits two separate models (one for females, one for males)
# as requested, rather than a single multi-group model.
#
# [fix 2] this version also removes the quadratic growth factor, as the
# full quadratic model was not identified (produced a "not positive definite"
# warning). we are now fitting a linear-only model.
#
# [fix 3] simplify growth structure: fit two correlated first-order
# linear growth models instead of a second-order model to aid identification.
#
# [fix 4] further simplify: fit completely separate first-order models
# for parent and youth reports to resolve persistent identification issues.
# -----------------------------------------------------------------------------
ordered_vars_p <- c(
  "peta_p_bl", "petb_p_bl", "petc_p_bl", "petd_p_bl",     # "pete_p_bl", 
  "peta_p_fu1", "petb_p_fu1", "petc_p_fu1", "petd_p_fu1", # "pete_p_fu1", 
  "peta_p_fu2", "petb_p_fu2", "petc_p_fu2", "petd_p_fu2", # "pete_p_fu2", 
  "peta_p_fu3", "petb_p_fu3", "petc_p_fu3", "petd_p_fu3", # "pete_p_fu3", 
  "peta_p_fu4", "petb_p_fu4", "petc_p_fu4", "petd_p_fu4", # "pete_p_fu4", 
  "peta_p_fu5", "petb_p_fu5", "petc_p_fu5", "petd_p_fu5", # "pete_p_fu5", 
  "peta_p_fu6", "petb_p_fu6", "petc_p_fu6", "petd_p_fu6" # ,"pete_p_fu6"
)
ordered_vars_y <- c(
  "peta_y_bl", "petb_y_bl", "petc_y_bl", "petd_y_bl", "pete_y_bl", 
  "peta_y_fu1", "petb_y_fu1", "petc_y_fu1", "petd_y_fu1", "pete_y_fu1", 
  "peta_y_fu2", "petb_y_fu2", "petc_y_fu2", "petd_y_fu2", "pete_y_fu2", 
  "peta_y_fu3", "petb_y_fu3", "petc_y_fu3", "petd_y_fu3", "pete_y_fu3", 
  "peta_y_fu4", "petb_y_fu4", "petc_y_fu4", "petd_y_fu4", "pete_y_fu4", 
  "peta_y_fu5", "petb_y_fu5", "petc_y_fu5", "petd_y_fu5", "pete_y_fu5", 
  "peta_y_fu6", "petb_y_fu6", "petc_y_fu6", "petd_y_fu6", "pete_y_fu6"
)


# define time scores for linear growth
time_scores_linear <- c(0, 1, 2, 3, 4, 5, 6)

# --- [fix 4] model syntax for females - parent report only ---
model_females_parent_lgcm <- paste0('
  # -------------------------------------------------
  # part 1a: parent-report measurement model (females)
  # -------------------------------------------------
  pub_p_bl  =~ 1*peta_p_bl + l2_p*petb_p_bl + l3_p*petc_p_bl + l4_p*petd_p_bl # + l5_p*pete_p_bl
  pub_p_fu1 =~ 1*peta_p_fu1 + l2_p*petb_p_fu1 + l3_p*petc_p_fu1 + l4_p*petd_p_fu1 # + l5_p*pete_p_fu1
  pub_p_fu2 =~ 1*peta_p_fu2 + l2_p*petb_p_fu2 + l3_p*petc_p_fu2 + l4_p*petd_p_fu2 # + l5_p*pete_p_fu2
  pub_p_fu3 =~ 1*peta_p_fu3 + l2_p*petb_p_fu3 + l3_p*petc_p_fu3 + l4_p*petd_p_fu3 # + l5_p*pete_p_fu3
  pub_p_fu4 =~ 1*peta_p_fu4 + l2_p*petb_p_fu4 + l3_p*petc_p_fu4 + l4_p*petd_p_fu4 # + l5_p*pete_p_fu4
  pub_p_fu5 =~ 1*peta_p_fu5 + l2_p*petb_p_fu5 + l3_p*petc_p_fu5 + l4_p*petd_p_fu5 # + l5_p*pete_p_fu5
  pub_p_fu6 =~ 1*peta_p_fu6 + l2_p*petb_p_fu6 + l3_p*petc_p_fu6 + l4_p*petd_p_fu6 # + l5_p*pete_p_fu6

  # -------------------------------------------------
  # part 2: first-order growth curves (parent only)
  # -------------------------------------------------
  i_parent =~ 1*pub_p_bl + 1*pub_p_fu1 + 1*pub_p_fu2 + 1*pub_p_fu3 + 1*pub_p_fu4 + 1*pub_p_fu5 + 1*pub_p_fu6
  s_parent =~ ', paste0(time_scores_linear, "*pub_p_", waves, collapse = " + "), '
  
  # estimate variances and covariance of growth factors
  i_parent ~~ i_parent + s_parent
  s_parent ~~ s_parent
  
  # fix residual variances of time-specific factors to 0
  pub_p_bl ~~ 0*pub_p_bl; pub_p_fu1 ~~ 0*pub_p_fu1; pub_p_fu2 ~~ 0*pub_p_fu2; 
  pub_p_fu3 ~~ 0*pub_p_fu3; pub_p_fu4 ~~ 0*pub_p_fu4; pub_p_fu5 ~~ 0*pub_p_fu5; pub_p_fu6 ~~ 0*pub_p_fu6;
')

# --- [fix 4] model syntax for females - youth report only ---
model_females_youth_lgcm <- paste0('
  # -------------------------------------------------
  # part 1b: youth-report measurement model (females)
  # -------------------------------------------------
  pub_y_bl  =~ 1*peta_y_bl + l2_y*petb_y_bl + l3_y*petc_y_bl + l4_y*petd_y_bl + l5_y*pete_y_bl
  pub_y_fu1 =~ 1*peta_y_fu1 + l2_y*petb_y_fu1 + l3_y*petc_y_fu1 + l4_y*petd_y_fu1 + l5_y*pete_y_fu1
  pub_y_fu2 =~ 1*peta_y_fu2 + l2_y*petb_y_fu2 + l3_y*petc_y_fu2 + l4_y*petd_y_fu2 + l5_y*pete_y_fu2
  pub_y_fu3 =~ 1*peta_y_fu3 + l2_y*petb_y_fu3 + l3_y*petc_y_fu3 + l4_y*petd_y_fu3 + l5_y*pete_y_fu3
  pub_y_fu4 =~ 1*peta_y_fu4 + l2_y*petb_y_fu4 + l3_y*petc_y_fu4 + l4_y*petd_y_fu4 + l5_y*pete_y_fu4
  pub_y_fu5 =~ 1*peta_y_fu5 + l2_y*petb_y_fu5 + l3_y*petc_y_fu5 + l4_y*petd_y_fu5 + l5_y*pete_y_fu5
  pub_y_fu6 =~ 1*peta_y_fu6 + l2_y*petb_y_fu6 + l3_y*petc_y_fu6 + l4_y*petd_y_fu6 + l5_y*pete_y_fu6
  
  # -------------------------------------------------
  # part 2: first-order growth curves (youth only)
  # -------------------------------------------------
  i_youth =~ 1*pub_y_bl + 1*pub_y_fu1 + 1*pub_y_fu2 + 1*pub_y_fu3 + 1*pub_y_fu4 + 1*pub_y_fu5 + 1*pub_y_fu6
  s_youth =~ ', paste0(time_scores_linear, "*pub_y_", waves, collapse = " + "), '

  # estimate variances and covariance of growth factors
  i_youth ~~ i_youth + s_youth
  s_youth ~~ s_youth

  # fix residual variances of time-specific factors to 0
  pub_y_bl ~~ 0*pub_y_bl; pub_y_fu1 ~~ 0*pub_y_fu1; pub_y_fu2 ~~ 0*pub_y_fu2;
  pub_y_fu3 ~~ 0*pub_y_fu3; pub_y_fu4 ~~ 0*pub_y_fu4; pub_y_fu5 ~~ 0*pub_y_fu5; pub_y_fu6 ~~ 0*pub_y_fu6;
')

# --- [fix 4] model syntax for males - parent report only ---
model_males_parent_lgcm <- paste0('
  # -------------------------------------------------
  # part 1a: parent-report measurement model (males)
  # -------------------------------------------------
  pub_p_bl  =~ 1*peta_p_bl + l2_p*petb_p_bl + l3_p*petc_p_bl + l4_p*petd_p_bl + l5_p*pete_p_bl
  pub_p_fu1 =~ 1*peta_p_fu1 + l2_p*petb_p_fu1 + l3_p*petc_p_fu1 + l4_p*petd_p_fu1 + l5_p*pete_p_fu1
  pub_p_fu2 =~ 1*peta_p_fu2 + l2_p*petb_p_fu2 + l3_p*petc_p_fu2 + l4_p*petd_p_fu2 + l5_p*pete_p_fu2
  pub_p_fu3 =~ 1*peta_p_fu3 + l2_p*petb_p_fu3 + l3_p*petc_p_fu3 + l4_p*petd_p_fu3 + l5_p*pete_p_fu3
  pub_p_fu4 =~ 1*peta_p_fu4 + l2_p*petb_p_fu4 + l3_p*petc_p_fu4 + l4_p*petd_p_fu4 + l5_p*pete_p_fu4
  pub_p_fu5 =~ 1*peta_p_fu5 + l2_p*petb_p_fu5 + l3_p*petc_p_fu5 + l4_p*petd_p_fu5 + l5_p*pete_p_fu5
  pub_p_fu6 =~ 1*peta_p_fu6 + l2_p*petb_p_fu6 + l3_p*petc_p_fu6 + l4_p*petd_p_fu6 + l5_p*pete_p_fu6

  # -------------------------------------------------
  # part 2: first-order growth curves (parent only)
  # -------------------------------------------------
  i_parent =~ 1*pub_p_bl + 1*pub_p_fu1 + 1*pub_p_fu2 + 1*pub_p_fu3 + 1*pub_p_fu4 + 1*pub_p_fu5 + 1*pub_p_fu6
  s_parent =~ ', paste0(time_scores_linear, "*pub_p_", waves, collapse = " + "), '
  
  # estimate variances and covariance of growth factors
  i_parent ~~ i_parent + s_parent
  s_parent ~~ s_parent
  
  # fix residual variances of time-specific factors to 0
  pub_p_bl ~~ 0*pub_p_bl; pub_p_fu1 ~~ 0*pub_p_fu1; pub_p_fu2 ~~ 0*pub_p_fu2; 
  pub_p_fu3 ~~ 0*pub_p_fu3; pub_p_fu4 ~~ 0*pub_p_fu4; pub_p_fu5 ~~ 0*pub_p_fu5; pub_p_fu6 ~~ 0*pub_p_fu6;
')

# --- [fix 4] model syntax for males - youth report only ---
model_males_youth_lgcm <- paste0('
  # -------------------------------------------------
  # part 1b: youth-report measurement model (males)
  # -------------------------------------------------
  pub_y_bl  =~ 1*peta_y_bl + l2_y*petb_y_bl + l3_y*petc_y_bl + l4_y*petd_y_bl + l5_y*pete_y_bl
  pub_y_fu1 =~ 1*peta_y_fu1 + l2_y*petb_y_fu1 + l3_y*petc_y_fu1 + l4_y*petd_y_fu1 + l5_y*pete_y_fu1
  pub_y_fu2 =~ 1*peta_y_fu2 + l2_y*petb_y_fu2 + l3_y*petc_y_fu2 + l4_y*petd_y_fu2 + l5_y*pete_y_fu2
  pub_y_fu3 =~ 1*peta_y_fu3 + l2_y*petb_y_fu3 + l3_y*petc_y_fu3 + l4_y*petd_y_fu3 + l5_y*pete_y_fu3
  pub_y_fu4 =~ 1*peta_y_fu4 + l2_y*petb_y_fu4 + l3_y*petc_y_fu4 + l4_y*petd_y_fu4 + l5_y*pete_y_fu4
  pub_y_fu5 =~ 1*peta_y_fu5 + l2_y*petb_y_fu5 + l3_y*petc_y_fu5 + l4_y*petd_y_fu5 + l5_y*pete_y_fu5
  pub_y_fu6 =~ 1*peta_y_fu6 + l2_y*petb_y_fu6 + l3_y*petc_y_fu6 + l4_y*petd_y_fu6 + l5_y*pete_y_fu6
  
  # -------------------------------------------------
  # part 2: first-order growth curves (youth only)
  # -------------------------------------------------
  i_youth =~ 1*pub_y_bl + 1*pub_y_fu1 + 1*pub_y_fu2 + 1*pub_y_fu3 + 1*pub_y_fu4 + 1*pub_y_fu5 + 1*pub_y_fu6
  s_youth =~ ', paste0(time_scores_linear, "*pub_y_", waves, collapse = " + "), '

  # estimate variances and covariance of growth factors
  i_youth ~~ i_youth + s_youth
  s_youth ~~ s_youth

  # fix residual variances of time-specific factors to 0
  pub_y_bl ~~ 0*pub_y_bl; pub_y_fu1 ~~ 0*pub_y_fu1; pub_y_fu2 ~~ 0*pub_y_fu2;
  pub_y_fu3 ~~ 0*pub_y_fu3; pub_y_fu4 ~~ 0*pub_y_fu4; pub_y_fu5 ~~ 0*pub_y_fu5; pub_y_fu6 ~~ 0*pub_y_fu6;
')


# -----------------------------------------------------------------------------
# step 3: fit the separate one-step models
# -----------------------------------------------------------------------------

# --- fit female models ---
cat("\n--- fitting final one-step lgcm for females - parent report (this may take time) ---\n")
fit_lgcm_females_parent <- growth(model_females_parent_lgcm, 
                                  data = females_wide, 
                                  estimator = "wlsmv", 
                                  missing = "pairwise",
                                  ordered = ordered_vars_p,
                                  se = "robust") # use parent-only list

cat("\n\n\n========================================================\n")
cat("  final linear growth model summary (females - parent) \n")
cat("========================================================\n")
print(summary(fit_lgcm_females_parent, fit.measures = TRUE, standardized = TRUE))

cat("\n--- fitting final one-step lgcm for females - youth report (this may take time) ---\n")
fit_lgcm_females_youth <- growth(model_females_youth_lgcm, 
                                 data = females_wide, 
                                 estimator = "wlsmv", 
                                 missing = "pairwise",
                                 ordered = ordered_vars_females_y) # use youth-only list

cat("\n\n\n========================================================\n")
cat("  final linear growth model summary (females - youth) \n")
cat("========================================================\n")
print(summary(fit_lgcm_females_youth, fit.measures = TRUE, standardized = TRUE))


# --- fit male models ---
cat("\n--- fitting final one-step lgcm for males - parent report (this may take time) ---\n")
fit_lgcm_males_parent <- growth(model_males_parent_lgcm, 
                                data = males_wide, 
                                estimator = "wlsmv", 
                                missing = "pairwise",
                                ordered = ordered_vars_males_p) # use parent-only list

cat("\n\n\n========================================================\n")
cat("  final linear growth model summary (males - parent) \n")
cat("========================================================\n")
print(summary(fit_lgcm_males_parent, fit.measures = TRUE, standardized = TRUE))


cat("\n--- fitting final one-step lgcm for males - youth report (this may take time) ---\n")
fit_lgcm_males_youth <- growth(model_males_youth_lgcm, 
                               data = males_wide, 
                               estimator = "wlsmv", 
                               missing = "pairwise",
                               ordered = ordered_vars_males_y) # use youth-only list

cat("\n\n\n========================================================\n")
cat("  final linear growth model summary (males - youth) \n")
cat("========================================================\n")
print(summary(fit_lgcm_males_youth, fit.measures = TRUE, standardized = TRUE))


cat("\n--- script complete. ---\n")

# next steps: 
# - extract factor scores (i_parent, s_parent, i_youth, s_youth) using lavPredict() 
#   from each of the four fitted models.
# - merge these scores into a single dataframe.
# - analyze correlations between parent and youth growth factors.
# - use the scores in further analyses (e.g., predicting outcomes).

long_data_f <- long_data %>% 
  filter(sex == 1) %>% 
  filter(!is.na(score)) %>%
  mutate(item_rep = paste0(item, "_", reporter)) %>%
  pivot_wider(
    id_cols    = c(id, wave),
    names_from = item_rep,
    values_from = score) 

items4  <- c("peta_p","petb_p","petc_p","petd_p")
item2   <- "pete_p"

df4 <- long_data_f %>%
  mutate(across(all_of(items4), ~ ordered(as.integer(.), levels = 1:4)))

df2 <- long_data_f %>%
  mutate(across(all_of(item2),  ~ ordered(as.integer(.), levels = 0:1)))

#%>% 
  mutate(
    peta_p = factor(as.numeric(as.character(peta_p)), levels = 1:4, ordered = TRUE),
    petb_p = factor(as.numeric(as.character(petb_p)), levels = 1:4, ordered = TRUE),
    petc_p = factor(as.numeric(as.character(petc_p)), levels = 1:4, ordered = TRUE),
    petd_p = factor(as.numeric(as.character(petd_p)), levels = 1:4, ordered = TRUE),
    pete_p = factor(as.numeric(as.character(pete_p)), levels = 0:1, ordered = TRUE)  # binary, keep it binary
  ) %>% 
  mutate(wave = factor(wave, levels = c("bl","fu1","fu2","fu3","fu4","fu5","fu6")))


str(long_data_f)

items <- c("peta_p","petb_p","petc_p","petd_p")


long_brms <- long_data_f %>%
  pivot_longer(all_of(items), names_to = "item", values_to = "resp") %>%
  # make sure responses are ordered 1..4 with no weird labels
  mutate(
    item = factor(item),
    resp = ordered(resp, levels = sort(unique(unlist(long_data_f[items]))))
  ) %>%
  drop_na(resp, age, wave)

# graded response model with smooth moderation of loadings and thresholds
# - latent factor is the person random effect
# - item random intercepts handle baseline difficulties
# - disc varies by item with its own smooth over age, plus wave shifts
f <- bf(
  resp ~ 1 + cs(s(age)) + cs(wave) + (1 | i | item) + (1 | p | id),
  disc ~ 0 + item + s(age, by = item) + wave,
  nl = TRUE
)

fit_mnlfa <- brm(
  formula = f,
  data = long_brms,
  family = brmsfamily("grm", link = "logit", threshold = "flexible"),
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed = 123
)

print(fit_mnlfa)









base_cfa <- '
  f_p =~ peta_p + petb_p + petc_p + petd_p 
'

fit_conf <- cfa(base_cfa, data=long_data_f, group="wave",
                estimator="wlsmv", missing = "pairwise", ordered = c("peta_y","petb_y","petc_y","petd_y","pete_y", "peta_p","petb_p","petc_p","petd_p"))
fit_met  <- cfa(base_cfa, data=long_data_f, group="wave",
                estimator="wlsmv", missing = "pairwise", ordered= c("peta_y","petb_y","petc_y","petd_y","pete_y", "peta_p","petb_p","petc_p","petd_p"),
                group.equal="loadings")
fit_th   <- cfa(base_cfa, data=long_data_f, group="wave",
                estimator="wlsmv", missing = "pairwise", ordered= c("peta_y","petb_y","petc_y","petd_y","pete_y", "peta_p","petb_p","petc_p","petd_p"),
                group.equal=c("loadings","thresholds"))

anova(fit_conf, fit_met)  # metric test (WLSMV uses DIFFTEST under the hood)
anova(fit_met, fit_th)  

summary(fit_met, fit.measures = TRUE)

score <- lavTestScore(fit_met)  # or fit_th if you're chasing thresholds
pt    <- parameterTable(fit_met)  # same model as score test


lookup <- pt %>%
  dplyr::mutate(
    pretty = dplyr::case_when(
      op == "=~" ~ paste0("[g", group, "] ", lhs, " =~ ", rhs),         # loading
      op == "|"  ~ paste0("[g", group, "] ", lhs, " | ", rhs),          # threshold (rhs like t1, t2)
      op == "~1" ~ paste0("[g", group, "] ", lhs, " ~ 1"),              # intercept (for continuous)
      op == "~~" ~ paste0("[g", group, "] ", lhs, " ~~ ", rhs),         # variance/cov
      TRUE       ~ paste0("[g", group, "] ", lhs, " ", op, " ", rhs)
    )
  ) %>%
  dplyr::select(plabel, pretty) %>%
  dplyr::distinct()

readable <- score$uni %>%
  dplyr::mutate(lhs_pl = lhs, rhs_pl = rhs) %>%
  dplyr::left_join(lookup, by = c("lhs_pl" = "plabel")) %>%
  dplyr::rename(lhs_pretty = pretty) %>%
  dplyr::left_join(lookup, by = c("rhs_pl" = "plabel")) %>%
  dplyr::rename(rhs_pretty = pretty) %>%
  dplyr::arrange(dplyr::desc(X2))

readable %>%
  dplyr::select(X2, df, p.value, lhs_pretty, rhs_pretty) %>%
  head(20)

fit_met_partial  <- cfa(base_cfa, data=long_data_f, group="wave",
                estimator="wlsmv", missing = "pairwise", ordered= c("peta_y","petb_y","petc_y","petd_y","pete_y", "peta_p","petb_p","petc_p","petd_p"),
                group.equal="loadings",
                group.partial = c("f=~peta_p", "f=~petc_y", "f=~petb_y", "f=~petc_p", "f=~pete_y", "f=~petd_y", "f=~petd_p"))

anova(fit_conf, fit_met_partial)

score <- lavTestScore(fit_met_partial)  # or fit_th if you're chasing thresholds
pt    <- parameterTable(fit_met_partial)  # same model as score test


lookup <- pt %>%
  dplyr::mutate(
    pretty = dplyr::case_when(
      op == "=~" ~ paste0("[g", group, "] ", lhs, " =~ ", rhs),         # loading
      op == "|"  ~ paste0("[g", group, "] ", lhs, " | ", rhs),          # threshold (rhs like t1, t2)
      op == "~1" ~ paste0("[g", group, "] ", lhs, " ~ 1"),              # intercept (for continuous)
      op == "~~" ~ paste0("[g", group, "] ", lhs, " ~~ ", rhs),         # variance/cov
      TRUE       ~ paste0("[g", group, "] ", lhs, " ", op, " ", rhs)
    )
  ) %>%
  dplyr::select(plabel, pretty) %>%
  dplyr::distinct()

readable <- score$uni %>%
  dplyr::mutate(lhs_pl = lhs, rhs_pl = rhs) %>%
  dplyr::left_join(lookup, by = c("lhs_pl" = "plabel")) %>%
  dplyr::rename(lhs_pretty = pretty) %>%
  dplyr::left_join(lookup, by = c("rhs_pl" = "plabel")) %>%
  dplyr::rename(rhs_pretty = pretty) %>%
  dplyr::arrange(dplyr::desc(X2))

readable %>%
  dplyr::select(X2, df, p.value, lhs_pretty, rhs_pretty) %>%
  head(20)
