## Cleaned Pipeline for Longitudinal Measurement Invariance Analysis
# This script implements the streamlined strategy of testing for partial metric
# and partial scalar invariance between adjacent waves for multiple groups.

# -----------------------------------------------------------------------------
# 1. SETUP
# -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(readr)
library(lavaan)
library(semTools)

set.seed(90025)

# Define file paths (please update if necessary)
pub_data_root <- file.path(
  "/Users/clarefmccann/University of Oregon Dropbox",
  "Clare McCann/mine/projects/abcd-projs/abcd-data-release-6.0/physical-health/puberty"
)

# -----------------------------------------------------------------------------
# 2. DATA LOADING AND PREPARATION
# -----------------------------------------------------------------------------

# Helper function to read, rename, and select columns
load_and_prep <- function(file_path, sex_code, reporter) {
  
  if (reporter == "parent") {
    name_map <- c(pete_p = "fpete", peta_p = "peta", petb_p = "petb", petc_p = "petc", petd_p = "petd", pdss_p = "PDSS")
    if (sex_code == 0) name_map["pete_p"] <- "mpete" # Use mpete for males
  } else {
    name_map <- c(pete_y = "fpete", peta_y = "peta", petb_y = "petb", petc_y = "petc", petd_y = "petd", pdss_y = "PDSS")
    if (sex_code == 0) name_map["pete_y"] <- "mpete"
  }
  
  read_csv(file_path, show_col_types = FALSE) %>%
    rename(any_of(name_map)) %>%
    select(id, wave, age, starts_with("pet"), starts_with("pdss")) %>%
    mutate(sex = sex_code)
}

# Load all four datasets
pub_f_p <- load_and_prep(file.path(pub_data_root, "filtered_parent_tannerstages_f.csv"), 1, "parent")
pub_m_p <- load_and_prep(file.path(pub_data_root, "filtered_parent_tannerstages_m.csv"), 0, "parent")
pub_f_y <- load_and_prep(file.path(pub_data_root, "filtered_youth_tannerstages_f.csv"), 1, "youth")
pub_m_y <- load_and_prep(file.path(pub_data_root, "filtered_youth_tannerstages_m.csv"), 0, "youth")

# Define wave mapping
map_wave <- c("ses-00A" = "bl", "ses-01A" = "fu1", "ses-02A" = "fu2", 
              "ses-03A" = "fu3", "ses-04A" = "fu4", "ses-05A" = "fu5", "ses-06A" = "fu6")

# Combine, join, and pivot to wide format in one pipeline
wide_data <- bind_rows(pub_f_p, pub_m_p) %>%
  full_join(bind_rows(pub_f_y, pub_m_y), by = c("id", "wave", "sex", "age")) %>%
  mutate(wave = recode(wave, !!!map_wave)) %>%
  pivot_wider(
    id_cols = c(id, sex),
    names_from = wave,
    values_from = c(starts_with(c("age", "pet")), starts_with("pdss")),
    names_glue = "{.value}_{wave}"
  ) %>% 
  select(-contains(c("m.x","f.x", "m.y", "f.y")))
  

# Create data subsets for males and females, converting all puberty items to numeric
females_wide <- wide_data %>%
  filter(sex == 1) %>%
  mutate(across(starts_with("peta") | starts_with("petb") | starts_with("petc") | 
                starts_with("petd") | starts_with("pete"), as.numeric))
males_wide <- wide_data %>%
  filter(sex == 0) %>%
  mutate(across(starts_with("peta") | starts_with("petb") | starts_with("petc") | 
                starts_with("petd") | starts_with("pete"), as.numeric))


data_to_analyze <- females_wide # Use females_wide, change if needed


# -----------------------------------------------------------------------------
# 1. DEFINE ITEM SETS AND ANALYSIS PARAMETERS
# -----------------------------------------------------------------------------

# Base item names (excluding menarche for scale analysis)
base_items_p <- c("peta_p", "petb_p", "petc_p", "petd_p")
base_items_y <- c("peta_y", "petb_y", "petc_y", "petd_y")

# Timepoints
timepoints <- c("bl", "fu1", "fu2", "fu3", "fu4", "fu5", "fu6")

# List to store results (optional, but good practice)
all_results <- list()

# -----------------------------------------------------------------------------
# 2. MAIN ANALYSIS LOOP
# -----------------------------------------------------------------------------

for (reporter in c("p", "y")) {
  
  base_items <- if (reporter == "p") base_items_p else base_items_y
  reporter_suffix <- paste0("_", reporter)
  
  for (tp in timepoints) {
    
    current_wave_id <- paste0(reporter_suffix, "_", tp)
    all_results[[current_wave_id]] <- list()
    
    cat(paste("\n\n\n========================================================"))
    cat(paste("\n  ANALYZING REPORTER:", reporter, " | TIMEPOINT:", tp))
    cat(paste("\n========================================================\n"))
    
    # --- Create list of current items ---
    current_items <- paste0(base_items, "_", tp)
    
    # --- Subset data for current items, ensuring they exist ---
    if (!all(current_items %in% names(data_to_analyze))) {
      cat("!!! Skipping: Not all items found in data for", current_wave_id, "\n")
      next
    }
    # Ensure the subset data is numeric for functions needing it
    current_data_numeric <- data_to_analyze[, current_items, drop = FALSE]
    # Explicitly convert potentially non-numeric columns to numeric
    current_data_numeric <- data.frame(lapply(current_data_numeric, function(x) {
      if(is.factor(x) || is.character(x)) as.numeric(as.character(x)) else as.numeric(x)
    }))
    
    # Check for constant columns after ensuring numeric
    variances <- apply(current_data_numeric, 2, var, na.rm = TRUE)
    if (any(is.na(variances)) || any(variances == 0, na.rm = TRUE)) {
      zero_var_items <- names(which(is.na(variances) | variances == 0))
      cat("!!! Skipping: Zero or NA variance detected in items:", paste(zero_var_items, collapse=", "), "for", current_wave_id, "\n")
      
      valid_items <- setdiff(current_items, zero_var_items)
      if(length(valid_items) < 2) {
        cat("!!! Not enough items with variance remaining. Skipping analyses for", current_wave_id, "\n")
        next
      }
      # Proceeding with valid items only
      current_items <- valid_items
      current_data_numeric <- current_data_numeric[, current_items, drop = FALSE]
      cat("Proceeding analysis with remaining items:", paste(current_items, collapse=", "), "\n")
    }
    
    # --- Analysis Section 1: Descriptives, Reliability, Clustering ---
    cat("\n--- Section 1: Descriptives, Reliability, Clustering ---\n")
    tryCatch({
      
      # Basic descriptives (use numeric data)
      descriptives <- psych::describe(current_data_numeric)
      print("Descriptive Statistics:")
      print(descriptives)
      all_results[[current_wave_id]][["descriptives"]] <- descriptives
      
      # Polychoric correlations
      poly_results <- polychoric(current_data_numeric, na.rm = TRUE, correct = 0)
      poly_cor <- poly_results$rho
      print("Polychoric Correlation Matrix:")
      print(round(poly_cor, 2))
      all_results[[current_wave_id]][["poly_cor"]] <- poly_cor
      
      # Reliability
      alpha_res <- psych::alpha(current_data_numeric, na.rm = TRUE, check.keys = TRUE)
      omega_res <- psych::omega(poly_cor, nfactors = 1, fm = "minres")
      print("Reliability Estimates:")
      print(alpha_res$total)
      print(paste("Omega Total:", round(omega_res$omega_t, 3)))
      all_results[[current_wave_id]][["alpha"]] <- alpha_res
      all_results[[current_wave_id]][["omega"]] <- omega_res
      
      # Hierarchical Clustering
      poly_cor_matrix <- as.matrix(poly_cor) # Ensure it's a basic matrix
      if(is.numeric(poly_cor_matrix) && !anyNA(poly_cor_matrix) && nrow(poly_cor_matrix) > 1){
        dist_input <- 1 - abs(poly_cor_matrix)
        if(is.matrix(dist_input) && is.numeric(dist_input) && !anyNA(dist_input)){
          dist_matrix <- as.dist(dist_input)
          hclust_res <- hclust(dist_matrix, method = "ward.D2")
          print("Hierarchical Clustering (Ward.D2 method on 1-|rho|):")
          # plot(hclust_res) # Optional plot command
          all_results[[current_wave_id]][["hclust"]] <- hclust_res
        } else {
          cat("!!! Skipping hclust: Calculated distance matrix is invalid (non-numeric or NAs).\n")
          print("Distance matrix summary:")
          print(summary(as.vector(dist_input)))
        }
      } else {
        cat("!!! Skipping hclust: Polychoric matrix is invalid (non-numeric, NAs, or < 2 items).\n")
        print("Polychoric matrix summary:")
        print(summary(as.vector(poly_cor_matrix)))
      }
      
    }, error = function(e) {
      cat("!!! ERROR in Section 1 for", current_wave_id, ":", e$message, "\n")
    }) # End Section 1 tryCatch
    
    # --- Analysis Section 2: Exploratory Factor Analysis (EFA) ---
    cat("\n--- Section 2: EFA (using polychoric) ---\n")
    tryCatch({
      
      # EFA - Unidimensional (using numeric data + cor = "poly")
      efa_uni <- fa(current_data_numeric, nfactors = 1, fm = "minres", rotate = "none", cor = "poly")
      print("EFA - Unidimensional:")
      print(efa_uni$loadings)
      print(efa_uni$Vaccounted)
      all_results[[current_wave_id]][["efa_uni"]] <- efa_uni
      
      # EFA - Correlated Factors (Placeholder)
      if (length(current_items) >= 4) {
        cat("Skipping Correlated Factors EFA - Need theory/PA for number of factors.\n")
      } else {
        cat("Skipping Correlated Factors EFA - Not enough items.\n")
      }
      
      # EFA - Bifactor (Placeholder)
      cat("Skipping Bifactor EFA - Requires theoretical specific factors.\n")
      
    }, error = function(e) {
      cat("!!! ERROR in Section 2 for", current_wave_id, ":", e$message, "\n")
    }) # End Section 2 tryCatch
    
    # --- Analysis Section 3: Confirmatory Factor Analysis (CFA) ---
    cat("\n--- Section 3: CFA (WLSMV estimator) ---\n")
    tryCatch({
      # CFA - Unidimensional
      cfa_uni_syntax <- paste("factor =~", paste(current_items, collapse = " + "))
      fit_cfa_uni <- cfa(cfa_uni_syntax, data = data_to_analyze,
                         ordered = current_items, estimator = "WLSMV",
                         missing = "pairwise")
      print("CFA - Unidimensional Fit:")
      print(fitMeasures(fit_cfa_uni, c("chisq.scaled", "df.scaled", "pvalue.scaled", "cfi.scaled", "rmsea.scaled", "srmr")))
      all_results[[current_wave_id]][["cfa_uni_fit"]] <- fit_cfa_uni
      
      # CFA - Correlated Factors (Placeholder)
      cat("Skipping Correlated Factors CFA - Need theoretical model.\n")
      
      # CFA - Bifactor (Placeholder)
      cat("Skipping Bifactor CFA - Need theoretical model.\n")
      
    }, error = function(e) {
      cat("!!! ERROR in Section 3 for", current_wave_id, ":", e$message, "\n")
    }) # End Section 3 tryCatch
    
    # --- Analysis Section 4: Item Response Theory (IRT) ---
    cat("\n--- Section 4: IRT (Unidimensional) ---\n")
    tryCatch({
      # Ensure data is numeric matrix
      irt_data_matrix <- as.matrix(current_data_numeric)
      
      # Check for non-finite values (Inf, -Inf, NaN) excluding NAs
      if (!is.numeric(irt_data_matrix) || any(!is.finite(irt_data_matrix[!is.na(irt_data_matrix)]))) {
        cat("!!! Skipping IRT: Data matrix is not purely numeric or contains non-finite values (excluding NA).\n")
      } else if (nrow(irt_data_matrix) < 2 || ncol(irt_data_matrix) < 2) {
        cat("!!! Skipping IRT: Data matrix has less than 2 rows or 2 columns after handling missing/variance.\n")
      } else {
        # IRT - Graded Response Model (GRM)
        n_items_irt <- ncol(irt_data_matrix)
        model_grm <- mirt.model(paste0('F = 1-', n_items_irt))
        
        # Use EM method, verbose = FALSE
        fit_grm <- mirt(irt_data_matrix, model_grm, itemtype = 'graded',
                        method = 'EM', verbose = FALSE)
        
        print("IRT - GRM Fit (using M2):")
        # M2 might fail if data is sparse, wrap in tryCatch
        m2_results <- tryCatch(M2(fit_grm), error = function(e) {
          cat("!!! M2 fit calculation failed:", e$message, "\n"); return(NULL)
        })
        if (!is.null(m2_results)) print(m2_results)
        
        print("IRT - GRM Item Parameters (a=discrimination, b=thresholds):")
        # Extract coefficients, handle potential errors
        coefs <- tryCatch(coef(fit_grm, simplify = TRUE, IRTpars = TRUE)$items, error = function(e) {
          cat("!!! Failed to extract GRM coefficients:", e$message, "\n"); return(NULL)
        })
        if (!is.null(coefs)) print(coefs)
        
        all_results[[current_wave_id]][["irt_grm_fit"]] <- fit_grm
      }
      
      # IRT - Nominal Response Model (Placeholder)
      cat("Skipping NRM - GRM is standard for ordered items.\n")
      
    }, error = function(e) {
      # Catch errors specifically from the mirt() function itself
      cat("!!! ERROR in Section 4 (mirt execution) for", current_wave_id, ":", e$message, "\n")
    }) # End Section 4 tryCatch
    
  } # End loop through timepoints
} # End loop through reporters

# --- 3. Explore Results ---
# names(all_results)
# all_results$`_p_bl`$descriptives
# summary(all_results$`_y_fu3`$cfa_uni_fit, fit.measures=TRUE)

cat("\n\n--- Analysis Complete --- \n")

### checking if adding another factor in a model with parent and youth items improves model fit

# define items
parent_items <- c("peta_p","petb_p","petc_p","petd_p")
youth_items  <- c("peta_y","petb_y","petc_y","petd_y","pete_y")
items <- c(parent_items, youth_items)

long_data <- data_to_analyze %>%
  pivot_longer(
    cols = -c(id, sex),
    names_to = c(".value", "wave"),
    names_pattern = "(.*)_(bl|fu[1-6])"
  ) %>%
  mutate(
    wave = factor(wave, levels = c("bl","fu1","fu2","fu3","fu4","fu5","fu6"))
  ) %>%
  arrange(id, wave)

# make sure they exist
items <- items[items %in% names(long_data)]
if (length(items) < 6) stop("not enough items found to run a meaningful efa.")

# pool across waves; coerce to integer; drop rows with any missing on these items
efa_dat <- long_data %>%
  select(all_of(items)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.integer(.x)))) %>%
  tidyr::drop_na()

if (nrow(efa_dat) < 200) message("tiny sample after drop_na; results will wobble.")

# polychoric correlations
pc <- psych::polychoric(efa_dat)
r_poly <- pc$rho

# 1-factor efa (minres); no rotation needed for 1 factor
efa_1f <- psych::fa(r_poly, nfactors = 1, fm = "minres")
cat("\n--- efa: 1 factor (minres) ---\n")
print(efa_1f$loadings, cutoff = 0.20, sort = TRUE)
cat("\nvariance accounted for (1f):\n")
print(efa_1f$Vaccounted)

# 2-factor efa (minres + oblimin) to allow correlated factors
efa_2f <- psych::fa(r_poly, nfactors = 2, fm = "minres", rotate = "oblimin")
cat("\n--- efa: 2 factors (minres, oblimin) ---\n")
print(efa_2f$loadings, cutoff = 0.20, sort = TRUE)
cat("\nphi (factor correlations):\n")
print(efa_2f$Phi)

# quick readout: which factor each item prefers in 2f
lf <- as.data.frame(unclass(efa_2f$loadings))
lf$item <- rownames(lf)
lf$primary <- ifelse(abs(lf$ML1) >= abs(lf$ML2), "F1", "F2")
cat("\nprimary loading by item (2f):\n")
print(lf[, c("item","ML1","ML2","primary")])

