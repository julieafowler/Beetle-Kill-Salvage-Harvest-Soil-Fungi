################################################################################
#                                                                              #
#     FULL COEFFICIENT TABLES FOR PUBLICATION                                  #
#                                                                              #
#     Generates comprehensive tables of ALL coefficients from the              #
#     Bayesian horseshoe analysis, not just selected variables                 #
#                                                                              #
#     Outputs:                                                                 #
#       - CSV tables with all coefficients, CIs, and probabilities             #
#       - Formatted tables ready for publication (main text + supplement)      #
#       - Effect sizes in original units (optional)                            #
#                                                                              #
################################################################################

## Created with the assistance of Claude AI Opus 4.5

library(tidyverse)

# =============================================================================
# SECTION 1: LOAD EXISTING RESULTS
# =============================================================================

#' This script assumes you have already run your main analysis and have:
#'   - model_results.csv (from Greenhouse_Bayesian_Analysis_6FactorLandUseHistory.R)
#'   - model_results_SH.csv (from Greenhouse_Bayesian_Analysis_2FactorSHStatus.R)
#'   - df: your original data (for computing original-scale effects)


# =============================================================================
# SECTION 2: STANDARDIZE COLUMN NAMES
# =============================================================================

standardize_column_names <- function(results_df) {
  #' Standardize column names from various possible formats
  #' R may convert special characters when reading CSV
  #' e.g., "2.5%" -> "X2.5.", "Pr(>0)" -> "Pr..0."
  #' 
  #' @param results_df Data frame with model results
  #' @return Data frame with standardized column names
  
  # Get current column names
  cols <- names(results_df)
  cat("  Original columns:", paste(cols, collapse = ", "), "\n")
  
  # Create new standardized names
  new_names <- cols
  
  # Map various possible column names to standard names
  for (i in seq_along(cols)) {
    col <- cols[i]
    
    # Outcome
    if (col %in% c("Outcome", "outcome")) new_names[i] <- "outcome"
    
    # Predictor/Parameter  
    if (col %in% c("Parameter", "Predictor", "parameter", "predictor")) new_names[i] <- "predictor"
    
    # Estimate
    if (col %in% c("Estimate", "estimate")) new_names[i] <- "estimate"
    
    # Family
    if (col %in% c("Family", "family")) new_names[i] <- "family"
    
    # SE
    if (col %in% c("SE", "se", "std.error")) new_names[i] <- "se"
    
    # CI lower - handle many variations
    if (col %in% c("2.5%", "X2.5.", "X2.5", "Q2.5", "CI_low", "ci_lower", "lower", "l-95% CI")) {
      new_names[i] <- "ci_lower"
    }
    
    # CI upper - handle many variations
    if (col %in% c("97.5%", "X97.5.", "X97.5", "Q97.5", "CI_high", "ci_upper", "upper", "u-95% CI")) {
      new_names[i] <- "ci_upper"
    }
    
    # Probability > 0
    if (col %in% c("Pr(>0)", "Pr..0.", "Pr...0.", "prob_positive", "prob_gt_zero", "P(>0)")) {
      new_names[i] <- "prob_positive"
    }
    
    # Probability nonzero
    if (col %in% c("Pr(nonzero)", "Pr.nonzero.", "prob_nonzero", "P(nonzero)")) {
      new_names[i] <- "prob_nonzero"
    }
    
    # Fold change columns
    if (col %in% c("Fold Change", "Fold.Change", "fold_change")) new_names[i] <- "fold_change"
    if (col %in% c("FC 2.5%", "FC.2.5.", "fc_lower")) new_names[i] <- "fc_lower"
    if (col %in% c("FC 97.5%", "FC.97.5.", "fc_upper")) new_names[i] <- "fc_upper"
  }
  
  # Apply new names
  names(results_df) <- new_names
  cat("  Standardized columns:", paste(new_names, collapse = ", "), "\n")
  
  return(results_df)
}


# =============================================================================
# SECTION 3: CREATE FULL COEFFICIENT TABLE
# =============================================================================

create_full_coefficient_table <- function(results_df, model_name = "Model") {
  #' Create a comprehensive coefficient table from model results
  #' 
  #' @param results_df Data frame with model results (from main analysis)
  #' @param model_name Name for labeling ("Land-Use History" or "SH")
  #' @return Formatted data frame ready for publication
  
  cat("Creating full coefficient table for:", model_name, "\n")
  
  # Standardize column names first
  results_df <- standardize_column_names(results_df)
  
  # Check required columns
  required_cols <- c("outcome", "predictor", "estimate", "ci_lower", "ci_upper", "prob_nonzero")
  missing_cols <- setdiff(required_cols, names(results_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Add prob_positive if missing (estimate sign-based approximation)
  if (!"prob_positive" %in% names(results_df)) {
    results_df$prob_positive <- ifelse(results_df$estimate > 0, 
                                        results_df$prob_nonzero, 
                                        1 - results_df$prob_nonzero)
  }
  
  # Create formatted table
  full_table <- results_df %>%
    mutate(
      # Clean predictor names for display
      predictor_clean = case_when(
        str_detect(predictor, "Treatment_LandUseHistory") ~ 
          str_replace(predictor, "Treatment_LandUseHistory", ""),
        str_detect(predictor, "Salvage_Harvest_Status") ~ 
          str_replace(predictor, "Salvage_Harvest_Status", "SH: "),
        str_detect(predictor, "Soil_PC") ~ predictor,
        str_detect(predictor, "vip_mean_ITS") ~ 
          str_replace_all(predictor, c(
            "vip_mean_ITS_rhizosphere_" = "VIP Rhizo ",
            "vip_mean_ITS_roottips_" = "VIP RootTips ",
            "_positiveassociations" = " (+)",
            "_negativeassociations" = " (-)"
          )),
        str_detect(predictor, "Fungal_") ~
          str_replace_all(predictor, c(
            "Fungal_" = "",
            "_Rhizo" = " (Rhizo)",
            "_RootTips" = " (RootTips)",
            "Shannons" = "Shannon's H",
            "ObservedFeatures" = "Richness",
            "EMF_RA" = "EMF Rel. Abund.",
            "EMF_ObservedFeatures" = "EMF Richness",
            "PlantPathogen_RA" = "Plant Pathogen RA",
            "FungalParasite_RA" = "Fungal Parasite RA"
          )),
        str_detect(predictor, "Prokayote") ~
          str_replace(predictor, "Prokayote_Shannons_Rhizo", "Prokaryote Shannon (Rhizo)"),
        predictor == "pH" ~ "pH",
        predictor == "Intercept" ~ "Intercept",
        TRUE ~ predictor
      ),
      
      # Format estimate with CI
      estimate_ci = sprintf("%.3f [%.3f, %.3f]", estimate, ci_lower, ci_upper),
      
      # Significance indicators
      sig_90 = prob_nonzero > 0.90,
      sig_95 = prob_nonzero > 0.95,
      sig_99 = prob_nonzero > 0.99,
      
      sig_stars = case_when(
        sig_99 ~ "***",
        sig_95 ~ "**",
        sig_90 ~ "*",
        TRUE ~ ""
      ),
      
      # Direction indicator
      direction = case_when(
        sig_90 & estimate > 0 ~ "+",
        sig_90 & estimate < 0 ~ "-",
        TRUE ~ ""
      ),
      
      # Probability formatted
      prob_fmt = sprintf("%.1f%%", prob_nonzero * 100)
    ) %>%
    arrange(outcome, desc(prob_nonzero))
  
  cat("  Created table with", nrow(full_table), "coefficients\n")
  
  return(full_table)
}


# =============================================================================
# SECTION 4: CREATE SUMMARY TABLE (SELECTED VARIABLES ONLY)
# =============================================================================

create_selected_variables_table <- function(full_table, threshold = 0.90) {
  #' Extract only the selected variables (prob_nonzero > threshold)
  #' 
  #' @param full_table Full coefficient table
  #' @param threshold Selection threshold (default 0.90)
  #' @return Filtered table with only selected variables
  
  selected <- full_table %>%
    filter(prob_nonzero > threshold) %>%
    filter(predictor_clean != "Intercept") %>%
    arrange(outcome, desc(prob_nonzero))
  
  cat("  Selected", nrow(selected), "variables with P(nonzero) >", threshold, "\n")
  
  return(selected)
}


# =============================================================================
# SECTION 5: CONVERT TO ORIGINAL UNITS
# =============================================================================

add_original_units <- function(coef_table, df, outcomes_info = NULL) {
  #' Add effect sizes in original units
  #' 
  #' @param coef_table Coefficient table
  #' @param df Original data frame (for computing means and SDs)
  #' @param outcomes_info Optional: pre-computed outcome statistics
  #' @return Table with original-unit columns added
  
  cat("Adding effect sizes in original units...\n")
  
 # Compute outcome statistics if not provided
  if (is.null(outcomes_info)) {
    # Check which outcome columns exist in df
    possible_outcomes <- c("Height", "MeanCanWidth", "AbovegroundBiomass", 
                           "RootDryMass", "RootToShootRatio")
    existing_outcomes <- possible_outcomes[possible_outcomes %in% names(df)]
    
    if (length(existing_outcomes) == 0) {
      cat("  Warning: No outcome columns found in df. Skipping original units.\n")
      return(coef_table)
    }
    
    outcomes_info <- data.frame(
      outcome = existing_outcomes,
      mean_val = sapply(existing_outcomes, function(x) mean(df[[x]], na.rm = TRUE)),
      sd_val = sapply(existing_outcomes, function(x) sd(df[[x]], na.rm = TRUE)),
      family = ifelse(existing_outcomes %in% c("Height", "MeanCanWidth"), 
                      "gaussian", "lognormal"),
      units = c("cm", "cm", "g", "g", "ratio")[match(existing_outcomes, possible_outcomes)],
      stringsAsFactors = FALSE
    )
    row.names(outcomes_info) <- NULL
  }
  
  # Join with coefficient table
  coef_with_units <- coef_table %>%
    left_join(outcomes_info, by = "outcome") %>%
    mutate(
      # For Gaussian outcomes: effect in original units = estimate * SD
      # For Lognormal outcomes: fold change = exp(estimate)
      effect_original = case_when(
        family == "gaussian" ~ estimate * sd_val,
        family == "lognormal" ~ exp(estimate) - 1,  # Percent change
        TRUE ~ NA_real_
      ),
      
      effect_original_fmt = case_when(
        family == "gaussian" ~ sprintf("%.2f %s", effect_original, units),
        family == "lognormal" ~ sprintf("%.1f%%", effect_original * 100),
        TRUE ~ NA_character_
      ),
      
      # Fold change for lognormal (easier to interpret)
      fold_change_calc = case_when(
        family == "lognormal" ~ exp(estimate),
        TRUE ~ NA_real_
      ),
      
      fold_change_fmt = case_when(
        family == "lognormal" ~ sprintf("%.2fx", fold_change_calc),
        TRUE ~ NA_character_
      )
    )
  
  cat("  Added original-unit effects for", nrow(coef_with_units), "coefficients\n")
  
  return(coef_with_units)
}


# =============================================================================
# SECTION 6: FORMAT FOR PUBLICATION
# =============================================================================

format_for_publication <- function(coef_table, table_type = "main") {
  #' Format table for publication (Word/PDF ready)
  #' 
  #' @param coef_table Coefficient table with all info
  #' @param table_type "main" for main text, "supplement" for full table
  #' @return Formatted table
  
  if (table_type == "main") {
    # Main text: selected variables only, key columns
    # Check if effect_original_fmt exists
    if ("effect_original_fmt" %in% names(coef_table)) {
      pub_table <- coef_table %>%
        filter(prob_nonzero > 0.90, predictor_clean != "Intercept") %>%
        select(
          Outcome = outcome,
          Predictor = predictor_clean,
          `Estimate [95% CI]` = estimate_ci,
          `P(nonzero)` = prob_fmt,
          `Direction` = direction,
          `Effect Size` = effect_original_fmt
        )
    } else {
      pub_table <- coef_table %>%
        filter(prob_nonzero > 0.90, predictor_clean != "Intercept") %>%
        select(
          Outcome = outcome,
          Predictor = predictor_clean,
          `Estimate [95% CI]` = estimate_ci,
          `P(nonzero)` = prob_fmt,
          `Direction` = direction
        )
    }
  } else {
    # Supplement: all predictors, full details
    pub_table <- coef_table %>%
      filter(predictor_clean != "Intercept") %>%
      select(
        Outcome = outcome,
        Predictor = predictor_clean,
        Estimate = estimate,
        `Lower 95%` = ci_lower,
        `Upper 95%` = ci_upper,
        `P(>0)` = prob_positive,
        `P(nonzero)` = prob_nonzero,
        Sig = sig_stars
      ) %>%
      mutate(
        Estimate = round(Estimate, 3),
        `Lower 95%` = round(`Lower 95%`, 3),
        `Upper 95%` = round(`Upper 95%`, 3),
        `P(>0)` = round(`P(>0)`, 3),
        `P(nonzero)` = round(`P(nonzero)`, 3)
      )
  }
  
  return(pub_table)
}


# =============================================================================
# SECTION 7: CREATE COMPARISON TABLE (6-LEVEL VS 2-LEVEL)
# =============================================================================

create_comparison_table <- function(table_landuse, table_mpb) {
  #' Create side-by-side comparison of Land-Use History vs SH models
  #' 
  #' @param table_landuse Full table from 6-level model
  #' @param table_mpb Full table from 2-level model
  #' @return Comparison table
  
  cat("Creating model comparison table...\n")
  
  # Get selected variables from each
  selected_landuse <- table_landuse %>%
    filter(prob_nonzero > 0.90, predictor_clean != "Intercept") %>%
    select(outcome, predictor_clean, estimate_landuse = estimate, 
           prob_landuse = prob_nonzero, sig_landuse = sig_stars)
  
  selected_mpb <- table_mpb %>%
    filter(prob_nonzero > 0.90, predictor_clean != "Intercept") %>%
    select(outcome, predictor_clean, estimate_mpb = estimate, 
           prob_mpb = prob_nonzero, sig_mpb = sig_stars)
  
  # Full join to see which variables selected in each
  comparison <- full_join(
    selected_landuse, 
    selected_mpb, 
    by = c("outcome", "predictor_clean")
  ) %>%
    mutate(
      selected_in = case_when(
        !is.na(estimate_landuse) & !is.na(estimate_mpb) ~ "Both",
        !is.na(estimate_landuse) ~ "Land-Use Only",
        !is.na(estimate_mpb) ~ "SH Only",
        TRUE ~ "Neither"
      )
    ) %>%
    arrange(outcome, predictor_clean)
  
  cat("  Variables selected in both models:", 
      sum(comparison$selected_in == "Both"), "\n")
  cat("  Variables selected only in Land-Use:", 
      sum(comparison$selected_in == "Land-Use Only"), "\n")
  cat("  Variables selected only in SH:", 
      sum(comparison$selected_in == "SH Only"), "\n")
  
  return(comparison)
}


# =============================================================================
# SECTION 8: MAIN FUNCTION
# =============================================================================

generate_coefficient_tables <- function(results_landuse, results_mpb, df = NULL,
                                        output_prefix = "coef_table") {
  #' Main function to generate all coefficient tables
  #' 
  #' @param results_landuse Results from 6-level treatment analysis
  #' @param results_mpb Results from 2-level SH analysis
  #' @param df Original data frame (optional, for original-unit effects)
  #' @param output_prefix Prefix for output file names
  #' @return List with all tables
  
  cat("\n")
  cat("================================================================\n")
  cat("  GENERATING FULL COEFFICIENT TABLES                            \n")
  cat("================================================================\n\n")
  
  # Create full tables
  cat("Step 1: Creating full coefficient tables...\n")
  full_landuse <- create_full_coefficient_table(results_landuse, "Land-Use History")
  full_mpb <- create_full_coefficient_table(results_mpb, "Salvage_Harvest_Status")
  
  # Add original units if df is provided
  if (!is.null(df)) {
    cat("\nStep 2: Adding original-unit effect sizes...\n")
    full_landuse <- add_original_units(full_landuse, df)
    full_mpb <- add_original_units(full_mpb, df)
  } else {
    cat("\nStep 2: Skipping original-unit effects (no data frame provided)\n")
  }
  
  # Create selected-only tables
  cat("\nStep 3: Creating selected variables tables...\n")
  selected_landuse <- create_selected_variables_table(full_landuse)
  selected_mpb <- create_selected_variables_table(full_mpb)
  
  # Create comparison table
  cat("\nStep 4: Creating model comparison table...\n")
  comparison <- create_comparison_table(full_landuse, full_mpb)
  
  # Format for publication
  cat("\nStep 5: Formatting for publication...\n")
  pub_main_landuse <- format_for_publication(full_landuse, "main")
  pub_main_mpb <- format_for_publication(full_mpb, "main")
  pub_supp_landuse <- format_for_publication(full_landuse, "supplement")
  pub_supp_mpb <- format_for_publication(full_mpb, "supplement")
  
  # Save all tables
  cat("\nStep 6: Saving tables...\n")
  
  write.csv(full_landuse, paste0(output_prefix, "_full_landuse.csv"), row.names = FALSE)
  write.csv(full_mpb, paste0(output_prefix, "_full_mpb.csv"), row.names = FALSE)
  write.csv(selected_landuse, paste0(output_prefix, "_selected_landuse.csv"), row.names = FALSE)
  write.csv(selected_mpb, paste0(output_prefix, "_selected_mpb.csv"), row.names = FALSE)
  write.csv(comparison, paste0(output_prefix, "_comparison.csv"), row.names = FALSE)
  write.csv(pub_main_landuse, paste0(output_prefix, "_pub_main_landuse.csv"), row.names = FALSE)
  write.csv(pub_main_mpb, paste0(output_prefix, "_pub_main_mpb.csv"), row.names = FALSE)
  write.csv(pub_supp_landuse, paste0(output_prefix, "_pub_supplement_landuse.csv"), row.names = FALSE)
  write.csv(pub_supp_mpb, paste0(output_prefix, "_pub_supplement_mpb.csv"), row.names = FALSE)
  
  cat("\nOutput files created:\n")
  cat("  ", paste0(output_prefix, "_full_landuse.csv"), "(all coefficients, 6-level)\n")
  cat("  ", paste0(output_prefix, "_full_mpb.csv"), "(all coefficients, 2-level)\n")
  cat("  ", paste0(output_prefix, "_selected_landuse.csv"), "(selected only, 6-level)\n")
  cat("  ", paste0(output_prefix, "_selected_mpb.csv"), "(selected only, 2-level)\n")
  cat("  ", paste0(output_prefix, "_comparison.csv"), "(model comparison)\n")
  cat("  ", paste0(output_prefix, "_pub_main_landuse.csv"), "(publication main text)\n")
  cat("  ", paste0(output_prefix, "_pub_main_mpb.csv"), "(publication main text)\n")
  cat("  ", paste0(output_prefix, "_pub_supplement_landuse.csv"), "(publication supplement)\n")
  cat("  ", paste0(output_prefix, "_pub_supplement_mpb.csv"), "(publication supplement)\n")
  
  cat("\n================================================================\n")
  cat("  COEFFICIENT TABLES COMPLETE                                    \n")
  cat("================================================================\n")
  
  return(list(
    full_landuse = full_landuse,
    full_mpb = full_mpb,
    selected_landuse = selected_landuse,
    selected_mpb = selected_mpb,
    comparison = comparison,
    pub_main_landuse = pub_main_landuse,
    pub_main_mpb = pub_main_mpb,
    pub_supplement_landuse = pub_supp_landuse,
    pub_supplement_mpb = pub_supp_mpb
  ))
}


# =============================================================================
# SECTION 9: EXECUTE
# =============================================================================

# To run:
# 
# # Load your results
results_landuse <- read.csv("model_results.csv")
results_mpb <- read.csv("model_results_SH.csv")
# 
# # Generate all tables (without original units)
coef_tables <- generate_coefficient_tables(results_landuse, results_mpb)
# 
# # OR with original units (requires df)
# coef_tables <- generate_coefficient_tables(results_landuse, results_mpb, df)

cat("Full coefficient tables script loaded.\n")
cat("Run: coef_tables <- generate_coefficient_tables(results_landuse, results_mpb, df)\n")
