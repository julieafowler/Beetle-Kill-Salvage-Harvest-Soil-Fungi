#### DIAGNOSTIC SCRIPT: To be run before any soil PCA script or any univariate or multivariate Bayesian runs ####             
## Checks distributions, missingness, and helps decide on transformations
## Created with the assistance of Claude AI Opus 4.5

#### Packages ####

library(tidyverse)
library(naniar) ## missing data exploration
library(moments) ## looks at skewness, kurtosis of data, etc. 
library(ggcorrplot) ## correlation plots 

#### Load Data ####

setwd()
df <- read.delim("Bayesian_AllData_ExcludedThoseWithOnlyNAsAcross.txt", row.names = 1)

#### Information on Missingness ####

run_missingness_diagnostics <- function(df) {
  
  cat("\n========================================\n")
  cat("MISSINGNESS DIAGNOSTICS\n")
  cat("========================================\n\n")
  
  # Basic summary
  cat("Dataset dimensions:", nrow(df), "rows x", ncol(df), "columns\n\n")
  
  # Missingness by variable
  miss_summary <- df %>%
    summarise(across(everything(), ~sum(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
    mutate(
      pct_missing = round(n_missing / nrow(df) * 100, 1),
      missingness_category = case_when(
        pct_missing == 0 ~ "Complete",
        pct_missing < 10 ~ "Low (<10%)",
        pct_missing < 30 ~ "Moderate (10-30%)",
        pct_missing < 50 ~ "High (30-50%)",
        pct_missing < 80 ~ "Very High (50-80%)",
        pct_missing < 100 ~ "Extreme (80-99%)",
        TRUE ~ "Fully Missing"
      )
    ) %>%
    arrange(desc(pct_missing))
  
  cat("Variables by missingness category:\n")
  print(table(miss_summary$missingness_category))
  
  cat("\n\nVariables with >50% missing:\n")
  print(miss_summary %>% filter(pct_missing > 50))
  
  # Check for 100% missing (must be removed)
  fully_missing <- miss_summary %>% filter(pct_missing == 100)
  if (nrow(fully_missing) > 0) {
    cat("\n\nWARNING: These variables are 100% missing and MUST be removed:\n")
    print(fully_missing$variable)
  }
  
  # Complete cases
  cat("\n\nComplete cases (no missing values):", sum(complete.cases(df)), "\n")
  
  # Visualization
  p1 <- miss_summary %>%
    ggplot(aes(x = reorder(variable, pct_missing), y = pct_missing, fill = missingness_category)) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    labs(title = "Missing Data by Variable", x = "", y = "% Missing")
  
  return(list(summary = miss_summary, plot = p1))
}

run_missingness_diagnostics(df)


#### Information on Distributions of Data #### 

run_distribution_diagnostics <- function(df, vars) {
  
  cat("\n========================================\n")
  cat("DISTRIBUTION DIAGNOSTICS\n")
  cat("========================================\n\n")
  
  results <- data.frame()
  plots <- list()
  
  for (var in vars) {
    if (!var %in% names(df)) {
      cat("  Skipping", var, "(not found)\n")
      next
    }
    
    x <- df[[var]]
    x <- x[!is.na(x)]
    
    if (length(x) < 5) {
      cat("  Skipping", var, "(too few non-NA values)\n")
      next
    }
    
    # Calculate statistics
    sk <- skewness(x)
    kt <- kurtosis(x) - 3  # Excess kurtosis
    
    result <- data.frame(
      variable = var,
      n = length(x),
      n_missing = sum(is.na(df[[var]])),
      pct_missing = round(sum(is.na(df[[var]])) / nrow(df) * 100, 1),
      min = min(x),
      max = max(x),
      mean = mean(x),
      median = median(x),
      sd = sd(x),
      skewness = round(sk, 2),
      excess_kurtosis = round(kt, 2),
      has_zeros = any(x == 0),
      has_negatives = any(x < 0),
      suggested_transform = case_when(
        any(x < 0) ~ "none (has negatives)",
        abs(sk) < 0.5 ~ "none (symmetric)",
        abs(sk) < 1.0 ~ "consider sqrt",
        abs(sk) >= 1.0 & any(x == 0) ~ "log1p",
        abs(sk) >= 1.0 ~ "log",
        TRUE ~ "check manually"
      ),
      suggested_family = case_when(
        any(x < 0) ~ "gaussian",
        abs(sk) < 0.5 ~ "gaussian",
        all(x > 0) & sk > 1 ~ "lognormal or gamma",
        TRUE ~ "gaussian (default)"
      )
    )
    
    results <- rbind(results, result)
    
    # Create histogram
    plots[[var]] <- ggplot(data.frame(x = x), aes(x = x)) +
      geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.7) +
      geom_vline(xintercept = mean(x), color = "red", linetype = "dashed") +
      theme_minimal() +
      labs(title = paste(var, "\n(skew =", round(sk, 2), ")"),
           x = "", y = "Count")
  }
  
  cat("Distribution summary:\n")
  print(results %>% select(variable, n, skewness, suggested_transform, suggested_family))
  
  return(list(summary = results, plots = plots))
}


# Define variable groups
outcomes <- c("Height", "MeanCanWidth", "AbovegroundBiomass", 
              "RootDryMass", "RootToShootRatio")

soil_chem_vars <- c("Total_C_percent", "Total_N_percent", 
                    "Na_waterextract_mgperL", "NH4_waterextract_mgperL",
                    "K_waterextract_mgperL", "Mg_waterextract_mgperL",
                    "Ca_waterextract_mgperL", "Cl_waterextract_mgperL",
                    "NO3_waterextract_mgperL", "PO4_waterextract_mgperL",
                    "SO4_waterextract_mgperL")

microbial_vars <- c("Fungal_Shannons_Rhizo", "Fungal_ObservedFeatures_Rhizo",
                    "Fungal_EMF_RA_Rhizo", "Fungal_EMF_ObservedFeatures_Rhizo",
                    "Fungal_PlantPathogen_RA_Rhizo", "Fungal_FungalParasite_RA_Rhizo",
                    "Fungal_Shannons_RootTips", "Fungal_ObservedFeatures_RootTips",
                    "Fungal_EMF_RA_RootTips", "Fungal_EMF_ObservedFeatures_RootTips",
                    "Fungal_PlantPathogen_RA_RootTips", "Fungal_FungalParasite_RA_RootTips",
                    "Prokayote_Shannons_Rhizo", "Prokayote_ObservedFeatures_Rhizo")

results <- list()

# 3. Predictor distributions - soil chemistry
cat("\n--- Soil Chemistry Distributions ---\n")
results$soil_dist <- run_distribution_diagnostics(df, soil_chem_vars)

# 4. Predictor distributions - microbial
cat("\n--- Microbial Variable Distributions ---\n")
results$micro_dist <- run_distribution_diagnostics(df, microbial_vars)

## Ultimately, not transforming the predictors as we later choose to use horseshoe priors for the Bayesian models
## Which are robust to non-normality, etc. Scaling/standardization happens in the Bayesian scripts for predictors. 
## VIP score predictors (my sPLS predictors, using mean) are continuous/bounded, but can be checked by:
vip_vars <- names(df)[grep("^vip_mean_", names(df))]
results$vip_dist <- run_distribution_diagnostics(df, vip_vars)


#### Correlations ####

run_correlation_diagnostics <- function(df, vars, save_plot = TRUE, 
                                        filename = "correlation_plot.pdf",
                                        width = 14, height = 14) {
  
  cat("\n========================================\n")
  cat("CORRELATION DIAGNOSTICS\n")
  cat("========================================\n\n")
  
  # Select only numeric variables that exist
  vars_exist <- vars[vars %in% names(df)]
  df_numeric <- df %>% 
    select(all_of(vars_exist)) %>%
    select(where(is.numeric))
  
  if (ncol(df_numeric) < 2) {
    cat("Not enough numeric variables for correlation analysis\n")
    return(NULL)
  }
  
  # Pairwise complete correlations
  cor_matrix <- cor(df_numeric, use = "pairwise.complete.obs")
  
  # Find highly correlated pairs (>0.8)
  high_cor <- which(abs(cor_matrix) > 0.8 & abs(cor_matrix) < 1, arr.ind = TRUE)
  
  high_cor_df <- NULL
  
  if (nrow(high_cor) > 0) {
    cat("WARNING: Highly correlated variable pairs (|r| > 0.8):\n\n")
    
    high_cor_df <- data.frame(
      var1 = rownames(cor_matrix)[high_cor[, 1]],
      var2 = colnames(cor_matrix)[high_cor[, 2]],
      correlation = cor_matrix[high_cor]
    ) %>%
      filter(var1 < var2) %>%  # Remove duplicates
      arrange(desc(abs(correlation)))
    
    print(high_cor_df)
    
    cat("\n\nCONSIDERATION: You may want to remove one variable from each highly correlated pair,")
    cat("\nor the horseshoe prior will handle this automatically by shrinking redundant effects.\n")
  } else {
    cat("No highly correlated pairs found (|r| > 0.8)\n")
  }
  
  # Correlation plot using ggcorrplot (better sizing and export)
  # Requires: install.packages("ggcorrplot") if not already installed
  if (!requireNamespace("ggcorrplot", quietly = TRUE)) {
    cat("\nInstalling ggcorrplot for better correlation visualization...\n")
    install.packages("ggcorrplot", quiet = TRUE)
  }
  library(ggcorrplot)
  
  # Create the plot
  cor_plot <- ggcorrplot(cor_matrix, 
                         type = "upper",
                         lab = FALSE,
                         colors = c("#2166AC", "white", "#B2182B"),
                         ggtheme = theme_minimal(),
                         title = "Predictor Correlations") +
    theme(
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  # Display the plot
  print(cor_plot)
  
  # Save the plot if requested
  if (save_plot) {
    ggsave(filename, cor_plot, width = width, height = height, device = "pdf")
    cat("\nCorrelation plot saved to:", filename, "\n")
  }
  
  return(list(cor_matrix = cor_matrix, high_correlations = high_cor_df, plot = cor_plot))
}


run_correlation_diagnostics(df, c(soil_chem_vars, microbial_vars))

## Ultimately, not worrying about highly correlated variables too too much 
## Because horseshoe priors shrink one of a highly correlated pair towards zero 


#### Checks for Outcome Variables ####

check_outcomes <- function(df) {
  
  cat("\n========================================\n")
  cat("OUTCOME VARIABLE CHECKS\n")
  cat("========================================\n\n")
  
  outcomes <- c("Height", "MeanCanWidth", "AbovegroundBiomass", 
                "RootDryMass", "RootToShootRatio")
  
  results <- list()
  
  for (out in outcomes) {
    if (!out %in% names(df)) {
      cat(out, ": NOT FOUND\n\n")
      next
    }
    
    x <- df[[out]]
    x <- x[!is.na(x)]
    
    cat(out, ":\n")
    cat("  n =", length(x), " (", sum(is.na(df[[out]])), " missing)\n")
    cat("  Range: [", round(min(x), 3), ",", round(max(x), 3), "]\n")
    cat("  Mean (SD):", round(mean(x), 3), "(", round(sd(x), 3), ")\n")
    cat("  Skewness:", round(skewness(x), 2), "\n")
    
    # Shapiro-Wilk test for normality (if n < 5000)
    if (length(x) < 5000 && length(x) >= 3) {
      sw_test <- shapiro.test(x)
      cat("  Shapiro-Wilk p-value:", round(sw_test$p.value, 4), "\n")
    }
    
    # Recommendation
    sk <- skewness(x)
    if (abs(sk) < 0.5) {
      cat("  RECOMMENDATION: Gaussian family likely appropriate\n")
    } else if (abs(sk) < 1.0) {
      cat("  RECOMMENDATION: Consider Gaussian or lognormal\n")
    } else {
      if (all(x > 0)) {
        cat("  RECOMMENDATION: Consider lognormal or Gamma family\n")
      } else {
        cat("  RECOMMENDATION: Highly skewed, consider transformation\n")
      }
    }
    cat("\n")
    
    results[[out]] <- list(
      n = length(x),
      mean = mean(x),
      sd = sd(x),
      skewness = sk,
      min = min(x),
      max = max(x)
    )
  }
  
  # Create histogram panel
  outcome_data <- df %>%
    select(any_of(outcomes)) %>%
    pivot_longer(everything(), names_to = "outcome", values_to = "value") %>%
    filter(!is.na(value))
  
  p <- ggplot(outcome_data, aes(x = value)) +
    geom_histogram(bins = 25, fill = "steelblue", color = "white", alpha = 0.7) +
    facet_wrap(~outcome, scales = "free") +
    theme_minimal() +
    labs(title = "Outcome Variable Distributions", x = "", y = "Count")
  
  return(list(summary = results, plot = p))
}


check_outcomes(df) 

# Result:
#   
#   outcome_families <- list(
#     Height = gaussian(),
#     MeanCanWidth = gaussian(),
#     AbovegroundBiomass = lognormal(),
#     RootDryMass = lognormal(),
#     RootToShootRatio = lognormal()
#   )


#### Checking Correlations Between VIP (sPLS) Variables - Moving Forward Just Moving Vip_Means ####

check_vip_correlations <- function(df) {
  
  cat("\n========================================\n")
  cat("VIP VARIABLE CORRELATION CHECK\n")
  cat("========================================\n\n")
  
  cat("VIP sum and mean variables are derived from the same data and will be\n")
  cat("highly correlated. Check correlations below to decide whether to use\n")
  cat("sum, mean, or both.\n\n")
  
  # Find VIP variables
  vip_vars <- names(df)[grep("^vip_", names(df))]
  
  if (length(vip_vars) == 0) {
    cat("No VIP variables found\n")
    return(NULL)
  }
  
  # Group by outcome
  outcomes <- c("height", "meancanwidth", "abovegroundbiomass", 
                "rootdrymass", "roottoshootratio")
  
  for (out in outcomes) {
    out_vips <- vip_vars[grep(out, vip_vars)]
    if (length(out_vips) < 2) next
    
    cat("\n", toupper(out), "VIP variables:\n")
    
    vip_data <- df %>% select(all_of(out_vips))
    
    # Check sum vs mean correlation for rhizosphere
    sum_rhizo <- out_vips[grep("sum.*rhizo", out_vips)]
    mean_rhizo <- out_vips[grep("mean.*rhizo", out_vips)]
    
    if (length(sum_rhizo) > 0 && length(mean_rhizo) > 0) {
      for (i in seq_along(sum_rhizo)) {
        if (i <= length(mean_rhizo)) {
          r <- cor(df[[sum_rhizo[i]]], df[[mean_rhizo[i]]], use = "complete.obs")
          if (!is.na(r)) {
            cat("  Rhizo sum vs mean (", 
                ifelse(grepl("positive", sum_rhizo[i]), "pos", "neg"), 
                "): r = ", round(r, 3), "\n", sep = "")
          }
        }
      }
    }
    
    # Check sum vs mean correlation for roottips
    sum_root <- out_vips[grep("sum.*roottips", out_vips)]
    mean_root <- out_vips[grep("mean.*roottips", out_vips)]
    
    if (length(sum_root) > 0 && length(mean_root) > 0) {
      for (i in seq_along(sum_root)) {
        if (i <= length(mean_root)) {
          r <- cor(df[[sum_root[i]]], df[[mean_root[i]]], use = "complete.obs")
          if (!is.na(r)) {
            cat("  RootTips sum vs mean (", 
                ifelse(grepl("positive", sum_root[i]), "pos", "neg"), 
                "): r = ", round(r, 3), "\n", sep = "")
          }
        }
      }
    }
  }
  
  cat("\n\nRECOMMENDATION: If sum and mean are highly correlated (r > 0.9),\n")
  cat("consider using only ONE of them to avoid multicollinearity.\n")
  cat("The 'mean' is often more interpretable.\n")
}


check_vip_correlations(df)


#### Checking Treatment Balance ####

check_treatment_balance <- function(df) {
  
  cat("\n========================================\n")
  cat("TREATMENT BALANCE CHECK\n")
  cat("========================================\n\n")
  
  if (!"Treatment_LandUseHistory" %in% names(df)) {
    cat("Treatment_LandUseHistory not found\n")
    return(NULL)
  }
  
  cat("Treatment_LandUseHistory levels:\n")
  print(table(df$Treatment_LandUseHistory, useNA = "ifany"))
  
  if ("MPB_Impacted" %in% names(df)) {
    cat("\nMPB_Impacted levels:\n")
    print(table(df$MPB_Impacted, useNA = "ifany"))
    
    cat("\nCross-tabulation:\n")
    print(table(df$Treatment_LandUseHistory, df$MPB_Impacted, useNA = "ifany"))
  }
  
  # Check outcome availability by treatment
  outcomes <- c("Height", "MeanCanWidth", "AbovegroundBiomass", 
                "RootDryMass", "RootToShootRatio")
  
  cat("\n\nOutcome availability by Treatment_LandUseHistory:\n")
  for (out in outcomes) {
    if (out %in% names(df)) {
      tab <- df %>%
        group_by(Treatment_LandUseHistory) %>%
        summarise(
          n_total = n(),
          n_available = sum(!is.na(.data[[out]])),
          pct_available = round(n_available / n_total * 100, 1)
        )
      cat("\n", out, ":\n")
      print(as.data.frame(tab))
    }
  }
}

check_treatment_balance(df)





## Conclusions 

# Based on the above tests, the following adjustments will be made to the final analysis script:
#   
#   Summary of Updates
#   Treatment variable: Using Treatment_LandUseHistory seperately in models from MPB_Impacted - perfectly nested, so testing seperately (two univariate runs)
#   Prokaryote variable: Now uses only Prokayote_Shannons_Rhizo (dropped ObservedFeatures due to r=0.89 correlation)
#   VIP variables: Now uses only vip_mean_* (dropped all vip_sum_* since r=1 correlation)
#   100% missing VIPs: Excluded the two roottips negative association variables
#   Outcome families: Height & MeanCanWidth: gaussian() / AbovegroundBiomass, RootDryMass, RootToShootRatio: lognormal()


# Interpreting Lognormal Results
# For the three outcomes using lognormal() family:
#   
#   Coefficients are on the log scale
# exp(β) gives the fold change (multiplicative effect)
# Example: β = 0.3 → exp(0.3) = 1.35 → 35% increase per 1 SD change in predictor
# Example: β = -0.2 → exp(-0.2) = 0.82 → 18% decrease per 1 SD change in predictor



# Updated Predictors per Outcome Model
# Treatment:   Treatment_LandUseHistory    1 (5 df for 6 levels)
# Soil PCs:    Soil_PC1, Soil_PC2, Soil_PC3    3
# Rhizo Fungal: Shannon, ObservedFeatures, EMF_RA, EMF_ObservedFeatures, PlantPathogen_RA, FungalParasite_RA    6
# RootTip Fungal:  Same structure as Rhizo.  6
# Prokaryote: Prokayote_Shannons_Rhizo only. 1
# VIP (outcome-specific): 3-4 depending on outcome.   3-4
# Total predictors~22-23
# 
# Whats Excluded and Why
# MPB_Impacted:   Perfectly nested in Treatment_LandUseHistory
# Texture_Percent*:    Constant within treatment (confounded with site)
# Prokayote_ObservedFeatures_Rhizo:   r=0.89 with Shannon
# vip_sum_*:    r=1 with vip_mean_*
# Two 100% missing VIPs:   No data








