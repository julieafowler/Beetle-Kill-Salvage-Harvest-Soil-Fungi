#### Bayesian Analysis - Univariate (5 Plant Outcomes) - w/ Multiple Imputation and Variable Selection for Later Multivariate Work ####
## This script uses LandUseHistory (6 variables) as the primary treatment in the models; there is another script that uses MPB-impact status (impacted by MPB or not) as 
## the dominant treatment - cannot do both in the same model as perfectly nested but both are relevant to the study 
## Assumes have already run Greenhouse_SoilChemistry_PCA_forBayesian.R and have three PCs for use here 
## Created with the assistance of Claude AI Opus 4.5


#### Setup and Packages ####

# Install packages 
# install.packages(c("brms", "mice", "tidyverse", "ggplot2", "patchwork", 
#                    "bayesplot", "loo", "tidybayes", "GGally", "naniar",
#                    "corrplot", "posterior",
#                    "projpred", "cmdstanr"))

# For brms, cmdstanr backend is recommended for complex models
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# cmdstanr::install_cmdstan()

library(tidyverse)
library(brms)
library(mice)
library(naniar)        # For missingness visualization
library(bayesplot)     # For MCMC diagnostics
library(loo)           # For model comparison
library(tidybayes)     # For tidy posterior extraction
library(patchwork)     # For combining plots
library(corrplot)      # For correlation matrices
library(posterior)     # For posterior summaries

# Set seed for reproducibility
set.seed(42)

# Set brms options
options(mc.cores = parallel::detectCores() - 1)  # Leave one core free
options(brms.backend = "cmdstanr")  # Use cmdstanr if available, otherwise "rstan"

#### Data Prep ####

# --- Define variable groups ---

# Outcome variables
outcomes <- c("Height", "MeanCanWidth", "AbovegroundBiomass", 
              "RootDryMass", "RootToShootRatio")

# Treatment variables (categorical)
# NOTE: MPB_Impacted is perfectly nested within Treatment_LandUseHistory
#       (each LandUseHistory level maps to exactly one MPB status)
#       Using only Treatment_LandUseHistory for more detailed analysis
treatments <- c("Treatment_LandUseHistory")

# Soil chemistry variables (pH listed for reference but handled separately - not in PCA)
soil_chem_vars <- c("pH", "Total_C_percent", "Total_N_percent", 
                    "Na_waterextract_mgperL", "NH4_waterextract_mgperL",
                    "K_waterextract_mgperL", "Mg_waterextract_mgperL",
                    "Ca_waterextract_mgperL", "Cl_waterextract_mgperL",
                    "NO3_waterextract_mgperL", "PO4_waterextract_mgperL",
                    "SO4_waterextract_mgperL")

# Texture variables (100% complete)
# NOTE: Texture is confounded with Treatment_LandUseHistory (each treatment = 1 site = 1 texture)
#       NOT included in regression models, but retained for reference
texture_vars <- c("Texture_PercentClay", "Texture_PercentSilt", "Texture_PercentSand")

# Rhizosphere fungal variables
rhizo_fungal_vars <- c("Fungal_Shannons_Rhizo", "Fungal_ObservedFeatures_Rhizo",
                       "Fungal_EMF_RA_Rhizo", "Fungal_EMF_ObservedFeatures_Rhizo",
                       "Fungal_PlantPathogen_RA_Rhizo", "Fungal_FungalParasite_RA_Rhizo")

# Root tip fungal variables
roottip_fungal_vars <- c("Fungal_Shannons_RootTips", "Fungal_ObservedFeatures_RootTips",
                         "Fungal_EMF_RA_RootTips", "Fungal_EMF_ObservedFeatures_RootTips",
                         "Fungal_PlantPathogen_RA_RootTips", "Fungal_FungalParasite_RA_RootTips")

# Prokaryote variables
# NOTE: Shannons and ObservedFeatures are highly correlated (r=0.89)
#       Using only Shannons for parsimony
prokaryote_vars <- c("Prokayote_Shannons_Rhizo")

# VIP variables by outcome (outcome-specific)
# NOTE: Using only vip_mean_* variables (not vip_sum_*) because sum and mean are perfectly correlated (r=1)
# NOTE: vip_mean_ITS_roottips_meancanwidth_negativeassociations and 
#       vip_mean_ITS_roottips_abovegroundbiomass_negativeassociations are 100% missing - EXCLUDED

vip_vars_height <- c("vip_mean_ITS_rhizosphere_height_positiveassociations",
                     "vip_mean_ITS_rhizosphere_height_negativeassociations",
                     "vip_mean_ITS_roottips_height_positiveassociations",
                     "vip_mean_ITS_roottips_height_negativeassociations")

vip_vars_meancanwidth <- c("vip_mean_ITS_rhizosphere_meancanwidth_positiveassociations",
                           "vip_mean_ITS_rhizosphere_meancanwidth_negativeassociations",
                           "vip_mean_ITS_roottips_meancanwidth_positiveassociations")
                           # EXCLUDED: vip_mean_ITS_roottips_meancanwidth_negativeassociations (100% missing)

vip_vars_abovegroundbiomass <- c("vip_mean_ITS_rhizosphere_abovegroundbiomass_positiveassociations",
                                  "vip_mean_ITS_rhizosphere_abovegroundbiomass_negativeassociations",
                                  "vip_mean_ITS_roottips_abovegroundbiomass_positiveassociations")
                                  # EXCLUDED: vip_mean_ITS_roottips_abovegroundbiomass_negativeassociations (100% missing)

vip_vars_rootdrymass <- c("vip_mean_ITS_rhizosphere_rootdrymass_positiveassociations",
                          "vip_mean_ITS_rhizosphere_rootdrymass_negativeassociations",
                          "vip_mean_ITS_roottips_rootdrymass_positiveassociations",
                          "vip_mean_ITS_roottips_rootdrymass_negativeassociations")

vip_vars_roottoshootratio <- c("vip_mean_ITS_rhizosphere_roottoshootratio_positiveassociations",
                               "vip_mean_ITS_rhizosphere_roottoshootratio_negativeassociations",
                               "vip_mean_ITS_roottips_roottoshootratio_positiveassociations",
                               "vip_mean_ITS_roottips_roottoshootratio_negativeassociations")

# --- Outcome family specifications (based on diagnostic results) ---
# Height: Gaussian (skew=0.14, Shapiro p=0.79, symmetric)
# MeanCanWidth: Gaussian (skew=0.05, Shapiro p=0.20, symmetric)
# AbovegroundBiomass: Lognormal (skew=1.05, right-skewed, positive values)
# RootDryMass: Lognormal (skew=1.97, right-skewed, positive values)
# RootToShootRatio: Lognormal (skew=2.10, right-skewed, positive values)

outcome_families <- list(
  Height = gaussian(),
  MeanCanWidth = gaussian(),
  AbovegroundBiomass = lognormal(),
  RootDryMass = lognormal(),
  RootToShootRatio = lognormal()
)

# For lognormal family, coefficients represent log-scale effects
# exp(beta) gives the multiplicative effect (fold-change)


#### Missingness ####

visualize_missingness <- function(df) {
  
  # Overall missingness summary
  cat("\n=== MISSINGNESS SUMMARY ===\n")
  miss_summary <- df %>%
    summarise(across(everything(), ~sum(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
    mutate(pct_missing = round(n_missing / nrow(df) * 100, 1)) %>%
    arrange(desc(pct_missing))
  
  print(miss_summary, n = 30)
  
  # Visualize missingness pattern
  p1 <- gg_miss_var(df, show_pct = TRUE) +
    theme_minimal() +
    labs(title = "Missing Values by Variable")
  
  # Missingness upset plot (shows combinations)
  p2 <- gg_miss_upset(df, nsets = 10) 
  
  # Heatmap of missingness
  p3 <- vis_miss(df, sort_miss = TRUE, cluster = TRUE) +
    theme(axis.text.x = element_text(angle = 90, size = 6))
  
  return(list(summary = miss_summary, 
              var_plot = p1, 
              upset_plot = p2, 
              heatmap = p3))
}

# Run missingness diagnostics 
miss_diagnostics <- visualize_missingness(df)
print(miss_diagnostics$var_plot)
print(miss_diagnostics$heatmap)


#### Data Preprocessing #### 

preprocess_data <- function(df, scale_vars = NULL) {
  
  df_processed <- df
  
  # --- Convert factors ---
  df_processed <- df_processed %>%
    mutate(
      Treatment_LandUseHistory = factor(Treatment_LandUseHistory),
      # Set Old_Forest as reference level (all comparisons relative to Old_Forest)
      Treatment_LandUseHistory = relevel(Treatment_LandUseHistory, ref = "Old_Forest")
    )
  
  # --- Scaling (after imputation is better, but prepare the function) ---
  # Note: Scale AFTER imputation for continuous predictors
  
  return(df_processed)
}

# Scale continuous predictors (run after imputation)
scale_predictors <- function(df, vars_to_scale) {
  df_scaled <- df
  
  for (var in vars_to_scale) {
    if (var %in% names(df_scaled)) {
      df_scaled[[var]] <- as.numeric(scale(df_scaled[[var]]))
    }
  }
  
  return(df_scaled)
}


#### Multiple Imputation with Mice ####

perform_imputation <- function(df, m = 50, maxit = 20, seed = 42) {
  
  cat("\n=== STARTING MULTIPLE IMPUTATION ===\n")
  cat("Number of imputations:", m, "\n")
  cat("Max iterations:", maxit, "\n\n")
  
  # Remove variables with 100% missing
  vars_to_remove <- names(df)[colSums(is.na(df)) == nrow(df)]
  if (length(vars_to_remove) > 0) {
    cat("Removing variables with 100% missing:", paste(vars_to_remove, collapse = ", "), "\n")
    df <- df %>% select(-all_of(vars_to_remove))
  }
  
  # Check predictor matrix - MICE handles this automatically but can customize
  # For block-missing patterns, predictive mean matching (pmm) works well
  
  # Set up imputation methods
  # pmm = predictive mean matching (good for continuous, preserves distribution)
  # logreg = logistic regression (for binary)
  # polyreg = polytomous regression (for factors with >2 levels)
  
  init <- mice(df, maxit = 0)
  meth <- init$method
  pred <- init$predictorMatrix
  
  # For continuous variables with zeros, pmm is usually fine
  # meth["variable_name"] <- "pmm"
  
  # Remove outcomes as predictors of each other during imputation
  outcomes_in_df <- intersect(outcomes, names(df))
  for (out in outcomes_in_df) {
    pred[out, outcomes_in_df] <- 0
  }
  
  # Run imputation
  imp <- mice(df, 
              m = m, 
              maxit = maxit, 
              method = meth,
              predictorMatrix = pred,
              seed = seed,
              printFlag = TRUE)
  
  # Diagnostics
  cat("\n=== IMPUTATION DIAGNOSTICS ===\n")
  
  # Convergence plots
  p_conv <- plot(imp)
  
  # Density plots comparing imputed to observed
  p_dens <- densityplot(imp)
  
  # Strip plots
  p_strip <- stripplot(imp)
  
  return(list(
    mids = imp,
    convergence_plot = p_conv,
    density_plot = p_dens,
    strip_plot = p_strip
  ))
}


#### Load and Merge Soil PCs from Greenhouse_SoilChemistry_PCA_forBayesian.R Script ####

# This section loads soil PCA results from Greenhouse_SoilChemistry_PCA_forBayesian.R
# and merges the PC scores with the data here 
#
# The soil PCA was computed using multiple imputation -> PCA -> averaging
# This properly handles the 83% missingness in soil chemistry, except pH which will be included separately 

load_soil_pcs <- function(soil_pca_result) {
  # Extract the averaged PC scores from the soil PCA result object
  # soil_pca_result should be the output from run_soil_pca_workflow()
  
  pc_scores <- soil_pca_result$primary$avg_scores
  
  cat("\n=== SOIL PC SCORES LOADED ===\n")
  cat("Number of samples:", nrow(pc_scores), "\n")
  cat("Number of PCs:", ncol(pc_scores), "\n")
  cat("PC names:", paste(names(pc_scores), collapse = ", "), "\n")
  
  return(pc_scores)
}

# Function to add soil PCs to a dataframe
add_soil_pcs_to_data <- function(df, soil_pc_scores) {
  # Check row counts match
  if (nrow(df) != nrow(soil_pc_scores)) {
    stop("Row counts don't match! df: ", nrow(df), " rows; soil_pc_scores: ", nrow(soil_pc_scores), " rows")
  }
  
  # Check row names match (if both have row names)
  if (!is.null(rownames(df)) && !is.null(rownames(soil_pc_scores))) {
    if (!identical(rownames(df), rownames(soil_pc_scores))) {
      # Check if they're the same but in different order
      if (setequal(rownames(df), rownames(soil_pc_scores))) {
        warning("Row names match but are in different order. Reordering soil_pc_scores to match df.")
        soil_pc_scores <- soil_pc_scores[rownames(df), , drop = FALSE]
      } else {
        stop("Row names don't match between df and soil_pc_scores!")
      }
    }
    cat("Row names verified: all", nrow(df), "rows match.\n")
  } else {
    warning("Row names not available on one or both objects. Assuming rows are in same order.")
  }
  
  df_with_pcs <- cbind(df, soil_pc_scores)
  return(df_with_pcs)
}


#### Specify Priors ####

# Regularized horseshoe prior for variable selection
# The horseshoe prior:
#   - Shrinks SMALL effects strongly toward zero (noise)
#   - Leaves LARGE effects relatively untouched (signal)
#   - Automatically identifies which is which

# PARAMETERS:
# p = total number of predictors
# p0 = expected number of relevant predictors (we guess ~1/3)
# tau_scale = global shrinkage parameter (computed from p0/p and n)

get_horseshoe_prior <- function(p, p0 = NULL, tau_scale = NULL) {
  if (is.null(p0)) p0 <- ceiling(p / 3)
  if (is.null(tau_scale)) tau_scale <- p0 / (p - p0) / sqrt(103)
  
  prior <- c(
    set_prior(paste0("horseshoe(df = 1, scale_global = ", tau_scale, ", df_global = 1)"), class = "b"),
    set_prior("student_t(3, 0, 2.5)", class = "Intercept"), # standard "weakly informative" prior, heavy tails allow for large values if data supports them
    set_prior("student_t(3, 0, 2.5)", class = "sigma")
  )
  
  return(prior)
}

# Standard weakly informative priors (for comparison or simpler models)
get_standard_prior <- function() {
  prior <- c(
    prior(normal(0, 1), class = "b"),           # Coefficients (for scaled predictors)
    prior(student_t(3, 0, 2.5), class = "Intercept"),
    prior(student_t(3, 0, 2.5), class = "sigma")
  )
  return(prior)
}


#### Model Fitting Functions ####

# --- Fit single outcome model with imputed data ---
fit_brms_with_imputation <- function(imp, formula, family = gaussian(), 
                                     prior = NULL, iter = 4000, warmup = 2000,
                                     chains = 4, cores = 4, seed = 42,
                                     control = list(adapt_delta = 0.95)) {
  
  # brms can handle mids objects directly with brm_multiple()
  fit <- brm_multiple(
    formula = formula,
    data = imp,
    family = family,
    prior = prior,
    iter = iter,
    warmup = warmup,
    chains = chains,
    cores = cores,
    seed = seed,
    control = control,
    combine = TRUE,  # Combine across imputations
    file = NULL      # Don't save to file
  )
  
  return(fit)
}

# --- Alternative: Loop through imputations manually ---

# --- Pool results manually (Rubin's rules) ---
pool_brms_results <- function(fits) {
  # Extract posterior samples from each fit
  # This is a simplified pooling - brm_multiple does this automatically
  
  all_draws <- map(fits, ~as_draws_df(.))
  
  # Combine draws
  combined_draws <- bind_rows(all_draws, .id = "imputation")
  
  # Summary
  summary_pooled <- combined_draws %>%
    select(-imputation, -.chain, -.iteration, -.draw) %>%
    summarise(across(everything(), list(
      mean = mean,
      sd = sd,
      q2.5 = ~quantile(., 0.025),
      q97.5 = ~quantile(., 0.975)
    )))
  
  return(list(draws = combined_draws, summary = summary_pooled))
}


#### Model Building for Individual Outcomes (5 Plant Outcomes) ####

# Build formula for each outcome
# Using horseshoe prior for automatic variable selection

build_outcome_formulas <- function(include_soil_pcs = TRUE, n_soil_pcs = 3) {
  
  # Base predictors (treatment only - texture is confounded with treatment)
  # NOTE: MPB_Impacted excluded because it's perfectly nested within Treatment_LandUseHistory
  base_preds <- c("Treatment_LandUseHistory")
  
  # NOTE: Texture NOT included - confounded with treatment (each treatment = 1 site = 1 texture)
  
  # Soil PCs (pre-computed from Greenhouse_SoilChemistry_PCA_forBayesian.R) + pH as separate predictor
  soil_pc_preds <- if (include_soil_pcs) c(paste0("Soil_PC", 1:n_soil_pcs), "pH") else NULL
  
  # Microbial predictors (shared structure)
  rhizo_preds <- c("Fungal_Shannons_Rhizo", "Fungal_ObservedFeatures_Rhizo",
                   "Fungal_EMF_RA_Rhizo", "Fungal_EMF_ObservedFeatures_Rhizo",
                   "Fungal_PlantPathogen_RA_Rhizo", "Fungal_FungalParasite_RA_Rhizo")
  
  roottip_preds <- c("Fungal_Shannons_RootTips", "Fungal_ObservedFeatures_RootTips",
                     "Fungal_EMF_RA_RootTips", "Fungal_EMF_ObservedFeatures_RootTips",
                     "Fungal_PlantPathogen_RA_RootTips", "Fungal_FungalParasite_RA_RootTips")
  
  # NOTE: Only Shannons for prokaryotes (ObservedFeatures dropped due to r=0.89 correlation)
  prokaryote_preds <- c("Prokayote_Shannons_Rhizo")
  
  # Outcome-specific VIP predictors
  # Using ONLY mean (not sum - they're perfectly correlated)
  
  formulas <- list()
  
  # HEIGHT
  height_vip <- c("vip_mean_ITS_rhizosphere_height_positiveassociations",
                  "vip_mean_ITS_rhizosphere_height_negativeassociations",
                  "vip_mean_ITS_roottips_height_positiveassociations",
                  "vip_mean_ITS_roottips_height_negativeassociations")
  
  height_preds <- c(base_preds, soil_pc_preds, 
                    rhizo_preds, roottip_preds, prokaryote_preds, height_vip)
  formulas$Height <- as.formula(paste("Height ~", paste(height_preds, collapse = " + ")))
  
  # MEAN CANOPY WIDTH
  canwidth_vip <- c("vip_mean_ITS_rhizosphere_meancanwidth_positiveassociations",
                    "vip_mean_ITS_rhizosphere_meancanwidth_negativeassociations",
                    "vip_mean_ITS_roottips_meancanwidth_positiveassociations")
                    # Note: vip_mean_ITS_roottips_meancanwidth_negativeassociations is 100% missing
  
  canwidth_preds <- c(base_preds, soil_pc_preds,
                      rhizo_preds, roottip_preds, prokaryote_preds, canwidth_vip)
  formulas$MeanCanWidth <- as.formula(paste("MeanCanWidth ~", paste(canwidth_preds, collapse = " + ")))
  
  # ABOVEGROUND BIOMASS
  biomass_vip <- c("vip_mean_ITS_rhizosphere_abovegroundbiomass_positiveassociations",
                   "vip_mean_ITS_rhizosphere_abovegroundbiomass_negativeassociations",
                   "vip_mean_ITS_roottips_abovegroundbiomass_positiveassociations")
                   # Note: vip_mean_ITS_roottips_abovegroundbiomass_negativeassociations is 100% missing
  
  biomass_preds <- c(base_preds, soil_pc_preds,
                     rhizo_preds, roottip_preds, prokaryote_preds, biomass_vip)
  formulas$AbovegroundBiomass <- as.formula(paste("AbovegroundBiomass ~", paste(biomass_preds, collapse = " + ")))
  
  # ROOT DRY MASS
  rootmass_vip <- c("vip_mean_ITS_rhizosphere_rootdrymass_positiveassociations",
                    "vip_mean_ITS_rhizosphere_rootdrymass_negativeassociations",
                    "vip_mean_ITS_roottips_rootdrymass_positiveassociations",
                    "vip_mean_ITS_roottips_rootdrymass_negativeassociations")
  
  rootmass_preds <- c(base_preds, soil_pc_preds,
                      rhizo_preds, roottip_preds, prokaryote_preds, rootmass_vip)
  formulas$RootDryMass <- as.formula(paste("RootDryMass ~", paste(rootmass_preds, collapse = " + ")))
  
  # ROOT TO SHOOT RATIO
  ratio_vip <- c("vip_mean_ITS_rhizosphere_roottoshootratio_positiveassociations",
                 "vip_mean_ITS_rhizosphere_roottoshootratio_negativeassociations",
                 "vip_mean_ITS_roottips_roottoshootratio_positiveassociations",
                 "vip_mean_ITS_roottips_roottoshootratio_negativeassociations")
  
  ratio_preds <- c(base_preds, soil_pc_preds,
                   rhizo_preds, roottip_preds, prokaryote_preds, ratio_vip)
  formulas$RootToShootRatio <- as.formula(paste("RootToShootRatio ~", paste(ratio_preds, collapse = " + ")))
  
  # Print predictor counts
  cat("\n=== PREDICTOR COUNTS PER OUTCOME ===\n")
  for (out in names(formulas)) {
    n_preds <- length(attr(terms(formulas[[out]]), "term.labels"))
    cat(out, ":", n_preds, "predictors\n")
  }
  
  return(formulas)
}


# FINAL PREDICTOR COUNTS:
#
# For each outcome:
#   - Treatment: 1 variable (creates 5 dummy variables for 6 levels)
#   - Soil PCs: 3 (PC1, PC2, PC3)
#   - Rhizo fungal: 6
#   - RootTip fungal: 6
#   - Prokaryote: 1
#   - VIPs: 3-4 (outcome-specific)
#
# Total: ~20-21 predictors per outcome model
# With n=103 samples, this is feasible with horseshoe regularization


#### Variable Selection ####

# After fitting with horseshoe prior, identify important variables
# Variables with credible intervals excluding zero are "selected"

# INTERPRETING prob_nonzero:
#
# prob_nonzero = 0.95 means:
#   - 95% of posterior samples are on ONE side of zero
#   - Strong evidence effect is non-zero
#   - Variable is "selected"
#
# prob_nonzero = 0.55 means:
#   - Only 55% on one side
#   - Nearly as likely to be positive as negative
#   - Effect is uncertain, likely noise
#   - Variable is NOT selected (horseshoe shrunk it toward 0)

# Version for pooled results across multiple imputations
extract_important_variables_pooled <- function(pooled_result, threshold = 0.9) {
  # Extract from pooled draws instead of single fit
  # This properly accounts for imputation uncertainty
  
  post_summary <- pooled_result$draws %>%
    select(starts_with("b_")) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    group_by(parameter) %>%
    summarise(
      mean = mean(value),
      sd = sd(value),
      q2.5 = quantile(value, 0.025),
      q10 = quantile(value, 0.10),
      q90 = quantile(value, 0.90),
      q97.5 = quantile(value, 0.975),
      prob_positive = mean(value > 0),
      prob_negative = mean(value < 0),
      prob_nonzero = pmax(prob_positive, prob_negative)
    ) %>%
    mutate(
      significant_90 = (q10 > 0 | q90 < 0),  # 80% CI excludes zero
      significant_95 = (q2.5 > 0 | q97.5 < 0),  # 95% CI excludes zero
      direction = case_when(
        prob_positive > threshold ~ "positive",
        prob_negative > threshold ~ "negative",
        TRUE ~ "uncertain"
      )
    ) %>%
    arrange(desc(prob_nonzero))
  
  # Remove intercept from selection
  post_summary <- post_summary %>%
    filter(parameter != "b_Intercept")
  
  return(post_summary)
}


#### Model Diagnostics ####

diagnose_model <- function(fit) {
  
  diagnostics <- list()
  
  # --- MCMC diagnostics ---
  cat("\n=== MCMC DIAGNOSTICS ===\n")
  
  # Rhat (should be < 1.01)
  rhat_vals <- rhat(fit)
  cat("Rhat range:", range(rhat_vals, na.rm = TRUE), "\n")
  if (any(rhat_vals > 1.01, na.rm = TRUE)) {
    warning("Some Rhat values > 1.01 - chains may not have converged!")
  }
  
  # ESS (effective sample size - should be > 400)
  ess_bulk <- neff_ratio(fit)
  cat("ESS ratio range:", range(ess_bulk, na.rm = TRUE), "\n")
  
  # Trace plots
  diagnostics$trace <- mcmc_trace(fit, regex_pars = "b_")
  
  # Rank plots (alternative to trace)
  diagnostics$rank <- mcmc_rank_overlay(fit, regex_pars = "b_")
  
  # --- Posterior predictive checks ---
  cat("\n=== POSTERIOR PREDICTIVE CHECKS ===\n")
  
  diagnostics$pp_check <- pp_check(fit, ndraws = 100)
  diagnostics$pp_stat <- pp_check(fit, type = "stat", stat = "mean")
  diagnostics$pp_stat_2d <- pp_check(fit, type = "stat_2d", stat = c("mean", "sd"))
  
  # --- LOO-CV ---
  cat("\n=== LOO-CV ===\n")
  diagnostics$loo <- loo(fit)
  print(diagnostics$loo)
  
  # Check for problematic observations
  if (any(diagnostics$loo$diagnostics$pareto_k > 0.7)) {
    warning("Some Pareto k values > 0.7 - LOO estimates may be unreliable for these observations")
    problematic <- which(diagnostics$loo$diagnostics$pareto_k > 0.7)
    cat("Problematic observations:", problematic, "\n")
  }
  
  return(diagnostics)
}


#### Figures ####

# --- Coefficient plot for POOLED results ---
plot_coefficients_pooled <- function(pooled_result, outcome_name = "", exclude_intercept = TRUE) {
  
  # Extract posterior summary from pooled draws
  post <- pooled_result$draws %>%
    select(starts_with("b_")) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value")
  
  # Calculate summaries
  post_summary <- post %>%
    group_by(parameter) %>%
    summarise(
      mean = mean(value),
      q2.5 = quantile(value, 0.025),
      q10 = quantile(value, 0.10),
      q90 = quantile(value, 0.90),
      q97.5 = quantile(value, 0.975),
      prob_positive = mean(value > 0)
    ) %>%
    mutate(
      # Clean parameter names
      parameter_clean = str_remove(parameter, "^b_"),
      # Determine significance
      significant_95 = (q2.5 > 0 | q97.5 < 0),
      significant_80 = (q10 > 0 | q90 < 0)
    )
  
  if (exclude_intercept) {
    post_summary <- post_summary %>% 
      filter(parameter_clean != "Intercept")
  }
  
  # Order by effect size
  post_summary <- post_summary %>%
    mutate(parameter_clean = fct_reorder(parameter_clean, mean))
  
  # Create plot
  p <- ggplot(post_summary, aes(x = mean, y = parameter_clean)) +
    # 95% CI
    geom_errorbarh(aes(xmin = q2.5, xmax = q97.5), height = 0, linewidth = 0.5, color = "gray40") +
    # 80% CI
    geom_errorbarh(aes(xmin = q10, xmax = q90), height = 0, linewidth = 1.2, color = "gray20") +
    # Point estimate
    geom_point(aes(color = significant_95), size = 2.5) +
    # Zero line
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
    # Colors
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "steelblue"),
                       labels = c("95% CI includes 0", "95% CI excludes 0"),
                       name = "Significance") +
    # Theme
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      title = paste("Coefficient Estimates:", outcome_name, "(Pooled)"),
      subtitle = "Pooled across multiple imputations",
      x = "Standardized Effect Size",
      y = ""
    )
  
  return(p)
}

# --- Variable importance plot ---
plot_variable_importance <- function(selection_results, outcome_name = "") {
  
  # Plot probability that effect is non-zero
  p <- selection_results %>%
    mutate(parameter = str_remove(parameter, "^b_")) %>%
    mutate(parameter = fct_reorder(parameter, prob_nonzero)) %>%
    ggplot(aes(x = prob_nonzero, y = parameter, fill = direction)) +
    geom_col() +
    geom_vline(xintercept = 0.9, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_fill_manual(values = c("positive" = "#2166AC", 
                                 "negative" = "#B2182B", 
                                 "uncertain" = "gray60")) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      title = paste("Variable Importance:", outcome_name),
      subtitle = "Probability that effect is non-zero (dashed line = 90% threshold)",
      x = "Probability",
      y = "",
      fill = "Direction"
    )
  
  return(p)
}


#### Results Tables ####

# Version for pooled results across multiple imputations
create_results_table_pooled <- function(pooled_result, outcome_name = "") {
  
  # Determine if this outcome uses lognormal family
  is_lognormal <- outcome_name %in% c("AbovegroundBiomass", "RootDryMass", "RootToShootRatio")
  
  # Extract and format results from pooled draws
  results <- pooled_result$draws %>%
    select(starts_with("b_")) %>%
    pivot_longer(everything(), names_to = "Parameter", values_to = "value") %>%
    group_by(Parameter) %>%
    summarise(
      Estimate = mean(value),
      SE = sd(value),
      `2.5%` = quantile(value, 0.025),
      `97.5%` = quantile(value, 0.975),
      `Pr(>0)` = mean(value > 0),
      `Pr(nonzero)` = max(mean(value > 0), mean(value < 0))
    ) %>%
    mutate(
      Parameter = str_remove(Parameter, "^b_"),
      across(where(is.numeric), ~round(., 3)),
      Outcome = outcome_name,
      Family = ifelse(is_lognormal, "lognormal", "gaussian")
    )
  
  # For lognormal outcomes, add fold-change interpretation
  if (is_lognormal) {
    results <- results %>%
      mutate(
        `Fold Change` = round(exp(Estimate), 3),
        `FC 2.5%` = round(exp(`2.5%`), 3),
        `FC 97.5%` = round(exp(`97.5%`), 3)
      ) %>%
      select(Outcome, Family, Parameter, Estimate, SE, `2.5%`, `97.5%`, `Pr(>0)`, `Pr(nonzero)`,
             `Fold Change`, `FC 2.5%`, `FC 97.5%`)
  } else {
    results <- results %>%
      select(Outcome, Family, Parameter, Estimate, SE, `2.5%`, `97.5%`, `Pr(>0)`, `Pr(nonzero)`)
  }
  
  return(results)
}


#### Main Workflow ####

run_full_analysis <- function(df, soil_pca_result, m_imputations = 50) {
  
  cat("\n")
  cat("================================================================\n")
  cat("       BAYESIAN ANALYSIS OF GREENHOUSE PLANT OUTCOMES          \n")
  cat("================================================================\n\n")
  
  # --- Step 1: Preprocessing ---
  cat("STEP 1: Preprocessing data...\n")
  
  # Check distributions and determine transformations
  # (Add log transforms as needed based on distribution checks)
  df_processed <- preprocess_data(df)
  
  # --- Step 2: Add pre-computed soil PCs ---
  cat("\nSTEP 2: Adding pre-computed soil PCs from Greenhouse_SoilChemistry_PCA_forBayesian.R...\n")
  
  soil_pc_scores <- load_soil_pcs(soil_pca_result)
  df_with_pcs <- add_soil_pcs_to_data(df_processed, soil_pc_scores)
  
  # --- Step 3: Imputation (for microbial and VIP variables) ---
  cat("\nSTEP 3: Performing multiple imputation (microbial/VIP variables)...\n")
  cat("        (Soil chemistry already imputed during PCA step)\n")
  
  imp_result <- perform_imputation(df_with_pcs, m = m_imputations)
  imp <- imp_result$mids
  
  # --- Step 4: Scale predictors in each imputed dataset ---
  cat("\nSTEP 4: Scaling predictors in imputed datasets...\n")
  
  imputed_dfs <- list()
  
  # Define continuous variables to scale (excluding texture - confounded with treatment)
  continuous_vars <- c(rhizo_fungal_vars, roottip_fungal_vars,
                       prokaryote_vars, paste0("Soil_PC", 1:3), "pH")
  # Add VIP vars
  all_vip <- unique(c(vip_vars_height, vip_vars_meancanwidth, vip_vars_abovegroundbiomass,
                      vip_vars_rootdrymass, vip_vars_roottoshootratio))
  continuous_vars <- c(continuous_vars, all_vip)
  
  for (i in 1:m_imputations) {
    df_i <- complete(imp, i)
    
    # Scale only variables that exist in the data
    continuous_vars_exist <- continuous_vars[continuous_vars %in% names(df_i)]
    df_i <- scale_predictors(df_i, continuous_vars_exist)
    
    imputed_dfs[[i]] <- df_i
  }
  
  # --- Step 5: Build formulas ---
  cat("\nSTEP 5: Building model formulas...\n")
  
  formulas <- build_outcome_formulas(include_soil_pcs = TRUE, n_soil_pcs = 3)
  
  # --- Step 6: Fit individual models with horseshoe priors ---
  cat("\nSTEP 6: Fitting individual outcome models with horseshoe priors...\n")
  cat("        (This may take a while...)\n\n")
  cat("        Family specifications:\n")
  cat("          Height, MeanCanWidth: Gaussian\n")
  cat("          AbovegroundBiomass, RootDryMass, RootToShootRatio: Lognormal\n\n")
  
  fits <- list()              # Store first fit for diagnostics
  pooled_results <- list()    # Store pooled results for inference
  selection_results <- list()
  
  for (outcome in names(formulas)) {
    cat("  Fitting model for:", outcome, "\n")
    
    # Get the appropriate family for this outcome
    outcome_family <- outcome_families[[outcome]]
    cat("    Using family:", outcome_family$family, "\n")
    
    # Count predictors for horseshoe prior
    n_pred <- length(attr(terms(formulas[[outcome]]), "term.labels"))
    prior_hs <- get_horseshoe_prior(p = n_pred)
    
    # Fit across all imputations
    fits_imp <- list()
    for (i in 1:m_imputations) {
      cat("      Imputation", i, "of", m_imputations, "\r")
      fits_imp[[i]] <- brm(
        formula = formulas[[outcome]],
        data = imputed_dfs[[i]],
        family = outcome_family,  # Use outcome-specific family
        prior = prior_hs,
        iter = 4000,
        warmup = 2000,
        chains = 4,
        cores = 4,
        seed = 42 + i,
        control = list(adapt_delta = 0.99, max_treedepth = 12),
        silent = 2
      )
    }
    cat("\n")
    
    # Pool results across all imputations
    pooled <- pool_brms_results(fits_imp)
    pooled_results[[outcome]] <- pooled
    fits[[outcome]] <- fits_imp[[1]]  # Store first fit for diagnostics only
    
    # Variable selection - USE POOLED RESULTS
    selection_results[[outcome]] <- extract_important_variables_pooled(pooled)
    
    cat("    Complete. Top predictors (pooled across", m_imputations, "imputations):\n")
    print(head(selection_results[[outcome]][, c("parameter", "mean", "prob_nonzero")], 5))
    cat("\n")
  }
  
  # --- Step 7: Identify selected variables ---
  cat("\nSTEP 7: Selecting variables for final model...\n")
  
  selected_vars <- map(selection_results, function(res) {
    # Select variables with >90% probability of non-zero effect
    selected <- res %>%
      filter(prob_nonzero > 0.9) %>%
      pull(parameter) %>%
      str_remove("^b_")
    
    # Always include treatments
    selected <- union(selected, treatments)
    
    return(selected)
  })
  
  cat("\nSelected variables by outcome:\n")
  for (out in names(selected_vars)) {
    cat(" ", out, ":", length(selected_vars[[out]]), "variables\n")
    cat("   ", paste(selected_vars[[out]], collapse = ", "), "\n")
  }
  
  # --- Step 8: Create figures ---
  cat("\nSTEP 8: Creating publication figures...\n")
  
  # Use pooled results for coefficient plots
  coef_plots <- map2(pooled_results, names(pooled_results), 
                     ~plot_coefficients_pooled(.x, .y))
  importance_plots <- map2(selection_results, names(selection_results), 
                           ~plot_variable_importance(.x, .y))
  
  # Save individual plots
  for (out in names(coef_plots)) {
    ggsave(paste0("coef_", out, ".pdf"), coef_plots[[out]], width = 10, height = 8)
    ggsave(paste0("importance_", out, ".pdf"), importance_plots[[out]], width = 10, height = 8)
  }
  cat("  Saved coefficient and importance plots for each outcome.\n")
  
  # Combined figure
  # pub_fig <- create_publication_figure(coef_plots, importance_plots, pca_result$plots)
  # ggsave("publication_figure.pdf", pub_fig, width = 16, height = 20, dpi = 300)
  
  # --- Step 9: Create results tables ---
  cat("\nSTEP 9: Creating results tables...\n")
  
  # Use pooled results for the final table
  results_table <- map2_dfr(pooled_results, names(pooled_results), 
                            ~create_results_table_pooled(.x, .y))
  write.csv(results_table, "model_results.csv", row.names = FALSE)
  cat("  Saved: model_results.csv\n")
  
  # --- Step 10: Run diagnostics (on first fit only) ---
  cat("\nSTEP 10: Running model diagnostics...\n")
  
  diagnostics <- map(fits, diagnose_model)
  
  # Save diagnostic plots
  for (out in names(diagnostics)) {
    ggsave(paste0("ppcheck_", out, ".pdf"), diagnostics[[out]]$pp_check, width = 8, height = 6)
    ggsave(paste0("trace_", out, ".pdf"), diagnostics[[out]]$trace, width = 12, height = 8)
  }
  
  cat("\n")
  cat("================================================================\n")
  cat("                    ANALYSIS COMPLETE                           \n")
  cat("================================================================\n")
  
  return(list(
    imputation = imp_result,
    fits = fits,                    # First imputation fits (for diagnostics)
    pooled_results = pooled_results, # Pooled across all imputations (for inference)
    selection = selection_results,
    selected_vars = selected_vars,
    diagnostics = diagnostics,
    results_table = results_table
  ))
}


#### Run Workflow ####

# 1. First, run the Greenhouse_SoilChemistry_PCA_forBayesian.R to create soil chemistry PCs:
#    
    # source("Greenhouse_SoilChemistry_PCA_forBayesian.R")
    setwd()
    df <- read.delim("Bayesian_AllData_ExcludedThoseWithOnlyNAsAcross_RemovedProblemColumns.txt", row.names = 1)
    soil_results <- soil_pca_result ## Already ran, just assigning to the name used in this script
#
# 2. Then, run this script with the soil PCA results:
#
    # source("greenhouse_bayesian_analysis_FINAL.R")
    results <- run_full_analysis(df, soil_results, m_imputations = 50)
#
# 3. Access results:
    results$selection$Height          # Variable selection for Height
    results$results_table             # Combined results table (pooled)
    results$pooled_results$Height     # Pooled posterior draws for Height
    results$fits$Height               # First imputation fit (for diagnostics)
    results$diagnostics$Height        # Diagnostics for Height


    # Check number of rows in pooled draws
    nrow(results$pooled_results$Height$draws)
    # Should be ~400,000 (50 imputations × 8000 draws)
    
    # Compare to a single fit
    nrow(as_draws_df(results$fits$Height))
    # Should be ~8,000
    
    View(results$results_table)
    
    
#### Notes ####
    
#
# 1. Workflow Order
#    - Run Greenhouse_SoilChemistry_PCA_forBayesian.R FIRST with m=50 imputations
#    - This creates stable PC scores averaged across imputations
#    - Then run this script with the soil_results object
#
# 2. Confounding:
#    - Treatment_LandUseHistory is confounded with Site and Texture
#    - Each treatment used soil from a specific site
#    - Cannot statistically separate land use effects from site effects
#    - Texture variables excluded from regression for this reason
#    - MPB_Impacted excluded (perfectly nested within Treatment)
#    - Discuss this limitation in your methods/discussion
#
# 3. Variables to Exclude:
#    - MPB_Impacted: Perfectly nested within Treatment_LandUseHistory
#    - Texture: Confounded with treatment (each treatment = 1 site)
#    - Prokayote_ObservedFeatures: r=0.89 with Shannons, redundant
#    - All vip_sum_*: r=1.0 with vip_mean_*, using mean only
#    - Two VIP variables: 100% missing, excluded
#
# 4. Outcome Distributions:
#    - Height, MeanCanWidth: Gaussian family (symmetric)
#    - AbovegroundBiomass, RootDryMass, RootToShootRatio: Lognormal family (right-skewed)
#
# 5. Lognormal Interpretation:
#    - Coefficients are on log scale
#    - exp(beta) gives fold-change per 1 SD predictor increase
#    - beta = 0.3 → exp(0.3) = 1.35 → 35% increase
#    - beta = -0.2 → exp(-0.2) = 0.82 → 18% decrease
#
# 6. Variable Selection:
#    - Horseshoe priors provide automatic variable selection
#    - Variables with >90% posterior probability of non-zero effect are "selected"
#    - This is more Bayesian than stepwise/projpred


    
    
    
    

    
    
    
    
    
