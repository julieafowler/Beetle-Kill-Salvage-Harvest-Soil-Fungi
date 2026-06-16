#### Multivariate Bayesian Models with Multiple Imputation ####
## Created with the assistance of Claude AI Opus 4.5

## Using the predictor selection from the univariate models 

## Two versions: Land-Use History (6-level treatment), Salvage_Harvest_Status (2-level treatment)

## Multiple imputations instead of single imputation
## Proper pooling of coefficients using Rubin's rules
## Pooling of residual correlations using Fisher-z transformation
## Hypothesis tests on pooled estimates

#### Packages & Seed #### 

library(tidyverse)
library(brms)
library(bayesplot)
library(tidybayes)
library(corrplot)
library(patchwork)
library(mice)

setwd()
set.seed(42)
options(mc.cores = parallel::detectCores() - 1)
options(brms.backend = "cmdstanr")


#### Set Number of Imputations ####

# Number of imputations for multivariate analysis
# Note: Univariate used m=50, but multivariate is more computationally expensive.
# m=20 provides good imputation uncertainty coverage (Rubin, 1987).
M_IMPUTATIONS <- 20  


#### Prepare Multiple Imputed Datasets ####

prepare_mv_datasets <- function(df, soil_pca_result, m = M_IMPUTATIONS) {
  #' Prepare m complete datasets for multivariate analysis
  #' 
  #' @param df Original data frame

  #' @param soil_pca_result Results from Greenhouse_SoilChemistry_PCA_forBayesian.R

  #' @param m Number of imputations
  #' @return List with mice object and list of complete datasets
  
  cat("Preparing", m, "imputed datasets for multivariate analysis...\n")
  
  # Outcomes
  outcomes <- c("Height", "MeanCanWidth", "AbovegroundBiomass", 
                "RootDryMass", "RootToShootRatio")
  
  # Variables to impute (same as main analysis)
  microbial_vars <- c(
    "Fungal_Shannons_Rhizo", "Fungal_ObservedFeatures_Rhizo",
    "Fungal_EMF_RA_Rhizo", "Fungal_EMF_ObservedFeatures_Rhizo",
    "Fungal_PlantPathogen_RA_Rhizo", "Fungal_FungalParasite_RA_Rhizo",
    "Fungal_Shannons_RootTips", "Fungal_ObservedFeatures_RootTips",
    "Fungal_EMF_RA_RootTips", "Fungal_EMF_ObservedFeatures_RootTips",
    "Fungal_PlantPathogen_RA_RootTips", "Fungal_FungalParasite_RA_RootTips",
    "Prokayote_Shannons_Rhizo"
  )
  
  vip_vars <- c(
    "vip_mean_ITS_rhizosphere_height_positiveassociations",
    "vip_mean_ITS_rhizosphere_height_negativeassociations",
    "vip_mean_ITS_roottips_height_positiveassociations",
    "vip_mean_ITS_roottips_height_negativeassociations",
    "vip_mean_ITS_rhizosphere_meancanwidth_positiveassociations",
    "vip_mean_ITS_rhizosphere_meancanwidth_negativeassociations",
    "vip_mean_ITS_roottips_meancanwidth_positiveassociations",
    "vip_mean_ITS_rhizosphere_abovegroundbiomass_positiveassociations",
    "vip_mean_ITS_rhizosphere_abovegroundbiomass_negativeassociations",
    "vip_mean_ITS_roottips_abovegroundbiomass_positiveassociations",
    "vip_mean_ITS_rhizosphere_rootdrymass_positiveassociations",
    "vip_mean_ITS_rhizosphere_rootdrymass_negativeassociations",
    "vip_mean_ITS_roottips_rootdrymass_positiveassociations",
    "vip_mean_ITS_roottips_rootdrymass_negativeassociations",
    "vip_mean_ITS_rhizosphere_roottoshootratio_positiveassociations",
    "vip_mean_ITS_rhizosphere_roottoshootratio_negativeassociations",
    "vip_mean_ITS_roottips_roottoshootratio_positiveassociations",
    "vip_mean_ITS_roottips_roottoshootratio_negativeassociations"
  )
  
  # Filter to existing vars
  microbial_vars <- microbial_vars[microbial_vars %in% names(df)]
  vip_vars <- vip_vars[vip_vars %in% names(df)]
  
  # Select variables for imputation
  vars_to_use <- c(outcomes, "Treatment_LandUseHistory", "Salvage_Harvest_Status",
                   "Texture_PercentClay", "Texture_PercentSilt", "Texture_PercentSand",
                   "pH", microbial_vars, vip_vars)
  vars_to_use <- vars_to_use[vars_to_use %in% names(df)]
  
  df_sub <- df[, vars_to_use]
  
  # Convert factors
  df_sub$Treatment_LandUseHistory <- factor(df_sub$Treatment_LandUseHistory)
  df_sub$Treatment_LandUseHistory <- relevel(df_sub$Treatment_LandUseHistory, ref = "Old_Forest")
  df_sub$Salvage_Harvest_Status <- factor(df_sub$Salvage_Harvest_Status)
  df_sub$Salvage_Harvest_Status <- relevel(df_sub$Salvage_Harvest_Status, ref = "Non_Salvage_Harvested")
  
  # Store original rownames
  original_rownames <- rownames(df_sub)
  
  # Set up imputation - outcomes don't predict each other
  pred_matrix <- mice::quickpred(df_sub, minpuc = 0.1, mincor = 0.1)
  for (out in outcomes) {
    for (other_out in setdiff(outcomes, out)) {
      if (other_out %in% colnames(pred_matrix)) {
        pred_matrix[out, other_out] <- 0
      }
    }
  }
  
  # Run imputation with m datasets
  cat("  Running mice with m =", m, "imputations...\n")
  imp <- mice::mice(df_sub, m = m, maxit = 20, predictorMatrix = pred_matrix, 
                    method = "pmm", seed = 42, printFlag = FALSE)
  
  # Get soil PC scores
  soil_pc_scores <- soil_pca_result$primary$avg_scores
  
  # Create list of complete datasets
  datasets <- list()
  
  for (i in 1:m) {
    df_complete <- mice::complete(imp, i)
    rownames(df_complete) <- original_rownames
    
    # Add soil PCs
    common_rows <- intersect(rownames(df_complete), rownames(soil_pc_scores))
    
    if (length(common_rows) > 0) {
      df_complete <- df_complete[common_rows, ]
      df_complete$Soil_PC1 <- soil_pc_scores[common_rows, "Soil_PC1"]
      df_complete$Soil_PC2 <- soil_pc_scores[common_rows, "Soil_PC2"]
      df_complete$Soil_PC3 <- soil_pc_scores[common_rows, "Soil_PC3"]
    } else if (nrow(df_complete) == nrow(soil_pc_scores)) {
      df_complete$Soil_PC1 <- soil_pc_scores[, "Soil_PC1"]
      df_complete$Soil_PC2 <- soil_pc_scores[, "Soil_PC2"]
      df_complete$Soil_PC3 <- soil_pc_scores[, "Soil_PC3"]
    } else {
      stop("Cannot match soil PC scores to imputed data!")
    }
    
    # Scale continuous predictors
    continuous_vars <- c("pH", "Soil_PC1", "Soil_PC2", "Soil_PC3", 
                         microbial_vars, vip_vars)
    continuous_vars <- continuous_vars[continuous_vars %in% names(df_complete)]
    
    for (var in continuous_vars) {
      df_complete[[var]] <- scale(df_complete[[var]])[, 1]
    }
    
    # Create log-transformed versions for multivariate model
    df_complete$log_AbovegroundBiomass <- log(df_complete$AbovegroundBiomass)
    df_complete$log_RootDryMass <- log(df_complete$RootDryMass)
    df_complete$log_RootToShootRatio <- log(df_complete$RootToShootRatio)
    
    datasets[[i]] <- df_complete
  }
  
  cat("  Created", m, "complete datasets with n =", nrow(datasets[[1]]), "each\n")
  
  return(list(
    mice_object = imp,
    datasets = datasets,
    m = m
  ))
}


#### Fit Multivariate Model to Single Dataset ####

fit_mv_model_landuse <- function(data, imp_num = 1, verbose = TRUE) {
  #' Fit Land-Use History multivariate model to one imputed dataset
  
  if (verbose) {
    cat("  Fitting Land-Use History model (imputation", imp_num, ")...\n")
  }
  
  # Define formulas - all Gaussian to allow rescor
  bf_height <- bf(Height ~ Treatment_LandUseHistory)  # No Soil_PC2
  bf_width <- bf(MeanCanWidth ~ Treatment_LandUseHistory)  # No Soil_PC2
  bf_biomass <- bf(log_AbovegroundBiomass ~ Treatment_LandUseHistory + Soil_PC2)  # Soil_PC2 selected
  bf_root <- bf(log_RootDryMass ~ Treatment_LandUseHistory)  # No extras
  bf_ratio <- bf(log_RootToShootRatio ~ Treatment_LandUseHistory + 
                   vip_mean_ITS_rhizosphere_roottoshootratio_positiveassociations)  # VIP selected
  
  mv_fit <- brm(
    bf_height + bf_width + bf_biomass + bf_root + bf_ratio + set_rescor(TRUE),
    data = data,
    prior = prior(lkj(2), class = "rescor"),
    chains = 4,
    iter = 4000,
    warmup = 2000,
    cores = 4,
    control = list(adapt_delta = 0.95),
    seed = 42 + imp_num,  # Different seed per imputation
    silent = 2,
    refresh = 0
  )
  
  return(mv_fit)
}


fit_mv_model_mpb <- function(data, imp_num = 1, verbose = TRUE) {
  # Fit Salvage_Harvest_Status multivariate model to one imputed dataset
  
  if (verbose) {
    cat("  Fitting Salvage_Harvest_Status model (imputation", imp_num, ")...\n")
  }
  
  # Define formulas - all Gaussian to allow rescor
  bf_height <- bf(Height ~ Salvage_Harvest_Status)  # No Soil_PC2
  bf_width <- bf(MeanCanWidth ~ Salvage_Harvest_Status)  # No Soil_PC2
  bf_biomass <- bf(log_AbovegroundBiomass ~ Salvage_Harvest_Status + Soil_PC2)  # Soil_PC2 selected
  bf_root <- bf(log_RootDryMass ~ Salvage_Harvest_Status)  # No extras
  bf_ratio <- bf(log_RootToShootRatio ~ Salvage_Harvest_Status + Soil_PC2 + 
                   vip_mean_ITS_rhizosphere_roottoshootratio_positiveassociations)  # Both selected
  
  mv_fit <- brm(
    bf_height + bf_width + bf_biomass + bf_root + bf_ratio + set_rescor(TRUE),
    data = data,
    prior = prior(lkj(2), class = "rescor"),
    chains = 4,
    iter = 4000,
    warmup = 2000,
    cores = 4,
    control = list(adapt_delta = 0.95),
    seed = 42 + imp_num,
    silent = 2,
    refresh = 0
  )
  
  return(mv_fit)
}


#### Fit Models Across All Imputations ####

fit_all_mv_models <- function(datasets, model_type = "landuse") {
  #' Fit multivariate model to all imputed datasets
  #' 
  #' @param datasets List of imputed datasets
  #' @param model_type "landuse" or "mpb"
  #' @return List of fitted brms models
  
  m <- length(datasets)
  
  cat("\n================================================================\n")
  cat("  FITTING", toupper(model_type), "MODELS ACROSS", m, "IMPUTATIONS\n")
  cat("================================================================\n\n")
  
  fits <- list()
  
  for (i in 1:m) {
    cat("Imputation", i, "of", m, ":\n")
    
    if (model_type == "landuse") {
      fits[[i]] <- fit_mv_model_landuse(datasets[[i]], imp_num = i)
    } else {
      fits[[i]] <- fit_mv_model_mpb(datasets[[i]], imp_num = i)
    }
    
    cat("  Complete.\n\n")
  }
  
  return(fits)
}


#### Pool Coefficient Estimates (Rubin's Rules) ####

pool_mv_coefficients <- function(fits_list) {
  #' Pool coefficient estimates across imputations using Rubin's rules
  #' 
  #' @param fits_list List of brms model fits
  #' @return Data frame with pooled estimates
  
  cat("Pooling coefficient estimates using Rubin's rules...\n")
  
  m <- length(fits_list)
  
  # Extract fixed effects from each model
  all_coefs <- lapply(fits_list, function(fit) {
    fe <- fixef(fit, summary = TRUE)
    data.frame(
      parameter = rownames(fe),
      estimate = fe[, "Estimate"],
      se = fe[, "Est.Error"],
      stringsAsFactors = FALSE
    )
  })
  
  # Get all parameter names
  all_params <- unique(unlist(lapply(all_coefs, function(x) x$parameter)))
  
  # Pool each parameter
  pooled <- data.frame()
  
  for (param in all_params) {
    # Extract estimates and SEs for this parameter across imputations
    estimates <- sapply(all_coefs, function(x) {
      idx <- which(x$parameter == param)
      if (length(idx) > 0) x$estimate[idx] else NA
    })
    
    ses <- sapply(all_coefs, function(x) {
      idx <- which(x$parameter == param)
      if (length(idx) > 0) x$se[idx] else NA
    })
    
    # Remove NAs
    valid <- !is.na(estimates)
    estimates <- estimates[valid]
    ses <- ses[valid]
    m_valid <- length(estimates)
    
    if (m_valid == 0) next
    
    # Rubin's rules
    # Pooled estimate: mean of estimates
    Q_bar <- mean(estimates)
    
    # Within-imputation variance: mean of squared SEs
    U_bar <- mean(ses^2)
    
    # Between-imputation variance
    B <- var(estimates)
    
    # Total variance
    T_var <- U_bar + (1 + 1/m_valid) * B
    
    # Pooled SE
    pooled_se <- sqrt(T_var)
    
    # Degrees of freedom (Barnard-Rubin adjustment)
    if (B > 0) {
      lambda <- (1 + 1/m_valid) * B / T_var
      df_old <- (m_valid - 1) / lambda^2
      # Simplified: use large df approximation for CIs
      df <- max(df_old, 100)
    } else {
      df <- Inf
    }
    
    # 95% CI
    ci_low <- Q_bar - 1.96 * pooled_se
    ci_high <- Q_bar + 1.96 * pooled_se
    
    # Probability calculations (approximate from pooled posterior)
    # Using normal approximation
    prob_positive <- pnorm(0, mean = Q_bar, sd = pooled_se, lower.tail = FALSE)
    prob_nonzero <- max(prob_positive, 1 - prob_positive)
    
    # Parse outcome and predictor from parameter name
    parts <- str_split(param, "_", n = 2)[[1]]
    outcome_raw <- parts[1]
    predictor <- ifelse(length(parts) > 1, parts[2], "Intercept")
    
    # Clean outcome name
    outcome <- case_when(
      outcome_raw == "logAbovegroundBiomass" ~ "AbovegroundBiomass",
      outcome_raw == "logRootDryMass" ~ "RootDryMass",
      outcome_raw == "logRootToShootRatio" ~ "RootToShootRatio",
      TRUE ~ outcome_raw
    )
    
    pooled <- bind_rows(pooled, data.frame(
      parameter = param,
      outcome = outcome,
      predictor = predictor,
      Estimate = Q_bar,
      SE = pooled_se,
      CI_low = ci_low,
      CI_high = ci_high,
      prob_positive = prob_positive,
      prob_nonzero = prob_nonzero,
      m_imputations = m_valid,
      within_var = U_bar,
      between_var = B,
      stringsAsFactors = FALSE
    ))
  }
  
  # Sort by outcome and effect size
  pooled <- pooled %>%
    arrange(outcome, desc(abs(Estimate)))
  
  cat("  Pooled", nrow(pooled), "parameters across", m, "imputations\n")
  
  return(pooled)
}


#### Pool Residuel Correlations (Fisher-Z Transformation) ####

pool_residual_correlations <- function(fits_list) {
  #' Pool residual correlations across imputations using Fisher-z transformation
  #' 
  #' @param fits_list List of brms model fits
  #' @return List with summary and correlation matrix
  
  cat("Pooling residual correlations using Fisher-z transformation...\n")
  
  m <- length(fits_list)
  
  # Extract rescor parameters from each model
  all_rescor <- lapply(fits_list, function(fit) {
    draws <- as_draws_df(fit) %>%
      select(starts_with("rescor__"))
    
    # Get posterior means and SDs
    data.frame(
      parameter = names(draws),
      mean = sapply(draws, mean),
      sd = sapply(draws, sd),
      stringsAsFactors = FALSE
    )
  })
  
  # Get all rescor parameter names
  all_params <- unique(unlist(lapply(all_rescor, function(x) x$parameter)))
  
  # Pool each correlation
  pooled <- data.frame()
  
  for (param in all_params) {
    # Extract correlations for this parameter
    correlations <- sapply(all_rescor, function(x) {
      idx <- which(x$parameter == param)
      if (length(idx) > 0) x$mean[idx] else NA
    })
    
    sds <- sapply(all_rescor, function(x) {
      idx <- which(x$parameter == param)
      if (length(idx) > 0) x$sd[idx] else NA
    })
    
    valid <- !is.na(correlations)
    correlations <- correlations[valid]
    sds <- sds[valid]
    m_valid <- length(correlations)
    
    if (m_valid == 0) next
    
    # Fisher-z transformation for pooling correlations
    # z = 0.5 * ln((1+r)/(1-r)) = atanh(r)
    z_values <- atanh(pmin(pmax(correlations, -0.999), 0.999))  # Bound to avoid Inf
    
    # Pool z values
    z_bar <- mean(z_values)
    z_var <- var(z_values)
    
    # Approximate within-imputation variance for z
    # SE(z) ≈ 1/sqrt(n-3), but we use the posterior SD transformed
    n <- 103  # Sample size
    within_var_z <- mean((sds / (1 - correlations^2))^2)  # Delta method approximation
    
    # Total variance
    total_var_z <- within_var_z + (1 + 1/m_valid) * z_var
    
    # Back-transform to correlation scale
    pooled_r <- tanh(z_bar)
    
    # CI on z scale, then back-transform
    z_ci_low <- z_bar - 1.96 * sqrt(total_var_z)
    z_ci_high <- z_bar + 1.96 * sqrt(total_var_z)
    ci_low <- tanh(z_ci_low)
    ci_high <- tanh(z_ci_high)
    
    # Significance: does CI exclude zero?
    significant <- (ci_low > 0 & ci_high > 0) | (ci_low < 0 & ci_high < 0)
    
    pooled <- bind_rows(pooled, data.frame(
      parameter = param,
      mean = pooled_r,
      sd = sqrt(total_var_z) * (1 - pooled_r^2),  # Approximate SD on r scale
      q2.5 = ci_low,
      q97.5 = ci_high,
      significant = significant,
      m_imputations = m_valid,
      stringsAsFactors = FALSE
    ))
  }
  
  # Create correlation matrix
  outcomes_display <- c("Height", "MeanCanWidth", "AbovegroundBiomass", 
                        "RootDryMass", "RootToShootRatio")
  
  cor_matrix <- matrix(1, nrow = 5, ncol = 5, 
                       dimnames = list(outcomes_display, outcomes_display))
  
  for (i in 1:nrow(pooled)) {
    param <- pooled$parameter[i]
    val <- pooled$mean[i]
    
    # Parse outcome names
    parts <- str_split(param, "__")[[1]]
    if (length(parts) >= 3) {
      o1_raw <- parts[2]
      o2_raw <- parts[3]
      
      o1 <- case_when(
        o1_raw == "logAbovegroundBiomass" ~ "AbovegroundBiomass",
        o1_raw == "logRootDryMass" ~ "RootDryMass",
        o1_raw == "logRootToShootRatio" ~ "RootToShootRatio",
        TRUE ~ o1_raw
      )
      o2 <- case_when(
        o2_raw == "logAbovegroundBiomass" ~ "AbovegroundBiomass",
        o2_raw == "logRootDryMass" ~ "RootDryMass",
        o2_raw == "logRootToShootRatio" ~ "RootToShootRatio",
        TRUE ~ o2_raw
      )
      
      if (o1 %in% outcomes_display & o2 %in% outcomes_display) {
        cor_matrix[o1, o2] <- val
        cor_matrix[o2, o1] <- val
      }
    }
  }
  
  cat("  Pooled", nrow(pooled), "correlations across", m, "imputations\n")
  
  return(list(
    summary = pooled,
    matrix = cor_matrix
  ))
}


#### Pooled Hypothesis Tests ####

test_pooled_hypotheses <- function(fits_landuse, fits_mpb) {
  #' Run hypothesis tests on pooled estimates across imputations
  #' 
  #' Strategy: Run hypothesis test on each imputed model, then pool results
  #' 
  #' @param fits_landuse List of Land-Use History model fits

  #' @param fits_mpb List of SH model fits
  #' @return List of hypothesis test results
  
  cat("\n================================================================\n")
  cat("  POOLED HYPOTHESIS TESTS ACROSS IMPUTATIONS\n")
  cat("================================================================\n\n")
  
  results <- list()
  m <- length(fits_mpb)
  
  # Helper function to pool hypothesis test results
  pool_hypothesis <- function(fits, hyp_string, hyp_name) {
    estimates <- numeric(length(fits))
    cis_low <- numeric(length(fits))
    cis_high <- numeric(length(fits))
    post_probs <- numeric(length(fits))
    
    for (i in seq_along(fits)) {
      tryCatch({
        h <- hypothesis(fits[[i]], hyp_string)
        estimates[i] <- h$hypothesis$Estimate
        cis_low[i] <- h$hypothesis$CI.Lower
        cis_high[i] <- h$hypothesis$CI.Upper
        post_probs[i] <- ifelse(!is.null(h$hypothesis$Post.Prob), 
                                 h$hypothesis$Post.Prob, NA)
      }, error = function(e) {
        estimates[i] <<- NA
        cis_low[i] <<- NA
        cis_high[i] <<- NA
        post_probs[i] <<- NA
      })
    }
    
    # Pool using Rubin's rules (simplified)
    valid <- !is.na(estimates)
    if (sum(valid) == 0) return(NULL)
    
    pooled_est <- mean(estimates[valid])
    between_var <- var(estimates[valid])
    within_var <- mean((cis_high[valid] - cis_low[valid])^2 / (2*1.96)^2)
    total_var <- within_var + (1 + 1/sum(valid)) * between_var
    pooled_se <- sqrt(total_var)
    
    pooled_ci_low <- pooled_est - 1.96 * pooled_se
    pooled_ci_high <- pooled_est + 1.96 * pooled_se
    pooled_prob <- mean(post_probs[valid], na.rm = TRUE)
    
    return(list(
      name = hyp_name,
      estimate = pooled_est,
      ci_low = pooled_ci_low,
      ci_high = pooled_ci_high,
      post_prob = pooled_prob,
      m_valid = sum(valid)
    ))
  }
  
  # --- Land-Use History Hypotheses ---
  cat("Land-Use History Model (pooled across", m, "imputations):\n")
  cat("-------------------------------------------------------\n")
  
  # H1: 1stPreMPB effect equal on Height vs Width?
  h1 <- pool_hypothesis(
    fits_landuse,
    "Height_Treatment_LandUseHistory1stPreMPB = MeanCanWidth_Treatment_LandUseHistory1stPreMPB",
    "H1: 1stPreMPB effect equal on Height vs Width"
  )
  if (!is.null(h1)) {
    cat("H1:", h1$name, "\n")
    cat("   Pooled difference:", round(h1$estimate, 3), 
        "95% CI [", round(h1$ci_low, 3), ",", round(h1$ci_high, 3), "]\n\n")
    results$h1_landuse <- h1
  }

  
  # --- SH Hypotheses ---
  cat("SH Model (pooled across", m, "imputations):\n")
  cat("--------------------------------------------\n")
  
  # H3: SH effect equal on Height vs Width?
  h3 <- pool_hypothesis(
    fits_mpb,
    "Height_Salvage_Harvest_StatusSH = MeanCanWidth_Salvage_Harvest_StatusSH",
    "H3: SH effect equal on Height vs Width"
  )
  if (!is.null(h3)) {
    cat("H3:", h3$name, "\n")
    cat("   Pooled difference:", round(h3$estimate, 3),
        "95% CI [", round(h3$ci_low, 3), ",", round(h3$ci_high, 3), "]\n\n")
    results$h3_mpb <- h3
  }
  
  # H4: SH effect stronger on AbovegroundBiomass than RootDryMass?
  h4 <- pool_hypothesis(
    fits_mpb,
    "logAbovegroundBiomass_Salvage_Harvest_StatusSH < logRootDryMass_Salvage_Harvest_StatusSH",
    "H4: SH reduces Aboveground more than Root"
  )
  if (!is.null(h4)) {
    cat("H4:", h4$name, "\n")
    cat("   Pooled difference:", round(h4$estimate, 3),
        "95% CI [", round(h4$ci_low, 3), ",", round(h4$ci_high, 3), "]\n")
    cat("   Pooled Prob(H4 true):", round(h4$post_prob, 3), "\n\n")
    results$h4_mpb <- h4
  }
  
  return(results)
}


#### Create Custom Hypothesis Tests ####

test_custom_hypothesis <- function(fits_list, hypothesis_string, hypothesis_name = "Custom") {
  #' Test a custom hypothesis across all imputed models and pool results
  #' 
  #' @param fits_list List of brms model fits
  #' @param hypothesis_string Hypothesis string in brms format
  #' @param hypothesis_name Descriptive name for the hypothesis
  #' @return Pooled hypothesis test results
  #' 
  #' @examples
  #' # Test if Height effect = Canopy Width effect for SH
  #' test_custom_hypothesis(
  #'   fits_mpb, 
  #'   "Height_Salvage_Harvest_StatusSH = MeanCanWidth_Salvage_Harvest_StatusSH",
  #'   "SH effect equal on Height vs Width"
  #' )
  #' 
  #' # Test if a correlation exceeds 0.5
  #' test_custom_hypothesis(
  #'   fits_mpb,
  #'   "rescor__Height__MeanCanWidth > 0.5",
  #'   "Height-Width correlation > 0.5"
  #' )
  
  cat("\nTesting custom hypothesis:", hypothesis_name, "\n")
  cat("Formula:", hypothesis_string, "\n\n")
  
  m <- length(fits_list)
  estimates <- numeric(m)
  cis_low <- numeric(m)
  cis_high <- numeric(m)
  post_probs <- numeric(m)
  evid_ratios <- numeric(m)
  
  for (i in 1:m) {
    tryCatch({
      h <- hypothesis(fits_list[[i]], hypothesis_string)
      estimates[i] <- h$hypothesis$Estimate
      cis_low[i] <- h$hypothesis$CI.Lower
      cis_high[i] <- h$hypothesis$CI.Upper
      post_probs[i] <- ifelse(!is.null(h$hypothesis$Post.Prob), 
                               h$hypothesis$Post.Prob, NA)
      evid_ratios[i] <- ifelse(!is.null(h$hypothesis$Evid.Ratio),
                                h$hypothesis$Evid.Ratio, NA)
    }, error = function(e) {
      cat("  Error in imputation", i, ":", e$message, "\n")
      estimates[i] <<- NA
    })
  }
  
  # Pool results
  valid <- !is.na(estimates)
  m_valid <- sum(valid)
  
  if (m_valid == 0) {
    cat("ERROR: Hypothesis test failed for all imputations\n")
    return(NULL)
  }
  
  # Rubin's rules
  pooled_est <- mean(estimates[valid])
  between_var <- var(estimates[valid])
  within_var <- mean((cis_high[valid] - cis_low[valid])^2 / (2*1.96)^2)
  total_var <- within_var + (1 + 1/m_valid) * between_var
  pooled_se <- sqrt(total_var)
  
  pooled_ci_low <- pooled_est - 1.96 * pooled_se
  pooled_ci_high <- pooled_est + 1.96 * pooled_se
  pooled_prob <- mean(post_probs[valid], na.rm = TRUE)
  
  cat("Results (pooled across", m_valid, "imputations):\n")
  cat("  Estimate:", round(pooled_est, 4), "\n")
  cat("  95% CI: [", round(pooled_ci_low, 4), ",", round(pooled_ci_high, 4), "]\n")
  if (!is.na(pooled_prob)) {
    cat("  Posterior probability:", round(pooled_prob, 4), "\n")
  }
  cat("  Between-imputation variance:", round(between_var, 6), "\n")
  cat("  Within-imputation variance:", round(within_var, 6), "\n")
  
  return(list(
    name = hypothesis_name,
    formula = hypothesis_string,
    estimate = pooled_est,
    se = pooled_se,
    ci_low = pooled_ci_low,
    ci_high = pooled_ci_high,
    post_prob = pooled_prob,
    between_var = between_var,
    within_var = within_var,
    m_valid = m_valid
  ))
}


#### Model Diagnostics ####

run_pooled_diagnostics <- function(fits_list, model_name = "Model") {
  #' Run diagnostics across all imputed models
  
  cat("\n================================================================\n")
  cat("  DIAGNOSTICS:", model_name, "(across", length(fits_list), "imputations)\n")
  cat("================================================================\n\n")
  
  m <- length(fits_list)
  
  # Collect diagnostics from each model
  all_rhat <- c()
  all_ess_bulk <- c()
  all_ess_tail <- c()
  all_divergences <- c()
  
  for (i in 1:m) {
    fit <- fits_list[[i]]
    
    fit_summary <- summary(fit)$fixed
    all_rhat <- c(all_rhat, fit_summary$Rhat)
    all_ess_bulk <- c(all_ess_bulk, fit_summary$Bulk_ESS)
    all_ess_tail <- c(all_ess_tail, fit_summary$Tail_ESS)
    
    np <- nuts_params(fit)
    n_div <- sum(subset(np, Parameter == "divergent__")$Value)
    all_divergences <- c(all_divergences, n_div)
  }
  
  cat("Convergence (Rhat) across all imputations:\n")
  cat("  Range:", round(min(all_rhat, na.rm = TRUE), 4), "-", 
      round(max(all_rhat, na.rm = TRUE), 4), "\n")
  cat("  All < 1.01?", ifelse(all(all_rhat < 1.01, na.rm = TRUE), "YES ✓", "NO ✗"), "\n\n")
  
  cat("Effective Sample Size across all imputations:\n")
  cat("  Bulk ESS range:", round(min(all_ess_bulk, na.rm = TRUE)), "-", 
      round(max(all_ess_bulk, na.rm = TRUE)), "\n")
  cat("  Tail ESS range:", round(min(all_ess_tail, na.rm = TRUE)), "-", 
      round(max(all_ess_tail, na.rm = TRUE)), "\n")
  cat("  Min ESS > 400?", ifelse(min(c(all_ess_bulk, all_ess_tail), na.rm = TRUE) > 400, 
                                  "YES ✓", "NO - consider more iterations"), "\n\n")
  
  cat("Divergent Transitions:\n")
  cat("  Per imputation:", paste(all_divergences, collapse = ", "), "\n")
  cat("  Total:", sum(all_divergences), "across", m * 8000, "total iterations\n")
  cat("  Acceptable (< 1%)?", 
      ifelse(sum(all_divergences) / (m * 8000) < 0.01, "YES ✓", "NO ✗"), "\n")
  
  return(list(
    rhat_range = range(all_rhat, na.rm = TRUE),
    ess_bulk_range = range(all_ess_bulk, na.rm = TRUE),
    ess_tail_range = range(all_ess_tail, na.rm = TRUE),
    divergences_per_imp = all_divergences,
    divergences_total = sum(all_divergences)
  ))
}


#### Create Figures ####

create_pooled_figures <- function(pooled_coef_landuse, pooled_coef_mpb,
                                  pooled_rescor_landuse, pooled_rescor_mpb,
                                  pooled_hypotheses) {
  #' Create publication-ready figures from pooled results
  
  cat("\nCreating publication figures from pooled results...\n")
  
  figures <- list()
  
  # -------------------------------------------------------------------------
  # FIGURE 1: Residual correlation heatmaps (side-by-side)
  # -------------------------------------------------------------------------
  
  prep_cor_data <- function(cor_matrix, model_name) {
    cor_df <- as.data.frame(as.table(cor_matrix))
    names(cor_df) <- c("Outcome1", "Outcome2", "Correlation")
    cor_df$Model <- model_name
    
    cor_df$Outcome1 <- factor(cor_df$Outcome1, 
                              levels = c("Height", "MeanCanWidth", "AbovegroundBiomass", 
                                         "RootDryMass", "RootToShootRatio"),
                              labels = c("Height", "Canopy\nWidth", "Aboveground\nBiomass",
                                         "Root\nMass", "Root:Shoot\nRatio"))
    cor_df$Outcome2 <- factor(cor_df$Outcome2,
                              levels = c("Height", "MeanCanWidth", "AbovegroundBiomass", 
                                         "RootDryMass", "RootToShootRatio"),
                              labels = c("Height", "Canopy\nWidth", "Aboveground\nBiomass",
                                         "Root\nMass", "Root:Shoot\nRatio"))
    return(cor_df)
  }
  
  cor_data_landuse <- prep_cor_data(pooled_rescor_landuse$matrix, "Land-Use History (6-level)")
  cor_data_mpb <- prep_cor_data(pooled_rescor_mpb$matrix, "Salvage_Harvest_Status (2-level)")
  cor_data_combined <- bind_rows(cor_data_landuse, cor_data_mpb)
  
  fig1 <- ggplot(cor_data_combined, aes(x = Outcome1, y = Outcome2, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3.5) +
    facet_wrap(~Model) +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                         midpoint = 0, limits = c(-0.5, 1),
                         name = "Residual\nCorrelation") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = "Residual Correlations Among Plant Outcomes (Pooled)",
      subtitle = paste0("Pooled across ", M_IMPUTATIONS, " imputations using Fisher-z transformation"),
      x = "", y = ""
    )
  
  figures$correlation_heatmaps <- fig1
  ggsave("mv_MI_fig1_correlation_heatmaps.pdf", fig1, width = 12, height = 6, device = "pdf")
  cat("  Saved: mv_MI_fig1_correlation_heatmaps.pdf\n")
  
  # -------------------------------------------------------------------------
  # FIGURE 2: SH effects forest plot
  # -------------------------------------------------------------------------
  
  mpb_effects <- pooled_coef_mpb %>%
    filter(str_detect(predictor, "Salvage_Harvest_StatusSH")) %>%
    mutate(
      outcome_clean = case_when(
        outcome == "Height" ~ "Height",
        outcome == "MeanCanWidth" ~ "Canopy Width",
        outcome == "AbovegroundBiomass" ~ "Aboveground Biomass",
        outcome == "RootDryMass" ~ "Root Dry Mass",
        outcome == "RootToShootRatio" ~ "Root:Shoot Ratio"
      ),
      outcome_clean = factor(outcome_clean, 
                             levels = rev(c("Height", "Canopy Width", "Aboveground Biomass",
                                            "Root Dry Mass", "Root:Shoot Ratio"))),
      significant = prob_nonzero > 0.90,
      sig_label = case_when(
        prob_nonzero > 0.99 ~ "***",
        prob_nonzero > 0.95 ~ "**",
        prob_nonzero > 0.90 ~ "*",
        TRUE ~ ""
      )
    )
  
  fig2 <- ggplot(mpb_effects, aes(x = Estimate, y = outcome_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.25, 
                   color = "gray40", linewidth = 0.8) +
    geom_point(aes(fill = significant), shape = 21, size = 4, color = "black") +
    geom_text(aes(label = sig_label, x = CI_high + 0.05), hjust = 0, size = 5) +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "steelblue"),
                      guide = "none") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 11)
    ) +
    labs(
      title = "Effect of SH-Affected Soils on Plant Outcomes (Pooled)",
      subtitle = paste0("Pooled across ", M_IMPUTATIONS, " imputations using Rubin's rules"),
      x = "Standardized Effect Size",
      y = "",
      caption = "Error bars = 95% credible intervals. * P>90%, ** P>95%, *** P>99%"
    ) +
    coord_cartesian(xlim = c(min(mpb_effects$CI_low, na.rm = TRUE) - 0.1, 
                             max(mpb_effects$CI_high, na.rm = TRUE) + 0.15))
  
  figures$mpb_forest_plot <- fig2
  ggsave("mv_MI_fig2_mpb_effects.pdf", fig2, width = 8, height = 5, device = "pdf")
  cat("  Saved: mv_MI_fig2_mpb_effects.pdf\n")
  
}


#### Main Workflow ####

run_multivariate_analysis_MI <- function(df, soil_pca_result, m = M_IMPUTATIONS) {
  #' Main workflow for multivariate analysis with multiple imputation
  #' 
  #' @param df Original data frame
  #' @param soil_pca_result Results from Greenhouse_SoilChemistry_PCA_forBayesian.R
  #' @param m Number of imputations (default: M_IMPUTATIONS)
  #' @return List with all results
  
  cat("\n")
  cat("****************************************************************\n")
  cat("*                                                              *\n")
  cat("*   MULTIVARIATE BAYESIAN ANALYSIS WITH MULTIPLE IMPUTATION   *\n")
  cat("*                                                              *\n")
  cat("*   m =", sprintf("%-3d", m), "imputations with proper pooling               *\n")
  cat("*                                                              *\n")
  cat("****************************************************************\n\n")
  
  start_time <- Sys.time()
  
  # --- Step 1: Prepare multiple imputed datasets ---
  cat("STEP 1: Preparing", m, "imputed datasets...\n")
  imp_data <- prepare_mv_datasets(df, soil_pca_result, m = m)
  
  # --- Step 2: Fit Land-Use History models ---
  cat("\nSTEP 2: Fitting Land-Use History models across imputations...\n")
  cat("  (This may take", m * 0.25, "-", m * 0.5, "minutes)\n")
  fits_landuse <- fit_all_mv_models(imp_data$datasets, model_type = "landuse")
  
  # --- Step 3: Fit SH models ---
  cat("\nSTEP 3: Fitting Salvage_Harvest_Status models across imputations...\n")
  fits_mpb <- fit_all_mv_models(imp_data$datasets, model_type = "mpb")
  
  # --- Step 4: Run diagnostics ---
  cat("\nSTEP 4: Running diagnostics...\n")
  diag_landuse <- run_pooled_diagnostics(fits_landuse, "Land-Use History")
  diag_mpb <- run_pooled_diagnostics(fits_mpb, "Salvage_Harvest_Status")
  
  # --- Step 5: Pool coefficient estimates ---
  cat("\nSTEP 5: Pooling coefficient estimates...\n")
  pooled_coef_landuse <- pool_mv_coefficients(fits_landuse)
  pooled_coef_mpb <- pool_mv_coefficients(fits_mpb)
  
  # --- Step 6: Pool residual correlations ---
  cat("\nSTEP 6: Pooling residual correlations...\n")
  pooled_rescor_landuse <- pool_residual_correlations(fits_landuse)
  pooled_rescor_mpb <- pool_residual_correlations(fits_mpb)
  
  # --- Step 7: Hypothesis tests ---
  cat("\nSTEP 7: Running pooled hypothesis tests...\n")
  pooled_hypotheses <- test_pooled_hypotheses(fits_landuse, fits_mpb)
  
  # --- Step 8: Create figures ---
  cat("\nSTEP 8: Creating publication figures...\n")
  figures <- create_pooled_figures(
    pooled_coef_landuse, pooled_coef_mpb,
    pooled_rescor_landuse, pooled_rescor_mpb,
    pooled_hypotheses
  )
  
  # --- Step 9: Save results ---
  cat("\nSTEP 9: Saving results...\n")
  write.csv(pooled_coef_landuse, "mv_MI_coefficients_landuse.csv", row.names = FALSE)
  write.csv(pooled_coef_mpb, "mv_MI_coefficients_mpb.csv", row.names = FALSE)
  write.csv(pooled_rescor_landuse$summary, "mv_MI_rescor_summary_landuse.csv", row.names = FALSE)
  write.csv(pooled_rescor_mpb$summary, "mv_MI_rescor_summary_mpb.csv", row.names = FALSE)
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  
  # --- Print summary ---
  cat("\n")
  cat("================================================================\n")
  cat("                    ANALYSIS COMPLETE                           \n")
  cat("================================================================\n")
  cat("\nTotal time:", round(as.numeric(elapsed), 1), "minutes\n")
  cat("Imputations:", m, "\n")
  
  cat("\n--- DIAGNOSTICS SUMMARY ---\n")
  cat("Land-Use History Model:\n")
  cat("  Rhat range:", round(diag_landuse$rhat_range[1], 4), "-", 
      round(diag_landuse$rhat_range[2], 4), "\n")
  cat("  Total divergences:", diag_landuse$divergences_total, "\n")
  
  cat("Salvage_Harvest_Status Model:\n")
  cat("  Rhat range:", round(diag_mpb$rhat_range[1], 4), "-", 
      round(diag_mpb$rhat_range[2], 4), "\n")
  cat("  Total divergences:", diag_mpb$divergences_total, "\n")
  
  cat("\n--- POOLED RESIDUAL CORRELATIONS ---\n")
  cat("\nLand-Use History Model:\n")
  print(round(pooled_rescor_landuse$matrix, 2))
  
  cat("\nSalvage_Harvest_Status Model:\n")
  print(round(pooled_rescor_mpb$matrix, 2))
  
  cat("\n--- OUTPUT FILES ---\n")
  cat("  mv_MI_coefficients_landuse.csv\n")
  cat("  mv_MI_coefficients_mpb.csv\n")
  cat("  mv_MI_rescor_summary_landuse.csv\n")
  cat("  mv_MI_rescor_summary_mpb.csv\n")
  cat("  mv_MI_fig1_correlation_heatmaps.pdf\n")
  cat("  mv_MI_fig2_mpb_effects.pdf\n")
  cat("  mv_MI_fig3_soil_effects.pdf\n")
  
  return(list(
    imputation_data = imp_data,
    fits = list(landuse = fits_landuse, mpb = fits_mpb),
    diagnostics = list(landuse = diag_landuse, mpb = diag_mpb),
    pooled_coefficients = list(landuse = pooled_coef_landuse, mpb = pooled_coef_mpb),
    pooled_correlations = list(landuse = pooled_rescor_landuse, mpb = pooled_rescor_mpb),
    pooled_hypotheses = pooled_hypotheses,
    figures = figures,
    m_imputations = m,
    elapsed_time = elapsed
  ))
}


#### Execute ####

# To run with default m=10 imputations:
#mv_results_MI <- run_multivariate_analysis_MI(df, soil_results)

# To run with m=20 imputations (more rigorous):
mv_results_MI <- run_multivariate_analysis_MI(df, soil_results, m = 20)

parnames(mv_results_MI$fits$landuse[[1]])

## Interested in 1stPreMPB vs. 2ndPostMPB as this seems to be the big divergence within SH statuses

# To test a custom hypothesis after running:
#   test_custom_hypothesis(
#     mv_results_MI$fits$mpb,
#     "Height_Salvage_Harvest_StatusSH = MeanCanWidth_Salvage_Harvest_StatusSH",
#     "SH effect equal on Height vs Width"
#   )


# In hypothesis_comparison_figure.R - these loops will give you new values
outcomes <- c("Height", "MeanCanWidth", "logAbovegroundBiomass", 
              "logRootDryMass", "logRootToShootRatio")

for (outcome in outcomes) {
  hyp_string <- paste0(outcome, "_Treatment_LandUseHistory1stPreMPB > ",
                       outcome, "_Treatment_LandUseHistory2ndPostMPB")
  test_custom_hypothesis(mv_results_MI$fits$landuse, hyp_string, "Impacts of 1stPreMPB and 2ndPostMPB on Plant Outcomes")
}


# SH vs Non-SH (testing if SH is worse, i.e., coefficient < 0)
outcomes <- c("Height", "MeanCanWidth", "logAbovegroundBiomass", 
              "logRootDryMass", "logRootToShootRatio")

for (outcome in outcomes) {
  hyp_string <- paste0(outcome, "_Salvage_Harvest_StatusSH < 0")
  test_custom_hypothesis(mv_results_MI$fits$mpb, hyp_string, 
                         paste("SH reduces", outcome))
}


## Management comparisons 
outcomes <- c("Height", "MeanCanWidth", "logAbovegroundBiomass", 
              "logRootDryMass", "logRootToShootRatio")

# Comparison 1: 1stPostMPB vs 2ndPostMPB
cat("\n=== COMPARISON 1: 1stPostMPB vs 2ndPostMPB ===\n")
for (outcome in outcomes) {
  hyp_string <- paste0(outcome, "_Treatment_LandUseHistory1stPostMPB > ",
                       outcome, "_Treatment_LandUseHistory2ndPostMPB")
  test_custom_hypothesis(mv_results_MI$fits$landuse, hyp_string,
                         paste("1stPostMPB > 2ndPostMPB for", outcome))
}

# Comparison 2: Recent_Cut vs 2ndPostMPB
cat("\n=== COMPARISON 2: Recent_Cut vs 2ndPostMPB ===\n")
for (outcome in outcomes) {
  hyp_string <- paste0(outcome, "_Treatment_LandUseHistoryRecent_Cut > ",
                       outcome, "_Treatment_LandUseHistory2ndPostMPB")
  test_custom_hypothesis(mv_results_MI$fits$landuse, hyp_string,
                         paste("Recent_Cut > 2ndPostMPB for", outcome))
}

# Comparison 3: Recent_Cut vs 1stPostMPB
cat("\n=== COMPARISON 3: Recent_Cut vs 1stPostMPB ===\n")
for (outcome in outcomes) {
  hyp_string <- paste0(outcome, "_Treatment_LandUseHistoryRecent_Cut > ",
                       outcome, "_Treatment_LandUseHistory1stPostMPB")
  test_custom_hypothesis(mv_results_MI$fits$landuse, hyp_string,
                         paste("Recent_Cut > 1stPostMPB for", outcome))
}

## Figure to represent these findings is in file hypothesis_comparison_figure.R or management_comparison_figure.R





cat("\nMultivariate analysis script with multiple imputation loaded.")
cat("\nDefault imputations: m =", M_IMPUTATIONS)
cat("\nRun with: mv_results_MI <- run_multivariate_analysis_MI(df, soil_results)\n")
