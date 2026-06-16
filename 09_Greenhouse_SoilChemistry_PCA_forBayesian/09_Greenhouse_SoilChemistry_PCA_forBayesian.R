#### Soil Chemistry PCA Module - with Imputation - Deals with Missingness/Heavy Patchiness Ahead of Bayesian ####
## Creates new soil chem variables (PCs) that are uncorrelated with one another, better interpretability, less overfitting, deals with multicollinearity
## Soil pH is excluded here as it seemed to throw everything off - adding as its own variable with mice imputation in the later Bayesian model 
## Created with the assistance of Claude AI Opus 4.5

#### Packages ####

library(tidyverse)
library(mice) ## Imputation 
library(FactoMineR)
library(factoextra)
library(missMDA) ## For PCA with missing data


#### Approach A: Imputate Soil Chemistry Values, then Perform PCA ####

## This approach:
## 1. Uses mice to impute missing soil chemistry values
## 2. Performs PCA on each imputed dataset
## 3. Averages PC scores across imputations (or uses first imputation)

setwd("/Users/juliefowler/OneDrive - Colostate/PhD Work/Studies/SFSP Clear Cuts/Paper Writing/GitHub/Files/11. Greenhouse_SoilChemistry_PCA_forBayesian")
df <- read.delim("Bayesian_AllData_ExcludedThoseWithOnlyNAsAcross_RemovedProblemColumns.txt", row.names = 1) ## Data with the two empty VIP columns removed per Diagnostic Script

impute_soil_then_pca <- function(df, soil_vars, m = 50, n_pcs = 3, seed = 42) {
  
  cat("\n========================================\n")
  cat("APPROACH A: IMPUTE SOIL CHEMISTRY, THEN PCA\n")
  cat("========================================\n\n")
  
  # Check which soil vars exist
  soil_vars_exist <- soil_vars[soil_vars %in% names(df)]
  cat("Soil variables found:", length(soil_vars_exist), "of", length(soil_vars), "\n")
  cat("Missing:", setdiff(soil_vars, soil_vars_exist), "\n\n")
  
  # Extract soil chemistry subset
  soil_data <- df %>% select(all_of(soil_vars_exist))
  soil_data
  
  # Check missingness
  n_complete <- sum(complete.cases(soil_data))
  cat("Complete cases for soil chemistry:", n_complete, "of", nrow(soil_data), "\n")
  cat("This represents", round(n_complete/nrow(soil_data)*100, 1), "% of samples\n\n")
  
  if (n_complete < 10) {
    cat("WARNING: Very few complete cases. Imputation will rely heavily on\n")
    cat("available predictors. Consider whether soil chemistry can be\n")
    cat("meaningfully imputed in your context.\n\n")
  }
  
  # --- Step 1: Impute soil chemistry ---
  cat("Step 1: Imputing soil chemistry with mice...\n")
  
  # Use all available variables as predictors for imputation
  # Including treatments, texture, outcomes, and microbial data
  
  # Create imputation dataset with helpful predictors
  # Salvage_Harvest_Status and Texture are totally nested in Treatment_LandUseHistory, so just using that
  imp_vars <- c(soil_vars_exist, 
              "Treatment_LandUseHistory",
              "Height", "MeanCanWidth", "AbovegroundBiomass", "RootDryMass") ## Root:Shoot Ratio made of RootDryMass and Aboveground Biomass so not necessary to include
  imp_vars <- imp_vars[imp_vars %in% names(df)]
  
  imp_data <- df %>% select(all_of(imp_vars))
  
  # Run mice
  set.seed(seed)
  imp <- mice(imp_data, m = m, maxit = 20, method = "pmm", printFlag = FALSE)
  
  cat("  Imputation complete.\n\n")
  
  # --- Step 2: PCA on each imputed dataset ---
  cat("Step 2: Performing PCA on each imputed dataset...\n")
  
  pca_results <- list()
  all_scores <- list()
  all_loadings <- list()
  
  for (i in 1:m) {
    df_imp <- complete(imp, i)
    soil_imp <- df_imp %>% select(all_of(soil_vars_exist))
    
    # Scale and run PCA
    pca_i <- prcomp(soil_imp, center = TRUE, scale. = TRUE)
    pca_results[[i]] <- pca_i
    
    # Extract scores
    scores_i <- as.data.frame(pca_i$x[, 1:min(n_pcs, ncol(pca_i$x))])
    names(scores_i) <- paste0("Soil_PC", 1:ncol(scores_i))
    all_scores[[i]] <- scores_i
    
    # Extract loadings
    loadings_i <- as.data.frame(pca_i$rotation[, 1:min(n_pcs, ncol(pca_i$rotation))])
    all_loadings[[i]] <- loadings_i
  } 
  
  cat("  PCA complete for all", m, "imputations.\n\n")
  
  # --- Step 3: Average/summarize across imputations ---
  cat("Step 3: Summarizing PCA results across imputations...\n")
  
  # Average PC scores
  avg_scores <- Reduce(`+`, all_scores) / m
  
  # Average loadings (for interpretation)
  avg_loadings <- Reduce(`+`, all_loadings) / m
  
  # Variance explained (from first imputation - should be similar across)
  var_explained <- pca_results[[1]]$sdev^2 / sum(pca_results[[1]]$sdev^2)
  
  cat("\n  Variance explained by each PC (averaged):\n")
  for (i in 1:n_pcs) {
    cat("    PC", i, ":", round(var_explained[i] * 100, 1), "%\n")
  }
  cat("    Cumulative:", round(sum(var_explained[1:n_pcs]) * 100, 1), "%\n")
  
  cat("\n  PC Loadings (averaged across imputations):\n")
  print(round(avg_loadings, 3))
  
  # --- Visualization ---
  # Use first imputation for plots
  p_scree <- fviz_eig(pca_results[[1]], addlabels = TRUE) +
    labs(title = "Soil Chemistry PCA - Scree Plot",
         subtitle = "(Based on first imputation)")
  
  p_var <- fviz_pca_var(pca_results[[1]], col.var = "contrib",
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = TRUE) +
    labs(title = "Soil Chemistry PCA - Variable Contributions")
  
  p_biplot <- fviz_pca_biplot(pca_results[[1]], repel = TRUE,
                              col.var = "#2E9FDF", col.ind = "#696969") +
    labs(title = "Soil Chemistry PCA - Biplot")
  
  return(list(
    imputation = imp,
    pca_results = pca_results,
    avg_scores = avg_scores,
    avg_loadings = avg_loadings,
    var_explained = var_explained,
    all_scores = all_scores,  # Keep individual scores for proper MI analysis
    plots = list(scree = p_scree, var = p_var, biplot = p_biplot)
  ))
}


#### Alternative Method (missMDA; regularized iterative PCA or PCA with missing data) - for confirmation ####

# missMDA can perform PCA directly with missing data using regularized 
# iterative PCA (useful as a comparison/robustness check)

pca_with_missing <- function(df, soil_vars, n_pcs = 3) {
  
  cat("\n========================================\n")
  cat("ALTERNATIVE: PCA WITH MISSING DATA (missMDA)\n")
  cat("========================================\n\n")
  
  soil_vars_exist <- soil_vars[soil_vars %in% names(df)]
  soil_data <- df %>% select(all_of(soil_vars_exist))
  
  cat("Using missMDA::imputePCA for regularized iterative PCA...\n")
  cat("This method handles missing data directly within the PCA algorithm.\n\n")
  
  # Estimate number of components
  cat("Estimating optimal number of components...\n")
  n_comp_est <- estim_ncpPCA(soil_data, ncp.min = 1, ncp.max = 5, method = "Regularized")
  cat("  Estimated optimal components:", n_comp_est$ncp, "\n\n")
  
  # Perform PCA with imputation
  pca_result <- imputePCA(soil_data, ncp = n_pcs, method = "Regularized")
  
  # The result contains:
  # - completeObs: the imputed data
  # - fittedX: the fitted values from PCA
  
  # Run standard PCA on completed data for visualization
  pca_final <- prcomp(pca_result$completeObs, center = TRUE, scale. = TRUE)
  
  # Extract scores
  scores <- as.data.frame(pca_final$x[, 1:n_pcs])
  names(scores) <- paste0("Soil_PC", 1:n_pcs)
  
  # Variance explained
  var_explained <- pca_final$sdev^2 / sum(pca_final$sdev^2)
  
  cat("Variance explained:\n")
  for (i in 1:n_pcs) {
    cat("  PC", i, ":", round(var_explained[i] * 100, 1), "%\n")
  }
  
  cat("\nLoadings:\n")
  print(round(pca_final$rotation[, 1:n_pcs], 3))
  
  # Visualization
  p_scree <- fviz_eig(pca_final, addlabels = TRUE) +
    labs(title = "Soil Chemistry PCA (missMDA) - Scree Plot")
  
  p_var <- fviz_pca_var(pca_final, col.var = "contrib",
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = TRUE) +
    labs(title = "Soil Chemistry PCA (missMDA) - Variable Contributions")
  
  return(list(
    pca = pca_final,
    imputed_data = pca_result$completeObs,
    scores = scores,
    loadings = pca_final$rotation[, 1:n_pcs],
    var_explained = var_explained,
    plots = list(scree = p_scree, var = p_var)
  ))
}


#### Sensitivity Analysis #### 

# Run models with and without soil chemistry to assess sensitivity

run_sensitivity_analysis <- function(df, imp_full, soil_pca_result, outcome, 
                                      other_predictors, iter = 4000, seed = 42) {
  
  cat("\n========================================\n")
  cat("SENSITIVITY ANALYSIS FOR:", outcome, "\n")
  cat("========================================\n\n")
  
  # Define two models:
  # Model 1: WITHOUT soil chemistry (uses only complete predictors)
  # Model 2: WITH soil chemistry PCs (requires imputation)
  
  # --- Model 1: Without soil chemistry ---
  cat("Model 1: Without soil chemistry PCs\n")
  cat("  Predictors:", paste(other_predictors, collapse = ", "), "\n\n")
  
  formula_no_soil <- as.formula(paste(outcome, "~", paste(other_predictors, collapse = " + ")))
  
  # This model can use all data (no imputation needed for soil)
  # But still needs imputation for microbial variables
  
  # --- Model 2: With soil chemistry PCs ---
  soil_pc_names <- names(soil_pca_result$avg_scores)
  all_predictors <- c(other_predictors, soil_pc_names)
  
  cat("Model 2: With soil chemistry PCs\n")
  cat("  Predictors:", paste(all_predictors, collapse = ", "), "\n\n")
  
  formula_with_soil <- as.formula(paste(outcome, "~", paste(all_predictors, collapse = " + ")))
  
  # Results comparison would be done in the main analysis
  # Here we just return the formulas and relevant info
  
  return(list(
    formula_no_soil = formula_no_soil,
    formula_with_soil = formula_with_soil,
    n_predictors_no_soil = length(other_predictors),
    n_predictors_with_soil = length(all_predictors)
  ))
}


#### Helper - Add Soil PCs to Multiple Imputed Datasets ####

add_soil_pcs_to_imputations <- function(imp, soil_pca_result) {
  
  # For proper multiple imputation analysis, we should use the 
  # PC scores from each corresponding imputation
  
  m <- imp$m
  
  # Create list of complete datasets with PCs added
  imputed_with_pcs <- list()
  
  for (i in 1:m) {
    df_i <- complete(imp, i)
    
    # Add soil PC scores from this imputation
    if (i <= length(soil_pca_result$all_scores)) {
      pc_scores <- soil_pca_result$all_scores[[i]]
    } else {
      # If fewer soil imputations than full imputations, use averaged
      pc_scores <- soil_pca_result$avg_scores
    }
    
    df_i <- cbind(df_i, pc_scores)
    imputed_with_pcs[[i]] <- df_i
  }
  
  return(imputed_with_pcs)
}


#### PC Interpretation Helper ####

## After PCA, you need to INTERPRET what each PC represents.
## This function helps by identifying high-loading variables.

## threshold = 0.3: variables with |loading| > 0.3 are considered important

interpret_soil_pcs <- function(loadings, var_explained, threshold = 0.3) {
  
  cat("\n========================================\n")
  cat("SOIL PC INTERPRETATION GUIDE\n")
  cat("========================================\n\n")
  
  n_pcs <- ncol(loadings)
  
  for (i in 1:n_pcs) {
    cat("PC", i, "(", round(var_explained[i] * 100, 1), "% variance):\n")
    
    # Find variables with high loadings
    high_pos <- rownames(loadings)[loadings[, i] > threshold]
    high_neg <- rownames(loadings)[loadings[, i] < -threshold]
    
    if (length(high_pos) > 0) {
      cat("  High positive loadings:\n")
      for (var in high_pos) {
        cat("    ", var, ": ", round(loadings[var, i], 3), "\n")
      }
    }
    
    if (length(high_neg) > 0) {
      cat("  High negative loadings:\n")
      for (var in high_neg) {
        cat("    ", var, ": ", round(loadings[var, i], 3), "\n")
      }
    }
    
    # Suggest interpretation
    cat("\n  Suggested interpretation:\n")
    
    # Common patterns in soil chemistry
    nutrients <- c("Total_N_percent", "Total_C_percent", "NH4_waterextract_mgperL",
                   "NO3_waterextract_mgperL", "PO4_waterextract_mgperL")
    cations <- c("Ca_waterextract_mgperL", "Mg_waterextract_mgperL", 
                 "K_waterextract_mgperL", "Na_waterextract_mgperL")
    
    high_all <- c(high_pos, high_neg)
    
    if (any(nutrients %in% high_all)) {
      cat("    - May represent soil fertility/nutrient availability\n")
    }
    if (any(cations %in% high_all)) {
      cat("    - May represent cation exchange/base saturation\n")
    }
    if ("Total_C_percent" %in% high_all && "Total_N_percent" %in% high_all) {
      cat("    - May represent organic matter content\n")
    }
    if ("SO4_waterextract_mgperL" %in% high_all) {
      cat("    - May be influenced by sulfur deposition or weathering\n")
    }
    
    cat("\n")
  }
}


#### Main PCA Workflow ####

run_soil_pca_workflow <- function(df, m_imputations = 50, n_pcs = 3, seed = 42) {
  
  cat("\n")
  cat("################################################################\n")
  cat("#           SOIL CHEMISTRY PCA WORKFLOW                        #\n")
  cat("################################################################\n")
  
  soil_vars <- c("Total_C_percent", "Total_N_percent", 
                 "Na_waterextract_mgperL", "NH4_waterextract_mgperL",
                 "K_waterextract_mgperL", "Mg_waterextract_mgperL",
                 "Ca_waterextract_mgperL", "Cl_waterextract_mgperL",
                 "NO3_waterextract_mgperL", "PO4_waterextract_mgperL",
                 "SO4_waterextract_mgperL")
  
  # Run main approach (impute then PCA)
  result_a <- impute_soil_then_pca(df, soil_vars, m = m_imputations, 
                                    n_pcs = n_pcs, seed = seed)
  
  # Run alternative approach (PCA with missing data) for comparison
  cat("\n--- Running alternative method for comparison ---\n")
  result_alt <- tryCatch(
    pca_with_missing(df, soil_vars, n_pcs = n_pcs),
    error = function(e) {
      cat("missMDA method failed:", e$message, "\n")
      cat("Continuing with mice-based approach only.\n")
      return(NULL)
    }
  )
  
  # Compare results if both succeeded
  if (!is.null(result_alt)) {
    cat("\n========================================\n")
    cat("COMPARISON OF METHODS\n")
    cat("========================================\n\n")
    
    # Compare loadings
    cat("Correlation of PC1 loadings between methods:\n")
    r <- cor(result_a$avg_loadings[, 1], result_alt$loadings[, 1])
    cat("  r =", round(r, 3), "\n")
    
    if (abs(r) > 0.9) {
      cat("  Methods are highly consistent.\n")
    } else if (abs(r) > 0.7) {
      cat("  Methods show moderate agreement.\n")
    } else {
      cat("  WARNING: Methods show substantial disagreement.\n")
      cat("  Consider which is more appropriate for your data.\n")
    }
  }
  
  # Interpretation guide
  interpret_soil_pcs(result_a$avg_loadings, result_a$var_explained)
  
  cat("\n################################################################\n")
  cat("#           SOIL PCA WORKFLOW COMPLETE                         #\n")
  cat("################################################################\n")
  
  return(list(
    primary = result_a,
    alternative = result_alt
  ))
}



#### Run ####

df <- read.delim("Bayesian_AllData_ExcludedThoseWithOnlyNAsAcross_RemovedProblemColumns.txt", row.names = 1)
soil_pca_result <- run_soil_pca_workflow(df, m_imputations = 50, n_pcs = 3)

## Save plots
ggsave("soil_pca_scree.pdf", soil_pca_result$primary$plots$scree, width = 8, height = 6)
ggsave("soil_pca_biplot.pdf", soil_pca_result$primary$plots$biplot, width = 10, height = 8)
ggsave("soil_pca_variables.pdf", soil_pca_result$primary$plots$var, width = 10, height = 8)

# cat("\nSoil PCA module loaded. Run: soil_results <- run_soil_pca_workflow(df)\n")




