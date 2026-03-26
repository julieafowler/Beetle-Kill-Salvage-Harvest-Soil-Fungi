#### Soil PCA Tables for Supplemental, etc. ####
## Created with the assistance of Claude AI Opus 4.5

library(tidyverse)

#### Extract PCA Loadings ####

extract_pca_loadings <- function(soil_results) {
  #' Extract PCA loadings from soil_results object
  #' 
  #' Uses first imputation's loadings rather than averaged loadings
  #' because PC signs can flip across imputations, causing averaged
  #' loadings to cancel out toward zero. The first imputation is
  #' representative and matches the visualization plots.
  #' 
  #' @param soil_results Results from soil_pca_module.R
  #' @return Data frame with loadings
  
  cat("Extracting PCA loadings (from first imputation)...\n")
  
  # Get loadings from first imputation (not averaged - sign flipping issue)
  loadings_raw <- soil_results$primary$pca_results[[1]]$rotation[, 1:3]
  
  # Convert to data frame
  loadings_df <- as.data.frame(loadings_raw)
  loadings_df$variable <- rownames(loadings_df)
  
  # Rename columns
  names(loadings_df)[1:3] <- c("Soil_PC1", "Soil_PC2", "Soil_PC3")
  
  # Move variable to first column
  loadings_df <- loadings_df %>%
    select(variable, everything())
  
  # Reset row names
  rownames(loadings_df) <- NULL
  
  cat("  Extracted loadings for", nrow(loadings_df), "variables\n")
  
  return(loadings_df)
}


#### Extract Variance Explained ####

extract_variance_explained <- function(soil_results) {
  #' Extract variance explained by each PC
  #' 
  #' @param soil_results Results from soil_pca_module.R
  #' @return Data frame with variance explained
  
  cat("Extracting variance explained...\n")
  
  var_exp <- soil_results$primary$var_explained
  
  # Create data frame
  var_df <- data.frame(
    PC = paste0("Soil_PC", 1:length(var_exp)),
    variance_explained = var_exp,
    percent = var_exp * 100,
    cumulative_percent = cumsum(var_exp) * 100
  )
  
  cat("  PC1:", round(var_df$percent[1], 1), "%, ",
      "PC2:", round(var_df$percent[2], 1), "%, ",
      "PC3:", round(var_df$percent[3], 1), "%\n")
  
  return(var_df)
}


#### Clean Variable Names ####

clean_soil_variable_names <- function(loadings_df) {
  #' Clean soil variable names for publication
  #' 
  #' @param loadings_df Data frame with loadings
  #' @return Data frame with cleaned names
  
  loadings_df <- loadings_df %>%
    mutate(
      variable_clean = case_when(
        # Actual variable names from your data
        variable == "Total_C_percent" ~ "Total C (%)",
        variable == "Total_N_percent" ~ "Total N (%)",
        variable == "Na_waterextract_mgperL" ~ "Sodium (Na)",
        variable == "NH4_waterextract_mgperL" ~ "Ammonium (NH4)",
        variable == "K_waterextract_mgperL" ~ "Potassium (K)",
        variable == "Mg_waterextract_mgperL" ~ "Magnesium (Mg)",
        variable == "Ca_waterextract_mgperL" ~ "Calcium (Ca)",
        variable == "Cl_waterextract_mgperL" ~ "Chloride (Cl)",
        variable == "NO3_waterextract_mgperL" ~ "Nitrate (NO3)",
        variable == "PO4_waterextract_mgperL" ~ "Phosphate (PO4)",
        variable == "SO4_waterextract_mgperL" ~ "Sulfate (SO4)",
        variable == "pH" ~ "pH",
        TRUE ~ variable  # Fallback: keep original name
      )
    ) %>%
    select(variable, variable_clean, everything())
  
  return(loadings_df)
}


#### Format Loadings Table ####

format_loadings_for_publication <- function(loadings_df, var_explained_df, 
                                             highlight_threshold = 0.3) {
  #' Format loadings table for publication
  #' 
  #' @param loadings_df Data frame with loadings
  #' @param var_explained_df Data frame with variance explained
  #' @param highlight_threshold Threshold for highlighting strong loadings
  #' @return Formatted table
  
  cat("Formatting loadings for publication...\n")
  
  # Clean names and format
  pub_table <- loadings_df %>%
    clean_soil_variable_names() %>%
    mutate(
      # Round loadings
      PC1 = round(Soil_PC1, 3),
      PC2 = round(Soil_PC2, 3),
      PC3 = round(Soil_PC3, 3),
      
      # Identify strong loadings (|loading| > threshold)
      PC1_strong = abs(Soil_PC1) > highlight_threshold,
      PC2_strong = abs(Soil_PC2) > highlight_threshold,
      PC3_strong = abs(Soil_PC3) > highlight_threshold,
      
      # Format with asterisks for strong loadings
      PC1_fmt = ifelse(PC1_strong, paste0(PC1, "*"), as.character(PC1)),
      PC2_fmt = ifelse(PC2_strong, paste0(PC2, "*"), as.character(PC2)),
      PC3_fmt = ifelse(PC3_strong, paste0(PC3, "*"), as.character(PC3))
    ) %>%
    arrange(desc(abs(Soil_PC1))) %>%  # Sort by PC1 loading magnitude
    select(
      Variable = variable_clean,
      `PC1` = PC1_fmt,
      `PC2` = PC2_fmt,
      `PC3` = PC3_fmt
    )
  
  # Add variance explained row
  var_row <- data.frame(
    Variable = "Variance Explained (%)",
    PC1 = sprintf("%.1f", var_explained_df$percent[1]),
    PC2 = sprintf("%.1f", var_explained_df$percent[2]),
    PC3 = sprintf("%.1f", var_explained_df$percent[3])
  )
  
  pub_table <- bind_rows(pub_table, var_row)
  
  cat("  Strong loadings (|r| >", highlight_threshold, ") marked with *\n")
  
  return(pub_table)
}


#### Create Interpretation Summary ####

create_interpretation_summary <- function(loadings_df, threshold = 0.3) {
  #' Create text summary of what each PC represents
  #' 
  #' @param loadings_df Data frame with loadings
  #' @param threshold Threshold for "strong" loadings
  #' @return Text summary
  
  cat("\nCreating PC interpretation summary...\n")
  
  loadings_clean <- loadings_df %>%
    clean_soil_variable_names()
  
  interpretation <- list()
  
  for (pc_num in 1:3) {
    pc_col <- paste0("Soil_PC", pc_num)
    
    # Get strong positive and negative loadings
    strong_pos <- loadings_clean %>%
      filter(.data[[pc_col]] > threshold) %>%
      arrange(desc(.data[[pc_col]])) %>%
      pull(variable_clean)
    
    strong_neg <- loadings_clean %>%
      filter(.data[[pc_col]] < -threshold) %>%
      arrange(.data[[pc_col]]) %>%
      pull(variable_clean)
    
    # Create interpretation
    if (length(strong_pos) > 0 | length(strong_neg) > 0) {
      pos_text <- if (length(strong_pos) > 0) {
        paste("High values:", paste(strong_pos, collapse = ", "))
      } else ""
      
      neg_text <- if (length(strong_neg) > 0) {
        paste("Low values:", paste(strong_neg, collapse = ", "))
      } else ""
      
      interpretation[[pc_col]] <- paste(pos_text, neg_text, sep = ". ")
    }
  }
  
  return(interpretation)
}


#### Create PC Interpretation Table ####

create_pc_interpretation_table <- function(loadings_df, var_explained_df, threshold = 0.3) {
  #' Create a table summarizing what each PC represents
  #' 
  #' @param loadings_df Data frame with loadings
  #' @param var_explained_df Data frame with variance explained
  #' @param threshold Threshold for strong loadings
  #' @return Interpretation table
  
  loadings_clean <- loadings_df %>%
    clean_soil_variable_names()
  
  interpretation_table <- data.frame(
    PC = character(),
    Variance_Explained = character(),
    Strong_Positive_Loadings = character(),
    Strong_Negative_Loadings = character(),
    Interpretation = character(),
    stringsAsFactors = FALSE
  )
  
  for (pc_num in 1:3) {
    pc_name <- paste0("Soil_PC", pc_num)
    pc_col <- paste0("Soil_PC", pc_num)
    
    strong_pos <- loadings_clean %>%
      filter(.data[[pc_col]] > threshold) %>%
      arrange(desc(.data[[pc_col]])) %>%
      pull(variable_clean)
    
    strong_neg <- loadings_clean %>%
      filter(.data[[pc_col]] < -threshold) %>%
      arrange(.data[[pc_col]]) %>%
      pull(variable_clean)
    
    # Auto-interpretation based on loadings
    interp <- case_when(
      pc_num == 1 & any(str_detect(strong_pos, "N|C|Organic")) ~ 
        "Overall soil fertility / organic matter",
      pc_num == 2 & any(str_detect(strong_pos, "Ca|Mg|Base")) ~ 
        "Base cation availability",
      pc_num == 3 ~ "Micronutrient / acidity gradient",
      TRUE ~ "Mixed gradient"
    )
    
    interpretation_table <- bind_rows(interpretation_table, data.frame(
      PC = pc_name,
      Variance_Explained = sprintf("%.1f%%", var_explained_df$percent[pc_num]),
      Strong_Positive_Loadings = paste(strong_pos, collapse = ", "),
      Strong_Negative_Loadings = paste(strong_neg, collapse = ", "),
      Interpretation = interp,
      stringsAsFactors = FALSE
    ))
  }
  
  return(interpretation_table)
}


#### Main Function ####

generate_soil_pca_tables <- function(soil_results, output_prefix = "soil_pca") {
  #' Main function to generate all soil PCA tables
  #' 
  #' @param soil_results Results from soil_pca_module.R
  #' @param output_prefix Prefix for output file names
  #' @return List with all tables
  
  cat("\n")
  cat("================================================================\n")
  cat("  GENERATING SOIL PCA TABLES FOR PUBLICATION                    \n")
  cat("================================================================\n\n")
  
  # Extract loadings
  loadings <- extract_pca_loadings(soil_results)
  
  # Extract variance explained
  var_explained <- extract_variance_explained(soil_results)
  
  # Format for publication
  pub_loadings <- format_loadings_for_publication(loadings, var_explained)
  
  # Create interpretation table
  interp_table <- create_pc_interpretation_table(loadings, var_explained)
  
  # Create text interpretation
  interp_text <- create_interpretation_summary(loadings)
  
  # Save tables
  cat("\nSaving tables...\n")
  
  # Raw loadings
  write.csv(loadings %>% clean_soil_variable_names() %>%
              select(variable, variable_clean, Soil_PC1, Soil_PC2, Soil_PC3),
            paste0(output_prefix, "_loadings_raw.csv"), row.names = FALSE)
  
  # Publication-formatted loadings
  write.csv(pub_loadings, paste0(output_prefix, "_loadings_publication.csv"), row.names = FALSE)
  
  # Variance explained
  write.csv(var_explained, paste0(output_prefix, "_variance_explained.csv"), row.names = FALSE)
  
  # Interpretation table
  write.csv(interp_table, paste0(output_prefix, "_interpretation.csv"), row.names = FALSE)
  
  cat("\nOutput files created:\n")
  cat("  ", paste0(output_prefix, "_loadings_raw.csv"), "\n")
  cat("  ", paste0(output_prefix, "_loadings_publication.csv"), "\n")
  cat("  ", paste0(output_prefix, "_variance_explained.csv"), "\n")
  cat("  ", paste0(output_prefix, "_interpretation.csv"), "\n")
  
  # Print interpretation summary
  cat("\n")
  cat("================================================================\n")
  cat("  PC INTERPRETATION SUMMARY                                      \n")
  cat("================================================================\n")
  for (i in 1:nrow(interp_table)) {
    cat("\n", interp_table$PC[i], "(", interp_table$Variance_Explained[i], "):\n")
    cat("  Interpretation:", interp_table$Interpretation[i], "\n")
    if (nchar(interp_table$Strong_Positive_Loadings[i]) > 0) {
      cat("  Positive loadings:", interp_table$Strong_Positive_Loadings[i], "\n")
    }
    if (nchar(interp_table$Strong_Negative_Loadings[i]) > 0) {
      cat("  Negative loadings:", interp_table$Strong_Negative_Loadings[i], "\n")
    }
  }
  
  return(list(
    loadings_raw = loadings %>% clean_soil_variable_names(),
    loadings_publication = pub_loadings,
    variance_explained = var_explained,
    interpretation = interp_table,
    interpretation_text = interp_text
  ))
}


#### Execute ####

# To run:
# 
# # Make sure soil_results is loaded from your previous analysis
pca_tables <- generate_soil_pca_tables(soil_pca_result)
# 
# # Access specific tables:
# pca_tables$loadings_publication     # Ready for publication
# pca_tables$interpretation           # PC interpretation table

cat("\nSoil PCA tables script loaded.")
cat("\nRun with: pca_tables <- generate_soil_pca_tables(soil_results)\n")
