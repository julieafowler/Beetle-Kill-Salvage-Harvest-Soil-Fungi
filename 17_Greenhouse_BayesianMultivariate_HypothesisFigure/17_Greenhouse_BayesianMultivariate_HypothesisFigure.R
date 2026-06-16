#### Hypothesis Test Comparison Figure ####
## Created with the assistance of Claude AI Opus 4.5

## Creates a combined forest plot showing:
## 1stPreMPB vs 2ndPostMPB contrasts (6-level model)
## Salvage_Harvested vs Non_Salvage_Harvested effects (2-level model)

## Uses manual y-positioning to ensure correct alignment

library(tidyverse)


## To test all outcomes between 1stPreMPB vs. 2ndPreMPB:

outcomes <- c("Height", "MeanCanWidth", "logAbovegroundBiomass", 
              "logRootDryMass", "logRootToShootRatio")

for (outcome in outcomes) {
  cat("\n========================================\n")
  
  # NO "b_" prefix - hypothesis() adds it automatically
  hyp_string <- paste0(outcome, "_Treatment_LandUseHistory1stPreMPB > ",
                       outcome, "_Treatment_LandUseHistory2ndPostMPB")
  
  test_custom_hypothesis(
    mv_results_MI$fits$landuse,
    hyp_string,
    paste("1stPreMPB > 2ndPostMPB for", outcome)
  )
}


## General Salvage_Harvested vs. Non_Salvage_Harvested  
## The coefficient Salvage_Harvest_Status represents the difference: Salvage_Harvested Non_Salvage_Harvested. So testing < 0 asks "Is Salvage_Harvested worse than Non_Salvage_Harvested?"

parnames(mv_results_MI$fits$mpb[[1]])

outcomes <- c("Height", "MeanCanWidth", "logAbovegroundBiomass", 
              "logRootDryMass", "logRootToShootRatio")

for (outcome in outcomes) {
  cat("\n========================================\n")
  
  # Test if SH effect is negative (i.e., Salvage_Harvested < Non_Salvage_Harvested)
  hyp_string <- paste0(outcome, "_Salvage_Harvest_Status < 0")
  
  test_custom_hypothesis(
    mv_results_MI$fits$mpb,
    hyp_string,
    paste("SH reduces", outcome)
  )
}


# =============================================================================
# CREATE HYPOTHESIS COMPARISON FIGURE
# =============================================================================

create_hypothesis_figure <- function(output_file = "mv_MI_fig_hypothesis_comparison.pdf") {
  
  cat("Creating hypothesis test comparison figure...\n")
  
  # -------------------------------------------------------------------------
  # Input results - create each comparison separately with explicit y_offset
  # -------------------------------------------------------------------------
  
  # Blue (1stPreMPB vs 2ndPostMPB) - TOP position (y_offset = +0.15)
  landuse_results <- data.frame(
    outcome = c("Height", "Canopy Width", "Aboveground\nBiomass", 
                "Root Mass", "Root:Shoot\nRatio"),
    outcome_num = c(5, 4, 3, 2, 1),
    estimate = c(3.030, 3.361, 1.514, 1.300, -0.359),
    ci_low = c(2.484, 2.808, 1.134, 0.879, -0.718),
    ci_high = c(3.575, 3.914, 1.894, 1.721, -0.0004),
    prob = c(1.000, 1.000, 1.000, 1.000, 0.044),
    comparison = "1stPreMPB vs 2ndPostMPB",
    y_offset = 0.15,
    stringsAsFactors = FALSE
  )
  
  # Red (Non_Salvage_Harvested vs Salvage_Harvested) - BOTTOM position (y_offset = -0.15)
  # Signs flipped so positive = Non_Salvage_Harvested advantage
  mpb_results <- data.frame(
    outcome = c("Height", "Canopy Width", "Aboveground\nBiomass", 
                "Root Mass", "Root:Shoot\nRatio"),
    outcome_num = c(5, 4, 3, 2, 1),
    estimate = c(1.525, 1.727, 0.669, 0.489, -0.121),
    ci_low = c(1.168, 1.341, 0.418, 0.217, -0.366),
    ci_high = c(1.883, 2.113, 0.920, 0.761, 0.125),
    prob = c(1.000, 1.000, 1.000, 0.999, 0.199),
    comparison = "Non_Salvage_Harvested vs Salvage_Harvested",
    y_offset = -0.15,
    stringsAsFactors = FALSE
  )
  
  # Combine datasets
  all_results <- bind_rows(landuse_results, mpb_results) %>%
    mutate(
      y_pos = outcome_num + y_offset,
      comparison = factor(comparison, levels = c("1stPreMPB vs 2ndPostMPB", "Non_Salvage_Harvested vs Salvage_Harvested")),
      significant = prob > 0.95,
      sig_label = case_when(
        prob > 0.999 ~ "***",
        prob > 0.99 ~ "**",
        prob > 0.95 ~ "*",
        TRUE ~ ""
      )
    )
  
  # -------------------------------------------------------------------------
  # Create the figure
  # -------------------------------------------------------------------------
  
  # Color palette
  colors <- c("1stPreMPB vs 2ndPostMPB" = "#2166AC", 
              "Non_Salvage_Harvested vs Salvage_Harvested" = "#B2182B")
  
  fig <- ggplot(all_results, aes(x = estimate, y = y_pos, color = comparison)) +
    # Reference line at zero
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    
    # Error bars (95% CI)
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), 
                   height = 0.15, linewidth = 0.8) +
    
    # Points (filled if significant)
    geom_point(aes(fill = ifelse(significant, as.character(comparison), "white")), 
               shape = 21, size = 4, stroke = 1.2) +
    
    # Significance stars
    geom_text(aes(label = sig_label, x = ci_high + 0.12), 
              hjust = 0, size = 5, fontface = "bold", show.legend = FALSE) +
    
    # Color scales
    scale_color_manual(values = colors, name = "Comparison") +
    scale_fill_manual(values = c(colors, "white" = "white"), guide = "none") +
    
    # Y-axis labels
    scale_y_continuous(
      breaks = c(1, 2, 3, 4, 5),
      labels = c("Root:Shoot\nRatio", "Root Mass", "Aboveground\nBiomass", 
                 "Canopy Width", "Height"),
      limits = c(0.5, 5.5)
    ) +
    
    # Theme
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      axis.text.y = element_text(size = 11),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray30"),
      plot.caption = element_text(size = 9, color = "gray50", hjust = 0)
    ) +
    
    # Labels
    labs(
      title = "Land-Use History Effects on Seedling Growth",
      subtitle = "Hypothesis tests from multivariate Bayesian models (m=20 imputations, pooled)",
      x = "Effect Size (standardized or log-scale)",
      y = "",
      caption = "Filled points indicate P > 95%. Significance: * P>95%, ** P>99%, *** P>99.9%\nPositive values indicate advantage for 1stPreMPB (blue) or Non_Salvage_Harvested (red) soils."
    ) +
    
    # Axis limits
    coord_cartesian(xlim = c(-0.8, 4.5)) +
    
    # Legend formatting
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 3, fill = NA)))
  
  # Save figure
  ggsave(output_file, fig, width = 10, height = 6, device = "pdf")
  cat("  Saved:", output_file, "\n")
  
  return(fig)
}


# =============================================================================
# ALTERNATIVE VERSION WITH ANNOTATION
# =============================================================================

create_hypothesis_figure_annotated <- function(output_file = "mv_MI_fig_hypothesis_annotated.pdf") {
  
  cat("Creating annotated hypothesis test figure...\n")
  
  # Blue (1stPreMPB vs 2ndPostMPB) - TOP position
  landuse_results <- data.frame(
    outcome = c("Height", "Canopy Width", "Aboveground\nBiomass", 
                "Root Mass", "Root:Shoot\nRatio"),
    outcome_num = c(5, 4, 3, 2, 1),
    estimate = c(2.778, 3.280, 1.238, 1.115, -0.116),
    ci_low = c(2.051, 2.588, 0.825, 0.624, -0.528),
    ci_high = c(3.505, 3.971, 1.651, 1.605, 0.296),
    prob = c(1.000, 1.000, 1.000, 0.9998, 0.316),
    comparison = "1stPreMPB vs 2ndPostMPB",
    y_offset = 0.15,
    stringsAsFactors = FALSE
  )
  
  # Red (Non_Salvage_Harvested vs Salvage_Harvested) - BOTTOM position (values already flipped)
  mpb_results <- data.frame(
    outcome = c("Height", "Canopy Width", "Aboveground\nBiomass", 
                "Root Mass", "Root:Shoot\nRatio"),
    outcome_num = c(5, 4, 3, 2, 1),
    estimate = c(1.145, 1.414, 0.545, 0.276, -0.165),
    ci_low = c(0.715, 0.946, 0.280, -0.037, -0.419),
    ci_high = c(1.575, 1.881, 0.809, 0.588, 0.089),
    prob = c(1.000, 1.000, 0.9997, 0.930, 0.134),
    comparison = "Non_Salvage_Harvested vs Salvage_Harvested",
    y_offset = -0.15,
    stringsAsFactors = FALSE
  )
  
  all_results <- bind_rows(landuse_results, mpb_results) %>%
    mutate(
      y_pos = outcome_num + y_offset,
      comparison = factor(comparison, levels = c("1stPreMPB vs 2ndPostMPB", "Non_Salvage_Harvested vs Salvage_Harvested")),
      significant = prob > 0.95,
      sig_label = case_when(
        prob > 0.999 ~ "***",
        prob > 0.99 ~ "**",
        prob > 0.95 ~ "*",
        TRUE ~ ""
      )
    )
  
  # Color palette
  colors <- c("1stPreMPB vs 2ndPostMPB" = "#2166AC", 
              "Non_Salvage_Harvested vs Salvage_Harvested" = "#B2182B")
  
  fig <- ggplot(all_results, aes(x = estimate, y = y_pos, color = comparison)) +
    # Background annotation for outcome types
    annotate("rect", xmin = -Inf, xmax = Inf, 
             ymin = 2.5, ymax = 5.5, 
             fill = "gray95", alpha = 0.5) +
    annotate("text", x = 4.2, y = 4, label = "Aboveground", 
             color = "gray50", size = 3, fontface = "italic") +
    annotate("text", x = 4.2, y = 1.5, label = "Belowground/\nAllocation", 
             color = "gray50", size = 3, fontface = "italic") +
    
    # Reference line
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    
    # Error bars
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), 
                   height = 0.15, linewidth = 0.8) +
    
    # Points
    geom_point(aes(fill = ifelse(significant, as.character(comparison), "white")), 
               shape = 21, size = 4, stroke = 1.2) +
    
    # Significance stars
    geom_text(aes(label = sig_label, x = ci_high + 0.1), 
              hjust = 0, size = 5, fontface = "bold", show.legend = FALSE) +
    
    # Colors
    scale_color_manual(values = colors, name = "Comparison") +
    scale_fill_manual(values = c(colors, "white" = "white"), guide = "none") +
    
    # Y-axis labels
    scale_y_continuous(
      breaks = c(1, 2, 3, 4, 5),
      labels = c("Root:Shoot\nRatio", "Root Mass", "Aboveground\nBiomass", 
                 "Canopy Width", "Height"),
      limits = c(0.5, 5.5)
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      axis.text.y = element_text(size = 11),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray30"),
      plot.caption = element_text(size = 9, color = "gray50", hjust = 0)
    ) +
    
    labs(
      title = "Land-Use History Effects on Seedling Growth",
      subtitle = "Hypothesis tests from multivariate Bayesian models (m=20 imputations, pooled)",
      x = "Effect Size (standardized or log-scale)",
      y = "",
      caption = "Filled points indicate P > 95%. Significance: * P>95%, ** P>99%, *** P>99.9%\nPositive values indicate advantage for 1stPreMPB (blue) or Non_Salvage_Harvested (red) soils.\nLog-scale effects for biomass outcomes: 1 unit ~ 2.7x difference."
    ) +
    
    coord_cartesian(xlim = c(-0.8, 4.7)) +
    
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 3, fill = NA)))
  
  ggsave(output_file, fig, width = 10, height = 6, device = "pdf")
  cat("  Saved:", output_file, "\n")
  
  return(fig)
}


# =============================================================================
# RUN
# =============================================================================

fig1 <- create_hypothesis_figure()
fig2 <- create_hypothesis_figure_annotated()

cat("\nHypothesis comparison figures created!\n")




