################################################################################
#                                                                              #
#     MANAGEMENT CONDITION COMPARISON FIGURE                                   #
#                                                                              #
#     Compares 1stPostMPB, 2ndPostMPB, and Recent_Cut conditions               #
#     Uses manual y-positioning to ensure correct alignment                    #
#                                                                              #
################################################################################

## Created with the assistance of Claude AI Opus 4.5

library(tidyverse)


# =============================================================================
# MANAGEMENT CONDITION COMPARISONS
# =============================================================================

outcomes <- c("Height", "MeanCanWidth", "logAbovegroundBiomass", 
              "logRootDryMass", "logRootToShootRatio")

# -----------------------------------------------------------------------------
# Comparison 1: 1stPostMPB vs 2ndPostMPB
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("COMPARISON 1: 1stPostMPB vs 2ndPostMPB\n")
cat("============================================================\n")

for (outcome in outcomes) {
  hyp_string <- paste0(outcome, "_Treatment_LandUseHistory1stPostMPB > ",
                       outcome, "_Treatment_LandUseHistory2ndPostMPB")
  test_custom_hypothesis(
    mv_results_MI$fits$landuse,
    hyp_string,
    paste("1stPostMPB > 2ndPostMPB for", outcome)
  )
}

# -----------------------------------------------------------------------------
# Comparison 2: Recent_Cut vs 2ndPostMPB  
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("COMPARISON 2: Recent_Cut vs 2ndPostMPB\n")
cat("============================================================\n")

for (outcome in outcomes) {
  hyp_string <- paste0(outcome, "_Treatment_LandUseHistoryRecent_Cut > ",
                       outcome, "_Treatment_LandUseHistory2ndPostMPB")
  test_custom_hypothesis(
    mv_results_MI$fits$landuse,
    hyp_string,
    paste("Recent_Cut > 2ndPostMPB for", outcome)
  )
}

# -----------------------------------------------------------------------------
# Comparison 3: Recent_Cut vs 1stPostMPB
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat("COMPARISON 3: Recent_Cut vs 1stPostMPB\n")
cat("============================================================\n")

for (outcome in outcomes) {
  hyp_string <- paste0(outcome, "_Treatment_LandUseHistoryRecent_Cut > ",
                       outcome, "_Treatment_LandUseHistory1stPostMPB")
  test_custom_hypothesis(
    mv_results_MI$fits$landuse,
    hyp_string,
    paste("Recent_Cut > 1stPostMPB for", outcome)
  )
}


# =============================================================================
# CREATE MANAGEMENT COMPARISON FIGURE
# =============================================================================

create_management_figure <- function(output_file = "mv_MI_fig_management_comparison.pdf") {
  
  cat("Creating management condition comparison figure...\n")
  
  # -------------------------------------------------------------------------
  # Input results from hypothesis tests
  # -------------------------------------------------------------------------
  
  # Create each comparison separately with explicit y_offset
  # Purple (1stPostMPB vs 2ndPostMPB) - TOP position (y_offset = +0.22)
  comp1 <- data.frame(
    outcome = c("Height", "Canopy Width", "Aboveground\nBiomass", "Root Mass", "Root:Shoot\nRatio"),
    outcome_num = c(5, 4, 3, 2, 1),
    estimate = c(0.715, 0.727, 0.186, 0.528, 0.261),
    ci_low = c(0.187, 0.184, -0.165, 0.103, -0.106),
    ci_high = c(1.244, 1.270, 0.536, 0.952, 0.627),
    prob = c(0.987, 0.986, 0.813, 0.980, 0.890),
    comparison = "1stPostMPB vs 2ndPostMPB",
    y_offset = 0.22,
    stringsAsFactors = FALSE
  )
  
  # Green (Recent_Cut vs 2ndPostMPB) - MIDDLE position (y_offset = 0)
  comp2 <- data.frame(
    outcome = c("Height", "Canopy Width", "Aboveground\nBiomass", "Root Mass", "Root:Shoot\nRatio"),
    outcome_num = c(5, 4, 3, 2, 1),
    estimate = c(1.229, 1.472, 0.841, 1.378, 0.421),
    ci_low = c(0.664, 0.881, 0.468, 0.931, 0.027),
    ci_high = c(1.794, 2.063, 1.214, 1.826, 0.815),
    prob = c(0.9998, 1.000, 0.9999, 1.000, 0.968),
    comparison = "Recent_Cut vs 2ndPostMPB",
    y_offset = 0,
    stringsAsFactors = FALSE
  )
  
  # Orange (Recent_Cut vs 1stPostMPB) - BOTTOM position (y_offset = -0.22)
  comp3 <- data.frame(
    outcome = c("Height", "Canopy Width", "Aboveground\nBiomass", "Root Mass", "Root:Shoot\nRatio"),
    outcome_num = c(5, 4, 3, 2, 1),
    estimate = c(0.514, 0.745, 0.656, 0.851, 0.161),
    ci_low = c(0.023, 0.229, 0.321, 0.461, -0.181),
    ci_high = c(1.005, 1.260, 0.990, 1.241, 0.502),
    prob = c(0.958, 0.992, 0.9996, 0.9998, 0.791),
    comparison = "Recent_Cut vs 1stPostMPB",
    y_offset = -0.22,
    stringsAsFactors = FALSE
  )
  
  # Combine
  all_results <- bind_rows(comp1, comp2, comp3) %>%
    mutate(
      y_pos = outcome_num + y_offset,
      comparison = factor(comparison, levels = c("1stPostMPB vs 2ndPostMPB",
                                                  "Recent_Cut vs 2ndPostMPB",
                                                  "Recent_Cut vs 1stPostMPB")),
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
  colors <- c("1stPostMPB vs 2ndPostMPB" = "#7570B3",
              "Recent_Cut vs 2ndPostMPB" = "#1B9E77",
              "Recent_Cut vs 1stPostMPB" = "#D95F02")
  
  fig <- ggplot(all_results, aes(x = estimate, y = y_pos, color = comparison)) +
    # Reference line at zero
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    
    # Error bars
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), 
                   height = 0.15, linewidth = 0.8) +
    
    # Points - filled if significant, color matches comparison
    geom_point(aes(fill = ifelse(significant, as.character(comparison), "white")), 
               shape = 21, size = 3.5, stroke = 1.1) +
    
    # Significance stars - positioned at ci_high + offset
    geom_text(aes(label = sig_label, x = ci_high + 0.06), 
              hjust = 0, size = 4, fontface = "bold", show.legend = FALSE) +
    
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
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(t = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      plot.caption = element_text(size = 8, color = "gray50", hjust = 0)
    ) +
    
    # Labels
    labs(
      title = "Management Condition Effects on Seedling Growth",
      subtitle = "Pairwise comparisons among disturbed land-use histories (m=20 imputations, pooled)",
      x = "Effect Size Difference (standardized or log-scale)",
      y = "",
      caption = paste0(
        "Filled points indicate P > 95%. Significance: * P>95%, ** P>99%, *** P>99.9%\n",
        "Positive values indicate first condition in comparison is greater.\n",
        "Log-scale for biomass outcomes: 0.7 ~ 2x difference, 1.0 ~ 2.7x difference."
      )
    ) +
    
    # Axis limits
    coord_cartesian(xlim = c(-0.5, 2.3)) +
    
    # Legend formatting
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 3, fill = NA)))
  
  # Save as PDF
  ggsave(output_file, fig, width = 9, height = 6, device = "pdf")
  cat("  Saved:", output_file, "\n")
  
  return(fig)
}


# =============================================================================
# TWO-PANEL VERSION
# =============================================================================

create_management_figure_twopanel <- function(output_file = "mv_MI_fig_management_twopanel.pdf") {
  
  cat("Creating two-panel management figure...\n")
  
  library(patchwork)
  
  # Panel A data: Comparisons to 2ndPostMPB
  # Purple (1stPostMPB vs 2ndPostMPB) - TOP position
  comp1_a <- data.frame(
    outcome = c("Height", "Canopy Width", "Aboveground\nBiomass", "Root Mass", "Root:Shoot\nRatio"),
    outcome_num = c(5, 4, 3, 2, 1),
    estimate = c(0.659, 0.709, 0.122, 0.483, 0.315),
    ci_low = c(0.124, 0.158, -0.225, 0.062, -0.045),
    ci_high = c(1.193, 1.260, 0.469, 0.904, 0.675),
    prob = c(0.979, 0.982, 0.725, 0.971, 0.932),
    comparison = "1stPostMPB vs 2ndPostMPB",
    y_offset = 0.15,
    stringsAsFactors = FALSE
  )
  
  # Green (Recent_Cut vs 2ndPostMPB) - BOTTOM position
  comp2_a <- data.frame(
    outcome = c("Height", "Canopy Width", "Aboveground\nBiomass", "Root Mass", "Root:Shoot\nRatio"),
    outcome_num = c(5, 4, 3, 2, 1),
    estimate = c(1.160, 1.450, 0.764, 1.283, 0.486),
    ci_low = c(0.584, 0.847, 0.389, 0.834, 0.089),
    ci_high = c(1.735, 2.054, 1.139, 1.733, 0.884),
    prob = c(0.9994, 1.000, 0.9997, 1.000, 0.982),
    comparison = "Recent_Cut vs 2ndPostMPB",
    y_offset = -0.15,
    stringsAsFactors = FALSE
  )
  
  panel_a_data <- bind_rows(comp1_a, comp2_a) %>%
    mutate(
      y_pos = outcome_num + y_offset,
      comparison = factor(comparison, levels = c("1stPostMPB vs 2ndPostMPB",
                                                  "Recent_Cut vs 2ndPostMPB")),
      significant = prob > 0.95,
      sig_label = case_when(
        prob > 0.999 ~ "***",
        prob > 0.99 ~ "**",
        prob > 0.95 ~ "*",
        TRUE ~ ""
      )
    )
  
  colors_a <- c("1stPostMPB vs 2ndPostMPB" = "#7570B3",
                "Recent_Cut vs 2ndPostMPB" = "#1B9E77")
  
  panel_a <- ggplot(panel_a_data, aes(x = estimate, y = y_pos, color = comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), 
                   height = 0.15, linewidth = 0.8) +
    geom_point(aes(fill = ifelse(significant, as.character(comparison), "white")), 
               shape = 21, size = 3.5, stroke = 1.1) +
    geom_text(aes(label = sig_label, x = ci_high + 0.08), 
              hjust = 0, size = 4, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = colors_a, name = "") +
    scale_fill_manual(values = c(colors_a, "white" = "white"), guide = "none") +
    scale_y_continuous(
      breaks = c(1, 2, 3, 4, 5),
      labels = c("Root:Shoot\nRatio", "Root Mass", "Aboveground\nBiomass", 
                 "Canopy Width", "Height"),
      limits = c(0.5, 5.5)
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      axis.text.y = element_text(size = 9)
    ) +
    labs(
      title = "A) Comparisons to 2ndPostMPB",
      subtitle = "How much better are other conditions?",
      x = "Effect Size Difference",
      y = ""
    ) +
    coord_cartesian(xlim = c(-0.5, 2.3))
  
  # Panel B data: Recent_Cut vs 1stPostMPB
  panel_b_data <- data.frame(
    outcome = c("Height", "Canopy Width", "Aboveground\nBiomass", 
                "Root Mass", "Root:Shoot\nRatio"),
    outcome_num = c(5, 4, 3, 2, 1),
    estimate = c(0.501, 0.742, 0.642, 0.800, 0.171),
    ci_low = c(0.010, 0.222, 0.313, 0.419, -0.166),
    ci_high = c(0.992, 1.261, 0.971, 1.181, 0.508),
    prob = c(0.955, 0.991, 0.9995, 0.9997, 0.809),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      significant = prob > 0.95,
      sig_label = case_when(
        prob > 0.999 ~ "***",
        prob > 0.99 ~ "**",
        prob > 0.95 ~ "*",
        TRUE ~ ""
      )
    )
  
  panel_b <- ggplot(panel_b_data, aes(x = estimate, y = outcome_num)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), 
                   height = 0.2, color = "#D95F02", linewidth = 0.8) +
    geom_point(aes(fill = significant), 
               shape = 21, size = 3.5, stroke = 1.1, color = "#D95F02") +
    geom_text(aes(label = sig_label, x = ci_high + 0.05), 
              hjust = 0, size = 4, fontface = "bold", color = "#D95F02") +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#D95F02"),
                      guide = "none") +
    scale_y_continuous(
      breaks = c(1, 2, 3, 4, 5),
      labels = c("Root:Shoot\nRatio", "Root Mass", "Aboveground\nBiomass", 
                 "Canopy Width", "Height"),
      limits = c(0.5, 5.5)
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 9)
    ) +
    labs(
      title = "B) Recent_Cut vs 1stPostMPB",
      subtitle = "",
      x = "Effect Size Difference",
      y = ""
    ) +
    coord_cartesian(xlim = c(-0.5, 1.5))
  
  # Combine
  fig <- panel_a / panel_b +
    plot_annotation(
      title = "Management Condition Effects on Seedling Growth",
      subtitle = "Pairwise comparisons from multivariate Bayesian models (m=20 imputations)",
      caption = "Filled points = P > 95%. * P>95%, ** P>99%, *** P>99.9%. Positive = first condition better.",
      theme = theme(
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray30"),
        plot.caption = element_text(size = 8, color = "gray50")
      )
    )
  
  ggsave(output_file, fig, width = 8, height = 8, device = "pdf")
  cat("  Saved:", output_file, "\n")
  
  return(fig)
}


# =============================================================================
# RUN
# =============================================================================

fig1 <- create_management_figure()
fig2 <- create_management_figure_twopanel()

cat("\nManagement comparison figures created!\n")


