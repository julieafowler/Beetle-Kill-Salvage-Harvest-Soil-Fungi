#### Soil Chemistry Summary Tables ####
# Statistical tests:
#   - Kruskal-Wallis across 6 Land Use Histories
#   - Wilcoxon rank-sum (Mann-Whitney) for Salvage vs. Non_Salvage_Harvested

## Created with the assistance of Claude AI Sonnet 4.6

# ── Libraries ──────────────────────────────────────────────────────────────────
library(readxl)
library(dplyr)
library(gtsummary)
library(gt)
library(broom)
library(rstatix)
library(tidyr)
library(tibble)
library(purrr)

setwd()

# ── Load & prep data ───────────────────────────────────────────────────────────
df <- read_excel("SoilChemistry_forTable.xlsx",
                 col_types = c(
                   "text",    # Sample
                   "text",    # Land_Use_History
                   "text",    # Salvage_Harvest_Status
                   "numeric", # Total C (%)
                   "numeric", # Total N (%)
                   "numeric", # Na (mg/L)
                   "numeric", # NH4 (mg/L)
                   "numeric", # K (mg/L)
                   "numeric", # Mg (mg/L)
                   "numeric", # Ca (mg/L)
                   "numeric", # Cl (mg/L)
                   "numeric", # NO3 (mg/L)
                   "numeric", # PO4 (mg/L)
                   "numeric", # SO4 (mg/L)
                   "numeric"  # pH
                 )) %>%
  rename(
    Total_C = `Total C (%)`,
    Total_N = `Total N (%)`,
    Na_     = `Na (mg/L)`,
    NH4     = `NH4  (mg/L)`,
    K_      = `K  (mg/L)`,
    Mg      = `Mg  (mg/L)`,
    Ca      = `Ca  (mg/L)`,
    Cl      = `Cl  (mg/L)`,
    NO3     = `NO3  (mg/L)`,
    PO4     = `PO4  (mg/L)`,
    SO4     = `SO4  (mg/L)`
  ) %>%
  mutate(
    across(c(pH, Total_C, Total_N, Na_, NH4, K_, Mg, Ca, Cl, NO3, PO4, SO4),
           as.double),
    Land_Use_History = factor(
      Land_Use_History,
      levels = c("Old Growth", "1st Pre-MPB", "2nd Pre-MPB",
                 "1st Post-MPB", "2nd Post-MPB", "Recent Cut")
    ),
    Salvage_Harvest_Status = factor(
      Salvage_Harvest_Status,
      levels = c("Non_Salvage_Harvested", "Salvage_Harvested")
    )
  )

# Sanity check
stopifnot(all(sapply(df[response_vars <- c("pH","Total_C","Total_N","Na_","NH4","K_",
                                            "Mg","Ca","Cl","NO3","PO4","SO4")],
                     is.numeric)))

var_labels <- list(
  pH      = "pH",
  Total_C = "Total C (%)",
  Total_N = "Total N (%)",
  Na_     = "Na (mg/L)",
  NH4     = "NH4 (mg/L)",
  K_      = "K (mg/L)",
  Mg      = "Mg (mg/L)",
  Ca      = "Ca (mg/L)",
  Cl      = "Cl (mg/L)",
  NO3     = "NO3 (mg/L)",
  PO4     = "PO4 (mg/L)",
  SO4     = "SO4 (mg/L)"
)

# ── Helper: compute median (Q1, Q3) using R's default quantile type 7 ──────────
# gtsummary internally uses a different quantile type for some group sizes,
# producing wrong IQRs. We compute correct values manually and inject them.
fmt_med_iqr <- function(x, digits = 3) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  qs <- quantile(x, probs = c(0.25, 0.5, 0.75), type = 7, na.rm = TRUE)
  sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
          qs[2], qs[1], qs[3])
}

# ── TABLE 1: Summary by Land Use History (6 groups) ───────────────────────────

# Compute ALL statistics manually: median (Q1,Q3) per group + KW p-value
luh_levels <- levels(df$Land_Use_History)

tbl1_stats <- map_dfr(response_vars, function(var) {
  row <- tibble(variable = var)
  for (luh in luh_levels) {
    vals <- df[[var]][df$Land_Use_History == luh]
    col_name <- paste0("new_stat_", which(luh_levels == luh))  # prefix avoids collision
    row[[col_name]] <- fmt_med_iqr(vals)
  }
  row$new_p.value <- kruskal.test(x = df[[var]], g = df$Land_Use_History)$p.value
  row
})

# Build a minimal tbl_summary for structure only, then overwrite stat columns
tbl1 <- df %>%
  select(Land_Use_History, all_of(response_vars)) %>%
  tbl_summary(
    by        = Land_Use_History,
    type      = all_continuous() ~ "continuous",
    statistic = all_continuous() ~ "{median}",  # placeholder; overwritten below
    digits    = all_continuous() ~ 3,
    label     = var_labels,
    missing   = "no"
  ) %>%
  modify_table_body(~ {
    .x %>%
      left_join(tbl1_stats, by = "variable") %>%
      mutate(
        stat_1 = ifelse(row_type == "label", new_stat_1, NA_character_),
        stat_2 = ifelse(row_type == "label", new_stat_2, NA_character_),
        stat_3 = ifelse(row_type == "label", new_stat_3, NA_character_),
        stat_4 = ifelse(row_type == "label", new_stat_4, NA_character_),
        stat_5 = ifelse(row_type == "label", new_stat_5, NA_character_),
        stat_6 = ifelse(row_type == "label", new_stat_6, NA_character_),
        p.value = ifelse(row_type == "label", new_p.value, NA_real_)
      ) %>%
      select(-starts_with("new_"))
  }) %>%
  modify_header(
    label   ~ "**Soil Chemistry Variable**",
    p.value ~ "**Kruskal-Wallis p**"
  ) %>%
  modify_fmt_fun(
    p.value ~ function(x) style_pvalue(x, digits = 3)
  ) %>%
  modify_spanning_header(
    all_stat_cols() ~ "**Land Use History** — Median (Q1, Q3)"
  ) %>%
  bold_labels() %>%
  modify_caption(
    "**Table 1.** Soil chemistry by Land Use History (wet soils; n=36).
     Values are median (Q1, Q3). Kruskal-Wallis test p-values (asymptotic);
     pairwise comparisons in Table 3."
  ) %>%
  modify_footnote(
    all_stat_cols() ~
      "n = 6 per group for all variables except Na, NH4, K, Mg, Ca, Cl, NO3, PO4, and SO4,
       for which 2nd Pre-MPB has n = 3 due to missing values."
  )

# ── TABLE 2: Summary by Salvage Harvest Status (2 groups) ─────────────────────

salvage_levels <- levels(df$Salvage_Harvest_Status)

tbl2_stats <- map_dfr(response_vars, function(var) {
  row <- tibble(variable = var)
  for (sl in salvage_levels) {
    vals <- df[[var]][df$Salvage_Harvest_Status == sl]
    col_name <- paste0("new_stat_", which(salvage_levels == sl))  # prefix avoids collision
    row[[col_name]] <- fmt_med_iqr(vals)
  }
  g1 <- df[[var]][df$Salvage_Harvest_Status == salvage_levels[1]]
  g2 <- df[[var]][df$Salvage_Harvest_Status == salvage_levels[2]]
  row$new_p.value <- wilcox.test(g1[!is.na(g1)], g2[!is.na(g2)],
                                  exact = FALSE)$p.value
  row
})

tbl2 <- df %>%
  select(Salvage_Harvest_Status, all_of(response_vars)) %>%
  tbl_summary(
    by        = Salvage_Harvest_Status,
    type      = all_continuous() ~ "continuous",
    statistic = all_continuous() ~ "{median}",  # placeholder; overwritten below
    digits    = all_continuous() ~ 3,
    label     = var_labels,
    missing   = "no"
  ) %>%
  modify_table_body(~ {
    .x %>%
      left_join(tbl2_stats, by = "variable") %>%
      mutate(
        stat_1 = ifelse(row_type == "label", new_stat_1, NA_character_),
        stat_2 = ifelse(row_type == "label", new_stat_2, NA_character_),
        p.value = ifelse(row_type == "label", new_p.value, NA_real_)
      ) %>%
      select(-starts_with("new_"))
  }) %>%
  modify_header(
    label   ~ "**Soil Chemistry Variable**",
    p.value ~ "**Wilcoxon p**"
  ) %>%
  modify_fmt_fun(
    p.value ~ function(x) style_pvalue(x, digits = 3)
  ) %>%
  modify_spanning_header(
    all_stat_cols() ~ "**Harvest Status** — Median (Q1, Q3)"
  ) %>%
  bold_labels() %>%
  modify_caption(
    "**Table 2.** Soil chemistry by Salvage Harvest Status (wet soils; n=36).
     Values are median (Q1, Q3). Wilcoxon rank-sum (Mann-Whitney U) test."
  ) %>%
  modify_footnote(
    all_stat_cols() ~
      "n = 18 for pH, Total C, and Total N; n = 15 (Non-Salvage Harvested) and
       n = 18 (Salvage Harvested) for all other variables due to missing values
       in the 2nd Pre-MPB group."
  )

# ── Export ─────────────────────────────────────────────────────────────────────
tbl1 %>% as_gt() %>% gtsave("Table1_LandUseHistory.docx")
tbl2 %>% as_gt() %>% gtsave("Table2_SalvageHarvest.docx")

tbl1 %>% as_gt() %>% gtsave("Table1_LandUseHistory.html")
tbl2 %>% as_gt() %>% gtsave("Table2_SalvageHarvest.html")

cat("\nDone! Tables saved.\n")
