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
    # Explicitly coerce all response columns to double — belt and suspenders
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

# Sanity check — all response cols should be numeric
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

# ── TABLE 1: Summary by Land Use History (6 groups) ───────────────────────────
# Kruskal-Wallis p-values computed manually (x/g form = asymptotic, no overflow)
# then injected via modify_table_body to avoid gtsummary's internal re-run.

kw_pvals <- map_dfr(response_vars, function(var) {
  tibble(
    variable = var,
    p.value  = kruskal.test(x = df[[var]], g = df$Land_Use_History)$p.value
  )
})

tbl1 <- df %>%
  select(Land_Use_History, all_of(response_vars)) %>%
  tbl_summary(
    by        = Land_Use_History,
    type      = all_continuous() ~ "continuous",  # force continuous, no auto-detect
    statistic = all_continuous() ~ "{median} ({p25}, {p75})",
    digits    = all_continuous() ~ 3,
    label     = var_labels,
    missing   = "no"
  ) %>%
  modify_table_body(
    ~ .x %>%
      left_join(kw_pvals, by = "variable") %>%
      mutate(p.value = ifelse(row_type == "label", p.value, NA_real_))
  ) %>%
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
tbl2 <- df %>%
  select(Salvage_Harvest_Status, all_of(response_vars)) %>%
  tbl_summary(
    by        = Salvage_Harvest_Status,
    type      = all_continuous() ~ "continuous",
    statistic = all_continuous() ~ "{median} ({p25}, {p75})",
    digits    = all_continuous() ~ 3,
    label     = var_labels,
    missing   = "no"
  ) %>%
  add_p(
    test       = all_continuous() ~ "wilcox.test",
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  modify_header(
    label   ~ "**Soil Chemistry Variable**",
    p.value ~ "**Wilcoxon p**"
  ) %>%
  modify_spanning_header(
    all_stat_cols() ~ "**Harvest Status** — Median (Q1, Q3)"
  ) %>%
  bold_labels() %>%
  modify_caption(
    "**Table 2.** Soil chemistry by Salvage Harvest Status (wet soils; n=36).
     Values are median (Q1, Q3). Wilcoxon rank-sum (Mann-Whitney U) test."
  )

# ── Export ─────────────────────────────────────────────────────────────────────
tbl1 %>% as_gt() %>% gtsave("Table1_LandUseHistory.docx")
tbl2 %>% as_gt() %>% gtsave("Table2_SalvageHarvest.docx")

tbl1 %>% as_gt() %>% gtsave("Table1_LandUseHistory.html")
tbl2 %>% as_gt() %>% gtsave("Table2_SalvageHarvest.html")

cat("\nDone! Tables saved.\n")
