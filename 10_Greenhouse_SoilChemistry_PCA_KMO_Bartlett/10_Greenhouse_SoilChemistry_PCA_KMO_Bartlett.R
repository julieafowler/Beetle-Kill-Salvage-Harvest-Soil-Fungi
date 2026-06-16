###############################################################################
## SOIL CHEMISTRY PCA — FACTORABILITY DIAGNOSTICS: KMO + BARTLETT'S
###############################################################################
## Created with the assistance of Claude AI Opus 4.8
## PURPOSE -------------------------------------------------------------------
## These two tests check ONE thing: is the soil chemistry correlated enough that
## reducing it to principal components is a sensible thing to do at all? 
##
## A plain-language reference describing both tests (and the R functions used
## below) is the easystats 'performance' documentation:
##   https://easystats.github.io/performance/reference/check_factorstructure.html
## The underlying methods papers are Bartlett (1951), Kaiser (1970), and
## Kaiser & Rice (1974).
##
## WHY n = 17 -----------------------------------------------------------------
## These tests are run on the OBSERVED soil chemistry only (complete cases).
###############################################################################


## ---- Packages --------------------------------------------------------------
## 'psych' provides both KMO() and cortest.bartlett(). Install once if needed:
##   install.packages("psych")
library(psych)


## ---- 1. Load data and isolate the soil chemistry block ---------------------
setwd()
df <- read.delim("Bayesian_AllData_ExcludedThoseWithOnlyNAsAcross_RemovedProblemColumns.txt",
                 row.names = 1)

soil_vars <- c("Total_C_percent", "Total_N_percent",
               "Na_waterextract_mgperL", "NH4_waterextract_mgperL",
               "K_waterextract_mgperL",  "Mg_waterextract_mgperL",
               "Ca_waterextract_mgperL", "Cl_waterextract_mgperL",
               "NO3_waterextract_mgperL","PO4_waterextract_mgperL",
               "SO4_waterextract_mgperL")

soil <- df[, soil_vars]


## ---- 2. Report missingness structure ---------------------------------------
## Important nuance: the variables do NOT all have the same n. Total C and N were
## measured on more samples (n = 20 vs. n = 17) than the water-extractable panel, so listwise
## complete cases (what the PCA actually uses) is smaller than either.
cat("\nObserved n per soil variable:\n")
print(colSums(!is.na(soil)))                       # expect C/N = 20, extractables = 17
cat("\nSamples with C/N measured but missing the extractable panel:",
    sum(!is.na(soil$Total_C_percent) & is.na(soil$Na_waterextract_mgperL)), "\n")


## ---- 3. Build the complete-case matrix the PCA runs on ---------------------
## complete.cases() keeps only rows with ALL 11 variables present, because PCA
## needs every variable in a row. This is the n adequacy tests use.
soil_cc <- soil[complete.cases(soil), ]
n <- nrow(soil_cc)
p <- ncol(soil_cc)

cat("\nListwise complete cases used for PCA: n =", n, "across p =", p, "variables\n")
cat("Observation-to-variable ratio:", round(n / p, 2), ": 1\n")
## A low ratio (here ~1.6:1, below the ~3:1 rule-of-thumb floor)


## ---- 4. Correlation matrix -------------------------------------------------
## Both tests operate on the correlation matrix R (correlation-based PCA, i.e.
## variables standardized — matches scale. = TRUE in the main PCA script).
R <- cor(soil_cc)


## ---- 5. BARTLETT'S TEST OF SPHERICITY --------------------------------------
## QUESTION: are the variables correlated at all, or is R indistinguishable from
##           an identity matrix (all correlations = 0)?
##   H0: R = I  -> variables unrelated -> PCA is pointless.
##
## HOW TO READ IT: a significant result (p < 0.05) is necessary but is a LOW bar
## — with real data it is almost always significant, so passing it is reassuring
## but not, by itself, strong evidence. (cortest.bartlett needs n when given R.)
bart <- cortest.bartlett(R, n = n)
cat("\n--- Bartlett's test of sphericity ---\n")
cat("chi-square =", round(bart$chisq, 1),
    " df =", bart$df,
    " p =", format(bart$p.value, scientific = TRUE, digits = 3), "\n")
cat(ifelse(bart$p.value < 0.05,
           "  -> Reject H0: variables are correlated; PCA is justified.\n",
           "  -> Fail to reject H0: PCA not justified.\n"))


## ---- 6. KAISER-MEYER-OLKIN (KMO) MEASURE OF SAMPLING ADEQUACY ---------------
## QUESTION: for each pair of variables, is their correlation "clean" (shared
##           with the whole set, good for PCA) or is it really a private partial
##           correlation between just those two (bad for PCA)? KMO summarizes the
##           ratio of plain correlations to partial correlations.
##
## RANGE 0-1. Kaiser & Rice (1974) labels:
##   >= .90 marvelous | .80s meritorious | .70s middling
##   .60s mediocre    | .50s miserable   | < .50 unacceptable (drop the variable)
##
## KMO() returns $MSA (overall) and $MSAi (one value per variable). Variables
## with individual MSA < .50 are the weak links — here, the base cations.
kmo <- KMO(R)
cat("\n--- Kaiser-Meyer-Olkin (KMO) ---\n")
cat("Overall KMO =", round(kmo$MSA, 3), "\n")
cat("Per-variable MSA (sorted, low = weakest):\n")
print(round(sort(kmo$MSAi), 3))


## ---- 7. Basis for retaining three components --------------------
## KMO/Bartlett do not choose the number of PCs. The three-component solution is
## justified by Kaiser's criterion: retain components with eigenvalue > 1 (i.e.
## that explain more than a single average variable's worth of variance).
ev <- eigen(R, only.values = TRUE)$values
cat("\n--- Eigenvalues (Kaiser's criterion: keep those > 1) ---\n")
print(round(ev, 3))
cat("Components with eigenvalue > 1:", sum(ev > 1), "\n")


