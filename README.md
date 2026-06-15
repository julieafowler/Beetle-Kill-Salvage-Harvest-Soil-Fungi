Scripts and data files accompanying the manuscript:

> **Salvage Harvest Timing Shapes Lodgepole Pine Seedling Performance Through Soil Legacy Effects**  
> Julie A. Fowler, Timothy S. Fegel, Emily K. Bechtold, Kya M. Sparks, Sarah J. Hart, Louise H. Comas, David M. Barnard, Charles C. Rhoades, Michael J. Wilkins 
> *Journal Name*, Year. DOI: [link]

---

## Overview

This repository contains analysis scripts for: (1) field surveys of soil fungal communities, and (2) a greenhouse experiment examining soil & root tip fungal communities and plant outcomes using soils from Field sites that span a gradient from beetle-impacted old growth forests to stands clear-cut before the outbreak to beetle-killed stands salvage harvested at varying intervals post-infestation.\
\
Note: Some language is inconsistent throughout the scripts/files. The land use history "Old Growth" may be refered to as "Old Forest" in places. "Early MPB" is referred to as "1st Post-MPB" and "Late MPB" is referred to as "2nd Post-MPB" in the resulting manuscript. 

---

## Repository Structure

Scripts are numbered in order of use in the analytical workflow:

| # | Folder | Description |
|---|--------|-------------|
| 1 | `Field_ITS_Analysis` | Field soil fungal (ITS) community analysis |
| 2 |
| 3 | `Greenhouse_Vegetation_Analysis` | Greenhouse plant response analysis |
| 4 | 
| 5 | `Greenhouse_16S_Analysis` | Greenhouse bacterial & archaeal (16S) community analysis |
| 6 | `Greenhouse_ITS_Analysis` | Greenhouse fungal (ITS) alpha/beta diversity, guild assignments, etc. |
| 7 | `Greenhouse_ITS_EMFOnly_Analysis` | Greenhouse ectomycorrhizal fungal community analyses |
| 8 | `Greenhouse_MAASLIN3_PlantOutcomes` | Greenhouse fungal MaAsLin3 differential abundance modeling |
| 9 | `Greenhouse_sPLS_PlantOutcomes` | Greenhouse fungal sparse PLS regression for plant outcomes |
| 10 | `Greenhouse_DiagnosticScript_forBayesian` | Data diagnostics & missingness assessment for downstream Bayesian work |
| 11 | `Greenhouse_SoilChemistry_PCA_forBayesian` | Soil PCA for Bayesian model inputs |
| 12 | `Greenhouse_SoilChemistry_PCA_Tables` | PCA tables for supplemental materials |
| 13 | `Greenhouse_Bayesian_Analysis_2FactorMPBStatus` | Univariate Bayesian models (non-MPB vs. MPB) |
| 14 | `Greenhouse_Bayesian_Analysis_6FactorLandUseHistory` | Univariate Bayesian models (6 land-use histories) |
| 15 | `Greenhouse_Bayesian_MultivariateAnalysis_MI` | Multivariate Bayesian models with selected predictors |
| 16 | `Greenhouse_Bayesian_CoefficientTables` | Coefficient tables for supplemental materials |
| 17 | `Greenhouse_BayesianMultivariate_HypothesisFigure` | Main non-MPB vs. MPB hypothesis figure generation |
| 18 | `Greenhouse_BayesianMultivariate_ManagementComparisonFigure` | Management comparison figure generation |

Each folder contains the relevant scripts and associated data files.

---

## Data Availability

Raw sequence data (field ITS; greenhouse 16S and ITS) are available at NCBI SRA under BioProject PRJNA1440534.

---

## Contact

Julie A. Fowler — julie.fowler@colostate.edu\
ORCID — https://orcid.org/0000-0003-2180-561X
