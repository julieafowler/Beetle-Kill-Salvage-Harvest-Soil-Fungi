#### Soil Fungi & Seedling Responses to Beetle Kill and Salvage Harvesting ####

Scripts and data files accompanying the manuscript:

> **[Paper Title]**  
> Author 1, Author 2, Author 3, et al.  
> *Journal Name*, Year. DOI: [link]

---

## Overview

This repository contains analysis scripts for: (1) field surveys of vegetation and soil fungal communities, and (2) a greenhouse experiment examining soil & root tip microbial communities (primarily fungal) and plant outcomes using soils from beetle-killed and salvage-harvested forest stands.

---

## Repository Structure

Scripts are numbered in order of use in the analytical workflow:

| # | Folder | Description |
|---|--------|-------------|
| 1 | `Field_ITS_Analysis` | Field soil fungal (ITS) community analysis |
| 2 | `Field_Vegetation_Boxplots` | Field vegetation data visualization |
| 3 | `Greenhouse_Vegetation_Analysis` | Greenhouse plant response analysis |
| 4 | `Greenhouse_SoilChemistryandTexture` | Greenhouse soil chemistry and texture characterization |
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

Julie A. Fowler — julie.fowler@colostate.edu
ORCID - https://orcid.org/0000-0003-2180-561X
