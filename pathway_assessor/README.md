## Usage
- [pathwayassessor.harmonic](#harmonic)
   - Pathway overrepresentation and underrepresentation scores derived from harmonic averaging of p-values
- [pathwayassessor.geometric](#geometric)
   - Pathway overrepresentation and underrepresentation scores derived from geometric averaging of p-values
- [pathwayassessor.minpval](#minpval)
   - Negative log of minimum p-values for each sample-pathway paid
- [pathwayassessor.all](#all)
   - Dictionary with harmonic, geometric, and min-p-val results

## Data Input
The expression matrix must be a dataframe with genes in rows and samples in columns. 
The rownames should be gene symbols. Rows with duplicate symbols will be averaged.

## pathwayassessor.harmonic
- [Description](#description)
- [Usage](#usage)
- [Arguments](#arguments)
- [Value](#value)
- [Notes](#notes)
- [Example](#example)


### Description

This function returns a dataframe of pathways x samples. 
Each value is an overrepresentation or underrepresentation score for a given pathway 
calculated by harmonically averaging all of the gene rank-based p-values for each gene 
that pathway. 

### Usage
```
DreamAI(data, k = 10, maxiter_MF = 10, ntree = 100,
  maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2),
  gamma_ADMIN = NA, gamma = 50, CV = FALSE,
  fillmethod = "row_mean", maxiter_RegImpute = 10,
  conv_nrmse = 1e-06, iter_SpectroFM = 40, method = c("KNN",
  "MissForest", "ADMIN", "Brinn", "SpectroFM", "RegImpute"),
  out = c("Ensemble"))
```