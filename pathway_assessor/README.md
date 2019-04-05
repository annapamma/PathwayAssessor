## Usage
- [pathwayassessor.harmonic](#harmonic)
   - Pathway overrepresentation and underrepresentation scores derived from harmonic averaging of p-values
- [pathwayassessor.geometric](#geometric)
   - Pathway overrepresentation and underrepresentation scores derived from geometric averaging of p-values
- [pathwayassessor.minpval](#minpval)
   - Negative log of minimum p-values for each sample-pathway paid
- [pathwayassessor.all](#all)
   - Dictionary with harmonic, geometric, and min-p-val results
- [Arguments](#arguments)
   - Arguments for pathway_assessor.harmonic, pathway_assessor.geometric, and pathway_assessor.min_p_val are the same.

## Data Input
The expression matrix must be a dataframe with genes in rows and samples in columns. 
The rownames should be gene symbols. Rows with duplicate symbols will be averaged.


## pathwayassessor.harmonic
- [Description](#description)
- [Usage](#usage)
- [Value](#value)
- [Notes](#notes)
- [Example](#example)


### Description

This function returns a dataframe of pathways x samples. 
Each value is an overrepresentation or underrepresentation score for a given pathway 
calculated by harmonically averaging all of the gene rank-based p-values for each gene 
that pathway. Values reported are the negative log of the harmonic average.

### Usage
```
pathwayassessor.harmonic(
        expression_table,
        pathways=None,
        db='kegg',
        ascending=True
)
```

### Arguments
For all, harmonic, geometric, and min_p_val.

| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|
| expression_table	       |	          | expression data frame with genes in rows and samples in columns.
| pathways         | None           |a dictionary with pathway names as keys and sets or lists of genes as values
| db 	       |	kegg	            |string indicating which of the included pathway databases to use. Options include: 'kegg', 'reactome', 'hmdb_smpdb', 'hallmark'
| ascending  		       | True	           | boolean for what direction to sort the expression table

Additional arguments for all:
| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|
| geometric	       |True	          | boolean of whether to calculate and include geometric average
| min_p_val         | True           |boolean of whether to calculate and include min p value