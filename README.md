# PathwayAssessor

- [Authors](#authors)
- [Overview](#Overview)
- [Installation](#installation)
- [Contributions](#contributions)

## Authors
Boris Reva<sup>1</sup>, Anna Calinawan<sup>1</sup>

<sup>1</sup>Icahn School of Medicine at Mount Sinai (USA)

## Overview

Comparison of tumors by activities of molecular pathways can help nominate disease driver pathways, 
identify targets for therapeutic intervention, and determine clinically relevant subtypes and diagnostic biomarkers. 
We introduce a new method for calculating overrepresentation and underrepresentation scores to assess a pathway’s 
activation and suppression, respectively. Similarly to the ssGSEA method, we assume that pathway activity can 
be compared across different samples based on ranking the expression of pathway genes.  
The novelty of this approach is in independent assessment of overrepresentation and underrepresentation of pathway’s 
genes from both top and bottom of the ranked gene list of a given sample.


## Installation

```
pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.python.org/pypi pathway-assessor-annapamma2==0.0.32

```

## Usage
```

pathway_assessor.all
pathway_assessor.harmonic
pathway_assessor.geometric
pathway_assessor.min_p_val

```

## References 
1. Barbie, D.A. , et al. “Systematic RNA interference reveals…”,  Nature. 2009 Nov 5;462(7269):108-12.
2. Abril-Rodrigues, G., Ribas A., “SnapShot: Immune Checkpoint Inhibitors”, Cancer Cell 31, June 12, 2017.
3. Aran, D., Hu Z., Butte, A. “xCell: digitally portraying the tissue…”, Genome Biology. (2017) 18:220. 

### Contributions

If you find small bugs, larger issues, or have suggestions, please email the maintainer at <anna.calinawan@mssm.edu>.
Contributions (via pull requests or otherwise) are welcome.