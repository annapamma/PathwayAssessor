# PathwayAssessor

- [Authors](#authors)
- [Overview](#Overview)
- [Installation](#installation)
- [Usage](#usage)
- [References](#references)
- [Contributions](#contributions)

## Authors
Boris Reva<sup>1</sup>, Anna Calinawan<sup>1</sup>

<sup>1</sup>Icahn School of Medicine at Mount Sinai (USA)

## Overview

Comparison of tumors by activities of molecular pathways can help nominate disease driver pathways, 
identify targets for therapeutic intervention, and determine clinically relevant subtypes and diagnostic biomarkers. 
We introduce a new method for calculating overrepresentation and underrepresentation scores to assess a pathway’s 
activation and suppression, respectively.
The novelty of this approach is in independent assessment of overrepresentation and underrepresentation of pathway’s 
genes from both top and bottom of the ranked gene list of a given sample.


## Installation
Requires: Python >= 3.6

To install the latest version of PathwayAssessor, run:

```
pip install pathway-assessor
```

## Usage
To see more information about the main functions, visit the [ReadMe included with the package code](https://github.com/annapamma/PathwayAssessor/tree/master/pathway_assessor).

```
import pathway_assessor as pa

pa.all()
pa.harmonic()
pa.geometric()
pa.min_p_val()
```


## References 
1. Barbie, D.A. , et al. “Systematic RNA interference reveals…”,  Nature. 2009 Nov 5;462(7269):108-12.
2. Abril-Rodrigues, G., Ribas A., “SnapShot: Immune Checkpoint Inhibitors”, Cancer Cell 31, June 12, 2017.
3. Aran, D., Hu Z., Butte, A. “xCell: digitally portraying the tissue…”, Genome Biology. (2017) 18:220. 

### Contributions

If you find small bugs, larger issues, or have suggestions, please email the maintainer at <anna.calinawan@mssm.edu>.
Contributions (via pull requests or otherwise) are welcome.
