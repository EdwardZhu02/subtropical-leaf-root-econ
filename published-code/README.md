## Code and data for: Non-correlated Variation of Leaf and Fine Root Traits in Subtropical Forest Plants

#### Datasets

- `rawdata/traits-individual.xlsx`: the original dataset containing raw trait values at individual-plant level
- `rawdata/Phylomaker_genus_family.csv`: a table matching plant family names to genus names, used for phylotree construction

#### Instructions for replicating the results

1. Install packages following `sessioninfo.txt`.

2. Run all analysis (`1-data-preprocess.R`, `2-phylotree-construction.R`, `3-calculate-PCA-nonphylo.R`, `4-calculate-PCA-phylo.R`) in order to generate necessary output data for plotting. The intermediate data are stored in `outdata` folder.

3. Run the plotting script (for example, `Fig1A-F.R`) to generate all the figures. The plots are stored in `outplts` folder.

