# Pituitary Gland Shiny App
Shiny app to visualize DE genes and miRNAs in the postnatal mouse pituitary gland.

## To run app locally
1. Clone git repository to local directory.
2. Navigate to directory containing cloned repository.
3. Install dependencies in R with `devtools::install()` (Set Working Directory to path containing `DESCRIPTION`).
4. Open `app.R` in RStudio. Click "Run App" to run locally.

## Data Browser
### Input
1. Navigate to the `Data Browser` tab.
2. Choose one option:
    1. `Type in genes/miRNAs`
        1. Enter in a comma-separated list of genes and/or miRNAs (can be a mixture of both)
    2. `Upload file`
        1. Upload newline-separated list of genes and/or miRNAs (can be a mixture of both)
        2. See `example_input_file.txt` as an example (located in `data/` directory)
    3. `Sex-biased genes and miRNAs`
        1. Choose a list of sex-biased genes and miRNAs pre-defined from Hou et al. 2022.
5. Press `Submit` to run.
### Output
#### Plot
1. Expression profile of user-inputted miRNAs and genes will be displayed for all profiled samples. Median profile is shown as filled points connected by lines. Individual replicates are shown as smaller unfilled points. Profiles are grouped by sex and age.
2. Gene and miRNA expression profiles displayed can be downloaded as `.pdf` or `.png` files.
3. For genes only, `StudyID` displayed on plots indicate study in which the given gene was previously described (see legend on the right side). All studies are cited in the `About` page.

#### DE gene/miRNA table
1. User can toggle between:
    1. `abs(FC) > 1.5, FDR < 0.05` (Only shows genes/miRNAs which pass this cutoff and are considered DE in Hou et al. 2022)
    2. `No cutoff` (Shows all genes/miRNAs)
2. DE table which is being displayed can be downloaded as `.txt` or `.csv` files.
3. For genes only, `StudyID` displayed on plots indicate study in which the given gene was previously described (see legend on the right side). All studies are cited in the `About` page.

#### Gene/miRNA expression values
1. User can toggle between:
    1. `Normalized counts` (Non-log2 counts; RUV-seq normalized to remove unwanted variation, see Hou et al. 2022 for more information)
    2. `log2(normCounts)` (log2-transformed RUV-seq normalized-counts)
2. Expression table which is being displayed can be downloaded as `.txt` or `.csv` files.

## miRNA-Gene Target Browser
### Input
1. Navigate to `miRNA-Gene Target Browser` tab.
2. Choose `miRNA` or `gene` and enter in a comma-separated list of miRNAs or genes respectively.
4. Choose filter:
    1. `DE`
        1. For sex-biased comparisons, miRNAs and genes which are sex-biased from any age are included
        2. For age-biased comparisons, miRNAs and genes which are age-biased from the SAME age comparison are included regardless of sex
    2. `All`
        1. All miRNAs and genes are considered (DE filter is ignored)
5. Press `Submit` to run.
    1. To try the different filter, switch the filter and press `Submit` to rerun.
### Output
1. Main panel displays correlation table filtered for user-inputted miRNAs/genes found in negatively correlated pairs (Spearman's rho < 0; FDR > 0.1)
2. `database`: miRNA-gene pair is computationally predicted (database: TargetScan) or experimentally validated (database: miRTarBase).
3. `gene/miRNA comparison`: Comparisons in which the given gene/miRNA is found to be DE.
4. Download:
    1. `Download edges` (Interaction table to be imported into Cytoscape `.txt`)
        1. `File > Import > Network from File...`
        2. Source Node: `mirna`
        3. Target Node: `gene`
        4. Edge Attribute: `rho,database`
    2. `Download nodes` (Nodes table to be imported into Cytoscape `.txt`)
        1. `File > Import > Table from File...`
        2. Key: `node`
        3. Attribute: `comparison,node_type`
    3. `Download gene/miRNA list` (List of all miRNAs and genes from correlation table to use as `File Input` for `Data Browser` .txt)
        1. See `Data Browser Input` example
