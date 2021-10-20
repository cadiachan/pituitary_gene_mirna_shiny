# Pituitary Gland Shiny App
Shiny app to visualize DE genes and miRNAs in the postnatal mouse pituitary gland.

## To run
1. Clone git repository to local directory.
2. Install dependencies in R with `devtools::install()` (Set Working Directory to path containing `DESCRIPTION`).
3. Install EDASeq using Bioconductor:  
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")    
BiocManager::install("EDASeq")
```
4. Open `app.R` in RStudio. Click "Run App" to run locally.

## Example
1. Navigate to the `Data Browser` tab. 
2. Upload `example_input_file.txt` (located in `data/` directory) or enter in a comma-separated list of genes and/or miRNAs.
3. Press `Submit` to run.
4. Press `Reset` to clear inputs.
