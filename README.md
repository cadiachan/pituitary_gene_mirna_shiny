# Pituitary Gland Shiny App
Shiny app to visualize DE genes and miRNAs in the postnatal mouse pituitary gland.

## To run:
1. Clone git repository to local directory.
2. Install dependencies in RStudio with `devtools::install()` (set Working Directory to where respository is cloned `DESCRIPTION` is).
3. Install EDASeq through using Bioconductor tools:  
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")    
BiocManager::install("EDASeq")
```
5. Open `app.R` in RStudio. Click "Run App" to run locally.

## Example
1. Navigate to the `Data Browser` tab. 
2. Upload `example_input_file.txt` (located in `data/` directory containing Shiny App files).
3. Press `Submit` to run.
