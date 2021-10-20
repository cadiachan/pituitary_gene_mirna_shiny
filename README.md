# Pituitary Gland Shiny App
Shiny app to visualize DE genes and miRNAs in the postnatal mouse pituitary gland.

## To run:
1. Clone git repository to local directory.
2. Run `devtools::install()` in your project R console (path where `DESCRIPTION` is) to install dependencies.
3. Install Bioconductor tools:  
`if (!requireNamespace("BiocManager", quietly = TRUE))`  
    `install.packages("BiocManager")`      
`BiocManager::install("EDASeq")`  
`BiocManager::install("edgeR")`  
`BiocManager::install("biomaRt")`  

5. Open `app.R` in RStudio. Click "Run App" to run locally.

## Example
1. Navigate to the `Data Browser` tab. 
2. Upload `example_input_file.txt` (located in `data/` directory containing Shiny App files).
3. Press `Submit` to run.
