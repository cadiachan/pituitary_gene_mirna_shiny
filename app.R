library(shiny)

source("scripts.R")
# Define UI ----
ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  # shinyUI(
  navbarPage("Mouse Pituitary Gland Data Hub",
             #### Home ####
             tabPanel("Home",icon = icon("home"),
                      h1("MOUSE PITUITARY GLAND DATA HUB"),
                      br(),
                      p("Shiny app in development to visualize gene and miRNA expression data generated",
                        "from 3'UTR-seq and small RNA-seq respectively."),
                      p("Raw data can be accessed at ",
                        a("E-MTAB-9460.", 
                          href = "https://www.ebi.ac.uk/arrayexpress/",
                          target = "_blank"),
                        "(Fix link once data is available)"),
                      br(),
                      h2("Features"),
                      h3("Data Browser"),
                      p("- Visualize gene and miRNA expression plots across postnatal ages and between sexes."),
                      p("- Quantify log2FC and FDR for DE genes and miRNAs for each comparison."),
                      p("- Intersect genes of interest with puberty- and pituitary disease-related gene list compendium curated from relevant published studies."),
                      p("- Output normalized and log2(normalized) counts for genes and miRNAs of interest."),
                      br(),
                      h3("miRNA-Gene Target Browser"),
                      p("- "),
                      br(),
                      br(),
                      br(),
                      h2("To-do <priority>"),
                      p("- Add in miRNA-gene correlation tab"),
                      p("- Add in default lists to show DE genes and miRNAs from comparisons in the manuscript"),
                      p("- Design graphical abstract for <Home> tab"),
                      br(),
                      h2("To-do <possible things to add>"),
                      p("- Add cell-type specificity to DE genes based on scMappR results?"),
                      p("- qPCR tab from Hou et al 2017"),
                      p("- Add filter for sex-bias or age-bias in DE tables (redundant feature of default lists?)")
             ),
             #### Data Browser ####
             tabPanel("Data Browser", icon = icon("chart-bar"),
                      sidebarLayout(
                        sidebarPanel(
                          shinyjs::useShinyjs(),
                          id = "side-panel",
                          width = 3,
                          h3("Select input method"),
                          br(),
                          
                          radioButtons("input_type",
                                       label = NULL,
                                       choices = list("Type in genes/miRNAs" = 1, "Upload file" = 2),
                                       selected = 1,
                                       inline = T),
                          uiOutput("add_helper_input"),
                          uiOutput("add_input_ui"),
                          div(style="display:inline-block",
                              actionButton("submit",
                                           label = "Submit",
                                           class = "btn-success")),
                          div(style="display:inline-block",
                              actionButton("reset",
                                           label = "Reset")),
                          br(),
                          br(),
                          # h3("Select studies to intersect"),
                          # checkboxGroupInput("pub_study",
                          #                    label = h4("Puberty-related gene lists"),
                          #                    choiceNames = list(
                          #                      HTML("Perry 2014 <a href = 'https://pubmed.ncbi.nlm.nih.gov/25231870/'>(PMID: 25231870)</a>"),
                          #                      HTML("Day 2015 <a href = 'https://pubmed.ncbi.nlm.nih.gov/26548314/'>(PMID: 26548314)</a>"),
                          #                      HTML("Day 2017 <a href = 'https://pubmed.ncbi.nlm.nih.gov/28436984/'>(PMID: 28436984)</a>"),
                          #                      HTML("Hollis 2020 <a href = 'https://pubmed.ncbi.nlm.nih.gov/32210231/'>(PMID: 32210231)</a>"),
                          #                      "IHH/Kallmann"
                          #                    ),
                          #                    choiceValues = list("Perry2014",
                          #                                        "Day2015_GWAS_VB",
                          #                                        "Day2017_nearest",
                          #                                        "Hollis2020_GWAS_VB_FH",
                          #                                        "IHH/Kallmann"
                          #                    ),
                          #                    selected = c("Perry2014",
                          #                                 "Day2015_GWAS_VB",
                          #                                 "Day2017_nearest",
                          #                                 "Hollis2020_GWAS_VB_FH",
                          #                                 "IHH/Kallmann")),
                          # checkboxGroupInput("pit_study",
                          #                    label = h4("Pituitary-related gene lists"),
                          #                    choiceNames = list(
                          #                      HTML("Ye 2015 <a href = 'https://pubmed.ncbi.nlm.nih.gov/26029870/'>(PMID: 26029870)</a>"),
                          #                      HTML("Fang 2016 <a href = 'https://pubmed.ncbi.nlm.nih.gov/27828722/'>(PMID: 27828722)</a>"),
                          #                      HTML("Hauser 2019 <a href = 'https://pubmed.ncbi.nlm.nih.gov/31139150/'>(PMID: 31139150)</a>"),
                          #                      HTML("Kurtoglu 2019 <a href = 'https://pubmed.ncbi.nlm.nih.gov/29739730/'>(PMID: 29739730)</a>")
                          #                    ),
                          #                    choiceValues = list("Ye2015_PA_GWAS",
                          #                                        "Fang2016_CPHD",
                          #                                        "Hauser2019_PA",
                          #                                        "Kurtoglu2019_hypopituitarism"
                          #                    ),
                          #                    selected = c("Ye2015_PA_GWAS",
                          #                                 "Fang2016_CPHD",
                          #                                 "Hauser2019_PA",
                          #                                 "Kurtoglu2019_hypopituitarism")),
                          h4("Text input examples:"),
                          h4("Genes: Lhb,ENSMUSG00000027120.7"),
                          h4("miRNAs: mmu-miR-224-5p,miR-383-5p"),
                          br(),
                          helpText(h5("It is recommended to input <20 genes/miRNAs for viewing on the browser.")),
                          helpText(h5("If more genes/miRNAs are inputted, expression values and DE table can be downloaded."))
                        ),
                        mainPanel(
                          shinyjs::useShinyjs(),
                          width = 9,
                          tabsetPanel(
                            type = "tabs",
                            tabPanel(
                              "Plot",
                              br(),
                              verbatimTextOutput("input_err"),
                              verbatimTextOutput("invalid_genes"),
                              
                              column(6,
                                     br(),
                                     tags$div(id = "save_gene_plot_text", h4("Save gene expression plots: ")),
                                     downloadButton("save_png_plot_genes", ".png"),
                                     downloadButton("save_pdf_plot_genes", ".pdf"),
                                     h3("Genes"),
                                     plotOutput("gene_plot",
                                                height = "auto",
                                                width = "70%")) ,
                              column(6,
                                     br(),
                                     tags$div(id = "save_mirna_plot_text", h4("Save miRNA expression plots: ")),
                                     downloadButton("save_png_plot_mirnas", ".png"),
                                     downloadButton("save_pdf_plot_mirnas", ".pdf"),
                                     h3("miRNAs"),
                                     plotOutput("mirna_plot",
                                                height = "auto",
                                                width = "70%"))
                            ),
                            tabPanel("DE gene table",
                                     br(),
                                     radioButtons("de_gene_filt",
                                                  label = h3("Cutoff for DE gene table"),
                                                  choices = list("abs(FC) > 1.5, FDR < 0.05" = 1, "No cutoff" = 2),
                                                  selected = 1,
                                                  inline = T),
                                     tags$div(id = "save_de_gene_text", h4("Save DE gene table: ")),
                                     downloadButton("save_txt_de_genes", ".txt"),
                                     downloadButton("save_csv_de_genes", ".csv"),
                                     tableOutput("de_gene_table"),
                            ),
                            tabPanel("DE miRNA table",
                                     br(),
                                     radioButtons("de_mirna_filt",
                                                  label = h3("Cutoff for DE miRNA table"),
                                                  choices = list("abs(FC) > 1.5, FDR < 0.05" = 1, "No cutoff" = 2),
                                                  selected = 1,
                                                  inline = T),
                                     tags$div(id = "save_de_mirna_text", h4("Save DE miRNA table: ")),
                                     downloadButton("save_txt_de_mirnas", ".txt"),
                                     downloadButton("save_csv_de_mirnas", ".csv"),
                                     tableOutput("de_mirna_table")),
                            tabPanel("Gene expression values",
                                     br(),
                                     radioButtons("count_type",
                                                  label = h3("Count type"),
                                                  choices = list("Normalized counts" = 1, "log2(normCounts)" = 2),
                                                  selected = 1,
                                                  inline = T),
                                     br(),
                                     tags$div(id = "save_gene_expr_text", h4("Save gene expression values: ")),
                                     downloadButton("save_txt_expr_genes", ".txt"),
                                     downloadButton("save_csv_expr_genes", ".csv"),
                                     tableOutput("gene_table"),
                            ),
                            tabPanel("miRNA expression values",
                                     br(),
                                     radioButtons("count_mirna_type",
                                                  label = h3("Count type"),
                                                  choices = list("Normalized counts" = 1, "log2(normCounts)" = 2),
                                                  selected = 1,
                                                  inline = T),
                                     br(),
                                     tags$div(id = "save_mirna_expr_text", h4("Save miRNA expression values: ")),
                                     downloadButton("save_txt_expr_mirnas", ".txt"),
                                     downloadButton("save_csv_expr_mirnas", ".csv"),
                                     tableOutput("mirna_table")
                            )
                          )
                        )) 
             ),
             #### miRNA-gene ####
             tabPanel(" miRNA-Gene Target Browser", icon = icon("project-diagram"),
                      sidebarLayout(
                        sidebarPanel(
                          shinyjs::useShinyjs(),
                          id = "side-panel_2",
                          width = 3,
                          h3("Input a list of genes OR miRNAs"),
                          # br(),
                          # h4("Gene examples: 'Lhb', 'ENSMUSG00000027120.7'"),
                          # h4("miRNA examples: 'mmu-miR-224-5p', 'miR-383-5p'"),
                          radioButtons("input_type_2",
                                       label = NULL,
                                       choices = list("miRNA" = 1, "gene" = 2),
                                       selected = 1,
                                       inline = T),
                          uiOutput("add_helper_input_2"),
                          uiOutput("add_input_ui_2"),
                          helpText(h5("Choose to show only differentially expressed (DE) miRNA-gene pairs or all miRNA-gene pairs.")),
                          radioButtons("filt_choice_2",
                                       label = NULL,
                                       choices = list("DE" = 1, "All" = 2),
                                       selected = 1,
                                       inline = T),
                          div(style="display:inline-block",
                              actionButton("submit_2",
                                           label = "Submit",
                                           class = "btn-success")),
                          div(style="display:inline-block",
                              actionButton("reset_2",
                                           label = "Reset")),
                          br(),
                          helpText(h5("Save and upload tables into Cytoscape to visualize as a network."))
                        ),
                        mainPanel(
                          width = 9,
                          tabsetPanel(
                            type = "tabs",
                            tabPanel(
                              "miRNA-target genes",
                              br(),
                              verbatimTextOutput("input_err_2"),
                              verbatimTextOutput("invalid_genes_2"),
                              div(style="display:inline-block",
                                  downloadButton("save_txt_corr_cyto", "Download table in Cytoscape input format")),
                              div(style="display:inline-block",
                                  downloadButton("save_txt_corr_input", "Download gene/miRNA list for input to Data Browser")),
                              br(),
                              tableOutput("corr_table")
                            )
                          )
                        )) 
             ),
             #### About ####
             tabPanel("About", icon=icon("quote-left"),
                      # h2("Credits"),
                      # p("Text about this app"),
                      # br(),
                      h2("Cite"),
                      h4("Hou H, Chan C, Yuki KE, et al. Postnatal developmental trajectory of sex-biased gene expression in the mouse pituitary gland.",
                         "Manuscript in preparation."),
                      br(),
                      h2("References"),
                      downloadButton("save_genelists", "Download gene list compendium"),
                      p("A compendium of puberty-related (#1-4,9) and pituitary gland disease-related (#5-8) gene lists were curated from the following studies:"),
                      tags$ol(
                        tags$li(h4("Perry JR, Day F, Elks CE, et al. Parent-of-origin-specific allelic associations among 106 genomic loci for age at menarche.",
                                   "Nature. 2014;514(7520):92-97. PMID: 25231870.",
                                   a("doi:10.1038/nature13545",
                                     href="https://pubmed.ncbi.nlm.nih.gov/25231870/",
                                     target = "_blank"))),
                        tags$li(h4("Day FR, Bulik-Sullivan B, Hinds DA, et al. Shared genetic aetiology of puberty timing between sexes and with",
                                   "health-related outcomes. Nat Commun. 2015;6:8842. Published 2015 Nov 9. PMID: 26548314.",
                                   a("doi:10.1038/ncomms9842",
                                     href="https://pubmed.ncbi.nlm.nih.gov/26548314/",
                                     target = "_blank"))),
                        tags$li(h4("Day FR, Thompson DJ, Helgason H, et al. Genomic analyses identify hundreds of variants associated with age at menarche",
                                   "and support a role for puberty timing in cancer risk. Nat Genet. 2017;49(6):834-841. PMID: 28436984.",
                                   a("doi:10.1038/ng.3841",
                                     href="https://pubmed.ncbi.nlm.nih.gov/28436984/",
                                     target = "_blank"))),
                        tags$li(h4("Hollis B, Day FR, Busch AS, et al. Genomic analysis of male puberty timing highlights shared genetic basis with hair colour",
                                   "and lifespan. Nat Commun. 2020;11(1):1536. Published 2020 Mar 24. PMID: 32210231.",
                                   a("doi:10.1038/s41467-020-14451-5",
                                     href="https://pubmed.ncbi.nlm.nih.gov/32210231/",
                                     target = "_blank"))),
                        tags$li(h4("Ye Z, Li Z, Wang Y, et al. Common variants at 10p12.31, 10q21.1 and 13q12.13 are associated with sporadic pituitary adenoma.",
                                   "Nat Genet. 2015;47(7):793-797. PMID: 26029870.",
                                   a("doi:10.1038/ng.3322",
                                     href="https://pubmed.ncbi.nlm.nih.gov/26029870/",
                                     target = "_blank"))),
                        tags$li(h4("Fang Q, George AS, Brinkmeier ML, et al. Genetics of Combined Pituitary Hormone Deficiency: Roadmap into the Genome Era.",
                                   "Endocr Rev. 2016;37(6):636-675. PMID: 27828722.",
                                   a("doi:10.1210/er.2016-1101",
                                     href="https://pubmed.ncbi.nlm.nih.gov/27828722/",
                                     target = "_blank"))),
                        tags$li( h4("Hauser BM, Lau A, Gupta S, Bi WL, Dunn IF. The Epigenomics of Pituitary Adenoma.",
                                    "Front Endocrinol (Lausanne). 2019;10:290. Published 2019 May 14. PMID: 31139150.",
                                    a("doi:10.3389/fendo.2019.00290",
                                      href="https://pubmed.ncbi.nlm.nih.gov/31139150/",
                                      target = "_blank"))),
                        tags$li(h4("Kurtoglu S, Ozdemir A, Hatipoglu N. Neonatal Hypopituitarism: Approaches to Diagnosis and Treatment.",
                                   "J Clin Res Pediatr Endocrinol. 2019;11(1):4-12. PMID: 29739730.",
                                   a("doi:10.4274/jcrpe.galenos.2018.2018.0036",
                                     href="https://pubmed.ncbi.nlm.nih.gov/32210231/",
                                     target = "_blank"))),
                        tags$li(h4("Isolated hypogonadotropic hypogonadism (IHH)/Kallmann syndrome.",
                                   "Gene list previously curated in Hou H, Uusküla-Reimand L, Makarem M, et al.",
                                   "Gene expression profiling of puberty-associated genes reveals abundant tissue and sex-specific changes across postnatal development.",
                                   "Hum Mol Genet. 2017;26(18):3585-3599. PMID: 28911201.",
                                   a("doi:10.1093/hmg/ddx246",
                                     href="https://pubmed.ncbi.nlm.nih.gov/28911201/",
                                     target = "_blank")))
                      ),
                      br(),
                      h2("Funding")
             )
  )
)



# Define server logic ----
server <- function(input, output) {
  #### Functions for miRNA-gene tab ####
  data_2 <- reactive({
    mirna_input <- F
    gene_input <- F
    # print(input$gene)
    if(!is.null(input$gene_2)) {
      # Check if genes are inputted by text box
      if(input$input_type_2 == 1 & nchar(input$gene_2) > 0) {
        mirna_input <- T
      }
      # Check if a file is uploaded
      else if(input$input_type_2 == 2) {
        gene_input <- T
      }else { mirna_input <- F; gene_input <- F}
    }
    
    print(mirna_input)
    print(gene_input)
    
    if(mirna_input == T & gene_input == F) {
      output$input_err_2 <- NULL
      data_2 <- parse_list(input$gene_2, type = "text")
      if(length(data_2[["genes"]]) > 0 & length(data_2[["mirnas"]]) > 0) {
        msg_2 <- paste0("Not a miRNA: ", paste0(data_2[["genes"]], collapse = ","))
        output$input_err_2 <- renderText({
          msg_2
        })
      }
      if(length(data_2[["mirnas"]]) == 0) {
        data_2 <- parse_list("", type = "text")
        msg_2 <- "Please input a valid list of miRNAs."
        output$input_err_2 <- renderText({
          msg_2
        })
      }
    }
    else if(mirna_input == F & gene_input == T) {
      data_2 <- parse_list(input$gene_2, type = "text")
      if(length(data_2[["mirnas"]]) > 0 & length(data_2[["genes"]]) > 0) {
        msg_2 <- paste0("Not a gene: ", paste0(data_2[["mirnas"]], collapse = ","))
        output$input_err_2 <- renderText({
          msg_2
        })
      }
      if(length(data_2[["genes"]]) == 0) {
        data_2 <- parse_list("", type = "text")
        msg_2 <- "Please input a valid list of genes."
        output$input_err_2 <- renderText({
          msg_2
        })
      }
    }
    else{ # Print error message if both fields are empty.
      data_2 <- parse_list("", type = "text")
      msg_2 <- "Please input a list of genes/miRNAs."
      output$input_err_2 <- renderText({
        msg_2
      })
    }
    return(data_2)
  })
  
  # Reset input values in the side-panel
  observeEvent(input$reset_2, {
    shinyjs::reset("side-panel_2")
  })
  
  # Add in UI based on input choice
  observeEvent(input$input_type_2, {
    if(input$input_type_2 == 1) {
      output$add_helper_input_2 <- renderUI({
        helpText(h5("Separate miRNAs with a comma"))
      })
      output$add_input_ui_2 <- renderUI({
        textInput("gene_2", label = NULL,
                  placeholder = "let-7a-5p,mmu-let-7e-5p,miR-224-5p"
        )
      })
    }
    if(input$input_type_2 == 2) {
      output$add_helper_input_2 <- renderUI({
        helpText(h5("Separate genes with a comma"))
      })
      output$add_input_ui_2 <- renderUI({
        textInput("gene_2", label = NULL,
                  placeholder = "Lhb,ENSMUSG00000027120.7,Gh,Prl"
        )
      })
    }
  })
  
  shinyjs::hide(id = "save_txt_corr_cyto")
  shinyjs::hide(id = "save_txt_corr_input")
  # Run functions in response to submit button
  # eventReactive events are delayed until the button is pressed
  press_submit_2 <- eventReactive(input$submit_2, {
    shinyjs::hide(id = "save_txt_corr_cyto")
    shinyjs::hide(id = "save_txt_corr_input")
    
    output$invalid_genes_2 <- NULL
    invalid_msg <- "Not found: "
    invalid_genes <- paste(data_2()[["invalid"]], collapse = ",")
    if(identical(data_2()[["invalid"]], character(0))) {
      output$invalid_genes_2 <- NULL
    }
    else if(length(data_2()[["invalid"]]) > 1) {
      output$invalid_genes_2 <- renderText({
        paste0(invalid_msg, invalid_genes)
      })
    }
    else if(length(data_2()[["invalid"]]) == 1 & data_2()[["invalid"]] != "") {
      output$invalid_genes_2 <- renderText({
        paste0(invalid_msg, invalid_genes)
      })
    } else{output$invalid_genes_2 <- NULL}
    # print(data_2())
    if(input$input_type_2 == 1) {
      if(input$filt_choice_2 == 1){
        corrtable <- print_corr_table(mirnalist = data_2()[["mirnas"]])
      }
      else {
        corrtable <- print_corr_table(mirnalist = data_2()[["mirnas"]],
                                      de_filt = F)
      }
    }
    
    if(input$input_type_2 == 2) {
      if(input$filt_choice_2 == 1){
        corrtable <- print_corr_table(genelist = data_2()[["genes"]])
      }
      else {
        corrtable <- print_corr_table(genelist = data_2()[["genes"]],
                                      de_filt = F)
      }
    }
    
    if(ncol(corrtable) > 1) { # Check that corr table is not just an error table with 1 column
      shinyjs::show(id = "save_txt_corr_cyto")
      shinyjs::show(id = "save_txt_corr_input")
    }
    
    output$save_txt_corr_cyto <- downloadHandler(
      filename <- function() {
        paste0(Sys.Date(), "_corr_table_cytoscape.txt")
      },
      content <- function(file) {
        write.table(corrtable, file,
                    quote = F, sep = "\t", col.names = T, row.names = F)
      }
    )
    
    output$save_txt_corr_input <- downloadHandler(
      filename <- function() {
        paste0(Sys.Date(), "_corr_gene_mirna_list.txt")
      },
      content <- function(file) {
        write.table(corrtable, file,
                    quote = F, sep = "\t", col.names = T, row.names = F)
      }
    )
    return(list("corrtable" = corrtable))
  })
  
  # Outputs correlation table
  output$corr_table <- renderTable({
    press_submit_2()[["corrtable"]]
  },
  striped = T,
  hover = T,
  digits = 5)
  

  
  #### Functions for data browser tab ####
  data <- reactive({
    text_input <- F
    file_input <- F
    # print(input$gene)
    if(!is.null(input$gene)) {
      # Check if genes are inputted by text box
      if(input$input_type == 1 & nchar(input$gene) > 0) {
        text_input <- T
      } 
      # Check if a file is uploaded
      else if(input$input_type == 2) {
        file_input <- T
      } else { file_input <- F; text_input <- F}
    }
    
    if(text_input == T & file_input == F) {
      output$input_err <- NULL
      data <- parse_list(input$gene, type = "text")
    }
    else if(text_input == F & file_input == T) {
      output$input_err <- NULL
      file <- input$gene
      read_file <- read.table(file$datapath, header = F)
      data <- parse_list(read_file, type = "file")
    }
    else{ # Print error message if both fields are empty.
      data <- parse_list("", type = "text")
      msg <- "Please input a list of genes/miRNAs."
      output$input_err <- renderText({
        msg
      })
    }
    return(data)
  })
  
  # Reset input values in the side-panel
  observeEvent(input$reset, {
    shinyjs::reset("side-panel")
  })
  
  # Add in UI based on input choice
  observeEvent(input$input_type, {
    if(input$input_type == 1) {
      output$add_helper_input <- renderUI({
        helpText(h5("Separate genes/miRNAs with a comma"))
      })
      output$add_input_ui <- renderUI({
        textInput("gene", label = NULL,
                  placeholder = "Lhb,let-7a-5p,mmu-let-7e-5p,ENSMUSG00000027120.7"
                  # value ="Lhb,let-7a-5p,mmu-let-7e-5p,ENSMUSG00000027120.7,miR-224-5p"
        )
      })
    }
    if(input$input_type == 2) {
      output$add_helper_input <- renderUI({
        helpText(h5("Upload a .txt file with genes/miRNAs on separate lines"))
      })
      output$add_input_ui <- renderUI({
        fileInput("gene", label = NULL, accept = ".txt")
      })
    }
  })
  
  #### Hide buttons during first initialization ####
  shinyjs::hide(id = "save_gene_plot_text")
  shinyjs::hide(id = "save_png_plot_genes")
  shinyjs::hide(id = "save_pdf_plot_genes")
  shinyjs::hide(id = "save_de_gene_text")
  shinyjs::hide(id = "save_txt_de_genes")
  shinyjs::hide(id = "save_csv_de_genes")
  shinyjs::hide(id = "save_gene_expr_text")
  shinyjs::hide(id = "save_txt_expr_genes")
  shinyjs::hide(id = "save_csv_expr_genes")
  
  shinyjs::hide(id = "save_mirna_plot_text")
  shinyjs::hide(id = "save_png_plot_mirnas")
  shinyjs::hide(id = "save_pdf_plot_mirnas")
  shinyjs::hide(id = "save_de_mirna_text")
  shinyjs::hide(id = "save_txt_de_mirnas")
  shinyjs::hide(id = "save_csv_de_mirnas")
  shinyjs::hide(id = "save_mirna_expr_text")
  shinyjs::hide(id = "save_txt_expr_mirnas")
  shinyjs::hide(id = "save_csv_expr_mirnas")
  
  # Run functions in response to submit button
  # eventReactive events are delayed until the button is pressed
  press_submit <- eventReactive(input$submit, {
    shinyjs::hide(id = "save_gene_plot_text")
    shinyjs::hide(id = "save_png_plot_genes")
    shinyjs::hide(id = "save_pdf_plot_genes")
    shinyjs::hide(id = "save_de_gene_text")
    shinyjs::hide(id = "save_txt_de_genes")
    shinyjs::hide(id = "save_csv_de_genes")
    shinyjs::hide(id = "save_gene_expr_text")
    shinyjs::hide(id = "save_txt_expr_genes")
    shinyjs::hide(id = "save_csv_expr_genes")
    
    shinyjs::hide(id = "save_mirna_plot_text")
    shinyjs::hide(id = "save_png_plot_mirnas")
    shinyjs::hide(id = "save_pdf_plot_mirnas")
    shinyjs::hide(id = "save_de_mirna_text")
    shinyjs::hide(id = "save_txt_de_mirnas")
    shinyjs::hide(id = "save_csv_de_mirnas")
    shinyjs::hide(id = "save_mirna_expr_text")
    shinyjs::hide(id = "save_txt_expr_mirnas")
    shinyjs::hide(id = "save_csv_expr_mirnas")
    
    output$invalid_genes <- NULL
    invalid_msg <- "Not found: "
    invalid_genes <- paste(data()[["invalid"]], collapse = ",")
    if(identical(data()[["invalid"]], character(0))) {
      output$invalid_genes <- NULL
    }
    else if(length(data()[["invalid"]]) > 1) {
      output$invalid_genes <- renderText({
        paste0(invalid_msg, invalid_genes)
      })
    }
    else if(length(data()[["invalid"]]) == 1 & data()[["invalid"]] != "") {
      output$invalid_genes <- renderText({
        paste0(invalid_msg, invalid_genes)
      })
    } else{output$invalid_genes <- NULL}
    
    gplot <- exprplot_hhtheme(genelist = data()[["genes"]],
                              count_data = utr_log2,
                              metadata = utr_meta,
                              counttype = "genes")
    num_genes <- length(gplot$plot_env$use_genelist)
    
    if(num_genes > 0) {
      shinyjs::show(id = "save_gene_plot_text")
      shinyjs::show(id = "save_png_plot_genes")
      shinyjs::show(id = "save_pdf_plot_genes")
      shinyjs::show(id = "save_gene_expr_text")
      shinyjs::show(id = "save_txt_expr_genes")
      shinyjs::show(id = "save_csv_expr_genes")
    }
    
    mplot <- exprplot_hhtheme(genelist = data()[["mirnas"]],
                              count_data = mirna_log2,
                              metadata = mirna_meta,
                              counttype = "mirnas")
    num_mirnas <- length(mplot$plot_env$use_genelist)
    
    if(num_mirnas > 0) {
      shinyjs::show(id = "save_mirna_plot_text")
      shinyjs::show(id = "save_png_plot_mirnas")
      shinyjs::show(id = "save_pdf_plot_mirnas")
      shinyjs::show(id = "save_mirna_expr_text")
      shinyjs::show(id = "save_txt_expr_mirnas")
      shinyjs::show(id = "save_csv_expr_mirnas")
    }
    
    gtable <- print_de_table(genelist = data()[["genes"]],
                             de_table = utr_de_table,
                             counttype = "genes")
    
    gtable_nofilt <- print_de_table(genelist = data()[["genes"]],
                                    de_table = utr_de_all_table,
                                    counttype = "genes")
    gtable <- add_study(gtable)
    gtable_nofilt <- add_study(gtable_nofilt)
    
    mtable <- print_de_table(genelist = data()[["mirnas"]],
                             de_table = mirna_de_table,
                             counttype = "mirnas")
    mtable_nofilt <- print_de_table(genelist = data()[["mirnas"]],
                                    de_table = mirna_de_all_table,
                                    counttype = "mirnas")
    
    norm_data_type <- expr_table(data()[["genes"]], utr_normcounts)
    log2_data_type <- expr_table(data()[["genes"]], utr_log2)
    
    norm_mirna_data_type <- expr_table(data()[["mirnas"]], mirna_normcounts)
    log2_mirna_data_type <- expr_table(data()[["mirnas"]], mirna_log2)
    return(list("gplot" = gplot,
                "num_genes" = num_genes,
                "mplot" = mplot,
                "num_mirnas" = num_mirnas,
                "gtable" = gtable,
                "gtable_nofilt" = gtable_nofilt,
                "mtable" = mtable,
                "mtable_nofilt" = mtable_nofilt,
                "norm_data_type" = norm_data_type,
                "log2_data_type" = log2_data_type,
                "norm_mirna_data_type" = norm_mirna_data_type,
                "log2_mirna_data_type" = log2_mirna_data_type))
  })
  
  # observed events occur the moment the button is pressed
  # Switches DE table between DE GENE table with cutoff and DE GENE table without cutoff
  # Also switches save buttons
  observeEvent(input$submit, {
    observeEvent(input$de_gene_filt, {
      gene_de <- press_submit()
      if(input$de_gene_filt == 1) {
        gene_de_table <- gene_de[["gtable"]]
        if(ncol(gene_de_table) > 1) {
          shinyjs::show(id = "save_de_gene_text")
          shinyjs::show(id = "save_txt_de_genes")
          shinyjs::show(id = "save_csv_de_genes")
        }
        # Save as txt file with DE filt
        output$save_txt_de_genes <- downloadHandler(
          filename <- function() {
            paste0(Sys.Date(), "_de_genes.txt")
          },
          content <- function(file) {
            write.table(gene_de_table, file,
                        quote = F, sep = "\t", col.names = T, row.names = F)
          }
        )
        # Save as csv file with DE filt
        output$save_csv_de_genes <- downloadHandler(
          filename <- function() {
            paste0(Sys.Date(), "_de_genes.csv")
          },
          content <- function(file) {
            write.csv(gene_de_table, file,
                      quote = F, row.names = F)
          }
        )
      }
      if(input$de_gene_filt == 2) {
        gene_de_table <- gene_de[["gtable_nofilt"]]
        if(ncol(gene_de_table) > 1) {
          shinyjs::show(id = "save_de_gene_text")
          shinyjs::show(id = "save_txt_de_genes")
          shinyjs::show(id = "save_csv_de_genes")
        }
        # Save as txt file with no DE filt
        output$save_txt_de_genes <- downloadHandler(
          filename <- function() {
            paste0(Sys.Date(), "_de_nofilt_genes.txt")
          },
          content <- function(file) {
            write.table(gene_de_table, file,
                        quote = F, sep = "\t", col.names = T, row.names = F)
          }
        )
        # Save as csv file with no DE filt
        output$save_csv_de_genes <- downloadHandler(
          filename <- function() {
            paste0(Sys.Date(), "_de_nofilt_genes.csv")
          },
          content <- function(file) {
            write.csv(gene_de_table, file,
                      quote = F, row.names = F)
          }
        )
      }
      output$de_gene_table <- renderTable({
        gene_de_table
      },
      striped = T,
      hover = T,
      digits = 5)
    })
  })
  
  # observed events occur the moment the button is pressed
  # Switches DE table between DE miRNA table with cutoff and DE miRNA table without cutoff
  # Also switches save buttons
  observeEvent(input$submit, {
    observeEvent(input$de_mirna_filt, {
      mirna_de <- press_submit()
      if(input$de_mirna_filt == 1) {
        mirna_de_table <- mirna_de[["mtable"]]
        if(ncol(mirna_de_table) > 1) {
          shinyjs::show(id = "save_de_mirna_text")
          shinyjs::show(id = "save_txt_de_mirnas")
          shinyjs::show(id = "save_csv_de_mirnas")
        }
        # Save as text file with DE filt
        output$save_txt_de_mirnas <- downloadHandler(
          filename <- function() {
            paste0(Sys.Date(), "_de_mirnas.txt")
          },
          content <- function(file) { 
            write.table(mirna_de_table, file,
                        quote = F, sep = "\t", col.names = T, row.names = F)
          }
        )
        # Save as csv file with DE filt
        output$save_csv_de_mirnas <- downloadHandler(
          filename <- function() {
            paste0(Sys.Date(), "_de_mirnas.csv")
          },
          content <- function(file) {
            write.csv(mirna_de_table[["mtable"]], file,
                      quote = F, row.names = F)
          }
        )
      }
      if(input$de_mirna_filt == 2) {
        mirna_de_table <- mirna_de[["mtable_nofilt"]]
        if(ncol(mirna_de_table) > 1) {
          shinyjs::show(id = "save_de_mirna_text")
          shinyjs::show(id = "save_txt_de_mirnas")
          shinyjs::show(id = "save_csv_de_mirnas")
        }
        # Save as txt file with no DE filt
        output$save_txt_de_mirnas <- downloadHandler(
          filename <- function() {
            paste0(Sys.Date(), "_de_nofilt_mirnas.txt")
          },
          content <- function(file) { 
            write.table(mirna_de_table, file,
                        quote = F, sep = "\t", col.names = T, row.names = F)
          }
        )
        # Save as csv file with no DE filt
        output$save_csv_de_mirnas <- downloadHandler(
          filename <- function() {
            paste0(Sys.Date(), "_de_nofilt_mirnas.csv")
          },
          content <- function(file) {
            write.csv(mirna_de_table, file,
                      quote = F, row.names = F)
          }
        )
      }
      output$de_mirna_table <- renderTable({
        mirna_de_table
      },
      striped = T,
      hover = T,
      digits = 5)
    })
  })
  
  # observed events occur the moment the button is pressed
  # Switches output table between normalized and log2norm gene counts
  observeEvent(input$submit, {
    observeEvent(input$count_type, {
      gene_expr <- press_submit()
      if(input$count_type == 1) {
        data_type <- gene_expr[["norm_data_type"]]
        err_gene_tab <- gene_expr[["norm_err_gene_tab"]]
        df_out_type <- (as.data.frame(cbind(sample=rownames(data_type), data_type)))
        output$save_txt_expr_genes <- downloadHandler(
          filename <- function() {
            paste0(Sys.Date(), "_norm_expression_genes.txt")
          },
          content <- function(file) {
            write.table(df_out_type, file,
                        quote = F, sep = "\t", col.names = T, row.names = F)
          }
        )
        output$save_csv_expr_genes <- downloadHandler(
          filename <- function() {
            paste0(Sys.Date(), "_norm_expression_genes.csv")
          },
          content <- function(file) {
            write.csv(df_out_type, file, quote = F, row.names = F)
          }
        )
      }
      if(input$count_type == 2) {
        data_type <- gene_expr[["log2_data_type"]]
        err_gene_tab <- gene_expr[["log2_err_gene_tab"]]
        df_out_type <- (as.data.frame(cbind(sample=rownames(data_type), data_type)))
        output$save_txt_expr_genes <- downloadHandler(
          filename <- function() {
            paste0(Sys.Date(), "_log2_expression_genes.txt")
          },
          content <- function(file) {
            write.table(df_out_type, file,
                        quote = F, sep = "\t", col.names = T, row.names = F)
          }
        )
        output$save_csv_expr_genes <- downloadHandler(
          filename <- function() {
            paste0(Sys.Date(), "_norm_expression_genes.csv")
          },
          content <- function(file) {
            write.csv(df_out_type, file, quote = F, row.names = F)
          }
        )
      }
      output$gene_table <- renderTable({
        data_type
      },
      striped = T,
      hover = T,
      rownames = T)
    })
  })
  
  
  # Switches output table between normalized and log2norm miRNA counts
  observeEvent(input$submit, {
    observeEvent(input$count_mirna_type, {
      mirna_expr <- press_submit()
      if(input$count_mirna_type == 1) {
        mirna_data_type <- mirna_expr[["norm_mirna_data_type"]]
        mirna_err_gene_tab <- mirna_expr[["norm_err_mirna_tab"]]
        df_out_type <- (as.data.frame(cbind(sample=rownames(mirna_data_type), mirna_data_type)))
        
        output$save_txt_expr_mirnas <- downloadHandler(
          filename <- function() {
            paste0(Sys.Date(), "_norm_expression_mirnas.txt")
          },
          content <- function(file) {
            write.table(df_out_type, file, quote = F, sep = "\t", col.names = T, row.names = F)
          }
        )
        output$save_csv_expr_mirnas <- downloadHandler(
          filename <- function() {
            paste0(Sys.Date(), "_norm_expression_mirnas.csv")
          },
          content <- function(file) {
            write.csv(df_out_type, file, quote = F, row.names = F)
          }
        )
      }
      if(input$count_mirna_type == 2) {
        mirna_data_type <- mirna_expr[["log2_mirna_data_type"]]
        mirna_err_gene_tab <- mirna_expr[["log2_err_mirna_tab"]]
        df_out_type <- (as.data.frame(cbind(sample=rownames(mirna_data_type), mirna_data_type)))
        
        output$save_txt_expr_mirnas <- downloadHandler(
          filename <- function() {
            paste0(Sys.Date(), "_log2_expression_mirnas.txt")
          },
          content <- function(file) {
            write.table(df_out_type, file, quote = F, sep = "\t", col.names = T, row.names = F)
          }
        )
        output$save_csv_expr_mirnas <- downloadHandler(
          filename <- function() {
            paste0(Sys.Date(), "_log2_expression_mirnas.csv")
          },
          content <- function(file) {
            write.csv(df_out_type, file, quote = F, row.names = F)
          }
        )
      }
      output$mirna_table <- renderTable({
        mirna_data_type
      },
      striped = T,
      hover = T,
      rownames = T)
    })
  })
  
  # Outputs gene plots
  # Plot heights are scaled by number of genes to plot
  output$gene_plot <- renderPlot({
    geneplot_list <- press_submit()["gplot"]
    return(geneplot_list)
  },
  height = function() {
    use_height <- press_submit()[["num_genes"]]
    if(use_height > 0) {
      return(use_height*300)
    }
    else {return(300)}# Height of 1 plot
  })
  
  #### miRNA plot outputs ####
  # Outputs miRNA plots
  # Plot heights are scaled by number of miRNAs to plot
  output$mirna_plot <- renderPlot({
    mirnaplot_list <- press_submit()["mplot"]
    return(mirnaplot_list)
  },
  height = function() {
    use_height <- press_submit()[["num_mirnas"]]
    if(use_height > 0) {
      return(use_height*300)
    }
    else { return(300)} # Height of 1 plot
  })
  
  
  output$save_png_plot_genes <- downloadHandler(
    filename <- function() {
      paste0(Sys.Date(), "_gene_plots.png")
    },
    content <- function(file) {
      gene_plot <- press_submit()[["gplot"]]
      num_gene <- press_submit()[["num_genes"]]
      ggsave(file, gene_plot,
             device = "png", height = num_gene*3,
             width = 4)
    },
  )  
  
  output$save_png_plot_mirnas <- downloadHandler(
    filename <- function() {
      paste0(Sys.Date(), "_mirna_plots.png")
    },
    content <- function(file) {
      mirna_plot <- press_submit()[["mplot"]]
      num_mirna <- press_submit()[["num_mirnas"]]
      ggsave(file, mirna_plot,
             device = "png", height = num_mirna*3,
             width = 4)
    },
  )  
  
  output$save_pdf_plot_genes <- downloadHandler(
    filename <- function() {
      paste0(Sys.Date(), "_gene_plots.pdf")
    },
    content <- function(file) {
      gene_plot <- press_submit()[["gplot"]]
      num_gene <- press_submit()[["num_genes"]]
      ggsave(file, gene_plot,
             device = "pdf", height = num_gene*3+1,
             width = 4)
    },
  )  
  
  output$save_pdf_plot_mirnas <- downloadHandler(
    filename <- function() {
      paste0(Sys.Date(), "_mirna_plots.pdf")
    },
    content <- function(file) {
      mirna_plot <- press_submit()[["mplot"]]
      num_mirna <- press_submit()[["num_mirnas"]]
      ggsave(file, mirna_plot,
             device = "pdf", height = num_mirna*3+1,
             width = 4)
    },
  )  
  
  # Save button is on About page
  # Curated gene list from studies listed on About page
  output$save_genelists <- downloadHandler(
    filename <- function() {
      paste0(Sys.Date(), "_study_gene_lists.txt")
    },
    content <- function(file) {
      write.table(pub_genes_split, file,
                  quote = F, sep = "\t", col.names = T, row.names = F)
    }
  )
  
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
