library(shiny)

source("scripts.R")
# Define UI ----
ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  # shinyUI(
  navbarPage("Mouse Pituitary Gland Data Hub",
             tabPanel("Home",icon = icon("home"),
                      h1("MOUSE PITUITARY GLAND DATA HUB"),
                      br(),
                      p("Shiny app in development to visualize gene and miRNA expression data generated",
                        "from 3'UTR-seq and small RNA-seq respectively."),
                      p("Raw data can be accessed at ",
                        a("E-MTAB-9460.", 
                          href = "https://www.ebi.ac.uk/arrayexpress/")),
                      br(),
                      h2("Features"),
                      p("- Visualize gene and miRNA expression plots across postnatal ages and between sexes."),
                      p("- Quantify log2FC and FDR for DE genes and miRNAs for each comparison."),
                      p("- Intersect genes of interest with gene lists from relevant published studies."),
                      p("- Output normalized and log2(normalized) counts for genes and miRNAs of interest.")
             ),
             tabPanel("Data Browser", icon = icon("chart-bar"),
                      sidebarLayout(
                        sidebarPanel(
                          shinyjs::useShinyjs(),
                          id = "side-panel",
                          width = 3,
                          h3("Input genes"),
                          helpText(h5("Separate genes/miRNAs with a comma")),
                          textInput("gene", label = NULL,
                                    placeholder = "Lhb,let-7a-5p,mmu-let-7e-5p,ENSMUSG00000027120.7"
                                    # value ="Lhb,let-7a-5p,mmu-let-7e-5p,ENSMUSG00000027120.7,miR-224-5p"
                          ),
                          h3("or"),
                          helpText(h5("Upload a .txt file with genes/miRNAs on separate lines")),
                          fileInput("file", label = NULL, accept = ".txt"),
                          div(style="display:inline-block",
                              actionButton("submit",
                                           label = "Submit",
                                           class = "btn-success")),
                          div(style="display:inline-block",
                              actionButton("reset",
                                           label = "Reset")),
                          br(),
                          br(),
                          h3("Select studies to intersect"),
                          checkboxGroupInput("pub_study",
                                             label = h4("Puberty-related"),
                                             choices = list("Perry 2014" = "Perry2014",
                                                            "Day 2015" = "Day2015_GWAS_VB",
                                                            "Day 2017" = "Day2017_nearest",
                                                            "Hollis 2020" = "Hollis2020_GWAS_VB_FH",
                                                            "IHH/Kallmann" = "IHH/Kallmann"
                                             ),
                                             selected = NULL),
                          checkboxGroupInput("pit_study",
                                             label = h4("Pituitary-related"),
                                             choices = list("Ye 2015" = "Ye2015_PA_GWAS",
                                                            "Fang 2016" = "Fang2016_CPHD",
                                                            "Hauser 2019" = "Hauser2019_PA",
                                                            "Kurtoglu 2019" = "Kurtoglu2019_hypopituitarism"
                                             ),
                                             selected = NULL),
                          helpText(h5("It is recommended to input ~20 genes/miRNAs for viewing on the browser.")),
                          helpText(h5("If more genes/miRNAs are inputted, expression values and DE table can be downloaded."))
                        ),
                        mainPanel(
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
                                     uiOutput("add_download_gene_plot"),
                                     h3("Genes"),
                                     plotOutput("gene_plot",
                                                height = "auto",
                                                width = "70%")),
                              column(6,
                                     br(),
                                     uiOutput("add_download_mirna_plot"),
                                     h3("miRNAs"),
                                     plotOutput("mirna_plot",
                                                height = "auto",
                                                width = "70%"))
                            ),
                            tabPanel("DE gene table",
                                     br(),
                                     uiOutput("add_download_de_gene"),
                                     tableOutput("de_gene_table"),
                            ),
                            tabPanel("DE miRNA table",
                                     br(),
                                     uiOutput("add_download_de_mirna"),
                                     tableOutput("de_mirna_table")),
                            tabPanel("Gene expression values",
                                     br(),
                                     radioButtons("count_type",
                                                  label = h3("Count type"),
                                                  choices = list("Normalized counts" = 1, "log2(normCounts)" = 2),
                                                  selected = 1,
                                                  inline = T),
                                     br(),
                                     uiOutput("add_download_expr_gene"),
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
                                     uiOutput("add_download_expr_mirna"),
                                     tableOutput("mirna_table")
                            )
                          )
                        )) 
             ),
             tabPanel("About", icon=icon("quote-left"),
                      h2("About"),
                      p("Text about this app"),
                      br(),
                      h2("Credits")
             )
  )
  # )
)



# Define server logic ----
server <- function(input, output) {
  
  # Reset input values in the side-panel
  observeEvent(input$reset, {
    shinyjs::reset("side-panel")
  })
  
  # Add in save buttons only when valid genes/miRNAs are inputted.
  add_ui <- function(type) {
    if(type == "gene") {
      output$add_download_gene_plot <- renderUI({
        tagList(
          h4("Save gene expression plots: "),
          downloadButton("save_png_plot_genes", ".png"),
          downloadButton("save_pdf_plot_genes", ".pdf")
        )
      })
      output$add_download_de_gene <- renderUI({
        tagList(
          h4("Save DE gene table: "),
          downloadButton("save_txt_de_genes", ".txt"),
          downloadButton("save_csv_de_genes", ".csv")
        )
      })
      output$add_download_expr_gene <- renderUI({
        tagList(
          h4("Save gene expression values: "),
          downloadButton("save_txt_expr_genes", ".txt"),
          downloadButton("save_csv_expr_genes", ".csv")
        )
      })
    }
    
    if(type == "mirna") {
      output$add_download_mirna_plot <- renderUI({
        tagList(
          h4("Save miRNA expression plots: "),
          downloadButton("save_png_plot_mirnas", ".png"),
          downloadButton("save_pdf_plot_mirnas", ".pdf")
        )
      })
      output$add_download_de_mirna <- renderUI({
        tagList(
          h4("Save DE miRNA table: "),
          downloadButton("save_txt_de_mirnas", ".txt"),
          downloadButton("save_csv_de_mirnas", ".csv")
        )
      })
      output$add_download_expr_mirna <- renderUI({
        tagList(
          h4("Save miRNA expression values: "),
          downloadButton("save_txt_expr_mirnas", ".txt"),
          downloadButton("save_csv_expr_mirnas", ".csv")
        )
      })
    }
  }
  
  # Run functions in response to submit button
  # eventReactive events are delayed until the button is pressed
  press_submit <- eventReactive(input$submit, {
    
    text_input <- F
    file_input <- F
    
    # Check if genes are inputted by text box
    if(nchar(input$gene) > 0) {
      text_input <- T
    }
    
    # Check if a file is uploaded
    if(!is.null(input$file)) {
      file_input <- T
    }
    
    if(text_input & file_input) { # Print error message if both fields are filled.
      msg <- "Please choose one input method."
      parse_input <- parse_list("", type = "text")
      output$input_err <- renderText({
        msg
      })
    }
    else if(text_input == T & file_input == F) {
      output$input_err <- NULL
      parse_input <- parse_list(input$gene, type = "text")
      output$add_download_plots <- renderUI({
        tagList(
          h5("Save plots: "),
          downloadButton("save_png_plots", ".png"),
          downloadButton("save_pdf_plots", ".pdf")
        )
      })
    }
    else if(text_input == F & file_input == T) {
      output$input_err <- NULL
      file <- input$file
      read_file <- read.table(file$datapath, header = F)
      parse_input <- parse_list(read_file, type = "file")
      output$add_download_plots <- renderUI({
        tagList(
          h5("Save plots: "),
          downloadButton("save_png_plots", ".png"),
          downloadButton("save_pdf_plots", ".pdf")
        )
      })
    }
    else{ # Print error message if both fields are empty.
      parse_input <- parse_list("", type = "text")
      msg <- "Please input a list of genes/miRNAs."
      output$input_err <- renderText({
        msg
      })
    }
    
    output$invalid_genes <- NULL
    if(length(parse_input[["invalid"]]) > 0) {
      if(parse_input[["invalid"]] != ""){
        invalid_msg <- "Not found: "
        invalid_genes <- paste(parse_input[["invalid"]], collapse = ",")
        output$invalid_genes <- renderText({
          paste0(invalid_msg, invalid_genes)
        })
      }
    } else { output$invalid_genes <- NULL}
    
    
    gplot <- exprplot_hhtheme(genelist = parse_input[["genes"]],
                              count_data = utr_log2,
                              metadata = pData(utr_obj),
                              counttype = "genes")
    num_genes <- length(gplot$plot_env$use_genelist)
    
    if(num_genes > 0) {
      add_ui(type = "gene")
    }
    
    mplot <- exprplot_hhtheme(genelist = parse_input[["mirnas"]],
                              count_data = mirna_log2,
                              metadata = pData(mirna_obj),
                              counttype = "mirnas")
    num_mirnas <- length(mplot$plot_env$use_genelist)
    
    if(num_mirnas > 0) {
      add_ui(type = "mirna")
    }
    
    gtable <- print_de_table(genelist = parse_input[["genes"]],
                             de_table = utr_de_table,
                             counttype = "genes")
    mtable <- print_de_table(genelist = parse_input[["mirnas"]],
                             de_table = mirna_de_table,
                             counttype = "mirnas")
    norm_data_type <- expr_table(parse_input[["genes"]], utr_normcounts)
    log2_data_type <- expr_table(parse_input[["genes"]], utr_log2)
    
    norm_mirna_data_type <- expr_table(parse_input[["mirnas"]], mirna_normcounts)
    log2_mirna_data_type <- expr_table(parse_input[["mirnas"]], mirna_log2)
    return(list("gplot" = gplot,
                "num_genes" = num_genes,
                "mplot" = mplot,
                "num_mirnas" = num_mirnas,
                "gtable" = gtable,
                "mtable" = mtable,
                "norm_data_type" = norm_data_type,
                "log2_data_type" = log2_data_type,
                "norm_mirna_data_type" = norm_mirna_data_type,
                "log2_mirna_data_type" = log2_mirna_data_type))
  })
  
  # observed events occur the moment the button is pressed
  # Switches output table between normalized and log2norm gene counts
  observeEvent(input$submit, {
    observeEvent(input$count_type, {
      if(input$count_type == 1) {
        data_type <- press_submit()[["norm_data_type"]]
        err_gene_tab <- press_submit()[["norm_err_gene_tab"]]
      }
      if(input$count_type == 2) {
        data_type <- press_submit()[["log2_data_type"]]
        err_gene_tab <- press_submit()[["log2_err_gene_tab"]]
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
      if(input$count_mirna_type == 1) {
        mirna_data_type <- press_submit()[["norm_mirna_data_type"]]
        mirna_err_gene_tab <- press_submit()[["norm_err_mirna_tab"]]
      }
      if(input$count_mirna_type == 2) {
        mirna_data_type <- press_submit()[["log2_mirna_data_type"]]
        mirna_err_gene_tab <- press_submit()[["log2_err_mirna_tab"]]
      }
      output$mirna_table <- renderTable({
        mirna_data_type
      },
      striped = T,
      hover = T,
      rownames = T)
    })
  })
  
  observeEvent(input$submit, {
    observeEvent(input$pub_study, {
      if(length(input$pub_study) > 0) {
        output$de_gene_table <- renderTable({
          de_genes <- press_submit()[["gtable"]]
          return(add_study(de_genes,
                           input$pub_study))
        },
        striped = T,
        hover = T,
        digits = 5)
      }
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
      return(use_height*375)
    }
    else {return(350)}# Height of 1 plot
  })
  
  # Outputs miRNA plots
  # Plot heights are scaled by number of miRNAs to plot
  output$mirna_plot <- renderPlot({
    mirnaplot_list <- press_submit()["mplot"]
    return(mirnaplot_list)
  },
  height = function() {
    use_height <- press_submit()[["num_mirnas"]]
    if(use_height > 0) {
      return(use_height*375)
    }
    else { return(350)} # Height of 1 plot
  })
  
  # Outputs DE gene table
  output$de_gene_table <- renderTable({
    press_submit()[["gtable"]]
  },
  striped = T,
  hover = T,
  digits = 5)
  
  # Outputs DE miRNA table
  output$de_mirna_table <- renderTable({
    press_submit()[["mtable"]]
  },
  striped = T,
  hover = T,
  digits = 5)
  
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
  
  
  # output$save_pdf_plots <- downloadHandler(
  #   filename <- function() {
  #     paste0(Sys.Date(), "_plots.pdf")
  #   },
  #   content <- function(file) {
  #     plot_gene <- F
  #     plot_mirna <- F
  #     
  #     gene_plot <- press_submit()[["gplot"]]
  #     num_gene <- press_submit()[["num_genes"]]
  #     
  #     mirna_plot <- press_submit()[["mplot"]]
  #     num_mirna <- press_submit()[["num_mirnas"]]
  #     
  #     if(num_gene > 0) {
  #       plot_gene <- T
  #     }
  #     if(num_mirna > 0) {
  #       plot_mirna <- T
  #     }
  #     
  #     if(plot_gene & plot_mirna) {
  #       use_height <- max(num_gene, num_mirna)
  #       
  #       ggsave(file,
  #              arrangeGrob(gene_plot,mirna_plot, ncol = 2),
  #              device = "pdf", height = use_height*3,
  #              width = 8)
  #     }
  #     else if(plot_gene & plot_mirna == F) {
  #       use_height <- num_gene
  #       ggsave(file, gene_plot,
  #              device = "pdf", height = use_height*3,
  #              width = 4)
  #     }
  #     else {
  #       use_height <- num_mirna
  #       ggsave(file, mirna_plot,
  #              device = "pdf", height = use_height*3,
  #              width = 4)
  #     }
  #   },
  # )  
  # 
  output$save_txt_de_genes <- downloadHandler(
    filename <- function() {
      paste0(Sys.Date(), "_de_genes.txt")
    },
    content <- function(file) {
      write.table(press_submit()[["gtable"]], file,
                  quote = F, sep = "\t", col.names = T, row.names = F)
    }
  )
  
  output$save_txt_de_mirnas <- downloadHandler(
    filename <- function() {
      paste0(Sys.Date(), "_de_mirnas.txt")
    },
    content <- function(file) { 
      write.table(press_submit()[["mtable"]], file,
                  quote = F, sep = "\t", col.names = T, row.names = F)
    }
  )
  
  output$save_csv_de_genes <- downloadHandler(
    filename <- function() {
      paste0(Sys.Date(), "_de_genes.csv")
    },
    content <- function(file) {
      write.csv(press_submit()[["gtable"]], file,
                quote = F, row.names = F)
    }
  )
  
  output$save_csv_de_mirnas <- downloadHandler(
    filename <- function() {
      paste0(Sys.Date(), "_de_mirnas.csv")
    },
    content <- function(file) {
      write.csv(press_submit()[["mtable"]], file,
                quote = F, row.names = F)
    }
  )
  
  output$save_txt_expr_genes <- downloadHandler(
    filename <- function() {
      paste0(Sys.Date(), "_expression_genes.txt")
    },
    content <- function(file) {
      df_out <- press_submit()
      if(input$count_type == 1) {
        df_out_type <- df_out[["norm_data_type"]]
        df_out_type <- (as.data.frame(cbind(sample=rownames(df_out_type), df_out_type)))
        
      }
      if(input$count_type == 2) {
        df_out_type <- df_out[["log2_data_type"]]
        df_out_type <- (as.data.frame(cbind(sample=rownames(df_out_type), df_out_type)))
      }
      write.table(df_out_type, file,
                  quote = F, sep = "\t", col.names = T, row.names = F)
    }
  )
  
  output$save_csv_expr_genes <- downloadHandler(
    filename <- function() {
      paste0(Sys.Date(), "_expression_genes.csv")
    },
    content <- function(file) {
      df_out <- press_submit()
      if(input$count_type == 1) {
        df_out_type <- df_out[["norm_data_type"]]
        df_out_type <- (as.data.frame(cbind(sample=rownames(df_out_type), df_out_type)))
        
      }
      if(input$count_type == 2) {
        df_out_type <- df_out[["log2_data_type"]]
        df_out_type <- (as.data.frame(cbind(sample=rownames(df_out_type), df_out_type)))
      }
      write.csv(df_out_type, file,
                quote = F, row.names = F)
    }
  )
  
  output$save_txt_expr_mirnas <- downloadHandler(
    filename <- function() {
      paste0(Sys.Date(), "_expression_mirnas.txt")
    },
    content <- function(file) {
      df_out <- press_submit()
      if(input$count_type == 1) {
        df_out_type <- df_out[["norm_mirna_data_type"]]
        df_out_type <- (as.data.frame(cbind(sample=rownames(df_out_type), df_out_type)))
        
      }
      if(input$count_type == 2) {
        df_out_type <- df_out[["log2_mirna_data_type"]]
        df_out_type <- (as.data.frame(cbind(sample=rownames(df_out_type), df_out_type)))
      }
      write.table(df_out_type, file,
                  quote = F, sep = "\t", col.names = T, row.names = F)
    }
  )
  
  output$save_csv_expr_mirnas <- downloadHandler(
    filename <- function() {
      paste0(Sys.Date(), "_expression_mirnas.csv")
    },
    content <- function(file) {
      df_out <- press_submit()
      if(input$count_type == 1) {
        df_out_type <- df_out[["norm_mirna_data_type"]]
        df_out_type <- (as.data.frame(cbind(sample=rownames(df_out_type), df_out_type)))
        
      }
      if(input$count_type == 2) {
        df_out_type <- df_out[["log2_mirna_data_type"]]
        df_out_type <- (as.data.frame(cbind(sample=rownames(df_out_type), df_out_type)))
      }
      write.csv(df_out_type, file,
                quote = F, row.names = F)
    }
  )
  
}

# Run the app ----
shinyApp(ui = ui, server = server)