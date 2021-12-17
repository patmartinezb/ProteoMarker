
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

library(BiocManager)


# Import libraries
lib <-  c("shiny", "shinydashboard", "DT", "shinyWidgets", "vroom", "DAtest", "markdown",
          "readr", "dplyr", "ggplot2", "stringr", "ggsci", "readr", "tidyr",
          "pheatmap", "reshape2", "imputeLCMD", "preprocessCore", "ggrepel", "shinycssloaders",
          "factoextra", "limma", "caret", "scales", "glmnet", "ROCR", "pvca", 
          "yardstick", "pcaMethods", "impute", "randomcoloR")


# Checking missing packages from list
new.packages <- setdiff(lib, installed.packages())

if (length(new.packages) != 0) invisible(lapply(new.packages, BiocManager::install, update = FALSE))

# DAtest

if (!require("DAtest", quietly = TRUE)){
  remotes::install_github("Russel88/DAtest", upgrade = "never", dependencies = TRUE, quiet = TRUE)
}

# Load
invisible(lapply(lib, require, character.only = TRUE))



## ui.R ##
sidebar <- dashboardSidebar(
    sidebarMenu(
        menuItem("Upload data", tabName = "upload", icon = icon("open-file", lib = "glyphicon")),
        menuItem("QC and filtering", icon = icon("wrench", lib = "glyphicon"),
                 menuItem("Missing values", tabName = "ms"),
                 menuItem("QC figures", tabName = "qc_figs"),
                 menuItem("Filtering", tabName = "filt")),
        menuItem("Imputation", tabName = "imp", icon = icon("edit", lib = "glyphicon")),
        menuItem("Normalization", tabName = "norm", icon = icon("balance-scale")),
        menuItem("Batch effect", tabName = "batch_all", icon = icon("search", lib = "glyphicon"),
                 menuItem("Check for batch effect", tabName = "batch"),
                 menuItem("Remove batch effect?", tabName = "rm_batch")),
        menuItem("DE analysis", tabName = "de", icon = icon("sort", lib = "glyphicon"),
                 menuItem("limma (2-group comparison)", tabName = "limma"),
                 menuItem("Student's T test", tabName = "ttest"),
                 menuItem("ANOVA", tabName = "aov")
                 ),
        menuItem("Biomarker selection", tabName = "bio_select", icon = icon("screenshot", lib = "glyphicon"),
                 menuItem("Data partition", tabName = "data"),
                 menuItem("Lasso regression", tabName = "lasso"),
                 menuItem("Evaluation", tabName = "eval")),
        menuItem("Help", tabName = "help", icon = icon("info-sign", lib = "glyphicon"))
        
    )
)

body <- dashboardBody(
    tabItems(
        tabItem(tabName = "upload",
                fluidRow(
                    
                    box(title = "Upload metadata", status = "primary", solidHeader = TRUE,
                        prettyCheckbox("ref_channel", label = HTML(paste0("<b>","Is there a reference channel?","</b>")), 
                                       icon = icon("check")),
                        fileInput("meta", "", buttonLabel = "Upload...", accept = c(".csv", ".txt")),
                        dataTableOutput("meta_table"),
                        "Supported formats are: .csv, .txt"),
                    
                    box(pickerInput(inputId = "origin", label = "Software used to process raw files:", width = "55%",
                                    choices = c("MaxQuant", "Proteome Discoverer or others"), 
                                    options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), 
                                    multiple = FALSE),
                        title = "Upload raw data matrix", status = "primary", solidHeader = TRUE,
                        fileInput("file", "", buttonLabel = "Upload...", accept = c(".csv", ".txt")),
                        dataTableOutput("file_table"),
                        "Supported formats are: .csv, .txt")
                ),
                fluidRow(
                    box(width = 12, title = "Trimmed data", "For MQ: Variables were selected, the normalization channel and", 
                        "contaminants were removed, 0 transformed to NAs and values were log2-transformed.", br(),
                        "For Proteome Discoverer and others: Variables were selected, 0 transformed to NAs and values were log2-transformed.",
                        status = "success", solidHeader = TRUE,
                        br(),
                        dataTableOutput("clean_raw"),
                        div(downloadButton("dl_raw", "Download .csv"), style="float:right"))
                )
        ),

        tabItem(tabName = "ms",
                fluidRow(
                    box(title = "Percentage of missing values per sample", status = "warning",
                        plotOutput("ms_pcts", height = "600px"),
                        selectInput("ext1", "Plot extension:", c("png", "jpeg", "svg"), width = "15%"),
                        downloadButton("dl_pct", label = "Download plot",  icon = icon("floppy-save", lib = "glyphicon"))),
                    
                    box(title = "Heatmap of missing values", status = "warning",
                        prettyCheckboxGroup("na_vars", "Select variables:", choices = NULL, inline = TRUE),
                        plotOutput("ms_heatmap", height = "650px"))
                ),
                fluidRow(
                    box(title = "Table: number of missing values per sample", status = "success",
                        dataTableOutput("ms_table"))
                    )
                ),
        
        tabItem(tabName = "qc_figs",
                fluidRow(
                    box(title = "Number of detected proteins per sample", status = "warning",
                        pickerInput(inputId = "n_prots_vars", label = "Select variables:", choices = NULL, 
                                    options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), 
                                    multiple = FALSE),
                        plotOutput("n_prots", height = "600px")),
                    
                    box(title = "Correlation among samples", status = "warning",
                        plotOutput("qc_corr", height = "600px"))
                    ),
                
                fluidRow(
                    box(title = "Boxplots of log2 values", status = "warning",
                        materialSwitch(inputId = "reorder", label = "Reorder by log2 abundance?", value = TRUE, status = "primary"),
                        pickerInput(inputId = "qc_fill_vars", label = "Select variables:", choices = NULL, 
                                    options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), 
                                    multiple = FALSE),
                        plotOutput("qc_boxplot", height = "600px")),
                    
                    box(title = "Density plot", status = "warning",
                        pickerInput(inputId = "den_vars", label = "Select variables:", choices = NULL, 
                                    options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), 
                                    multiple = FALSE),
                        
                        plotOutput("qc_den", height = "600px"))
                    )
                ),
        
        tabItem(tabName = "filt",
                fluidRow(
                    box(width = 12, title = "Filter data", solidHeader = TRUE, status = "success", "The threshold ranges from 0 to 1,",
                        "where 0 allows a 100% of missing values of a protein across all samples, and 1",
                        "allows only proteins detected across all samples.",
                        sliderInput("threshold", "", value = 0.7, min = 0, max = 1, round = TRUE, width = "50%"),
                        dataTableOutput("filt_table")),

                ),
                fluidRow(
                    box(title = "Heatmap after filtering", status = "warning",
                        prettyCheckboxGroup("filt_vars", "Select variables:", choices = NULL, inline = TRUE),
                        plotOutput("filt_heatmap", height = "600px")),
                    
                    box(title = "Number of detected proteins per sample after filtering", status = "warning",
                        plotOutput("filt_n_prots", height = "600px"))
                )
        ),
        
        tabItem(tabName = "imp",
                fluidRow(
                    box(width = 12, title = "Table with imputed data", status = "success",
                        materialSwitch(inputId = "imp_yn", label = "Do you want to impute?", value = TRUE, status = "success"),
                        pickerInput(inputId = "type_imp", label = "Select type of imputation:", 
                                    choices = c("MinDet", "MinProb", "Min", "KNN"), 
                                    options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), 
                                    multiple = FALSE, width = "30%"),
                        dataTableOutput("imp_table"),
                        div(downloadButton("dl_imp", "Download .csv"), style="float:right"))),
                
                fluidRow(
                  box(title = "Density plot after imputation", status = "warning",
                      pickerInput(inputId = "imp_den_vars", label = "Select variables:", choices = NULL, 
                                         options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), 
                                  multiple = FALSE),
                      plotOutput("imp_den", height = "700px")),
                  
                  box(title = "Boxplots after imputation", status = "warning",
                      plotOutput("imp_boxplot", height = "600px")))

                ),
        
        tabItem(tabName = "norm",
                fluidRow(
                  box(width = 12, title = "Table with normalized data", status = "success",
                      materialSwitch(inputId = "norm_yn", label = "Do you want to normalize?", value = TRUE, status = "success"),
                      pickerInput(inputId = "type_norm", label = "Select type of normalization:", 
                                  choices = c("quantNorm", "meanNorm", "medianNorm"), 
                                  options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), 
                                  multiple = FALSE, width = "30%"),
                      dataTableOutput("norm_table"),
                      div(downloadButton("download3", "Download .csv"), style="float:right"))
                ),
                
                fluidRow(
                  box(width = 12, status = "warning",
                    plotOutput("norm_boxplot", height = "650px", width = "98%"))
                  )
                ),
        
        tabItem(tabName = "batch",
                
                fluidRow(
                  box(title = "PCA", status = "warning",
                      pickerInput(inputId = "batch_pca_vars", label = "Select variable:", choices = NULL, 
                                  options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), 
                                  multiple = FALSE, width = "45%"),
                      plotOutput("batch_pca", height = "550px")),
                  
                  box(title = "Heatmap", status = "warning",
                      prettyCheckboxGroup("batch_vars", "Select variables:", choices = NULL, inline = TRUE),
                      plotOutput("batch_heatmap", height = "600px"))
                  ),
                
                fluidRow(
                  box(width = 12, title = "PVCA", status = "warning",
                      plotOutput("pvca", height = "650px"))
                  )
                ),
        
        tabItem(tabName = "rm_batch",
                
                fluidRow(
                  box(width = 12, title = "Remove batch effect", status = "primary",
                      pickerInput(inputId = "batch_rm_vars1", label = "Select variable 1:", choices = NULL, 
                                  options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), 
                                  multiple = FALSE, width = "45%"),
                      pickerInput(inputId = "batch_rm_vars2", label = "Select variable 2:", choices = NULL, 
                                  options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), 
                                  multiple = FALSE, width = "45%"),
                      actionButton("reset", "Reset"),
                      dataTableOutput("rm_batch_table"),
                      downloadButton("download1", "Download .csv"))
                  ),
                
                fluidRow(
                  box(title = "PCA after batch effect removal", status = "warning",
                      pickerInput(inputId = "batch_pca_vars2", label = "Select variable:", choices = NULL, 
                                  options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), 
                                  multiple = FALSE, width = "45%"),
                      plotOutput("batch_pca2", height = "550px")),
                  
                  box(title = "Heatmap after batch effect removal", status = "warning",
                      prettyCheckboxGroup("batch_vars2", "Select variables:", choices = NULL, inline = TRUE),
                      plotOutput("batch_heatmap2", height = "550px"))
                  )
                ),
        
        tabItem(tabName = "limma",

                fluidRow(
                  box(width = 9, title = "Differential expression analysis - limma (2-group comparison)", status = "primary",
                      pickerInput(inputId = "control_trt_limma", label = "Which one is the control group?", choices = NULL,
                                  options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                                  multiple = FALSE, width = "45%"),
                      materialSwitch(inputId = "df_use_limma", label = "Do you want to use the batch-corrected dataset?",
                                     value = FALSE, status = "primary"),
                      prettyCheckboxGroup("covars_limma", "Select covariables for modeling (if any):",
                                          choices = NULL, inline = TRUE),
                      dataTableOutput("res_limma"),
                      div(downloadButton("download_limma", "Download .csv"), style="float:right")
                      ),
                  box(width = 3, title = "Results table", status = "primary",
                      dataTableOutput("de_limma"))
                ),
                
                fluidRow(
                  box(title = "Summary graphics", status = "warning",
                      plotOutput("figs_limma", height = "600px")),
                  box(title = "Volcano plot", status = "warning",
                      plotOutput("volcano_limma", height = "600px"))
                )

        ),
        
        tabItem(tabName = "ttest",
                
                fluidRow(
                  box(width = 9, title = "Differential expression analysis - T test (2-group comparison)", status = "primary",
                      pickerInput(inputId = "control_trt", label = "Which one is the control group?", choices = NULL, 
                                options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), 
                                multiple = FALSE, width = "45%"),
                      materialSwitch(inputId = "df_use", label = "Do you want to use the batch-corrected dataset?", 
                                    value = TRUE, status = "primary"),
                      dataTableOutput("res_ttest"),
                      div(downloadButton("download2", "Download .csv"), style="float:right")),
                  box(width = 3, title = "Results table", status = "primary",
                      dataTableOutput("de_ttest"))
                  ),
                
                fluidRow(
                  box(title = "Summary graphics", status = "warning",
                      plotOutput("figs_res", height = "600px")),
                  box(title = "Volcano plot", status = "warning",
                      plotOutput("volcano_ttest", height = "600px"))
                  )
                ),
        
        
        tabItem(tabName = "aov",
                
                fluidRow(
                  box(width = 9, title = "Differential expression analysis - ANOVA (2+ comparison)", status = "primary",
                      materialSwitch(inputId = "df_use_aov", label = "Do you want to use the batch-corrected dataset?", 
                                     value = FALSE, status = "primary"),
                      pickerInput(inputId = "control_aov", label = "Choose a group to launch the analysis", choices = NULL, 
                                  options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), 
                                  multiple = FALSE, width = "45%"),
                      uiOutput("cityControls"),
                      withSpinner(dataTableOutput("res_aov")),
                      div(downloadButton("download_aov", "Download .csv"), style="float:right")),
                  box(width = 3, title = "Results table", status = "primary",
                      dataTableOutput("de_aov"))
                ),
                
                fluidRow(
                  box(title = "Summary graphics", status = "warning",
                      plotOutput("figs_aov", height = "600px")),
                  box(title = "Volcano plot", status = "warning",
                      plotOutput("volcano_aov", height = "600px"))
                )
        ),
        
        
        
        tabItem(tabName = "data",
                
                fluidRow(
                  box(width = 9, title = "Data partition into train and test subsets, for Lasso regression", status = "primary",
                      pickerInput(inputId = "control_trt2", label = "Which one is the control group?", choices = NULL, 
                                  options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), 
                                  multiple = FALSE, width = "50%"),
                      pickerInput(inputId = "res_use", label = "Which result do you want to use?", choices = c("limma", "T test", "ANOVA"),
                                  options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                                  multiple = FALSE, width = "30%"),
                      "The partition porcentage ranges from 0 to 1. If 0.7 is chosen, it means that the train dataset",
                      "will receive 70% of all samples - always maintaining the original proportion among groups.",
                      sliderInput("part", "", value = 0.7, min = 0, max = 1, round = TRUE, width = "50%"),
                      materialSwitch(inputId = "df_use2", label = "Do you want to use the batch-corrected dataset?", 
                                     value = FALSE, status = "primary"),
                      plotOutput("bio_dt", width = "80%", height = "600px"))
                  ),
                
                fluidRow(
                  box(dataTableOutput("dat_train"), title = "Train data", status = "warning"),
                  box(dataTableOutput("dat_test"), title = "Test data", status = "warning")
                  )
                ),
        
        tabItem(tabName = "lasso",
                
                fluidRow(
                  box(title = "Evolution of coefficients with respect to lambda", status = "warning",
                      plotOutput("lambda", height = "600px")),
                  box(title = "Crossvalidation to get the best lambda", status = "warning",
                      plotOutput("error", height = "600px"))
                  ),
                
                fluidRow(
                  box(title = "Representation of the model coefficients chosen my the model", status = "warning",
                      plotOutput("coefs_plot", height = "600px")),
                  box(dataTableOutput("coefs_dt"), title = "Predictors chosen by the model:", status = "warning")
                  )
                ),
        
        tabItem(tabName = "eval",
                
                fluidRow(
                  box(title = "AUC", status = "warning",
                      pickerInput(inputId = "auc_choice", label = "What AUC plot do you want to visualize?", 
                                  choices = c("Predictor by predictor", "As a whole"), 
                                  options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), 
                                  multiple = FALSE, width = "50%"),
                      plotOutput("auc", height = "600px"), width = 6),
                  box(plotOutput("conf_m_train"), title = "Confussion matrix for the train data:", status = "primary",
                      width = 3),
                  box(plotOutput("conf_m_test"), title = "Confussion matrix for the test data:", status = "primary",
                      width = 3)
                  ),
                
                fluidRow(
                  box(title = "PCA with Lasso predictors", status = "warning",
                      plotOutput("pca_lasso", height = "600px")),
                  box(title = "Heatmap with Lasso predictors", status = "warning",
                      plotOutput("heatmap_lasso", height = "600px"))
                  )
                ),
        
        tabItem(tabName = "help",
                
                fluidRow(
                  box(width = 9,
                      includeMarkdown("www/include.md")
                      )
                )
                
                )
        
        )
)

# Put them together into a dashboardPage
ui <- dashboardPage(
    dashboardHeader(title = "ProteoMarker"),
    sidebar,
    body
)



server <- function(input, output, session) {
    
  # Create reactives for starting datasets
  
    data <- reactive({
        
      req(input$file)
      
      ext <- tools::file_ext(input$file$name)
      switch(ext,
             csv = vroom::vroom(input$file$datapath, delim = ","),
             txt = vroom::vroom(input$file$datapath, delim = "\t"),
             validate("Invalid file; Please upload a .csv or .txt file")
      )
    })
    
    annot <- reactive({
        
      req(input$meta)
      
      ext <- tools::file_ext(input$meta$name)
      switch(ext,
             csv = vroom::vroom(input$meta$datapath, delim = ","),
             txt = vroom::vroom(input$meta$datapath, delim = "\t"),
             validate("Invalid file; Please upload a .csv or .txt file")
      )
    })
    
    meta <- reactive({
        
      req(input$meta)
      
      if (input$ref_channel == TRUE){
      
      filter_reference(annot())
        
      } else {
        
        clean_key(annot())
        
      }
    })
    
    
  # Dataset cleaning
    
    prot_mbr <- reactive({
        
      req(input$meta)
      
      if (input$origin == "MaxQuant"){
        
        select_mq_vars(data(), meta())
        
      } else {
        
        select_no_mq(data(), meta())
        
      }
      
      
    })
    
    dataf <- reactive({
      
      req(input$meta)
      
      if (input$origin == "MaxQuant"){
      
      clean_mq_data(prot_mbr())
        
      } else {
        
        clean_no_mq(prot_mbr())
        
      }
    })
    
  # Missing values tab
  ## Update variables for heatmap
    
    select_na_vars <- reactive({
        req(input$meta)
        
        meta() %>% 
            select(-Run,
                   -TechRepMixture,
                   -Fraction,
                   -Channel,
                   -BioReplicate,
                   -key) %>% 
            select_if(is.character) %>% 
            colnames(.)
        
    })
    
    observeEvent(input$meta,{
        updatePrettyCheckboxGroup(session, 
                                          "na_vars", "Select variables:", choices = select_na_vars(),
                                  inline = TRUE, selected = "Condition")
    })
    
  ## Create MS dataset for plotting
    
    missing.values <- reactive({
        
      req(input$meta)
      
      mis_vals(dataf())
        
    })
  
  # QC tab
  ## Update variables for number of prots barplot
    
    select_nprot_vars <- reactive({
      req(input$meta)
      
      meta() %>% 
        select(-Run,
               -TechRepMixture,
               -Fraction,
               -Channel,
               -key,
               -BioReplicate) %>% 
        select_if(is.character) %>% 
        colnames(.)
      
    })
    
    observeEvent(input$meta,{
      updatePickerInput(session, 
                        "n_prots_vars", "Select variables:", choices = select_nprot_vars(),
                        options = list(
                          `actions-box` = TRUE, 
                          size = 10,
                          `selected-text-format` = "count > 3"
                        ))
    })
    
  ## Update variables for density plot
    
    select_den_vars <- reactive({
        req(input$meta)
        
        meta() %>% 
            select(-Run,
                   -TechRepMixture,
                   -Fraction,
                   -Channel,
                   -key) %>% 
            select_if(is.character) %>% 
            colnames(.)
        
    })
    
    observeEvent(input$meta,{
      updatePickerInput(session, 
                        "den_vars", "Select variables:", choices = select_den_vars(),
                        options = list(
                          `actions-box` = TRUE, 
                          size = 10,
                          `selected-text-format` = "count > 3"
                        ))
    })
  
  ## Update variables for boxplot filling
    
    select_qc_fill_vars <- reactive({
      req(input$meta)
      
      meta() %>% 
        select(-Run,
               -TechRepMixture,
               -Fraction,
               -Channel,
               -BioReplicate,
               -key) %>% 
        select_if(is.character) %>% 
        colnames(.)
      
    })
    
    observeEvent(input$meta,{
      updatePickerInput(session, 
                        "qc_fill_vars", "Select variables:", choices = select_qc_fill_vars(),
                        options = list(
                          `actions-box` = TRUE, 
                          size = 10,
                          `selected-text-format` = "count > 3"
                        ))
    })
    
  # Filtering tab
    
    filt <- reactive({
      
      req(input$meta)  
      
      filt_nas(dataf(), threshold = input$threshold)
    })
    
  ## Update variables for heatmap
    
    select_filt_vars <- reactive({
        req(input$meta)
        
        meta() %>% 
            select(-Run,
                   -TechRepMixture,
                   -Fraction,
                   -Channel,
                   -BioReplicate,
                   -key) %>% 
            select_if(is.character) %>% 
            colnames(.)
        
    })
    
    observeEvent(input$meta,{
        updatePrettyCheckboxGroup(session, 
                                  "filt_vars", "Select variables:", choices = select_na_vars(),
                                  inline = TRUE, selected = "Condition")
    })
    
  # Imputation tab
    
    
    df.imp <- reactive({
      
      req(input$meta)
      
      if (input$imp_yn == TRUE){
        
        req(input$type_imp)
        
        impute_nas(filt(), input$type_imp)
        
      } else {
        
        no_impute(filt())
        
      }
        
    })
    
  ## Update variables for double density plot
    
    select_imp_den_vars <- reactive({
      req(input$meta)
      
      meta() %>% 
        select(-Run,
               -TechRepMixture,
               -Fraction,
               -Channel,
               -key) %>% 
        select_if(is.character) %>% 
        colnames(.)
      
    })
    
    observeEvent(input$meta,{
      updatePickerInput(session, 
                         "imp_den_vars", "Select variables:", choices = select_imp_den_vars(),
                        options = list(
                          `actions-box` = TRUE, 
                          size = 10,
                          `selected-text-format` = "count > 3"
                        ))
    })
    
  # Normalization tab
    
    df_norm <- reactive({
      
      req(input$meta)
      
      if (input$norm_yn == TRUE){
        
        req(input$type_norm)
        
        norm_prot(df.imp(), input$type_norm)
        
      } else {
        
        df.imp()
        
      }
      
      
      
    })
    
  # Check batch effect tab
  ## Update variables for PCA
    
    select_batch_pca_vars <- reactive({
      req(input$meta)
      
      meta() %>% 
        select(-Run,
               -TechRepMixture,
               -Fraction,
               -Channel,
               -key,
               -BioReplicate) %>% 
        select_if(is.character) %>% 
        colnames(.)
      
    })
    
    observeEvent(input$meta,{
      updatePickerInput(session, 
                        "batch_pca_vars", "Select variable:", choices = select_batch_pca_vars(),
                        options = list(
                          `actions-box` = TRUE, 
                          size = 10,
                          `selected-text-format` = "count > 3"
                        ))
    })
    
  ## Update variables for heatmap
    
    select_batch_vars <- reactive({
      req(input$meta)
      
      meta() %>% 
        select(-Run,
               -TechRepMixture,
               -Fraction,
               -Channel,
               -BioReplicate,
               -key) %>% 
        select_if(is.character) %>% 
        colnames(.)
      
    })
    
    observeEvent(input$meta,{
      updatePrettyCheckboxGroup(session, 
                                "batch_vars", "Select variables:", choices = select_batch_vars(),
                                inline = TRUE, selected = "Condition")
    })
    
    
  # Remove batch effect tab
  ## Update variables for removeBatchEffect function
    
    select_rmbatch_vars1 <- reactive({
      req(input$meta)
      
      meta() %>% 
        select(-Run,
               -TechRepMixture,
               -Fraction,
               -Channel,
               -BioReplicate,
               -key) %>% 
        colnames(.)
      
    })
    
    observeEvent(input$meta,{
      updatePickerInput(session, 
                                "batch_rm_vars1", "Select variable 1:", choices = select_rmbatch_vars1(),
                                selected = NA)
    })
    
    select_rmbatch_vars2 <- reactive({
      req(input$meta)
      
      meta() %>% 
        select(-Run,
               -TechRepMixture,
               -Fraction,
               -Channel,
               -BioReplicate,
               -key) %>% 
        colnames(.)
      
    })
    
    observeEvent(input$meta,{
      updatePickerInput(session, 
                        "batch_rm_vars2", "Select variable 2:", choices = select_rmbatch_vars2(),
                        selected = NA)
    })
  
  ## Creat new dataset without batch effect
    
    df_norm_wb <- reactive({
      
      req(input$batch_rm_vars1)
      
      if (is.null(input$batch_rm_vars1)){
        
        df_norm()
        
      }

      
      if (is.null(input$batch_rm_vars2)){
        
        batch_rm(df_norm(), meta(), batch1 = input$batch_rm_vars1)
        
      } else {
        
        req(input$batch_rm_vars2)
        
        batch_rm(df_norm(), meta(), batch1 = input$batch_rm_vars1, batch2 = input$batch_rm_vars2)
      }
      
    })
    
  ## Reset button
    
    observeEvent(input$reset, {
      updatePickerInput(session, 
                        "batch_rm_vars1", "Select variable 1:", choices = select_rmbatch_vars2(),
                        selected = NA)
      
      updatePickerInput(session, 
                        "batch_rm_vars2", "Select variable 2:", choices = select_rmbatch_vars2(),
                        selected = NA)
    })
    
  ## Update variables for plots after batch effect removal
    
    select_batch_pca_vars2 <- reactive({
      req(input$meta)
      
      meta() %>% 
        select(-Run,
               -TechRepMixture,
               -Fraction,
               -Channel,
               -key,
               -BioReplicate) %>% 
        select_if(is.character) %>% 
        colnames(.)
      
    })
    
    observeEvent(input$meta,{
      updatePickerInput(session, 
                        "batch_pca_vars2", "Select variable:", choices = select_batch_pca_vars2(),
                        options = list(
                          `actions-box` = TRUE, 
                          size = 10,
                          `selected-text-format` = "count > 3"
                        ))
    })
    
    select_batch_vars2 <- reactive({
      req(input$meta)
      
      meta() %>% 
        select(-Run,
               -TechRepMixture,
               -Fraction,
               -Channel,
               -BioReplicate,
               -key) %>% 
        select_if(is.character) %>% 
        colnames(.)
      
    })
    
    observeEvent(input$meta,{
      updatePrettyCheckboxGroup(session, 
                                "batch_vars2", "Select variables:", choices = select_batch_vars2(),
                                inline = TRUE, selected = "Condition")
    })
    
    
  # limma tab
  ## Select which one is the control group
    select_control_limma <- reactive({
      req(input$meta)
      
      meta() %>% 
        select(Condition) %>% 
        pull() %>% 
        as.factor() %>% 
        levels()
      
    })
    
    observeEvent(input$meta,{
      updatePickerInput(session, 
                        "control_trt_limma", "Which one is the control group?", choices = select_control_limma(),
                        selected = NA)
    })
    
    ## Update covariables
    
    select_covars <- reactive({
      req(input$meta)
      
      meta() %>% 
        select(-Run,
               -TechRepMixture,
               -Fraction,
               -Channel,
               -BioReplicate,
               -key) %>% 
        select_if(is.character) %>% 
        colnames(.)
      
    })
    
    observeEvent(input$meta,{
      updatePrettyCheckboxGroup(session, 
                                "covars_limma", "Select covariables for modeling (if any):", 
                                choices = select_covars(), inline = TRUE, selected = NULL)
    })
    
  ## Create res dataframe
    dat_limma <- reactive({
      
      req(input$control_trt_limma)
      
      if (input$df_use_limma == TRUE){
        
        validate(
          need(input$batch_rm_vars1, 'There is not a batch-corrected data matrix; consider modeling batch.'))
        
        de_limma(df_norm_wb(), meta(), control.name = input$control_trt_limma, covars = NULL)
        
      } else {
        
        if (is.null(input$covars_limma)){
          
          de_limma(df_norm(), meta(), control.name = input$control_trt_limma, covars = NULL)
          
        } else {
          
          req(input$covars_limma)
          
          df <- df_norm()
          meta <- meta()
          
          rownames(df) <- df$Protein.IDs
          
          df <- select_if(df, is.numeric)
          
          fac <- as.factor(meta$Condition)
          
          name.treat = levels(fac)[levels(fac) != input$control_trt_limma]
          
          pred <- factor(meta$Condition, levels = c(name.treat, input$control_trt_limma))
          
          DA.lim(data = as.matrix(df), predictor = pred, relative = FALSE, covars = get.covars(input$covars_limma, meta)) %>%
            select(- ordering,
                   - Method) %>%
            mutate(res = case_when((logFC > 0) & (pval.adj < 0.05) ~ "Up",
                                   (logFC < 0) & (pval.adj < 0.05) ~ "Down",
                                   TRUE ~ "no DE")) %>%
            rename(Protein.IDs = Feature,
                   log_fc = logFC,
                   p.val.adj = pval.adj,
                   p_val = pval) %>%
            mutate(log_pval_adj = -1 * log10(p.val.adj))
        }
        }
      })
  
    
  # DE T test tab
  ## Select which one is the control group
    select_control <- reactive({
      req(input$meta)
      
      meta() %>% 
        select(Condition) %>% 
        pull() %>% 
        as.factor() %>% 
        levels()
      
    })
    
    observeEvent(input$meta,{
      updatePickerInput(session, 
                        "control_trt", "Which one is the control group?", choices = select_control(),
                        selected = NA)
    })
    
  ## Create res dataframe
    
    dat_ttest <- reactive({
      
      req(input$control_trt)
      
      if (input$df_use == TRUE){
        
        validate(
          need(input$batch_rm_vars1, 'There is not a batch-corrected data matrix; choose the non-corrected one.'))
        
        ttest_prot(df_norm_wb(), df_norm(), meta(), name.con = input$control_trt)
        
      } else {
        
        ttest_prot(df_norm(), df_norm(), meta(), name.con = input$control_trt)
      }
      
    })
    
    
  # DE ANOVA tab
    
    select_control_aov <- reactive({
      req(input$meta)
      
      meta() %>% 
        select(Condition) %>% 
        pull() %>% 
        as.factor() %>% 
        levels()
      
    })
    
    observeEvent(input$meta,{
      updatePickerInput(session, 
                        "control_aov", "Choose a group to launch the analysis", choices = select_control_aov(),
                        selected = NA)
    })
    
    dat_aov <- eventReactive(input$control_aov, {
      
      req(input$control_aov)
      
      if (input$df_use_aov == TRUE){
        
        validate(
          need(input$batch_rm_vars1, 'There is not a batch-corrected data matrix; choose the non-corrected one.'))
        
        anova_prot(df_norm_wb(), meta(), control.name = input$control_aov)
        
      } else {
        
        anova_prot(df_norm(), meta(), control.name = input$control_aov)
      }
      
    })
    
    
    
  # Biomarker selection tab
    
    dat_fc <- reactive({
      
      req(input$control_trt2)
      
      if (input$res_use == "limma"){
        
        validate(
          need(input$control_trt_limma, 'No result has yet been created with limma'))
        
        dat_limma()
        
      } else if (input$res_use == "T test"){
        
        validate(
          need(input$control_trt, 'No result has yet been created with a T test'))
        
        dat_ttest()
        
      } else {
        
        validate(
          need(input$control_aov, 'No result has yet been created with ANOVA'))
        
        dat_aov()[[input$cities]]
        
      }
      
    })
    
    
    select_control2 <- reactive({
      req(input$meta)
      
      meta() %>% 
        select(Condition) %>% 
        pull() %>% 
        as.factor() %>% 
        levels()
      
    })
    
    observeEvent(input$meta,{
      updatePickerInput(session, 
                        "control_trt2", "Which one is the control group?", choices = select_control2(),
                        selected = NA)
    })
    
    de_prot_df <- reactive({
      
      req(input$control_trt2)
      
      validate(
        need(input$res_use, 'Select a result matrix.'))
      
      if (input$df_use2 == TRUE){
      
        validate(
          need(input$batch_rm_vars1, 'There is not a batch-corrected data matrix; choose the non-corrected one.'))
        
      prep_data_bio(df1 = df_norm_wb(), df2 = dat_fc(), meta = meta(), name.control = input$control_trt2)
        
      } else {

        prep_data_bio(df_norm(), dat_fc(), meta(), name.control = input$control_trt2)
      }
      
    })
    
    datos_train <- reactive({
      
      req(input$part)
      
      create_train(de_prot_df(), group = "Condition", p = input$part)
    })
    
    datos_test <- reactive({
      
      req(input$meta)
      
      create_test(de_prot_df(), group = "Condition", p = input$part)
    })
    
    ## Convert to matrix
    x_train <- reactive({
      req(input$meta)
      
      model.matrix(Condition ~ ., data = datos_train())[, -1]})
    y_train <- reactive({
      
      req(input$meta)
      datos_train()$Condition})

    x_test <- reactive({
      
      req(input$meta)
      model.matrix(Condition ~ ., data = datos_test())[, -1]})
    y_test <- reactive({
      req(input$meta)
      
      datos_test()$Condition})
    
    ## Run Lasso
    lasso <- reactive({
      
      req(input$control_trt2)
      
      all_lasso(x_train(), y_train(), x_test(), y_test(), meta(), name.control = input$control_trt2)
      
    })
    
    
    
    
  ########## OUTPUTS
    
  # Upload files tab
    output$file_table <- renderDataTable(
        
        DT::datatable(prot_mbr(), options = list(autoWidth = FALSE, scrollX = TRUE))
    )
    
    output$meta_table <- renderDataTable(
        DT::datatable(meta(), options = list(autoWidth = FALSE, scrollX = TRUE))
    )
    
    output$clean_raw <- renderDataTable(
        DT::datatable(dataf(), options = list(autoWidth = FALSE, scrollX = TRUE))
    )
    
    output$dl_raw <- downloadHandler(
      
      filename = function() {
        paste0("raw", ".csv")
      },
      content = function(file) {
        vroom::vroom_write(dataf(), file)
      }
    )
    
  # Missing values tab
    output$ms_table <- renderDataTable(
        DT::datatable(missing.values(), options = list(autoWidth = FALSE, scrollX = TRUE))
    )
    
    output$ms_pcts <- renderPlot(
        mis_vals_pcts(dataf(), meta())
    )
    
    unitSelect <- reactive({
      switch(input$ext1, 
             "png" = ".png", 
             "jpeg" = ".jpeg", 
             "svg" = ".svg")
    })
    
    output$dl_pct <- downloadHandler(
      filename = function() { paste('missing.vals.1', unitSelect(), sep='') },
      content = function(file) {
        ggsave(file, plot = mis_vals_pcts(dataf(), meta()), device = input$ext1)
      }
    )
    
    
    output$ms_heatmap <- renderPlot(
        na_heatmap(dataf(), meta(), vars = input$na_vars)
    )
    
  # QC tab
    
    output$n_prots <- renderPlot(
        plot_n_prots_sample(dataf(), meta(), fill = input$n_prots_vars)
    )
    
    output$qc_corr <- renderPlot(
        cor_plot(dataf())
    )
    
    output$qc_boxplot <- renderPlot(
      boxplot_ab(dataf(), meta(), fill = input$qc_fill_vars, reorder = input$reorder)
    )
    
    output$qc_den <- renderPlot(
        den_plot(dataf(), meta(), color = input$den_vars)
    )
    
  # Filtering tab
    
    output$filt_table <- renderDataTable(
        DT::datatable(filt(), options = list(autoWidth = FALSE, scrollX = TRUE))
    )
    
    output$filt_n_prots <- renderPlot(
        plot_n_prots_sample(filt(), meta())
    )
    
    output$filt_heatmap <- renderPlot(
        na_heatmap(filt(), meta(), vars = input$filt_vars)
    )
    
  # Imputation tab
    
    output$imp_table <- renderDataTable(
      DT::datatable(df.imp(), options = list(autoWidth = FALSE, scrollX = TRUE))
    )
    
    output$dl_imp <- downloadHandler(
      
      filename = function() {
        paste0("imputed.matrix", ".csv")
      },
      content = function(file) {
        vroom::vroom_write(df.imp(), file, delim = ",")
      }
    )
    
    output$imp_den <- renderPlot(
      compare_den(dataf(), df.imp(), meta(), color = input$imp_den_vars)
    )
    
    output$imp_boxplot <- renderPlot(
      boxplot_ab(df.imp(), meta(), fill = "Condition")
    )
    
  # Normalization tab
    
    output$norm_table <- renderDataTable(
      DT::datatable(df_norm(), options = list(autoWidth = FALSE, scrollX = TRUE))
    )
    
    output$norm_boxplot <- renderPlot(
      norm_boxplot(df.imp(), df_norm(), meta())
    )
    
    output$download3 <- downloadHandler(
      filename = function() {
        paste0("df_norm", ".csv")
      },
      content = function(file) {
        vroom::vroom_write(df_norm(), file, delim = ",")
      }
    )
    
  # Check batch effect tab
    
    output$batch_pca <- renderPlot({
      
      req(input$batch_pca_vars)
      
      pca_prot(df_norm(), meta(), var_p = input$batch_pca_vars)
    })
    
    output$pvca <- renderPlot({
        req(input$meta)
        pvca_plot(df_norm(), meta())
        
    })
    
    output$batch_heatmap <- renderPlot({
      
      req(input$batch_vars)
      
      heatmap_prot(df_norm(), meta(), vars = input$batch_vars)
    })
    
  # Remove batch effect tab
    
    output$rm_batch_table <- renderDataTable(
      DT::datatable(df_norm_wb(), options = list(autoWidth = FALSE, scrollX = TRUE))
    )
    
    output$download1 <- downloadHandler(
      filename = function() {
        paste0("df_norm_wb", ".csv")
      },
      content = function(file) {
        vroom::vroom_write(df_norm_wb(), file)
      }
    )
    
    output$batch_pca2 <- renderPlot({
      
      req(input$batch_rm_vars1)
      
      pca_prot(df_norm_wb(), meta(), var_p = input$batch_pca_vars2)
    })
    
    output$batch_heatmap2 <- renderPlot({
      
      req(input$batch_rm_vars1)
      
      heatmap_prot(df_norm_wb(), meta(), vars = input$batch_vars2)
    })
    
  # DE analysis
  ## limma
    output$res_limma <- renderDataTable(
      DT::datatable(dat_limma(), options = list(autoWidth = FALSE, scrollX = TRUE))
    )
    
    output$de_limma <- renderDataTable(
      table_des(dat_limma())
    )
    
    output$download_limma <- downloadHandler(
      filename = function() {
        paste0("dat_limma", ".csv")
      },
      content = function(file) {
        vroom::vroom_write(dat_limma(), file, delim = ",")
      }
    )
    
    output$figs_limma <- renderPlot(
      plot_res(dat_limma())
    )
    
    output$volcano_limma <- renderPlot(
      volcano_ttest(dat_limma())
    )
    
    
  ## t test
    output$res_ttest <- renderDataTable(
      DT::datatable(dat_ttest(), options = list(autoWidth = FALSE, scrollX = TRUE))
    )
    
    output$de_ttest <- renderDataTable(
      table_des(dat_ttest())
    )
    
    output$download2 <- downloadHandler(
      filename = function() {
        paste0("dat_ttest", ".csv")
      },
      content = function(file) {
        vroom::vroom_write(dat_ttest(), file)
      }
    )
    
    output$figs_res <- renderPlot(
      plot_res(dat_ttest())
    )

    output$volcano_ttest <- renderPlot(
      volcano_ttest(dat_ttest())
    )
    
    
  ## ANOVA
    output$cityControls <- renderUI({
      cities <- names(dat_aov())
      pickerInput(inputId = "cities", label = "Choose a comparison", cities, 
                  options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), 
                  multiple = FALSE, width = "45%")
    })
    
    
    output$res_aov <- renderDataTable(
      DT::datatable(dat_aov()[[input$cities]], options = list(autoWidth = FALSE, scrollX = TRUE))
    )
    
    
    output$download_aov <- downloadHandler(
      filename = function() {
        paste0("dat_aov_", input$cities, ".csv")
      },
      content = function(file) {
        vroom::vroom_write(dat_aov()[[input$cities]], file)
      }
    )
    
    output$de_aov <- renderDataTable(
      table_des(dat_aov()[[input$cities]])
    )
    
    output$figs_aov <- renderPlot({
      req(input$cities)
      
      plot_res_aov(dat_aov()[[input$cities]])
    })
    
    output$volcano_aov <- renderPlot(
      volcano_ttest(dat_aov()[[input$cities]])
    )
    
    
    # Biomarker selection
    ## Create data partition
    
    output$bio_dt <- renderPlot({
      
      req(input$control_trt2)
      
      data_part_plot(datos_train(), datos_test(), meta(), name.control = input$control_trt2)
    })
    
    output$dat_train <- renderDataTable(
      DT::datatable(datos_train(), options = list(autoWidth = FALSE, scrollX = TRUE))
    )
    
    output$dat_test <- renderDataTable(
      DT::datatable(datos_test(), options = list(autoWidth = FALSE, scrollX = TRUE))
    )
    
    ## Lasso results
    
    output$lambda <- renderPlot(
      lasso()$plot_first_model
    )
    
    output$error <- renderPlot(
      lasso()$fig_cv_error
    )
    
    output$coefs_plot <- renderPlot(
      lasso()$fig_coefs
    )
    
    output$coefs_dt <- renderDataTable(
      DT::datatable(lasso()$df_lasso, options = list(autoWidth = FALSE, scrollX = TRUE))
    )
    
    ## Evaluation
    output$conf_m_train <- renderPlot(
      lasso()$conf_matrix_train
    )
    
    output$conf_m_test <- renderPlot(
      lasso()$conf_matrix_test
    )
    
    output$auc <- renderPlot(
      
      if (input$auc_choice == "As a whole"){
        
        lasso()$roc_auc
        
      } else {
        
        plot_auc_all(lasso()$prot_names, datos_train(), x_test(), y_test())
        
      }
      
    )
    
    
    output$pca_lasso <- renderPlot(
      
      pca_lasso(de_prot_df(), meta(), lasso()$prot_names)
    )
    
    output$heatmap_lasso <- renderPlot(
      
      heatmap_lasso(de_prot_df(), meta(), lasso()$prot_names)
    )
    
}



shinyApp(ui, server)