# Author: Etienne CAMENEN
# Date: 2021
# Contact: etienne.camenen@gmail.com
# Key-words: omics, RGCCA, multi-block
# EDAM operation: analysis, correlation, visualisation
#
# Abstract: Performs multi-variate analysis (PCA, CCA, PLS, R/SGCCA, etc.)
# and produces textual and graphical outputs (e.g. variables and samples
# plots).
#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`.
#' @noRd

rm(list = ls())
options(shiny.maxRequestSize = 30 * 1024 ^ 2)

BSPLUS <- R.Version()$minor >= 3
setInfo <- function(., text) {
    shinyInput_label_embed(
        icon("question") %>%
            bs_embed_tooltip(title = text))
}

load_libraries <- function(librairies) {
    for (l in librairies) {
        if (!(l %in% installed.packages()[, "Package"]))
            utils::install.packages(l, repos = "cran.us.r-project.org")
    suppressPackageStartupMessages(
        library(
            l,
            character.only = TRUE,
            warn.conflicts = FALSE,
            quietly = TRUE
        ))
    }
}

load_libraries(c(
    "ggplot2",
    "scales",
    "igraph",
    "plotly",
    "visNetwork",
    "shiny",
    "shinyjs",
    "MASS",
    "DT",
    "devtools"
))

if (!("RGCCA" %in% installed.packages()[, "Package"]) ||
    as.double(paste(unlist(packageVersion("RGCCA"))[seq(2)], collapse = ".")) < 3.0) {
    devtools::install_github("rgcca-factory/RGCCA", ref = "3.0.0")
}

all_funcs <<- unclass(lsf.str(envir = asNamespace("RGCCA"), all = TRUE))
for (i in all_funcs)
    eval(parse(text = paste0(i, "<<-RGCCA:::", i)))

if (BSPLUS) {
    if (!("bsplus" %in% installed.packages()[, "Package"]))
        devtools::install_github("ijlyttle/bsplus", upgrade = "never")
    library("bsplus", warn.conflicts = FALSE, quietly = TRUE)
}

customDownloadbutton <- function(outputId, label = "Download"){
    tags$a(id = outputId, class = "btn btn-default shiny-download-link", href = "", 
           target = "_blank", download = NA, icon("camera"), label)
}

ui <- fluidPage(
    titlePanel("R/SGCCA - The Shiny graphical interface"),
    tags$div(
        tags$p(
                tags$a(
                    href = "mailto:etienne.camenen@gmail.com",
                    "Etienne CAMENEN,"
                ),
            "Ivan MOSZER, Arthur TENENHAUS"
        ),
        tags$i("Multi-block data analysis concerns the analysis of several sets of variables (blocks) observed on the same group of samples. The main aims of the RGCCA package are: to study the relationships between blocks and to identify subsets of variables of each block which are active in their relationships with the other blocks."),
        tags$br(), tags$br()
    ),
    tags$a(href = "https://github.com/ecamenen/rgcca-ui/blob/master/inst/shiny/tutorialShiny.md#description", "Go to the tutorial"),
    tags$strong("|"),
    tags$a(href = "https://www.youtube.com/watch?v=QCkEBsoP-tc", "Watch a demo", target = "_blank"),
    tags$br(), tags$br(),
    tags$style(".fa-camera {color:#c7c7c7}"),
    tags$style(".fa-camera:hover {color:#7c7c7c}"),
    tags$style("#connection_save, #ave_save, #bootstrap_t_save, #permutation_t_save {border-color:white; left: 0%}"),
    tags$style("#connection_save:hover, #ave_save:hover, #bootstrap_t_save:hover, #permutation_t_save:hover {background-color:white}"),
    tags$style("#connection_save:focus, #ave_save:focus, #bootstrap_t_save:focus, #permutation_t_save:focus {outline:none; background-color:white}"),
    tags$style("#connection_save:active, #ave_save:active, #bootstrap_t_save:active, #permutation_t_save:active {box-shadow:none}"),
    tags$style(".js-plotly-plot .plotly .modebar {left: 0%}"),
    # tags$style(
    #     HTML(".shiny-notification {
    #         position:fixed;
    #         top: calc(90%);
    #         left: calc(50%);
    #         }")
    #     ),
    useShinyjs(),
    sidebarLayout(sidebarPanel(
        tabsetPanel(
            id = "tabset",
            tabPanel(
                "Data",
                uiOutput("file_custom"),
                uiOutput("sep_custom"),
                checkboxInput(
                    inputId = "header",
                    label = "Consider first row as header",
                    value = TRUE
                ),
                actionButton(
                    inputId = "run_data",
                    label = "Use a default dataset")
            ),

            # Analysis parameters

            tabPanel(
                "RGCCA",
                uiOutput("analysis_type_custom"),
                uiOutput("scale_custom"),
                radioButtons(
                    "init",
                    label = "Mode of initialization",
                    choices = c(SVD = "svd",
                        Random = "random"),
                    selected = "svd"
                ),
                uiOutput("superblock_custom"),
                checkboxInput(
                    inputId = "supervised",
                    label = "Supervised analysis",
                    value = FALSE
                ),
                conditionalPanel(
                    condition = "input.supervised || input.analysis_type == 'RA'",
                uiOutput("blocks_names_response")),
                uiOutput("connection_custom"),
                checkboxInput(
                    inputId = "each_ncomp",
                    label = "Number of components for each block",
                    value = FALSE
                ),
                uiOutput("nb_compcustom"),
                uiOutput("tau_opt_custom"),
                uiOutput("each_tau_custom"),
                uiOutput("tau_custom"),
                uiOutput("tune_type_custom"),
                uiOutput("val_custom"),
                sliderInput(
                    inputId = "ncv",
                    label = "Number of cross-validation",
                    min = 1,
                    max = 100,
                    value = 1,
                    step = 1
                ),
                sliderInput(
                    inputId = "kfold",
                    label = "Number of folds",
                    min = 2,
                    max = 10,
                    value = 5,
                    step = 1
                ),
                actionButton(
                    inputId = "run_crossval",
                    label = "Run cross-validation"),
                uiOutput("nperm_custom"),
                actionButton(inputId = "run_perm",
                        label = "Run permutation"),
                # sliderInput(
                #     inputId = "power",
                #     label = "Power of the factorial",
                #     min = 2,
                #     max = 6,
                #     value = 2,
                #     step = 1
                # ),
                uiOutput("scheme_custom"),
                actionButton(
                    inputId = "run_analysis",
                    label = "Run analysis"),
                uiOutput("nboot_custom"),
                actionButton(inputId = "run_boot",
                        label = "Run bootstrap"),
                actionButton(
                    inputId = "run_crossval_single",
                    label = "Evaluate the model")
            ),

            # Graphical parameters

            tabPanel(
                "Graphic",
                radioButtons(
                    "format",
                    label = "Output image format",
                    choices = c(
                        `jpeg` = "jpeg",
                        `png` = "png"
                        #`svg` = "svg"
                        # `tiff` = "tiff",
                        # `pdf` = "pdf"
                    ),
                    selected = "png"
                ),
                checkboxInput(
                    inputId = "text",
                    label = "Display names",
                    value = TRUE
                ),
                uiOutput("blocks_names_custom_x"),
                uiOutput("blocks_names_custom_y"),
                uiOutput("compx_custom"),
                uiOutput("compy_custom"),
                uiOutput("nb_mark_custom"),
                uiOutput("response_custom"),
                checkboxInput(
                    inputId = "show_crossval",
                    label = "Display cross-validation",
                    value = TRUE
                ),
                radioButtons(
                    "indexes",
                    label = "Type of indexes",
                    choices = c(
                        Correlation = "loadings",
                        Weights = "weight")
                ),
                uiOutput("b_x_custom"),
                uiOutput("b_y_custom"),
                customDownloadbutton(outputId = "save_all", label = "Save all")
            )
        )
    ),

    mainPanel(
        tabsetPanel(
            type = "tabs",
            id = "navbar",
            tabPanel(
                "Connection",
                customDownloadbutton("connection_save", ""),
                visNetworkOutput("connectionPlot")
            ),
            tabPanel(
                "AVE",
                customDownloadbutton("ave_save", ""),
                plotOutput("AVEPlot")
            ),
            tabPanel(
                "Samples",
                plotlyOutput("samplesPlot", height = 500),
                customDownloadbutton("samples_save", "Save")
            ),
            tabPanel(
                "Corcircle",
                plotlyOutput("corcirclePlot", height = 500),
                customDownloadbutton("corcircle_save", "Save")
            ),
            tabPanel(
                "Fingerprint",
                plotlyOutput("fingerprintPlot", height = 700),
                customDownloadbutton("fingerprint_save", "Save")
            ),
            tabPanel(
                "Bootstrap",
                plotlyOutput("bootstrapPlot", height = 700),
                customDownloadbutton("bootstrap_save", "Save")
            ),
            tabPanel(
                "Bootstrap Summary",
                customDownloadbutton("bootstrap_t_save", ""),
                DT::dataTableOutput("bootstrapTable")
            ),
            tabPanel(
                "Permutation",
                plotlyOutput("permutationPlot", height = 700)
            ),
            tabPanel(
                "Permutation Summary",
                customDownloadbutton("permutation_t_save", ""),
                dataTableOutput("permutationTable")
            ),
            tabPanel(
                "Cross-validation",
                plotlyOutput("cvPlot", height = 700)
            )
        )
    ))
)
