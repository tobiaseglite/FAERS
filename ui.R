library(shiny)
library(shinydashboard)

data <- read.table(paste0("data/results_memory_sorted.txt"),
        stringsAsFactors = F,
        sep = "\t",
        h = T)

data <- data[!is.na(data$results.ROR),]

all_drugs <- unique(data$results.drug)

tags$style(".skin-blue .sidebar a { color: #444; }")

dashboardPage(
    
    # Application title
    dashboardHeader(title = "FAERS results"),
    
    # Sidebar with a slider input for number of bins 
    dashboardSidebar(
            selectInput(inputId = "result_name",
                label = "Type of adverse reaction",
                choices = c("anaemia",
                        "anorexia",
                	"anxiety",
                        "arthralgia",
                	"attention",
                	"dementia",
                	"depression",
                	"emotional",
                	"mania",
                	"memory",
                        "nasopharyngitis",
                	"panic_attack",
                	"paranoia",
                        "pneumonia",
                	"psychotic",
                	"suicide"),
                selected = "memory",
                multiple = FALSE,
                selectize = TRUE, width = NULL, size = NULL),
            selectInput(inputId = "type",
                label = "Results across",
                choices = c("all","indication","age","female","male"),
                selected = "all",
                multiple = FALSE,
                selectize = TRUE, width = NULL, size = NULL),
            selectInput(inputId = "drug",
                label = "Active ingredient",
                choices = all_drugs,
                selected = all_drugs[1],
                multiple = FALSE,
                selectize = TRUE, width = NULL, size = NULL),
            downloadButton("downloadData","Download results"),
            selectInput(inputId = "plotReference",
                label = "Illustrations (p/OR) order",
                choices = c("Indications",
                    "Order p all",
                    "Order p indication",
                    "Order p age",
                    "Order p female",
                    "Order p male"),
                selected = "Indications",
                multiple = FALSE,
                selectize = TRUE, width = NULL, size = NULL),
            sliderInput("range", "Range x-axis (ingredient):",
                  min = 1, max = nrow(data),
                  value = c(1,nrow(data))),
            sliderInput("rangeY", "Range y-axis (-log10(p)):",
                  min = -700, max = 700,
                  value = c(-300,300)),
            tags$style(".skin-blue .sidebar a { color: #444; margin-left: 14px;")
        ),
    dashboardBody(
        fluidRow(
            plotOutput("ORPlot",
                height = "200px",
                click = "ORPlot_click",
                brush = brushOpts(id = "ORPlot_brush",
                    resetOnNew = T),
                dblclick = "ORPlot_dblclick"),
            plotOutput("pvalPlot",
                height = "200px",
                click = "pvalPlot_click",
                brush = brushOpts(id = "pvalPlot_brush",
                    resetOnNew = T),
                dblclick = "pvalPlot_dblclick")
        ),
        br(),
        fluidRow(
            box(htmlOutput("textResults"), width = 12)
        ),
        fluidRow(
            column(6,
                tabBox(
                    title = "",
                    id = "tabset1", height = "300px",
                    tabPanel("Cross-table",
                        tableOutput("crossTable")),
                    tabPanel("Fisher's Exact Test",
                        tableOutput("fishersTable")),
                    tabPanel("Logistic reg.",
                        tableOutput("logregResults"),
                        downloadButton("downloadDataLog","Download results LogReg"),
                        tags$style(".skin-blue .sidebar a { color: #444; margin-left: 14px;"),
                        background = "black"),
                width = 12)
            ),
            column(6,box(title = "log(Odds ratios)",
                width = "500px",
                height = "300px",
                plotOutput("drugPlot",
                    height = "200px",
                    width = "400px")))
        ),
        fluidRow(
            tabBox(
                id = "tabset2",
                title = "DrugBank",
                tabPanel("Drug targets",
                    tableOutput("infoDrug")),
                tabPanel("GO annotations targets",
                    dataTableOutput("GOannot")),
                tabPanel("Description",
                    htmlOutput("infoDrugDescription"), width = 12),
                tabPanel("Pharmacodynamics",
                    htmlOutput("infoDrugPharmacodynamic"), width = 12),
                tabPanel("Indication",
                    htmlOutput("infoDrugIndication"), width = 12),
            width = 12)
        ),
        fluidRow(
            box(title = "The drugnames assigned to the active ingredient",
                dataTableOutput("drugActiveIngredients"),
                width = 12
            )
        ),
        fluidRow(
            box(title = "The 10 most frequent indications for the active ingredient",
                dataTableOutput("drugIndications"),
                width = 12
            )
        ),
        fluidRow(
            box(title = "The 10 most frequent active ingredients administered for the most frequent indication",
                dataTableOutput("topIndicationsDrug"),
                width = 12
            )
        )
    )
)

