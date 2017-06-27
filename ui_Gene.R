tabPanel("Gene", 
         fluidRow(column(4, 
                  wellPanel(
                  tags$textarea(id="NAME", rows=5, cols=10, "Drosha\nDgcr8\nDicer1\nAgo1\nAgo2"),
                  helpText("Note: Input should be Official Gene Symbol, e.g., Dicer1. Multiple names should be seperated with ENTER"),
                  checkboxInput(inputId = "Dgcr8", label = strong("Include Dgcr8 KO"), value = TRUE),
                  checkboxInput(inputId = "Drosha", label = strong("Include Drosha KO"), value = TRUE),
                  checkboxInput(inputId = "Dicer", label = strong("Include Dicer KO"), value = TRUE),
                  checkboxInput(inputId = "Ago12", label = strong("Include Ago12"), value = TRUE),
                  checkboxInput(inputId = "WT", label = strong("Include WT"), value = TRUE),
                  submitButton("Submit"),
                  helpText("Following Genes cannot be found; It is either not expressed or the gene name cannot be matched"),
                  textOutput("MissGene")
  )),
  column(8,
    tabsetPanel(
      tabPanel("Data Table", dataTableOutput("Data")),
      tabPanel('Scatter Plot', plotOutput("E14_Exp",height="600px",width="640px")),
      tabPanel('Heatmap', plotlyOutput("Heatmap")),
      tabPanel("Gene <-> miRNA", plotlyOutput("Dotsplot_Gene"),
               verbatimTextOutput("brush_Gene"))
    ))
))