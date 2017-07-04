tabPanel("miRNA", 
         fluidRow(column(4, 
                         wellPanel(
                           tags$textarea(id="NAME_MIR", rows=5, cols=10, "miR-92a-3p\nmiR-1983\nmiR-877-3p"),
                           helpText("Note: Input should be miRNA name, e.g., miR-92a-3p. Multiple names should be seperated with ENTER"),
                           checkboxInput(inputId = "Dgcr8_MIR", label = strong("Include Dgcr8 KO"), value = TRUE),
                           checkboxInput(inputId = "Drosha_MIR", label = strong("Include Drosha KO"), value = TRUE),
                           checkboxInput(inputId = "Dicer_MIR", label = strong("Include Dicer KO"), value = TRUE),
                           checkboxInput(inputId = "Ago12_MIR", label = strong("Include Ago12"), value = TRUE),
                           checkboxInput(inputId = "WT_MIR", label = strong("Include WT"), value = TRUE),
                           submitButton("Submit"),
                           helpText("Following miRNAs cannot be found; It is either not expressed or the gene name cannot be matched"),
                           textOutput("MissMIR")
                         )),
                  column(8,
                         tabsetPanel(
                           tabPanel("Data Table", dataTableOutput("Data_MIR")),
                           tabPanel('Scatter Plot', plotOutput("E14_Exp_MIR",height="600px",width="640px")),
                           tabPanel('Heatmap', plotlyOutput("Heatmap_MIR", width="100%")),
                           tabPanel("RIP-seq", plotOutput("RIP_MIR", height="600")),
                           tabPanel("miRNA <-> Gene", dataTableOutput("Data_Dots"),
                                    plotlyOutput("Dotsplot_MIR"),
                                    verbatimTextOutput("brush_MIR")),
                           tabPanel("Biplot", plotlyOutput("Dotsplot_Biplot"), verbatimTextOutput("brush_Biplot"))
                         ))
         ))