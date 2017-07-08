tabPanel("Network", 
         fluidRow(column(4, 
                         wellPanel(
                           tags$textarea(id="NAME_NET", rows=5, cols=10, "Drosha\nDgcr8\nDicer1\nAgo1\nAgo2"),
                           helpText("Note: Input should be Official Gene Symbol, e.g., Dicer1. Multiple names should be seperated with ENTER"),
                           submitButton("Submit"),
                           helpText("Following Genes cannot be found; It is either not expressed or the gene name cannot be matched"),
                           textOutput("MissGene_NET")
                         )),
                  column(8,
                         tabsetPanel(
                           tabPanel("Network", visNetworkOutput("network"))
                           
                         ))
         ))