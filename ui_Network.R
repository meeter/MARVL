tabPanel("Network", 
         fluidRow(column(4, 
                         wellPanel(
                           tags$textarea(id="NAME_NET", rows=5, cols=10, "Ivns1abp\nSirt6"),
                           checkboxInput("order", "Include Down-stream", FALSE),
                           helpText("Note: click this checkbox to include also down-stream targets of the input; otherwise only up-stream
                                    regulators are displayed"),
                           submitButton("Submit"),
                           helpText("Following Genes cannot be found; It is either not expressed or the gene name cannot be matched"),
                           textOutput("MissGene_NET")
                         )),
                  column(8,
                         tabsetPanel(
                           tabPanel("Network", visNetworkOutput("network"))
                           
                         ))
         ))