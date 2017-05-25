tabPanel("Help",
         fluidRow(
           column(4,wellPanel(
             h4("Instructions"),
             a("Gene", href = "#Gene"), br(),
             a("miRNA", href = "#miRNA"), br(),
             a("FAQ", href="#FAQ"), br()
           )
           ),#column
           column(8,
                  includeMarkdown("instructions/Instructions.md"))
         ))