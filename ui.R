library(shiny)

renderInputs <- function(prefix) {
   absolutePanel(top = 60, left = 0, right = 0,
    fixed = TRUE, style = "opacity: 0.9",
    div(
      style="padding: 8px; border-bottom: 1px solid #CCC; background: #FFFFEE;",	
     fluidRow(
       column(6,textInput("NAME", "Input Gene Symbol: e.g., Dgcr8", value = "Dgcr8"))
       ),
       submitButton("Submit")
         #tags$a(href="https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&position=chr15%3A12824815-12935291&hgsid=454727657_A6ynhaMAO3PirpGrkHbn1uZM0pxx", "UCSC session Link for all data", align="center")
 			 )
 								)       
}

shinyUI(fluidPage(theme="simplex.min.css",
  								tags$style(type="text/css",
    							"label {font-size: 12px;}",
    							".recalculating {opacity: 1.0;}"
),
##Application title
	tags$h2("CCLAB  RNAi-Mutant Sequencing Data",align="center"),
	fluidRow(
		tags$br(" "),
    tags$br(""),	
    tags$br(""),
    tags$br(""),
		helpText("Following Genes cannot be found; It is either not expressed or the gene name cannot be matched"),
		textOutput("MissGene"), 
    column(6, renderInputs("Gene Abundance in E14"))
    ),
    column(8, align="center", 
    plotOutput("E14_Exp",height="800px",width="800px")
    )
))


