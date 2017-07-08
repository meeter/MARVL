library(shiny)
library(plotly)

customHeaderPanel <- function(title,windowTitle=title){
  tagList(
    tags$head(
      tags$title(windowTitle),
      tags$link(rel="stylesheet", type="text/css",
                href="app.css"),
      tags$h1(a(href="http://www.ohsu.edu/xd/health/services/heart-vascular/"))
    )
  )
}

tagList(
  tags$head(
    tags$style(HTML(" .shiny-output-error-validation {color: darkred; } ")),
    tags$style(".mybuttonclass{background-color:#CD0000;} .mybuttonclass{color: #fff;} .mybuttonclass{border-color: #9E0000;}")
  ),
  titlePanel(div("MARVL: a Mirna And Rna-seq integrated VisuaLization web-service", align="center", style="color:darkred")),
  navbarPage(
    
    theme = "bootstrap.min.united.updated.css",
    title = "",
   
    ## =========================================================================== ##
    ## Gene TABS
    ## =========================================================================== ##
    source("ui_Gene.R",local=TRUE)$value,
    ## =========================================================================== ##
    ## miRNA TABS
    ## =========================================================================== ##
    source("ui_miRNA.R",local=TRUE)$value,
    ## =========================================================================== ##
    ## Network TABS
    ## =========================================================================== ##
    source("ui_Network.R",local=TRUE)$value,
    
    ## ============================================================================ ##
    ## HELP TAB
    ## ============================================================================ ##   
    source("Help.R",local=TRUE)$value,
    
    ## ==================================================================================== ##
    ## FOOTER
    ## ==================================================================================== ##              
    footer=p(hr(),p("Created by ", a("CCLAB", href="http://www.mhs.biol.ethz.ch/research/ciaudo.html"),
                                     " of ",align="center",width=4),
             p(("Institute of Molecular Health Science, ETH, Zurich"),align="center",width=4),
             p(("Copyright (C) 2017, code licensed under GPLv3"),align="center",width=4),
             p(("Code available on Github:"),a("https://github.com",href="https://github.com/"),align="center",width=4)
    ),
    ## ==================================================================================== ##
    ## end
    ## ==================================================================================== ## 
    tags$head(includeScript("google-analytics.js"))
  ) #end navbarpage
) #end taglist


	




