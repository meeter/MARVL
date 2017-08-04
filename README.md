---
output: html_document
---

# Shiny host
MARVL is hosted on shinyapps here:
<https://cclab.shinyapps.io/MARVL/>

# Local Installation
To run this app locally on your machine, download R or RStudio and run the following commands once to set up the environment:

```
install.packages(c("shiny", "plotly", "ggplot2", "ggrepel", "reshape2", "heatmaply", "gplots", 
		   "visNetwork", "htmlwidgets", "shinythemes","markdown","plyr","scales"))
```

You may now run the shiny app with just one command in R:

```
shiny::runGitHub("MARVL", "meeter")
```

# Instructions
Instructions can be found here: <https://github.com/meeter/MARVL/blob/master/instructions/Instructions.md> 


# Licensing
This shiny code is licensed under the GPLv3.
