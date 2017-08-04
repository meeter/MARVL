---
output:
  html_document:
    theme: united
---

# Help

The app is hosted on the website: https://cclab.shinyapps.io/MARVL/

Code can be found on github: https://github.com/meeter/MARVL

To run this app locally on your machine, download R or RStudio and run the following commands once to set up the environment:

```
install.packages(c("shiny", "plotly", "ggplot2", "ggrepel", "reshape2", "heatmaply", "gplots", "visNetwork", "htmlwidgets", "shinythemes","markdown","plyr","scales"))
```

You may now run the shiny app with just one command in R:

```
shiny::runGitHub("MARVL", "meeter")
```

<a name="Gene"></a> 

## Gene Level Analysis 

In this level you can

1. Exploring the abudance for interested genes across different RNAi-mutants. The abundance is indicated as Log2-TPM (transcripts per million reads)
2. Identification of miRNAs with strong correlation to the interested genes. The correlation is indicated as spearman correlation, with only those p-value < 0.05 and abs(spearman r) > 0.85. Negative values indicate negative correlation between miRNA and Gene, and vice versa.

<a name="miRNA"></a> 

## miRNA Level Analysis 

1. Exploring the abudance for interested miRNA across different RNAi-mutants. The abundance is indicated as Log2-normalized counts by DESeq2 [1].
2. Exploring the enrichment level of interested miRNA in Ago1/2 IP experiment. The enrichment level is indicated as Log2 fold change between IP and Knockout library.
3. Identification of miRNAs with strong correlation to the interested genes. This part is the same as the gene level. 

<a name="FAQ"></a>

## FAQ



