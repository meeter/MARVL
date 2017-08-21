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
install.packages(c("shiny", "plotly", "ggplot2", "ggrepel", "reshape2", "heatmaply", "gplots", 
		   "visNetwork", "htmlwidgets", "shinythemes","markdown","plyr","scales"))
```

You may now run the shiny app with just one command in R:

```
shiny::runGitHub("MARVL", "meeter")
```

<a name="Gene"></a> 

## Gene Level Analysis 

1. __Data Table__: Exploring the abudance for interested genes across different RNAi-mutants. The abundance is indicated as Log2-TPM (transcripts per million reads).
2. __Scatter Plot__: Abundance of interested genes in the two wild-type replicates. The abundance is normalized as above.
3. __Heatmap__: An interactive heatmap for interested genes. The abundance is indicated as above.
4. __Gene <-> miRNA__: An interactive scatterplot of miRNAs with strong correlation to the interested genes. To see the detailed information of Gene <-> miRNA, select the dots, and the corresponding information will be displayed below, including the abundance of genes (Log2-TPM) and miRNAs (Log2-normalized by RLE [1] using the median of snRNAs, snoRNAs and tRNAs). The correlation is indicated as spearman correlation, with only those p-value < 0.05 and abs(spearman r) > 0.70. Negative values indicate negative correlation between miRNA and Gene, and vice versa.

<a name="miRNA"></a> 

## miRNA Level Analysis 

1. __Data Table__: Exploring the abudance for interested miRNA across different RNAi-mutants. The abundance is indicated as Log2-normalized counts by RLE using the median of snRNAs, snoRNAs and tRNAs[1].
2. __Scatter Plot__: Abundance of interested miRNAs in the two wild-type replicates. The abundance is normalized as above.
3. __Heatmap__: An interactive heatmap for interested miRNAs. The abundance is indicated as above. 
4. __RIP-seq__: Exploring the enrichment level of interested miRNA in AGO1 and AGO2 IP experiment, using corresponding knockout as control. The enrichment level is indicated as Log2 fold change between IP and the knockout library.
5. __miRNA <-> Gene__: An interactive scatterplot of miRNAs with strong correlation to the interested genes. To see the detailed information of Gene <-> miRNA, select the dots, and the corresponding information will be displayed below, including the abundance of genes (Log2-TPM) and miRNAs (normalized as above). The correlation is indicated as spearman correlation, with only those p-value < 0.05 and abs(spearman r) > 0.70. Negative values indicate negative correlation between miRNA and Gene, and vice versa.
6. __Biplot__: An interactive MDS biplot of all expressed miRNAs, where dots stand for miRNAs and axes for features, including the abundances in WT, and Log2-FC between each knockout vs. WT. Reported miRTrons were colored and labeled in blue; predicted MI miRTrons in green; others in red. 95% confidence ellipses were drawn for each miRNA category. Detailed information will be displayed below upon cursor hovering, including the miRNA abundance in all samples.

<a name="Network"></a> 

## Combinatory Regulatory Network 

1. __Network__: A subnetwork of how interested genes are regulated across all samples. Potential regulators include miRNAs, transcription factors, and other epigentic regulators like histone modifiers and DNA methylation. The color of the vertex indicated the averaged FC between each KO versus WT, blue for down-regulation, and red for up-regulation. 




