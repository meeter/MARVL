library(ggplot2)
library(ggrepel)
library(gplots)
library(reshape2)
library(plotly)
load("WT.RData")

source("server_GeneralFun.R", local=TRUE)
source("server_Gene.R", local=TRUE)
source("server_miRNA.R", local=TRUE)


shinyServer(function(input, output) {
#################################################################
#####For Gene Level
#################################################################
  #####Get Missing Genes
  output$MissGene <- renderText({
    getMissGenes(input, GeneTE_1.tpm)
  })  
  #####Get ScatterPlot
  output$E14_Exp <- renderPlot({
    ScatterplotGene(input, GeneTE_1.tpm)
  })
  #####Get DataTable
  output$Data <- renderDataTable({
    TableGene(input, GeneTE_1.tpm)
  }) 
  #####Get Heatmap
  output$Heatmap <- renderPlotly({
    HeatmapGene(input, GeneTE_1.tpm)
  })
  #####Get Gene2miR, dotsplot and select function
 output$Dotsplot_Gene <- renderPlotly({
    tmp <- Gene2miR(input, GeneTE_1.tpm, WT_MIR, miR_Gene.MW.sel)
    Gene2miR_plot(tmp)
  })
  output$brush_Gene <- renderPrint({
    tmp <- Gene2miR(input, GeneTE_1.tpm, WT_MIR, miR_Gene.MW.sel)
    d <- event_data("plotly_selected")
    if (is.null(d)) "Selected events appear here (double-click to clear); Negative values indicate negative correlation" 
      else {
        tmp[d[["pointNumber"]]+1, c("ID","gene_name","spearman","spearman.p", "WT_1_Gene", "WT_2_Gene", "WT_1_miR","WT_2_miR")]
      }
  })
  
  
#################################################################
#####For miRNA Level
#################################################################  
  #####Get Missing Genes
  output$MissMIR <- renderText({
    getMissMIR(input, WT_MIR)
  })  
  #####Get ScatterPlot
  output$E14_Exp_MIR <- renderPlot({
    ScatterplotMIR(input, WT_MIR)
  })
  #####Get DataTable
  output$Data_MIR <- renderDataTable({
    TableMIR(input, WT_MIR)
  }) 
  #####Get Heatmap
  output$Heatmap_MIR <- renderPlotly({
    HeatmapMIR(input, WT_MIR)
  })
  #####Get Barplot for RIP-seq
  output$RIP_MIR <- renderPlot({
    BarplotMIR(input, RIP)
  })
  #####Get miR2Gene, dotsplot and select function
  output$Dotsplot_MIR <- renderPlotly({
    tmp <- miR2Gene(input, GeneTE_1.tpm, WT_MIR, miR_Gene.MW.sel)    
    miR2Gene_plot(tmp)
  })
  output$brush_MIR <- renderPrint({
    tmp <- miR2Gene(input, GeneTE_1.tpm, WT_MIR, miR_Gene.MW.sel)
    d <- event_data("plotly_selected")
    if (is.null(d)) "Selected events appear here (double-click to clear); Negative values indicate negative correlation" 
    else {
      tmp[d[["pointNumber"]]+1, c("ID","gene_name","spearman","spearman.p", "WT_1_Gene", "WT_2_Gene", "WT_1_miR","WT_2_miR")]
    }
  })
})
