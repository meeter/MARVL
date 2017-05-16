library(ggplot2)
library(ggrepel)
library(gplots)
library(reshape2)
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
  output$Heatmap <- renderImage({
    HeatmapGene(input, GeneTE_1.tpm)
  }, deleteFile = TRUE)
  
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
  output$Heatmap_MIR <- renderImage({
    HeatmapMIR(input, WT_MIR)
  }, deleteFile = TRUE)
  #####Get Barplot for RIP-seq
  output$RIP_MIR <- renderPlot({
    BarplotMIR(input, RIP)
  }) 
})
