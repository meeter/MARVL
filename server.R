library(ggplot2)
library(ggrepel)
library(gplots)
library(reshape2)
library(plotly)
library(heatmaply)
library(scales)
library(igraph)
library(plyr)
#library(grid)
library(visNetwork)
library(htmlwidgets)

load("WT.RData")

source("server_GeneralFun.R", local=TRUE)
source("server_Gene.R", local=TRUE)
source("server_miRNA.R", local=TRUE)
source("server_Network.R", local=TRUE)
source("ggbiplot.R", local=TRUE)


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
        tmp[d[["pointNumber"]]+1, c("ID","gene_name","spearman","spearman.p", "WT_1_Gene", "WT_2_Gene", "WT_1_miR","WT_2_miR", "targetScan")]
      }
  })

#################################################################
#####For Network
#################################################################  
  output$MissGene_NET <- renderText({
    getMissGenes(input, GeneTE_1.tpm)
  })  
  output$network <- renderVisNetwork({
    plotNetwork(input)
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
  output$Dotsplot_Biplot <- renderPlotly({
    NAME <- GetName(input$NAME_MIR)
    tmp <- prcomp(all.res.mds[, 2:7], scale=T, center=F)
    p.biplot <- ggbiplot(tmp, obs.scale = 1, var.scale = 1, alpha=0.5, 
                         #labels.size = 2, labels = all.res.mds$ID, 
                         groups = all.res.mds$miRTron, ellipse = TRUE, ellipse.prob = 0.95, circle = F) +
      scale_color_discrete(name = '') +
      #geom_text_repel(aes(label=all.res.mds$ID, col = all.res.mds$miRTron)) + 
      xlim(-5,3)  + ylim(-2,2)
      #theme(legend.direction = 'horizontal', legend.position = 'bottom') 
    ggplotly(p.biplot, source="A", height=400, width=800) %>% 
             layout(dragmode = "select")
            
  })
  output$brush_Biplot <- renderPrint({
    d <- event_data("plotly_hover", source="A")
    tmp <- prcomp(all.res.mds[, 2:7], scale=T, center=F)
    if (is.null(d)) "Selected miRNA appear here (double-click to clear);" 
    else {
      row.names(tmp$x)[grep(d[["x"]], tmp$x[, "PC1"])]
      #tmp$x[intersect(grep(d[["xvar"]], tmp$x[, "PC1"]), grep(d[["yvar"]], tmp$x[, "PC1"])), c("PC1", "PC2")]
      #all.res.mds[d[["pointNumber"]]+1, c("E14_P25","E14_P33","Dgcr8_Log2FC","Drosha_Log2FC","Dicer_Log2FC","Ago12_Log2FC")]
    }
  })
})
