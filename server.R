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
library(markdown)
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
  #####Intersection -log(p-value) of fisher test by intersecting input gene list with our 6 gene categories
  output$Fisher_Heatmap <- renderPlotly({
    tmp <- CateFisher(input, Dgcr8_DEG, Drosha_DEG, Dicer_DEG, Ago12_DEG, #DEGs from pairwise comparison
                      name.overlap.up, name.overlap.down, #consistently miss-regulated
                      name.MIG, name.MDG, name.Dicer, name.Ago12, #Specifically up-regulated, one vs rest
                      GeneTE_1.tpm) #GeneTE_1.tpm to provide expressed gene names
    CateFisher.plot(tmp)
    })
  output$brush_Fisher <- renderPrint({
    tmp <- CateFisher(input, Dgcr8_DEG, Drosha_DEG, Dicer_DEG, Ago12_DEG,
                      name.overlap.up, name.overlap.down, 
                      name.MIG, name.MDG, name.Dicer, name.Ago12, 
                      GeneTE_1.tpm)
    d.Fisher <- event_data("plotly_click")
    if (is.null(d.Fisher)) "Single-click the bar to display overlapped genes (double-click to clear)" 
    else {
      #d.Fisher
      print(data.frame(Overlapped_Gene = 
                unlist(strsplit(as.character(tmp[d.Fisher[["x"]], "Overlapped_Gene"]), split=";"))
                        ), row.names = FALSE)
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
  #####Biplot
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
    ggplotly(p.biplot, source="Biplot", height=400, width=800) %>% 
             layout(dragmode = "select")
            
  })
  output$brush_Biplot <- renderPrint({
    d <- event_data("plotly_hover", source="Biplot")
    tmp <- prcomp(all.res.mds[, 2:7], scale=T, center=F)
    if (is.null(d)) "Selected miRNA appear here (double-click to clear);" 
    else {
      WT_MIR[match(row.names(tmp$x)[grep(d[["x"]], tmp$x[, "PC1"])], gsub("mmu-", "", WT_MIR$ID)), 1:11]
      #tmp$x[intersect(grep(d[["xvar"]], tmp$x[, "PC1"]), grep(d[["yvar"]], tmp$x[, "PC1"])), c("PC1", "PC2")]
      #all.res.mds[d[["pointNumber"]]+1, c("E14_P25","E14_P33","Dgcr8_Log2FC","Drosha_Log2FC","Dicer_Log2FC","Ago12_Log2FC")]
    }
  })
})
