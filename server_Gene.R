source("server_GeneralFun.R", local=TRUE)

############################################################
#####Function 1: Get missing Genes
############################################################
getMissGenes <- function(input, data) {
      NAME <- GetName(input$NAME)
      MissGene <- NAME[is.na(match(tolower(NAME), tolower(data$"gene_name")))]
      print("Name of Missing Gene:")
      MissGene
}

############################################################
#####Function 2: ScatterPlot
############################################################
ScatterplotGene <- function(input, data) {
  NAME <- GetName(input$NAME)
  if (is.na(match(tolower(NAME), tolower(data$gene_name)))) {message("Check Gene Name")
  } else {
    NAME <- intersect(NAME, data$gene_name)
    data$Col <- ifelse(!is.na(match(data$gene_name, NAME)), "Hit", "Others")
    #data$Size <- ifelse(!is.na(match(data$gene_name, NAME)), 1.1, 1)
    p <- ggplot(data[1:2000,], aes(WT_1, WT_2)) +
      geom_point(size=2, aes(color=Col)) +
      geom_point(data=subset(data, Col == "Hit"), aes(WT_1, WT_2, color=Col), size=2) +
      scale_color_manual(values=c("brown", alpha("lightgrey", 0.5))) +
      geom_text_repel(data=subset(data, Col == "Hit"), aes(label = gene_name), 
                      size=4, color='brown', box.padding = unit(0.15, "lines"),
                      point.padding = unit(1, "lines")) +
      xlab(paste0("WT_1")) +
      ylab(paste0("WT_2")) + 
      ggtitle("TPM of All Expressed Genes in WT") + 
      theme_bw(base_size = 16)
    print(p)
  }
}

############################################################
#####Function 3: DataTable
############################################################
TableGene <- function(input, data) {
  NAME <- GetName(input$NAME)
  data[match(tolower(NAME), tolower(data$gene_name)), c(7,GetIndex_Gene(input))]
}

############################################################
#####Function 3: Heatmap
############################################################
HeatmapGene <- function(input, data) {
  NAME <- GetName(input$NAME)
  outfile <- tempfile(fileext='.png')
  png(outfile, width=700, height=632.44 + 11.66 * nrow(data[match(GetName(tolower(input$NAME)), tolower(data$"gene_name")),]))
  heatmap.data <- data[match(tolower(NAME), tolower(data$gene_name)), c(GetIndex_Gene(input),7)]
  heatmap.data <- heatmap.data[complete.cases(heatmap.data),] #Remove DataTable with unmatched Gene Symbol 
  dend.col<-as.dendrogram(
    hclust(as.dist(1-cor(scale(as.matrix(heatmap.data[,1:(ncol(heatmap.data)-1)])),method="pearson")),method = "complete", members=NULL))
  dend.row<-as.dendrogram(
    hclust(as.dist(1-cor(scale(t(as.matrix(heatmap.data[,1:(ncol(heatmap.data)-1)]))),method="pearson")),method = "complete", members=NULL))
  ColSideColors <- GetColor_Gene(input)$col
  heatmap.2(as.matrix(heatmap.data[,1:(ncol(heatmap.data)-1)]),dendrogram="row", Colv=F, Rowv=dend.row,
            col=bluered(1000),scale="row",ColSideColors=ColSideColors,
            key=T,keysize=0.8,symkey=FALSE, lhei=c(0.55,4),
            density.info="none", trace="none",cexRow=1.6,cexCol=1.6,margin=c(14,24),
            labRow = heatmap.data[,"gene_name"],main="Heatmap of Interested Genes"
  )
  legend("topright", fill=unique(ColSideColors), cex=1.2, bty="n", legend=unique(GetColor_Gene(input)$leg))
  dev.off()
  list(src = outfile,
       contentType = 'image/png',
       width = 800,
       height = 632.44 + 11.66 * nrow(data[match(tolower(NAME), tolower(data$"gene_name")),]),
       alt = "This is alternate text")
}

  
   
