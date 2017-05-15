library(ggplot2)
library(ggrepel)
library(gplots)
load("WT.RData")

GetName <- function(x) {return(unique(unlist(strsplit(x, split="\\n"))))}
GetIndex <- function(x) {
  if (x$Dgcr8) {include.idx_Dgcr8=8:9} else {include.idx_Dgcr8=NA}
  if (x$Drosha) {include.idx_Drosha=10:11} else {include.idx_Drosha=NA}
  if (x$Dicer) {include.idx_Dicer=12:13} else {include.idx_Dicer=NA}
  if (x$Ago12) {include.idx_Ago12=14:15} else {include.idx_Ago12=NA}
  if (x$WT) {include.idx_WT=16:17} else {include.idx_WT=NA}
  include.idx <- c(include.idx_Dgcr8, include.idx_Drosha, include.idx_Dicer, 
                   include.idx_Ago12, include.idx_WT)
  include.idx <- include.idx[!is.na(include.idx)]
  return(include.idx)
}
GetColor <- function(x) {
  if (x$Dgcr8) {include.col_Dgcr8=rep("Dgcr8",2)} else {include.col_Dgcr8=NA}
  if (x$Drosha) {include.col_Drosha=rep("Drosha",2)} else {include.col_Drosha=NA}
  if (x$Dicer) {include.col_Dicer=rep("Dicer",2)} else {include.col_Dicer=NA}
  if (x$Ago12) {include.col_Ago12=rep("Ago12",2)} else {include.col_Ago12=NA}
  if (x$WT) {include.col_WT=rep("WT",2)} else {include.col_WT=NA}
  include.leg <- c(include.col_Dgcr8, include.col_Drosha, include.col_Dicer, include.col_Ago12, include.col_WT)
  include.leg <- include.leg[!is.na(include.leg)]
  include.col <- gsub("Dgcr8",rep(rgb(215,14,14, maxColorValue=255),4), 
                      gsub("Drosha", rep(rgb(253,181,14, maxColorValue=255),2),
                          gsub("Dicer",rep(rgb(253,110,14, maxColorValue=255),2), 
                              gsub("Ago12", rep(rgb(242,242,73, maxColorValue=255),2),
                                  gsub("WT", rep(rgb(42,165,218, maxColorValue=255),2), 
                      include.leg)))))
  include <- list(col=include.col,leg=include.leg)
  return(include)
}


shinyServer(function(input, output) {
####Scatterplot of Abundance by RPKM in E14
    #NAME <- GetName(input$NAME)
    output$MissGene <- renderText({
      NAME <- GetName(input$NAME)
      MissGene <- NAME[is.na(match(tolower(NAME), tolower(GeneTE_1.tpm$"gene_name")))]
      print("Name of Missing Gene:")
      MissGene
     })  
    output$E14_Exp <- renderPlot({
      NAME <- GetName(input$NAME)
      if (is.na(match(tolower(NAME), tolower(GeneTE_1.tpm$gene_name)))) {message("Check Gene Name")
      } else {
        NAME <- intersect(NAME, GeneTE_1.tpm$gene_name)
        GeneTE_1.tpm$Col <- ifelse(!is.na(match(GeneTE_1.tpm$gene_name, NAME)), "Hit", "Others")
        #GeneTE_1.tpm$Size <- ifelse(!is.na(match(GeneTE_1.tpm$gene_name, NAME)), 1.1, 1)
        p <- ggplot(GeneTE_1.tpm[1:2000,], aes(WT_1, WT_2)) +
          geom_point(size=2, aes(color=Col)) +
          geom_point(data=subset(GeneTE_1.tpm, Col == "Hit"), aes(WT_1, WT_2, color=Col), size=2) +
          scale_color_manual(values=c("brown", alpha("lightgrey", 0.5))) +
          geom_text_repel(data=subset(GeneTE_1.tpm, Col == "Hit"), aes(label = gene_name), 
                          size=4, color='brown', box.padding = unit(0.15, "lines"),
                          point.padding = unit(1, "lines")) +
          xlab(paste0("WT_1")) +
          ylab(paste0("WT_2")) + 
          ggtitle("TPM of All Expressed Genes in WT") + 
          theme_bw(base_size = 16)
        print(p)
		  }
    })
    output$Data <- renderDataTable({
      NAME <- GetName(input$NAME)
      GeneTE_1.tpm[match(tolower(NAME), tolower(GeneTE_1.tpm$gene_name)), c(7,GetIndex(input))]
    }) 
    output$Heatmap <- renderImage({
      NAME <- GetName(input$NAME)
      outfile <- tempfile(fileext='.png')
      png(outfile, width=800, height=632.44 + 11.66 * nrow(GeneTE_1.tpm[match(GetName(tolower(input$NAME)), tolower(GeneTE_1.tpm$"Gene_Symbol")),]))
      heatmap.data <- GeneTE_1.tpm[match(NAME, GeneTE_1.tpm$gene_name), 7:11]
      heatmap.data <- GeneTE_1.tpm[match(tolower(NAME), tolower(GeneTE_1.tpm$gene_name)), c(GetIndex(input),7)]
      heatmap.data <- heatmap.data[complete.cases(heatmap.data),] #Remove DataTable with unmatched Gene Symbol 
      dend.col<-as.dendrogram(
        hclust(as.dist(1-cor(scale(as.matrix(heatmap.data[,1:(ncol(heatmap.data)-1)])),method="pearson")),method = "complete", members=NULL))
      dend.row<-as.dendrogram(
        hclust(as.dist(1-cor(scale(t(as.matrix(heatmap.data[,1:(ncol(heatmap.data)-1)]))),method="pearson")),method = "complete", members=NULL))
      ColSideColors <- GetColor(input)$col
      heatmap.2(as.matrix(heatmap.data[,1:(ncol(heatmap.data)-1)]),dendrogram="row", Colv=F, Rowv=dend.row,
                col=bluered(1000),scale="row",ColSideColors=ColSideColors,
                key=T,keysize=1,symkey=FALSE, lhei=c(0.55,4),
                density.info="none", trace="none",cexRow=1.6,cexCol=1.6,margin=c(14,24),
                labRow = heatmap.data[,"gene_name"],main="Heatmap of Interested Genes"
      )
      legend("topright", fill=unique(ColSideColors), cex=1.2, bty="n", legend=unique(GetColor(input)$leg))
      dev.off()
      list(src = outfile,
           contentType = 'image/png',
           width = 800,
           height = 632.44 + 11.66 * nrow(GeneTE_1.tpm[match(tolower(NAME), tolower(GeneTE_1.tpm$"gene)name")),]),
           alt = "This is alternate text")
    }, deleteFile = TRUE)
})
