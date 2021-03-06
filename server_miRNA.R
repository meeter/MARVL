source("server_GeneralFun.R", local=TRUE)

############################################################
#####Function 1: Get missing MIR
############################################################
getMissMIR <- function(input, data) {
  NAME <- paste("mmu-", GetName(input$NAME_MIR), sep="")
  MissMIR <- NAME[is.na(match(tolower(NAME), tolower(data$"ID")))]
  print("Name of Missing miRNAs :")
  MissMIR
}

############################################################
#####Function 2: ScatterPlot
############################################################
ScatterplotMIR <- function(input, data) {
  NAME <- paste("mmu-", GetName(input$NAME_MIR), sep="")
  if (is.na(match(tolower(NAME), tolower(data$ID)))) {message("Check Gene Name")
  } else {
    NAME <- intersect(NAME, data$ID)
    data$Col <- ifelse(!is.na(match(data$ID, NAME)), "Hit", "Others")
    #data$Size <- ifelse(!is.na(match(data$ID, NAME)), 1.1, 1)
    p <- ggplot(data, aes(WT_1, WT_2)) +
      geom_point(size=2, aes(color=Col)) +
      geom_point(data=subset(data, Col == "Hit"), aes(WT_1, WT_2, color=Col), size=2) +
      scale_color_manual(values=c("brown", alpha("lightgrey", 0.5))) +
      geom_label_repel(data=subset(data, Col == "Hit"), aes(label = ID), 
                       label.size=0.5, color='brown', box.padding = unit(0.15, "lines"),
                      point.padding = unit(1, "lines")) +
      xlab(paste0("WT_1")) +
      ylab(paste0("WT_2")) + 
      ggtitle("Log2-Normalized Counts of miRNAs in WT") + 
      theme_bw(base_size = 16)
    print(p)
  }
}

############################################################
#####Function 3: DataTable
############################################################
TableMIR <- function(input, data) {
  NAME <- paste("mmu-", GetName(input$NAME_MIR), sep="")
  data[match(tolower(NAME), tolower(data$ID)), c(1,GetIndex_MIR(input))]
}

############################################################
#####Function 4: Heatmap
############################################################
HeatmapMIR <- function(input, WT_MIR) {
  NAME <- paste("mmu-", GetName(input$NAME_MIR), sep="")
  #outfile <- tempfile(fileext='.png')
  #png(outfile, width=800, height=700 + 11.66 * nrow(data[match(tolower(NAME), tolower(data$"ID")),])) #632.44
  #heatmap.data <- data[match(NAME, data$ID), 7:11]
  heatmap.data <- WT_MIR[match(tolower(NAME), tolower(WT_MIR$ID)), c(GetIndex_MIR(input),1)]
  heatmap.data <- heatmap.data[complete.cases(heatmap.data),] #Remove DataTable with unmatched Gene Symbol 
  row.names(heatmap.data) <- gsub("mmu-", "", heatmap.data[, ncol(heatmap.data)])
  heatmap.data <- heatmap.data[, -ncol(heatmap.data)]
  #dend.col<-as.dendrogram(
  #  hclust(as.dist(1-cor(scale(as.matrix(heatmap.data[,1:(ncol(heatmap.data)-1)])),method="pearson")),method = "complete", members=NULL))
  #dend.row<-as.dendrogram(
  #  hclust(as.dist(1-cor(scale(t(as.matrix(heatmap.data[,1:(ncol(heatmap.data)-1)]))),method="pearson")),method = "complete", members=NULL))
  #ColSideColors <- GetColor_MIR(input)$col
  #heatmap.2(as.matrix(heatmap.data[,1:(ncol(heatmap.data)-1)]),dendrogram="none", Colv=F, Rowv=F,
  #          col=bluered(1000),scale="row",ColSideColors=ColSideColors,
  #          key=T,keysize=0.8,symkey=FALSE, lhei=c(0.55,4),
  #          density.info="none", trace="none",cexRow=1.6,cexCol=1.6, #margin=c(14,24),
  #          labRow = heatmap.data[,"ID"],main="Heatmap of Interested miRNAs"
  #)
  #legend("topright", fill=unique(ColSideColors), cex=1.2, bty="n", legend=unique(GetColor_MIR(input)$leg))
  data <- as.matrix(heatmap.data); data <- t(apply(data, 1, function(x){2*(x-mean(x))/(max(x)-min(x))}))
  plot_ly(y = row.names(heatmap.data), x = colnames(heatmap.data)[1:ncol(heatmap.data)], 
          z = data, colors = colorRamp(c("skyblue", "white", "red")), 
          type = "heatmap", colorbar = list(title = "Log2-Normalized Count")) %>%
    layout(xaxis = list(title = ""),  yaxis = list(title = ""), margin = list(l = 120, b = 100))
  #heatmaply(heatmap.data, scale='none', Rowv=F, Colv=F, subplot_widths=840,
  #          scale_fill_gradient_fun = scale_fill_gradient2(low = "skyblue", high = "red", midpoint = 5),
  #          key.title="Log2-Normalized Count", margins=c(120,100)
  #          )
  #dev.off()
  #return(list(src = outfile,
  #     contentType = 'image/png',
  #     width = 800,
  #     height = 700 + 11.66 * nrow(data[match(tolower(NAME), tolower(data$"ID")),]),
  #     alt = "This is alternate text"))
}

############################################################
#####Function 5: Barplot of RIP-seq
############################################################
BarplotMIR <- function(input, data) {
  NAME <- paste("mmu-", GetName(input$NAME_MIR), sep="")
  RIP.long <- melt(RIP[match(NAME, RIP$ID), ])
  RIP.long <- RIP.long[!is.na(RIP.long[,1]),]
  RIP.long$ID <- gsub("mmu-", "", as.character(RIP.long$ID))
  p1 <- ggplot(data=RIP.long, aes(x=ID,y=value, fill=variable)) +
    #scale_y_discrete(limits = value) + 
    geom_bar(position="dodge",stat="identity", width=0.5) + 
    #geom_line(data=BioType.long, aes(x=Var1, y=Cutoff)) +
    coord_flip() + xlab("") + ylab("Log2FC")
    ggtitle("Log2FC Comparing Ago1/2-RIP in E14 with KO")
  print(p1)
}

############################################################
#####Function 6: Dots plot for miR2Gene
############################################################
miR2Gene <- function(input, GeneTE_1.tpm, WT_MIR, miR_Gene.MW.sel){
  NAME <- GetName(input$NAME_MIR)
  tmp <- miR_Gene.MW.sel[!is.na(match(tolower(gsub("mmu-", "", miR_Gene.MW.sel[,1])), tolower(NAME))),]
  tmp <- subset(tmp, spearman < -0.7)
  tmp <- merge(tmp, GeneTE_1.tpm[, 7:17], by="gene_name", all.x=T)
  tmp <- merge(tmp, WT_MIR[, 1:11], by="ID", all.x=T)
  colnames(tmp) <- gsub(".y", "_miR", gsub(".x", "_Gene", colnames(tmp)))
  tmp$ID <- gsub("mmu-", "", tmp$ID)
  tmp$ID <- factor(tmp$ID); tmp$gene_name <- as.character(tmp$gene_name)
  tmp <- tmp[order(tmp$gene_name),]
  tmp
}

miR2Gene_plot <- function(tmp){
  p <- ggplot(tmp, aes(ID, spearman))
  p3 <- p + geom_point(position = position_jitter(width = 0.2))
  p3 <- p3 + coord_flip() + xlab("")
  ggplotly(p3, width=700, height = 300 + 11.66 * length(unique(tmp$ID)) ) %>% layout(dragmode = "select")
}

############################################################
#####Function 7: Biplot
############################################################
Biplot <- function(input) {
  NAME <- GetName(input$NAME_MIR)
  colnames(all.res.mds)[2:7] <- gsub("log2FoldChange","Log2FC", colnames(all.res.mds)[2:7])
  colnames(all.res.mds)[2:3] <- c("WT_1","WT_2")
  p.biplot <- ggbiplot(prcomp(all.res.mds[, 2:7], scale=T, center=F), obs.scale = 1, var.scale = 1, alpha=0.5, 
                       #labels.size = 2, labels = all.res.mds$ID, 
                       groups = all.res.mds$miRTron, ellipse = TRUE, ellipse.prob = 0.95, circle = F) +
              scale_color_discrete(name = '') +
              geom_text_repel(aes(label=all.res.mds$ID, col = all.res.mds$miRTron)) + 
              xlim(-5,4) + 
              theme(legend.direction = 'horizontal', legend.position = 'bottom')

  ggplotly(p.biplot,  source="biplot") %>% layout(dragmode = "select") %>%
          layout(autosize = T, width = 800, height=400)
}



