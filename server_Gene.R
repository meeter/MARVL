source("server_GeneralFun.R", local=TRUE)
library(visNetwork)

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
  NAME <- tolower(GetName(input$NAME))
  if (is.na(match(NAME, tolower(data$gene_name)))) {message("Check Gene Name")
  } else {
    NAME <- intersect(NAME, tolower(data$gene_name))
    data$Col <- ifelse(!is.na(match(tolower(data$gene_name), NAME)), "Hit", "Others")
    #data$Size <- ifelse(!is.na(match(data$gene_name, NAME)), 1.1, 1)
    p <- ggplot(data[1:2000,], aes(WT_1, WT_2)) +
      geom_point(size=2, aes(color=Col)) +
      geom_point(data=subset(data, Col == "Hit"), aes(WT_1, WT_2, color=Col), size=2) +
      scale_color_manual(values=c("brown", alpha("lightgrey", 0.5))) +
      geom_label_repel(data=subset(data, Col == "Hit"), aes(label = gene_name), 
                      label.size=0.5, color='brown', box.padding = unit(0.15, "lines"),
                      point.padding = unit(1, "lines")) +
      xlab(paste0("WT_1")) +
      ylab(paste0("WT_2")) + 
      ggtitle("Log2-TPM of All Expressed Genes in WT") + 
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
#####Function 4: Heatmap
############################################################
HeatmapGene <- function(input, GeneTE_1.tpm) {
  NAME <- GetName(input$NAME)
  #outfile <- tempfile(fileext='.png')
  #png(outfile, width=700, height=700 + 11.66 * nrow(data[match(GetName(tolower(input$NAME)), tolower(data$"gene_name")),]))
  heatmap.data <- GeneTE_1.tpm[match(tolower(NAME), tolower(GeneTE_1.tpm$gene_name)), c(GetIndex_Gene(input),7)]
  heatmap.data <- heatmap.data[complete.cases(heatmap.data),] #Remove DataTable with unmatched Gene Symbol 
  row.names(heatmap.data) <- gsub("mmu-", "", heatmap.data[, ncol(heatmap.data)])
  heatmap.data <- heatmap.data[, -ncol(heatmap.data)]
  #dend.col<-as.dendrogram(
    #hclust(as.dist(1-cor(scale(as.matrix(heatmap.data[,1:(ncol(heatmap.data)-1)])),method="pearson")),method = "complete", members=NULL))
  #dend.row<-as.dendrogram(
    #hclust(as.dist(1-cor(scale(t(as.matrix(heatmap.data[,1:(ncol(heatmap.data)-1)]))),method="pearson")),method = "complete", members=NULL))
  #ColSideColors <- GetColor_Gene(input)$col
  #heatmap.2(as.matrix(heatmap.data[,1:(ncol(heatmap.data)-1)]),dendrogram="none", Colv=F, Rowv=F,
  #          col=bluered(1000),scale="row",ColSideColors=ColSideColors,
  #          key=T,keysize=0.8,symkey=FALSE, lhei=c(0.55,4),
  #          density.info="none", trace="none",cexRow=1.6,cexCol=1.6, margin=c(2,2),
  #          labRow = heatmap.data[,"gene_name"],main="Heatmap of Interested Genes"
  #legend("topright", fill=unique(ColSideColors), cex=1.2, bty="n", legend=unique(GetColor_Gene(input)$leg))
  data <- as.matrix(heatmap.data); data <- t(apply(data, 1, function(x){2*(x-mean(x))/(max(x)-min(x))}))
  plot_ly(y = row.names(heatmap.data), x = colnames(heatmap.data)[1:ncol(heatmap.data)], 
          z = data, colors = colorRamp(c("skyblue", "white", "red")),
          type = "heatmap", colorbar = list(title = "Log2-TPM"), height = 300+11.66*nrow(heatmap.data)
          ) %>%
    layout(xaxis = list(title = ""),  yaxis = list(title = ""), margin = list(l = 100, b = 100))
  #heatmaply(heatmap.data, scale='none', Rowv=F, Colv=F, 
  #          scale_fill_gradient_fun = scale_fill_gradient2(low = "skyblue", high = "red", midpoint = 5),
  #          key.title="Log2-TPM", margins=c(120,100)
  #)
  #dev.off()
  #list(src = outfile,
  #     contentType = 'image/png',
  #     width = 800,
  #     height = 700 + 11.66 * nrow(data[match(tolower(NAME), tolower(data$"gene_name")),]),
  #     alt = "This is alternate text")
}

############################################################
#####Function 5: Dots plot for Gene2miR
############################################################
Gene2miR <- function(input, GeneTE_1.tpm, WT_MIR, miR_Gene.MW.sel){
  NAME <- GetName(input$NAME)
  tmp <- miR_Gene.MW.sel[!is.na(match(tolower(miR_Gene.MW.sel[,2]), tolower(NAME))),]
  tmp <- subset(tmp, spearman < -0.7)
  tmp <- merge(tmp, GeneTE_1.tpm[, 7:17], by="gene_name", all.x=T)
  tmp <- merge(tmp, WT_MIR[, 1:11], by="ID", all.x=T)
  colnames(tmp) <- gsub(".y", "_miR", gsub(".x", "_Gene", colnames(tmp)))
  tmp$gene_name <- factor(tmp$gene_name)
  tmp$ID <- gsub("mmu-", "", tmp$ID); tmp$gene_name <- as.character(tmp$gene_name)
  tmp <- tmp[order(tmp$gene_name),]
  tmp$targetScan <- paste("http://www.targetscan.org/cgi-bin/targetscan/mmu_71/targetscan.cgi?species=Mouse&gid=", tmp$gene_name, "&mir_sc=&mir_c=&mir_nc=&mir_vnc=&mirg=", sep="")
  tmp
}

Gene2miR_plot <- function(tmp){
  p <- ggplot(tmp, aes(gene_name, spearman))
  p3 <- p + geom_point(position = position_jitter(width = 0.2))
  p3 <- p3 + coord_flip() + xlab("")
  ggplotly(p3, width=600, height = 270 + 2.66 * length(unique(tmp$ID)) ) %>% 
  layout(dragmode = "select") 
}

############################################################
#####Function 6: Category Analysis
############################################################
CateFisher <- function(input, Dgcr8_DEG, Drosha_DEG, Dicer_DEG, Ago12_DEG,
                       name.overlap.up, name.overlap.down, 
                       name.MIG, name.MDG, name.Dicer, name.Ago12, 
                       GeneTE_1.tpm) {
  NAME <- GetName(input$NAME)
  expressed <- GeneTE_1.tpm$gene_name
  name.cate <- c("Dgcr8_DEG", "Drosha_DEG", "Dicer_DEG", "Ago12_DEG",
                 "name.overlap.up", "name.overlap.down",
                 "name.MIG", "name.MDG", "name.Dicer", "name.Ago12", "expressed")
  Fisher <- data.frame(matrix(data=NA, ncol = 4, nrow = 11)) 
  row.names(Fisher) <- c("Dgcr8", "Drosha", "Dicer", "Ago12", 
                         "Up", "Down",  "MIG", "MDG", "Dicer_Spec", "Ago12_Spec", "Expressed")
  colnames(Fisher) <- c("Enrichment_Score", "OR", "Number", "Overlapped_Gene")
  for (i in 1 : nrow(Fisher)) {
      tmp <- intersect(tolower(NAME), tolower(get(name.cate[i])))
      Fisher[i, 3] <- a <- length(tmp)
      Fisher[i, 4] <- paste(GeneTE_1.tpm[match(tmp, tolower(GeneTE_1.tpm$gene_name)), "gene_name"], collapse = ";")
      b <- length(get(name.cate[i])) - a
      c <- length(NAME) - a
      d <- 43629  - a - b - c ##Total Gene: 43629
      Fisher[i, 1] <- fisher.test(matrix(data = c(a,b,c,d), ncol=2))$p.value
      Fisher[i, 2] <- fisher.test(matrix(data = c(a,b,c,d), ncol=2))$estimate
  }
  Fisher[, 1] <- p.adjust(Fisher[, 1], method='BH')
  Fisher[, 1] <- -log(Fisher[, 1] + 1e-30, 10)
  for (i in 1 : nrow(Fisher)) { ##Enrichment score, enriched: -log(FDR, 10); depleted: log(FDR, 10)
    if (Fisher[i, 2] > 1) {Fisher[i, 1] <- Fisher[i, 1]
      } else {Fisher[i, 1] <- -1 * Fisher[i, 1]}
  }
  #Heatmap
  #plot_ly(y = row.names(Fisher.p), x = colnames(Fisher.p)[1:ncol(Fisher.p)], 
  #        z = as.matrix(Fisher.p[,1:ncol(Fisher.p)]), colors = colorRamp(c("skyblue", "white", "red")),
  #        type = "heatmap", colorbar = list(title = "Enrichment Score"), height = 300
  #) 
  #%>%
  #  layout(xaxis = list(title = ""),  yaxis = list(title = ""), margin = list(l = 100, b = 100))
  
  #Barplot by plot_ly
  #p <- plot_ly(x = as.numeric(Fisher.p), #y ~ reorder(colnames(Fisher.p), as.numeric(Fisher.p)),
  #             y = colnames(Fisher.p),  
  #             type = 'bar', orientation = 'h')  %>%
  #     layout(title = "Enrichment Score: > 2, enriched; < -2, depleted", 
  #            margin = list(l = 120, b = 100, t=50, r=50))
  
  #Barplot by ggplot
  Fisher$Gene_Cate <- row.names(Fisher)
  #Fisher.long <- melt(Fisher[, c(1,4,5)])
  Fisher$Gene_Cate <- factor(Fisher$Gene_Cate, 
                            levels=c("Dgcr8", "Drosha", "Dicer", "Ago12", 
                                     "Up", "Down",  
                                     "MIG", "MDG", "Dicer_Spec", "Ago12_Spec", 
                                     "Expressed"))
  Fisher
}
CateFisher.plot <- function(Fisher)
{
  p <- ggplot(Fisher, aes(y = Enrichment_Score, x = Gene_Cate, fill = Gene_Cate, text = paste("Gene:", Overlapped_Gene))) + 
       geom_bar(stat = "identity") + 
       theme(legend.position="none") +
       xlab("") + ylab("")
  ggplotly(p, tooltip = c("y", "x")) %>% 
    layout(dragmode = "click", title = "Enrichment Score: > 2, enriched; < -2, depleted", 
           margin = list(l = 50, b = 50, t=50, r=50)) 
       #coord_flip()  + 
  #ggplotly(tooltip = c("text"))
}
  

  

  
   
