source("server_GeneralFun.R", local=TRUE)
source("igraph_shape.R", local=TRUE)

############################################################
#####Function 1: Get missing Genes
############################################################
getMissGenes <- function(input, data) {
  NAME <- GetName(input$NAME_NET)
  MissGene <- NAME[is.na(match(tolower(NAME), tolower(data$"gene_name")))]
  print("Name of Missing Gene:")
  MissGene
}

############################################################
#####Function 2: Network Output
############################################################
plotNetwork <- function(input) {
  ######Order=0: direct neighbor of miss-regulated genes;
  ######Order=1: direct + 1st indirect neighbor of miss-regulated genes
  NAME <- GetName(input$NAME_NET)
  input_Network <- rbind(Gene.lasso.all.filt[!is.na(match(Gene.lasso.all.filt$Regulator, NAME)),],
                         Gene.lasso.all.filt[!is.na(match(Gene.lasso.all.filt$Gene, NAME)),])
  net <- graph.edgelist(as.matrix(unique(input_Network[!is.na(input_Network[,1]),])[, 1:2]), directed=T)
  ######Vertex Shape and color
  V(net)$Type <- sapply(V(net)$name, function(x) {
    if (length(grep("mir-|let-7", x)) > 0) return ("miR")
    else if (!is.na(match(x, TF.short.filt)) > 0) return ("TF")
    else if (!is.na(match(x, Hist.gene_name)) > 0 | !is.na(match(x, x5mc.gene_name)) > 0) return("HM_5mc")
    else return ("Gene")
  })
  vertex.shape <- ifelse(V(net)$Type == "TF", "circle", 
                         ifelse(V(net)$Type == "miR", "triangle", #
                                ifelse(V(net)$Type == "HM_5mc", "star", "square"))) 
  vertex.color  <- sapply(V(net)$name, function(x) {
    ifelse(length(grep("mir|let", x) > 0), 
           colorRampPalette(c("skyblue", "white", "pink"))(nrow(miR.res))[match(x, miR.res[order(miR.res$AvgLog2FC), ]$ID)],
           colorRampPalette(c("skyblue", "white", "pink"))(nrow(all.res))[match(x, all.res[order(all.res$AvgLog2FC), ]$gene_name)]
    )})
  vertex.frame.color<- ifelse(!is.na(match(V(net)$name, name.overlap.up)), "red", 
                              ifelse(!is.na(match(V(net)$name, name.MIG)), "orange",
                                     ifelse(!is.na(match(V(net)$name, name.MDG)), "yellow",		
                                            ifelse(!is.na(match(V(net)$name, name.Dicer)), "lightblue",
                                                   ifelse(!is.na(match(V(net)$name, name.Ago12)), "lightgreen",					
                                                          ifelse(!is.na(match(V(net)$name, name.overlap.down)), "blue", 
                                                                 "lightgrey"
                                                                 #rgb(211,211,211, max=255, alpha=50)
                                                          ))))))     
  #E(net)$color <- 'lightgrey'
  #E(net)[to(grep("Sirt6", V(net)$name))]$color <- 'steelblue'
  #E(net)[from(grep("Sirt6|Jund|Klf2", V(net)$name))]$color <- 'pink'
  
  ######Vertex Layout
  coords <- layout_on_grid(net); 
  coords.gene <- layout.norm(layout.fruchterman.reingold(induced_subgraph(net, which(V(net)$Type == "Gene"))), 
                             xmin=-6, xmax=6, ymin=-4, ymax=4)   
  coords.miR <- layout.norm(layout.fruchterman.reingold(induced_subgraph(net, which(V(net)$Type == "miR"))), 
                            xmin=-1.5, xmax=1.5, ymin=-8, ymax=-6)
  coords.TF <- layout.norm(layout.fruchterman.reingold(induced_subgraph(net, which(V(net)$Type == "TF"))),  
                           xmin=-8, xmax=-5, ymin=5, ymax=8)
  coords.HM_5mc <- layout.norm(layout.fruchterman.reingold(induced_subgraph(net, which(V(net)$Type == "HM_5mc"))), 
                               xmin=-1, xmax=8, ymin=4.5, ymax=8)
  
  layout <- coords; idx.miR <- idx.TF <- idx.HM_5mc <- idx.gene <- 0
  for (i in 1:length(V(net)$name)) {
    if (V(net)$Type[i] == "miR") {idx.miR <- idx.miR + 1; layout[i,] <- coords.miR[idx.miR,] * 0.2
    } else if (V(net)$Type[i] == "TF") {idx.TF <- idx.TF + 1; layout[i,] <- coords.TF[idx.TF,] * 0.25
    } else if (V(net)$Type[i] == "HM_5mc") {idx.HM_5mc <- idx.HM_5mc + 1; layout[i,] <- coords.HM_5mc[idx.HM_5mc,] * 0.3
    } else {idx.gene <- idx.gene + 1; layout[i,] <- coords.gene[idx.gene,] * 0.2}
  }
  
  ######Plot using VisNetwork
  net.Vis <- toVisNetworkData(net)
  net.Vis$nodes$shape <- vertex.shape
  net.Vis$nodes$color <- vertex.color
  lnodes <- data.frame(label = c("TF", "Epigenetic Regulator", "miRNA", "Gene"),  ###For legend
                       shape = c( "circle", "star", "triangle", "square"), 
                       color = "lightblue")
  visNetwork(nodes = net.Vis$nodes, edges = net.Vis$edges, height = "600px", width="100%") %>% 
  visInteraction(dragNodes = T, dragView = T, zoomView = F) %>%
  visLegend(addNodes = lnodes, useGroups = FALSE, width=0.10, main="Legend")
  #visOptions(manipulation = TRUE)  

  
  ######Plot using igraph																	
  #pdf("test.pdf",width=6, height=6, useDingbats=F)
  #plot(net, layout=layout, rescale=T,
  #     vertex.size = 8, vertex.shape=vertex.shape,
  #     vertex.color=vertex.color,
  #     vertex.label = V(net)$name,
  #     vertex.label.color = "black",
  #     vertex.frame.color = vertex.frame.color, vertex.frame.width = 1,
  #     vertex.label.cex=1,
  #     vertex.label.family = "sans",
  #     vertex.label.font = 0.5,  #vertex.label.dist=2, 
  #     edge.arrow.mode=2, edge.arrow.width=1
       #edge.color=E(net)$color
  #)
  #legend(bty="n","bottomleft",legend=c("Up","MIG", "MDG", "Dicer","Ago12", "Down","Others",  "TF", "HM_5mc", "miRNA", "Genes"), 
  #       col=c("red", "orange", "yellow", "lightblue", "lightgreen", "blue", "lightgrey",  rep("grey", 4)), 
  #       cex=5, pt.cex=c(rep(10,7), rep(10, 4)), pch=c(rep(19, 7), 16, 11, 17, 15))    

}
