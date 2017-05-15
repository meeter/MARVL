library(DEXSeq)
library(ggplot2)
load("WT.RData")
shinyServer(function(input, output) {
####Scatterplot of Abundance by RPKM in E14
    #NAME <- GetName(input$NAME)
    output$MissGene <- renderText({
      NAME <- input$NAME
      MissGene <- NAME[is.na(match(tolower(NAME), tolower(GeneTE_1.tpm$"gene_name")))]
      print("Name of Missing Gene:")
      MissGene
     })  
    output$E14_Exp <- renderPlot({
      NAME <- input$NAME  
      if (is.na(match(tolower(NAME), tolower(GeneTE_1.tpm$gene_name)))) {message("Check Gene Name")
      } else {plot(GeneTE_1.tpm$WT_1, GeneTE_1.tpm$WT_2, cex=0.4, cex.main=2.5, cex.lab=1.5, type="p", pch=3, col='grey', 
					xlab="WT_1", ylab="WT_2", main="Gene Abundance (Log2-Transformed RPKM) in E14", xlim=c(0,15), ylim=c(0,15))
      points(GeneTE_1.tpm[match(tolower(NAME), tolower(GeneTE_1.tpm$gene_name)), ]$WT_1,
	       		GeneTE_1.tpm[match(tolower(NAME), tolower(GeneTE_1.tpm$gene_name)), ]$WT_2,
		  col="brown3",pch=19,cex=1.4) 
		  }
    })
})
