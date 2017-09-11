GetName <- function(x) {return(gsub("( )+", "", unique(unlist(strsplit(x, split="\\n")))))}

GetIndex_Gene <- function(x) {
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

GetIndex_MIR <- function(x) {
  if (x$Dgcr8_MIR) {include.idx_Dgcr8=2:3} else {include.idx_Dgcr8=NA}
  if (x$Drosha_MIR) {include.idx_Drosha=4:5} else {include.idx_Drosha=NA}
  if (x$Dicer_MIR) {include.idx_Dicer=6:7} else {include.idx_Dicer=NA}
  if (x$Ago12_MIR) {include.idx_Ago12=8:9} else {include.idx_Ago12=NA}
  if (x$WT_MIR) {include.idx_WT=10:11} else {include.idx_WT=NA}
  include.idx <- c(include.idx_Dgcr8, include.idx_Drosha, include.idx_Dicer, 
                   include.idx_Ago12, include.idx_WT)
  include.idx <- include.idx[!is.na(include.idx)]
  return(include.idx)
}

GetColor_Gene <- function(x) {
  if (x$Dgcr8) {include.col_Dgcr8=rep("Dgcr8",2)} else {include.col_Dgcr8=NA}
  if (x$Drosha) {include.col_Drosha=rep("Drosha",2)} else {include.col_Drosha=NA}
  if (x$Dicer) {include.col_Dicer=rep("Dicer",2)} else {include.col_Dicer=NA}
  if (x$Ago12) {include.col_Ago12=rep("Ago12",2)} else {include.col_Ago12=NA}
  if (x$WT) {include.col_WT=rep("WT",2)} else {include.col_WT=NA}
  include.leg <- c(include.col_Dgcr8, include.col_Drosha, include.col_Dicer, include.col_Ago12, include.col_WT)
  include.leg <- include.leg[!is.na(include.leg)]
  include.col <- gsub("Dgcr8", rgb(215,14,14, maxColorValue=255), 
                      gsub("Drosha", rgb(253,181,14, maxColorValue=255),
                           gsub("Dicer",rgb(253,110,14, maxColorValue=255), 
                                gsub("Ago12", rgb(242,242,73, maxColorValue=255),
                                     gsub("WT", rgb(42,165,218, maxColorValue=255), 
                                          include.leg)))))
  include <- list(col=include.col,leg=include.leg)
  return(include)
}

GetColor_MIR <- function(x) {
  if (x$Dgcr8_MIR) {include.col_Dgcr8=rep("Dgcr8",2)} else {include.col_Dgcr8=NA}
  if (x$Drosha_MIR) {include.col_Drosha=rep("Drosha",2)} else {include.col_Drosha=NA}
  if (x$Dicer_MIR) {include.col_Dicer=rep("Dicer",2)} else {include.col_Dicer=NA}
  if (x$Ago12_MIR) {include.col_Ago12=rep("Ago12",2)} else {include.col_Ago12=NA}
  if (x$WT_MIR) {include.col_WT=rep("WT",2)} else {include.col_WT=NA}
  include.leg <- c(include.col_Dgcr8, include.col_Drosha, include.col_Dicer, include.col_Ago12, include.col_WT)
  include.leg <- include.leg[!is.na(include.leg)]
  include.col <- gsub("Dgcr8", rgb(215,14,14, maxColorValue=255), 
                      gsub("Drosha", rgb(253,181,14, maxColorValue=255),
                           gsub("Dicer",rgb(253,110,14, maxColorValue=255), 
                                gsub("Ago12", rgb(242,242,73, maxColorValue=255),
                                     gsub("WT", rgb(42,165,218, maxColorValue=255), 
                                          include.leg)))))
  include <- list(col=include.col,leg=include.leg)
  return(include)
}