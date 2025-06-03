
library(edgeR)
library(openxlsx)
library(stringr)

# Experiments.
EXP <- c(2,4,8,12,16)
PATH <- "alignments.sulf_clock."

dir.create("plots/volcano", recursive=TRUE, showWarnings=FALSE)
dir.create("reports", recursive=TRUE, showWarnings=FALSE)

# Gene annotations.
genes <- readRDS("r_data/gene_annotations_Cel.Rdata")

wb <- createWorkbook()

for(i in 1:length(EXP)) {
  
  # Read in the feature counts.
  fc <- read.table(paste0("counts/sulf_N2.featureCounts"), header=TRUE)
  
  row.names(fc) <- fc$Geneid
  colnames(fc) <- sub(PATH,"",colnames(fc))
  colnames(fc) <- sub(".bam","",colnames(fc))
  fc <- as.matrix(fc[,-c(1:6)])
  
  stopifnot(all(row.names(fc)==genes$gene_id))
  
  fc <- fc[,colnames(fc) %in% colnames(fc)[grepl(paste0("_",EXP[i],"_"),colnames(fc))]]
  
  show(colnames(fc))
  
  fc.col <- colnames(fc)
  col.n <- c()
  for(j in 1:length(fc.col)) {
    col.s <- str_split(fc.col[j], "_")[[1]]
    col.n <- append(col.n, paste0(col.s[4],"_d",col.s[3],"r",col.s[5]))
  }
  
  colnames(fc) <- col.n
  
  show(colnames(fc))
  
  # Grouping factor.
  group <- factor(c("T","T","T","T","C","C","C","C"))
  
  # Create the data object.
  y <- DGEList(counts=fc, group=group, genes=genes$gene_name)
  
  # Filter out lowly expressed genes.
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  # Normalize the library sizes.
  y <- calcNormFactors(y)
  
  # Save MDS plot.
  pdf(file=paste0("plots/mds_day",EXP[i],".pdf"), width=8, height=8)
  plotMDS(y)
  dev.off()
  
  # Create design matrix.
  design <- model.matrix(~group)
  
  show(design)
  
  # Estimate common dispersion and tagwise dispersions.
  y <- estimateDisp(y, design)
  
  # Perform quasi-likelihood F-tests for DE.
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit)
  
  # Save table.
  tags <- as.data.frame(topTags(qlf, n=nrow(qlf)))
  tags <- cbind(row.names(tags), tags)
  colnames(tags) <- c("geneWbid","gene","logFC","logCPM","F","pValue","FDR")
  
  pf <- c()
  for(j in 1:nrow(tags)) {
    pf <- append(pf, ifelse(tags$FDR[j]<=.05 && abs(tags$logFC[j])>=1, TRUE, FALSE))
  }
  tags$passFilter <- pf
  
  addWorksheet(wb=wb, sheetName=paste0("de_day",EXP[i]))
  modifyBaseFont(wb, fontSize=16)
  writeDataTable(wb=wb, sheet=i, x=tags)
  
}

saveWorkbook(wb, paste0("reports/two_group_comparisons.xlsx"), overwrite=TRUE)
