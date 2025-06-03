
library(edgeR)
library(reshape2)
library(stringr)

PATH <- "alignments.sulf_clock."

# Gene annotations.
genes <- readRDS("r_data/gene_annotations_Cel.Rdata")

# Read in the feature counts.
fc <- read.table(paste0("counts/sulf_N2.featureCounts"), header=TRUE)

row.names(fc) <- fc$Geneid
colnames(fc) <- sub(PATH,"",colnames(fc))
colnames(fc) <- sub(".bam","",colnames(fc))
fc <- as.matrix(fc[,-c(1:6)])

fc <- fc[,-c(21,28)] # remove bad samples

stopifnot(all(row.names(fc)==genes$gene_id))

show(colnames(fc))

# Grouping factor.
group <- factor(c(rep("T",20), rep("C",18)))
age <- as.numeric(sapply(colnames(fc), function(x) unlist(strsplit(x, "_"))[3], USE.NAMES=FALSE))

# Create the data object.
y <- DGEList(counts=fc, group=group, genes=genes$gene_name)

# Filter out lowly expressed genes.
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Normalize the library sizes.
y <- calcNormFactors(y)

cpm.y <- as.data.frame(cpm(y))
cpm.y$gene_id <- rownames(cpm.y)
cpm.y <- merge(cpm.y, genes, by="gene_id", all.x=TRUE, all.y=FALSE)

gene.name <- c()
for(i in 1:nrow(cpm.y)) {
  gene.name <- append(gene.name, ifelse(is.na(cpm.y$gene_name[i])==TRUE, cpm.y$gene_id[i], cpm.y$gene_name[i]))
}

rownames(cpm.y) <- gene.name
cpm.y <- as.matrix(cpm.y[c(2:39)])

saveRDS(cpm.y, file=paste0("r_data/cpm_all_wt_samples.Rdata"))

sample.metadata <- y$samples[c(1)]
sample.metadata$age <- age

# Save MDS plot.
pch <- sapply(group, function(x) ifelse(x=="C", 21, 23), USE.NAMES=FALSE)
colors <- sapply(group, function(x) ifelse(x=="C", "darkgreen", "magenta"), USE.NAMES=FALSE)

pdf(file=paste0("plots/mds_time_course_comparison.pdf"), width=5, height=5)
plotMDS(y, pch=pch, bg=colors, col=colors)
dev.off()
