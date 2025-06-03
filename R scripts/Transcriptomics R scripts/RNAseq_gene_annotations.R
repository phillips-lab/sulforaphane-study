
library(refGenome)
library(R.utils)

# C. elegans annotation file.
GTF_FN <- "../custom_ref/108/Caenorhabditis_elegans.WBcel235.108.gtf.gz"

dir.create("r_data", recursive=TRUE, showWarnings=FALSE)

# Extract WormBase gene IDs and gene names.
gunzip(GTF_FN, overwrite=TRUE, remove=FALSE)
gtf <- ensemblGenome()
read.gtf(gtf, filename=sub("[.]gz$","",GTF_FN))
genes <- gtf@ev$genes[,c("gene_id","gene_name")]

saveRDS(genes, file="r_data/gene_annotations_Cel.Rdata")
