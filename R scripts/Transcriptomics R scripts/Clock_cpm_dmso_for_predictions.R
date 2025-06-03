
library(edgeR)
library(openxlsx)
library(refGenome)
library(R.utils)
library(splines)

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

# Only the controls.
fc <- fc[,c(21:38)]
stopifnot(all(row.names(fc)==genes$gene_id))

show(colnames(fc))

# Grouping factor.
group <- sapply(colnames(fc), function(x) paste0("d",unlist(strsplit(x, "_"))[3]), USE.NAMES=FALSE)
days <- as.numeric(sub("d","",group))

# Create the data object.
y <- DGEList(counts=fc, group=group, genes=genes$gene_name)

# Filter out lowly expressed genes.
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Normalize the library sizes.
y <- calcNormFactors(y)

# Create design matrix.
X <- ns(days, df=3)
design <- model.matrix(~X)

show(design)

# Estimate common dispersion and tagwise dispersions.
y <- estimateDisp(y, design)

# Perform quasi-likelihood F-tests.
fit <- glmQLFit(y, design, robust=TRUE)

qlf <- glmQLFTest(fit, coef=2:4)

tags <- as.data.frame(topTags(qlf, n=nrow(qlf$table)))

cpm <- cpm(y, log=TRUE, prior.count=fit$prior.count)

saveRDS(cpm, file="r_data/cpm_dmso_for_predictions.Rdata")
