
library(dr4pl)
library(edgeR)
library(ggplot2)
library(openxlsx)
library(refGenome)
library(R.utils)
library(splines)

options(scipen=999)

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

coeff.tab <- data.frame(
  GeneWbid=character(),
  Gene=character(),
  UpperLimit.Estimate=numeric(),
  UpperLimit.StdErr=numeric(),
  UpperLimit.ConfLimLower=numeric(),
  UpperLimit.ConfLimUpper=numeric(),
  LowerLimit.Estimate=numeric(),
  LowerLimit.StdErr=numeric(),
  LowerLimit.ConfLimLower=numeric(),
  LowerLimit.ConfLimUpper=numeric(),
  Slope.Estimate=numeric(),
  Slope.StdErr=numeric(),
  Slope.ConfLimLower=numeric(),
  Slope.ConfLimUpper=numeric(),
  Log10EC50.Estimate=numeric(),
  Log10EC50.StdErr=numeric(),
  Log10EC50.ConfLimLower=numeric(),
  Log10EC50.ConfLimUpper=numeric(),
  GoF.F=numeric(),
  GoF.pVal=numeric(),
  stringsAsFactors=FALSE
)

# For each gene in the expression dataset, fit a logistic curve. Save the parameters for upper and lower bounds, slope, and inflection point along with their confidence intervals in a table.

fit.logistic <- function(dat, row){
  tryCatch(
    {
      model <- dr4pl(log.exp ~ days, data=dat, method.init="logistic")
      return(model)
    },
    error=function(e) {
      message("An error occurred")
      print(e)
      return(NA)
    },
    warning=function(w) {
      message("A warning occurred")
      print(w)
      return(NA)
    },
    finally={
      message(paste("Row",row,"done"))
    }
  )
}

for(i in 1:nrow(tags)) {
  
  log.exp <- cpm[rownames(tags)[i],]
  dat <- data.frame(days, log.exp)
  model <- fit.logistic(dat, i)
  
  if(!(is.na(model[1]))) {
    
    summ.model <- summary(model)
    gof.model <- dr4pl::gof(model)
    
    coeff <- as.data.frame(summ.model$coefficients)
    
    coeff.tab[nrow(coeff.tab)+1,] <- list(
      rownames(tags)[i],
      tags$genes[i],
      coeff[1,1],
      coeff[1,2],
      coeff[1,3],
      coeff[1,4],
      coeff[4,1],
      coeff[4,2],
      coeff[4,3],
      coeff[4,4],
      coeff[3,1],
      coeff[3,2],
      coeff[3,3],
      coeff[3,4],
      coeff[2,1],
      coeff[2,2],
      coeff[2,3],
      coeff[2,4],
      as.numeric(gof.model[[7]]),
      as.numeric(gof.model[[9]])
      )
    
  }
  
}

saveRDS(coeff.tab, file="r_data/logistic_estimates_controls.Rdata")

wb <- createWorkbook()
addWorksheet(wb=wb, sheetName="logistic_estimates")
modifyBaseFont(wb, fontSize=14)
writeDataTable(wb=wb, sheet=1, x=coeff.tab)
setColWidths(wb, sheet=1, cols=1:length(coeff.tab), widths="auto")
saveWorkbook(wb, "reports/logistic_estimates_controls.xlsx", overwrite=TRUE)
