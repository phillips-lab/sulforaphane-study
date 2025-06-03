
library(ggplot2)
library(introdataviz)
library(openxlsx)


genes <- readRDS("r_data/gene_annotations_Cel.Rdata")


GROUPS <- c("dmso","sulf")
R2_CUTOFF <- .8


# Select genes with significant temporal effect.

est.tab <- readRDS("r_data/logistic_estimates_controls.Rdata")

est.tab.sub <- subset(est.tab, ((Slope.ConfLimLower > 0 & Slope.ConfLimUpper > 0) | (Slope.ConfLimLower < 0 & Slope.ConfLimUpper < 0)) & LowerLimit.ConfLimUpper < UpperLimit.ConfLimLower)

cat("Total models initially selected -", nrow(est.tab.sub),"\n")

# Calculate R-squared.

r2 <- function(pred,actual) { 
  return(1-sum((pred-actual)^2)/sum((actual-mean(actual))^2))
}

cpm.dmso <- readRDS(paste0("r_data/cpm_dmso_for_predictions.Rdata"))

r2.val <- c()

for(i in 1:nrow(est.tab.sub)) {
  
  log.exp <- cpm.dmso[est.tab.sub$GeneWbid[i],]
  dat <- data.frame(day=c(rep(2,3),rep(4,3),rep(8,4),rep(12,4),rep(16,4)), log.exp)
  
  u <- est.tab.sub$UpperLimit.Estimate[i]
  l <- est.tab.sub$LowerLimit.Estimate[i]
  s <- est.tab.sub$Slope.Estimate[i]
  e <- 10^(est.tab.sub$Log10EC50.Estimate[i])
  
  pred <- u+((l-u)/(1+(dat$day/e)^s))
  r2.val <- append(r2.val, r2(pred, dat$log.exp))
  
}

est.tab.sub$R2 <- r2.val
est.tab.sub <- est.tab.sub[order(-est.tab.sub$R2),]

#r.2 <- hist(est.tab.sub$R2, breaks=100)
#plot(r.2, main="R-squared values, n=368", xlab="R-squared", col="grey")

est.tab.sub <- subset(est.tab.sub, R2>R2_CUTOFF)

cat("Total models kept -", nrow(est.tab.sub),"\n")

# Save a table of the final clock genes.

clock.genes <- est.tab.sub[c("GeneWbid","Gene")]

label <- c()
for(i in 1:nrow(clock.genes)) {
  if(is.na(clock.genes$Gene[i])) {
    label <- append(label, clock.genes$GeneWbid[i])
  } else {
    label <- append(label, clock.genes$Gene[i])
  }
}

clock.genes$Gene <- label

wb <- createWorkbook()
addWorksheet(wb=wb, sheetName="clock_genes")
writeDataTable(wb=wb, sheet=1, x=clock.genes)
saveWorkbook(wb, "reports/clock_genes.xlsx", overwrite=TRUE)


# Predict biological age of each gene from their expression values and logistic parameters.

AGES <- c(2,4,8,12,16)

pred.age <- data.frame(
  Gene=character(),
  x=numeric(),
  y=numeric(),
  samp_group=character(),
  stringsAsFactors=FALSE
)

solve.for.x <- function(y,u,l,s,i) {
  x = i*(((l-u)/(y-u)-1)^(1/s))
}

# controls

cpm.tab <- readRDS(paste0("r_data/cpm_dmso_for_predictions.Rdata"))
cpm.tab.sub <- subset(cpm.tab, rownames(cpm.tab) %in% est.tab.sub$GeneWbid)
cat("Total genes tested (control) -", nrow(cpm.tab.sub),"\n")

for(i in 1:length(AGES)) {
  for(j in 1:nrow(cpm.tab.sub)) {
    
    cpm.gene <- cpm.tab.sub[rownames(cpm.tab.sub)[j],]
    if(i==1) { cpm <- mean(c(cpm.gene[c(1:3)])) }
    if(i==2) { cpm <- mean(c(cpm.gene[c(4:6)])) }
    if(i==3) { cpm <- mean(c(cpm.gene[c(7:10)])) }
    if(i==4) { cpm <- mean(c(cpm.gene[c(11:14)])) }
    if(i==5) { cpm <- mean(c(cpm.gene[c(15:18)])) }
    
    est <- subset(est.tab.sub, GeneWbid==rownames(cpm.tab.sub)[j], select=c("GeneWbid","UpperLimit.Estimate","LowerLimit.Estimate","Slope.Estimate","Log10EC50.Estimate"))
    
    pred.x <- solve.for.x(cpm,
                          est$UpperLimit.Estimate,
                          est$LowerLimit.Estimate,
                          est$Slope.Estimate,
                          10^(est$Log10EC50.Estimate)
                          )
    
    pred.age[nrow(pred.age)+1,] <- list(est$GeneWbid, AGES[i], pred.x, "dmso")
    
  }
}

# treated

cpm.tab <- readRDS(paste0("r_data/cpm_sulf_for_predictions.Rdata"))
cpm.tab.sub <- subset(cpm.tab, rownames(cpm.tab) %in% est.tab.sub$GeneWbid)
cat("Total genes tested (treated) -", nrow(cpm.tab.sub),"\n")

for(i in 1:length(AGES)) {
  for(j in 1:nrow(cpm.tab.sub)) {
    
    cpm.gene <- cpm.tab.sub[rownames(cpm.tab.sub)[j],]
    cpm <- mean(c(cpm.gene[c((i*4-3):(i*4))]))
    
    est <- subset(est.tab.sub, GeneWbid==rownames(cpm.tab.sub)[j], select=c("GeneWbid","UpperLimit.Estimate","LowerLimit.Estimate","Slope.Estimate","Log10EC50.Estimate"))
    
    pred.x <- solve.for.x(cpm,
                          est$UpperLimit.Estimate,
                          est$LowerLimit.Estimate,
                          est$Slope.Estimate,
                          10^(est$Log10EC50.Estimate)
    )
    
    pred.age[nrow(pred.age)+1,] <- list(est$GeneWbid, AGES[i], pred.x, "sulf")
    
  }
}

pred.out <- merge(pred.age, genes, by.x="Gene", by.y="gene_id", all.x=TRUE, all.y=FALSE)
pred.out$samp_group <- replace(pred.out$samp_group, pred.out$samp_group=="dmso", "control")
pred.out$samp_group <- replace(pred.out$samp_group, pred.out$samp_group=="sulf", "treated")
pred.out <- pred.out[c(1,5,4,2,3)]
pred.out <- pred.out[order(pred.out$gene_name, pred.out$samp_group, pred.out$x), ]
pred.out[pred.out=="NaN"] <- NA
colnames(pred.out) <- c("geneWbid","gene","group","actual","predicted")

wb <- createWorkbook()
addWorksheet(wb=wb, sheetName="predicted_ages")
writeDataTable(wb=wb, sheet=1, x=pred.out)
saveWorkbook(wb, "reports/predicted_ages.xlsx", overwrite=TRUE)

pred.age <- subset(pred.age, x > 2)


# Kruskal–Wallis test at each age.
for(i in 2:length(AGES)) {
  cat("Age",AGES[i])
  kt <- subset(pred.out, actual==AGES[i], select=c("group","predicted"))
  kt$group <- as.factor(kt$group)
  res <- kruskal.test(predicted ~ group, data=kt)
  show(res)
}


# Violin plot of predicted ages for all clock genes in control and treated samples.
xy <- ggplot(pred.age, aes(x=factor(x), y=y, fill=samp_group)) +
  geom_split_violin(scale="width") +
  geom_point(pch=16, size=.5, position=position_jitterdodge(jitter.width=0.2)) +
  scale_x_discrete(expand = expand_scale(add = c(1,0.5))) +
  scale_y_continuous(limits=c(0,22), breaks=c(4,8,12,16)) +
  scale_fill_manual(values=c("#CCCCCC","#48884B")) +
  xlab("Actual age (days) of sample") + ylab("Predicted age (days)") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.line=element_line(color="black"),
        legend.position="none")

ggsave(paste0("plots/violin_plot_aging_clock.png"), width=7, height=8, xy, device="png")
ggsave(paste0("plots/violin_plot_aging_clock.svg"), width=7, height=8, xy, device="svg")

pred.age.s <- subset(pred.age, samp_group=="sulf", select=c("x","y"))
pred.age.s[nrow(pred.age.s)+1,] <- c(0,0)

# Fit a quadratic line to the predicted ages from the treated samples.
xy.quad <- ggplot(pred.age.s, aes(x=x, y=y)) +
  geom_point(pch=16, size=.5) +
  stat_smooth(method="lm", formula = y ~ x + I(x^2) + 0, se=FALSE, color="darkgreen", linewidth=1) +
  scale_x_continuous(limits=c(0,16.3), breaks=c(4,8,12,16)) +
  scale_y_continuous(limits=c(0,21.8), breaks=c(4,8,12,16)) +
  xlab("Actual age (days) of sample") + ylab("Predicted age (days)") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.line=element_line(color="black"),
        legend.position="none")

ggsave(paste0("plots/violin_plot_aging_clock_quadfit.png"), width=7, height=8, xy.quad, device="png")
ggsave(paste0("plots/violin_plot_aging_clock_quadfit.svg"), width=7, height=8, xy.quad, device="svg")
