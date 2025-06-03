
library(dr4pl)
library(ggplot2)
library(openxlsx)
library(R.utils)

options(scipen=999)

# Plot the expression of each clock gene individually (both control and treated).

cpm.dmso <- readRDS(paste0("r_data/cpm_dmso_for_predictions.Rdata"))
cpm.sulf <- readRDS(paste0("r_data/cpm_sulf_for_predictions.Rdata"))

dir.create("plots/fitted_expression_curves/individual", recursive=TRUE, showWarnings=FALSE)
dir.create("plots/fitted_expression_curves/for_grid", recursive=TRUE, showWarnings=FALSE)
clock.genes <- read.xlsx("reports/clock_genes.xlsx", sheet=1)

pred.y <- function(x,u,l,s,e) {
  u+((l-u)/(1+(x/e)^s))
}

for(i in 1:nrow(clock.genes)) {
  
  log.exp <- cpm.dmso[clock.genes$GeneWbid[i],]
  dat.dmso <- data.frame(day=c(rep(2,3),rep(4,3),rep(8,4),rep(12,4),rep(16,4)), log.exp)
  
  log.exp <- cpm.sulf[clock.genes$GeneWbid[i],]
  dat.sulf <- data.frame(day=c(rep(2,4),rep(4,4),rep(8,4),rep(12,4),rep(16,4)), log.exp)
  
  dr4pl.d <- dr4pl(log.exp ~ day, data=dat.dmso, method.init="logistic")
  
  # Plot logistic curve of the controls.
  p <- plot(dr4pl.d) +
    theme_classic() +
    geom_point(color="black", fill="black", shape=23, size=5.6) +
    scale_x_continuous(breaks=c(2,4,8,12,16))

  # Plot the sulf data.
  p <- p +
    geom_point(data=dat.sulf, aes(x=day, y=log.exp), color="#48884B", fill="#48884B", size=5.6)
  
  # Add a logistic curve to the sulf data.
  
  g.exclude <- c(
    "acs-1",
    "C32E12.1",
    "del-6",
    "F40F9.3",
    "F40G12.10",
    "F52H3.6",
    "irid-1",
    "kynu-1",
    "nep-20",
    "syp-1",
    "trxr-1",
    "utx-1",
    "W03F11.4",
    "WBGene00009651",
    "Y57G7A.6"
    ) # omit curve for genes with straight-line fit
  
  tryCatch(
    {
      if(!(clock.genes$Gene[i] %in% g.exclude)) {
        dr4pl.s <- dr4pl(log.exp ~ day, data=dat.sulf, method.init="logistic")
        summ.model <- summary(dr4pl.s)
        coeff <- as.data.frame(summ.model$coefficients)
        u <- coeff[1,1]; l <- coeff[4,1]; s <- coeff[3,1]; e <- 10^(coeff[2,1])
        p <- p + geom_function(fun = ~ pred.y(.x,u,l,s,e), color="#48884B", linewidth=1.5)
      }
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
      message(paste("Plot",i,"done"))
    }
  )
  
  p <- p +
    labs(title=clock.genes$Gene[i], x="Age in days of adulthood", y=expression(log[2]~CPM)) + 
    theme(plot.title=element_text(size=28),
          axis.title=element_text(size=20),
          axis.text=element_text(size=20))
  
  ggsave(filename=paste0("plots/fitted_expression_curves/individual/gene_",clock.genes$Gene[i],".png"), p, width=4.8, height=4.8, device="png")
  ggsave(filename=paste0("plots/fitted_expression_curves/individual/gene_",clock.genes$Gene[i],".svg"), p, width=4.8, height=4.8, device="svg")
  
  ylabs <- ggplot_build(p)$layout$panel_params[[1]]$y$get_labels()
  text.w <- max(strwidth(ylabs, units="inches", font=20, family="Arial"))*120
  
  pg <- p +
    labs(title=clock.genes$Gene[i], x=NULL, y=NULL) + 
    theme(plot.margin=unit(c(2,1,1,60-text.w),"pt"),
          plot.title=element_text(size=28),
          axis.title=element_text(size=20),
          axis.text=element_text(size=20))
  
  ggsave(filename=paste0("plots/fitted_expression_curves/for_grid/gene_",clock.genes$Gene[i],".png"), pg, width=4.8, height=4.5, device="png")
  ggsave(filename=paste0("plots/fitted_expression_curves/for_grid/gene_",clock.genes$Gene[i],".svg"), pg, width=dim, height=dim, device="svg")
  
}
