# Set Up ------------------------------------------------------------------

# load packages
library(readxl); library(coxme); library(multcomp); 
library(lme4); library(lubridate); 

# set active directory
activedir <- file.choose()
setwd(dirname(activedir))

# load data
comps <- read_excel(".xlsx")
date <- format(today(), "%m%d%y")

# backup dataset
origcomps <- comps

comps$Lab <- as.factor(comps$Lab)
comps$StartDate <- as.factor(comps$StartDate)
comps$Tech <- as.factor(comps$Tech)
comps$Compound <- as.factor(comps$Compound)
comps$Concentration <- as.factor(comps$Concentration)
comps$Rep <- as.factor(comps$Rep)
comps$death_id <- as.factor(comps$death_id)
comps$observation_id <- as.factor(comps$observation_id)
comps$plate_id <- as.factor(comps$plate_id)
comps$experiment_id <- as.factor(comps$experiment_id)

comps$Comp_Conc <- droplevels(with(comps, interaction(Compound,Concentration, sep="x")))


# set up multiple comparisons
comps_cc <- levels(comps$Comp_Conc)
assign(paste("comps","cc","pairs",sep="_"), c())
assign(paste("comps","cc","len",sep="_"), length(get(paste("comps","cc",sep="_"))))
for (k in 1:(comps_cc_len-1)) {
  comps_cc_pairs[k] <- c(paste(comps_cc[comps_cc_len+1-k],"-",comps_cc[1], " = 0"))
}


# Cox ---------------------------------------------------------------------

# main model, use to get the variance of the random effects
# multiple comparisons (estimates with p-values), used in figures


assign(paste("n2pd1073","cox",sep="_"),coxme(Surv(DeathAge,Dead) ~ Comp_Conc + (1|StartDate/Rep/plate_id), data=comps))
assign(paste("n2pd1073","cox","pairs",sep="_"),summary(glht(get(paste("n2pd1073","cox",sep="_")), linfct = mcp(Comp_Conc = comps_cc_pairs))))


# GLM ---------------------------------------------------------------------

# main model, use to get the standard deviation of the random effects
# confidence interval, of the standard deviation of the random effects
# (turn standard deviation into variance)
# multiple comparisons (estimates with p-values)


assign(paste("n2pd1073","lm",sep="_"),lmer(DeathAge ~ Comp_Conc + (1|StartDate/Rep/plate_id), data=comps))
assign(paste("n2pd1073","lm","prof",sep="_"),profile(get(paste("n2pd1073","lm",sep="_"))))
assign(paste("n2pd1073","lm","prof","CI",sep="_"),confint(get(paste("n2pd1073","lm","prof",sep="_"))))
assign(paste("n2pd1073","lm","pairs",sep="_"),summary(glht(get(paste("n2pd1073","lm",sep="_")), linfct = mcp(Comp_Conc = comps_cc_pairs))))



# Complete Output via Text ------------------------------------------------

# separated by strain
{
  sink(paste("n2pd1073","_",date,".txt",sep=""), append = FALSE)
  writeLines(paste("N2",sep=""))  
  writeLines("\nCox ####################\n")
  print(get(paste("n2pd1073","cox",sep="_")))
  writeLines("\nCox Planned Comparisons\n")
  print(get(paste("n2pd1073","cox","pairs",sep="_")))
  writeLines("\nLinear Model ####################\n")
  print(summary(get(paste("n2pd1073","lm",sep="_"))))
  print(get(paste("n2pd1073","lm","prof","CI",sep="_")))
  writeLines("\nLM Planned Comparisons\n")
  print(get(paste("n2pd1073","lm","pairs",sep="_")))
  sink()
}


# Quantiles and summary --------------------------------------------------------------

summary_list <- list()

survfit_obj <- survfit(Surv(DeathAge, Dead) ~ Comp_Conc, data = comps)

# Calculate quantiles on the survfit object
quantiles <- summary(survfit_obj)$table

# Calculate the 90th quantile
quantile_90_list <- sapply(levels(comps$Comp_Conc), function(comp_conc) {
  survfit_obj_comp <- survfit(Surv(DeathAge, Dead) ~ 1, data =comps[comps$Comp_Conc == comp_conc, ])
  quantile(survfit_obj_comp, probs = 0.9)$quantile[1]
})


summary_df <- data.frame(
  Strain = "n2_pd1073",
  Comp_Conc = rownames(quantiles),
  n_dead = quantiles[,"events"],
  n_censored = quantiles[,"records"] - quantiles[,"events"],
  n_total = quantiles[,"records"],
  median = quantiles[,"median"],
  LCL_95 = quantiles[,"0.95LCL"],
  UCL_95 = quantiles[,"0.95UCL"],
  rmean = quantiles[,"rmean"],
  se_rmean = quantiles[,"se(rmean)"],
  quantile_90 = quantile_90_list
)

summary_list <- summary_df

# Export to a CSV file
write.csv(summary_list, "summary_results.csv", row.names = FALSE)


