# Set Up ------------------------------------------------------------------

# load packages
library(readxl); library(coxme); library(multcomp); 
library(lme4); library(lubridate); 

# set active directory
activedir <- file.choose()
setwd(dirname(activedir))

# load data
comps <- read_excel("SFN_latelife.xlsx")
date <- format(today(), "%m%d%y")

# backup dataset
origcomps <- comps

comps$Compound <- as.factor(comps$Compound)
comps$Concentration <- as.factor(comps$Concentration)
comps$Condition <- as.factor(comps$Condition)
comps$StartDate <- as.factor(comps$StartDate)
comps$Tech <- as.factor(comps$Tech)
comps$Comp_Cond <- droplevels(with(comps, interaction(Compound,Condition, sep="x")))

levels(comps$Strain)
levels(comps$Comp_Cond)

comps_cc_pairs <- c("sulforaphanexstandard - CTRL_DMSOxstandard = 0",
                    "sulforaphanexD4_initiation - CTRL_DMSOxstandard = 0",
                    "sulforaphanexD8_initiation - CTRL_DMSOxstandard = 0",
                    "sulforaphanexD12_initiation - CTRL_DMSOxstandard = 0")


# Cox ---------------------------------------------------------------------

cox <- coxme(Surv(DeathAge,Dead) ~ Comp_Cond + (1|StartDate/Tech/plate_id), data=comps)
cox_pairs <- summary(glht(cox, linfct = mcp(Comp_Cond = get(paste("comps","cc","pairs",sep="_")))))


# GLM ---------------------------------------------------------------------

lm <- lmer(DeathAge ~ Comp_Cond + (1|StartDate/Tech/plate_id), data=comps)
lm_prof <- profile(lm)
lm_prof_ci <- confint(lm_prof)
lm_pairs <- summary(glht(lm, linfct = mcp(Comp_Cond = get(paste("comps","cc","pairs",sep="_")))))


# text output -------------------------------------------------------------
{
  sink(paste("SFN_latelife",date,".txt",sep="_"), append = FALSE)
  writeLines("N2_PD1073")  
  writeLines("\nCox\n")
  print(get(paste("cox",sep="_")))
  print(get(paste("cox","pairs",sep="_")))
  writeLines("\nLinear Model\n")
  print(summary(get(paste("lm",sep="_"))))
  print(get(paste("lm","prof","ci",sep="_")))
  print(get(paste("lm","pairs",sep="_")))
  sink()
} 



# summary table -----------------------------------------------------------

summary_list <- list()
{
 # Create the survfit object and assign it to a variable
survfit_obj <- survfit(Surv(DeathAge, Dead) ~ Comp_Cond, data = comps)
assign(paste("N2_late", "survfit", sep = "_"), survfit_obj)
  
  # Calculate quantiles on the survfit object
  quantiles <- summary(survfit_obj)$table
  
  # Print the strain and the quantiles for debugging
  print(paste("N2_late_survfit"))
  print(quantiles)
  
  quantile_90 <- quantile(survfit_obj, probs = 0.9)$quantile[1]
  
  summary_df <- data.frame(
    Strain = "N2",
    Comp_comc = rownames(quantiles),
    n_dead = quantiles[,"events"],
    n_censored = quantiles[,"records"] - quantiles[,"events"],
    n_total = quantiles[,"records"],
    median = quantiles[,"median"],
    LCL_95 = quantiles[,"0.95LCL"],
    UCL_95 = quantiles[,"0.95UCL"],
    rmean = quantiles[,"rmean"],
    se_rmean = quantiles[,"se(rmean)"],
    quantile_90 = quantile_90
  )
  
  # Append the data frame to the list
  summary_list[[length(summary_list) + 1]] <- summary_df
}

# Combine the list into a single data frame
summary_results <- do.call(rbind, summary_list)

# Write the data frame to a CSV file
write.csv(summary_results, "latelife_summary_results.csv", row.names = FALSE)

