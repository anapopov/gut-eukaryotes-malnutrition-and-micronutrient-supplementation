######## alpha diversity stats - mixed-effects model and GLM


###-----  load libraries
library(phyloseq)
library(ape)
library(tidyr)
library(lme4)
library(lmerTest)
library(sjPlot) 




###----- set directories and load data

setwd("LME_GLM/")


### import phyloseq object with bacterial or eukaryotic OTU table, taxonomic information and metadata
load("microbiome.rdata")



### define microbial set, diversity metrics and rarefaction read depths
MEASURE <- c("Observed","Shannon")

MICROBES <- "protozoa"  # bacteria, eukaryotes, protozoa, fungi

if(MICROBES=="bacteria") {
  RAREDEPTH <- 25000
  microbe.set <- bacteria
} else if(MICROBES=="eukaryotes") {
  RAREDEPTH <- 1000
  microbe.set <- eukaryotes
} else if(MICROBES=="protozoa") {
  RAREDEPTH <- 1000
  microbe.set <- protozoa
} else if(MICROBES=="fungi") {
  RAREDEPTH <- 1000
  microbe.set <- fungi
}



###----- calculate alpha diversities (averaged over 100 rarefactions)
for(j in 1:length(MEASURE)) {

  print(MEASURE[j])

  # create results matrix
  df.mat <- as.data.frame(matrix(ncol = 100, nrow = 160))
  rownames(df.mat) <- sample_names(microbe.set)
  
  # rarefy 100X and record diversity values
  for(i in 1:100) {
    print(i)

    df <- rarefy_even_depth(physeq = microbe.set, sample.size = RAREDEPTH)
    x <- plot_richness(df, measures = MEASURE[j], x = "group2")$data[,c("samples","value")]
    rownames(x) <- x$samples
    df.mat[rownames(df.mat),i] <- x[rownames(df.mat),"value"]

    rm(x, df)
  }
  
  # calculate mean
  df.mat$value <- rowMeans(df.mat, na.rm = TRUE)
  
  # calculate standard deviation
  df.mat$sd <- apply(df.mat[,1:100], MARGIN = 1, sd, na.rm = TRUE)

  # remove original values
  df.mat[,1:100] <- NULL
  
  # add metadata
  df.mat[rownames(df.mat),"nutritionalstatus"] <- sample_data(microbe.set)[rownames(df.mat),"group2"]
  df.mat[rownames(df.mat),"arm"] <- sample_data(microbe.set)[rownames(df.mat),"arm"]
  df.mat[rownames(df.mat),"age"] <- sample_data(microbe.set)[rownames(df.mat),"timepoint"]
  df.mat[rownames(df.mat),"site"] <- sample_data(microbe.set)[rownames(df.mat),"center"]
  df.mat[rownames(df.mat),"whz"] <- sample_data(microbe.set)[rownames(df.mat),"whz"]
  df.mat$subject <- gsub(pattern = "A|B", replacement = "", rownames(df.mat))

  # create object
  assign(paste0("df.mat.", MEASURE[j]), df.mat, envir = globalenv())

  rm(df.mat)  
}




###----- calculate stats

# choose diversity output
mat <- df.mat.Observed


# set factor levels
mat$site <- factor(mat$site, levels = c("Urban","Rural"))
mat$arm <- factor(mat$arm, levels = c("control","MNP","MNP+Zinc"))
mat$nutritionalstatus <- factor(mat$nutritionalstatus, levels = c("normal","malnourished"))
mat$whz <- as.numeric(mat$whz)


##- view preliminary four-way anova
anova.obs <- aov(formula =  value ~ arm * site * nutritionalstatus * age, data = mat)
summary(anova.obs)


##- fit Linear Mixed-Effects Models
# try various combinations

# fit site only
fit0 <- lmer(value ~ site + (1 | subject), data= mat, na.action=na.omit, REML= FALSE)
summary(fit0)

# fit arm only
fit1 <- lmer(value ~ arm + (1 | subject), data= mat, na.action=na.omit, REML= FALSE)
summary(fit1)

# fit arm and site
fit2 <- lmer(value ~ arm + site + (1 | subject), data= mat, na.action=na.omit, REML= FALSE)
summary(fit2)

# fit arm, site and nutritional status
fit3 <- lmer(value ~ arm + site + nutritionalstatus + (1 | subject), data= mat, na.action=na.omit, REML= FALSE)
summary(fit3)

fit4 <- lmer(value ~ age + arm + site + (1 | subject), data= mat, na.action=na.omit, REML= FALSE)
summary(fit4)

fit5 <- lmer(value ~ age + arm + site + nutritionalstatus + (1 | subject), data= mat, na.action=na.omit, REML= FALSE)
summary(fit5)

fit6 <- lmer(value ~ site + age*arm*nutritionalstatus + (1 | subject), data= mat, na.action=na.omit, REML= FALSE)
summary(fit6)

fit7 <- lmer(value ~ age + arm + site + whz + (1 | subject), data= mat, na.action=na.omit, REML= FALSE)
summary(fit7)


# compare models
anova(fit0, fit2)
anova(fit1, fit2, fit3, fit4, fit5, fit6)

tab_model(fit1, fit2, fit3, fit4, fit5, fit6, fit7,
          title = "LME for alpha diversity, model comparison", 
          file = paste0("LME_alpha_diversity_model_comparison", RAREDEPTH,".html"))

# low ICC values - samples for same child have low correlation


##----------------------------------------------------------------------------------------
##- use general linear models (due to low intraperson correlation, ICC)

# fit all variables
fit1 <- glm(value ~ age*nutritionalstatus*site*arm, data = mat, na.action=na.omit)
summary(fit1)

# export full model
tab_model(fit1, title = "Full GLM for microbial richness", 
          file = paste0("GLM_full_microbial_richness_rarefied", RAREDEPTH,".html"))


# perform stepwise model selection with function 'stepAIC'
library(tidyverse)
library(caret)
library(leaps)
library(MASS)


# Fit the full model 
full.model <- glm(value ~ site*arm*age*nutritionalstatus, data = mat, na.action=na.omit)
summary(full.model)

# Stepwise regression model
step.model <- stepAIC(full.model, direction = "both", trace = TRUE)
summary(step.model)


# verify final model and export

### End of document.