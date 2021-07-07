########## Beta-diversity analysis

library("phyloseq")
library("ape")
library("ggplot2")
library("reshape2")
library("vegan")
library("dplyr")
library("tidyr")
library("ggpubr")
library("DESeq2")


###-----  Define paths and import data 
setwd("beta_diversity/")
load("microbiome.rdata")


# filter low abundance 'noisy' taxa, that therefore have high variance, and transform
MINSUBJECT <- 0.05 #5%
MINCOUNT <- 5

DATATRANSFORM <- "rarefy" #deseq, rarefy, filter, log
AGEGROUP <- "all" #all, 12mo, 24mo


# choose dataset: bacteria, eukaryotes, protozoa, fungi

MICROBES <- "bacteria"


if(MICROBES=="bacteria") {
  RAREDEPTH <- 25000
  df <- bacteria
} else if(MICROBES=="eukaryotes") {
  RAREDEPTH <- 1000
  df <- eukaryotes
} else if(MICROBES=="protozoa") {
  RAREDEPTH <- 1000
  df <- protozoa
} else if(MICROBES=="fungi") {
  RAREDEPTH <- 1000
  df <- fungi
}



# rarefy OR transform
if(DATATRANSFORM=="rarefy") {
  df.filt <- rarefy_even_depth(physeq = df, sample.size = RAREDEPTH, rngseed = 1234)
  
} else if(DATATRANSFORM=="deseq") {
  df.filt <- filter_taxa(df, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)
  
  df.deseq <- phyloseq_to_deseq2(df.filt, ~timepoint)
  dds <- DESeq(df.deseq)
  dds <- estimateSizeFactors(df.deseq) #normalize
  NORMCOUNTS <- counts(dds, normalized=TRUE) #extract normalized counts
  NORMCOUNTS <- otu_table(NORMCOUNTS, taxa_are_rows = TRUE) # add normalized table to phyloseq object
  otu_table(df.filt) <- NORMCOUNTS
  rm(NORMCOUNTS,df.deseq,dds)
  
} else if(DATATRANSFORM=="relabund") {
  df.filt <- transform_sample_counts(physeq = df, function(x) x/sum(x))
  
} else if(DATATRANSFORM=="log") {
  df.filt <- filter_taxa(df, function(x) sum(x > MINCOUNT) > (MINSUBJECT*length(x)), TRUE)
  otu_table(df.filt) <- otu_table(df.filt) + 1
  df.filt <- transform_sample_counts(df.filt, log)
}



# subset if needed
if(AGEGROUP=="12mo"){
  df.filt <- subset_samples(physeq = df.filt, timepoint=="12mo")
} else if(AGEGROUP=="24mo"){
  df.filt <- subset_samples(physeq = df.filt, timepoint=="12mo")
}



###----- calculate distances and ordinate
set.seed(1234)


DIST <- "bray" # unifrac, wunifrac, bray
METHOD <- "NMDS" # NMDS, PCoA

p <- plot_ordination(df.filt, ordination = ordinate(df.filt, method = METHOD, distance=DIST), 
                     color="timepoint", axes = c(1,2)) +
        theme_bw() + geom_point(size=3) + ggtitle(paste0(METHOD, ", ", DIST)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + stat_ellipse(geom = "polygon", type="norm", alpha=0.05, aes(fill=arm), linetype = 2) 

print(p)



###----- check dimensionality and fit NMDS with vegan
DIST <- "bray"

# scree plot for NMDS - identify number of dimensions to minimize stress
library(goeveg)
dimcheckMDS(t(otu_table(df.filt)), k = 10, distance = DIST, autotransform = TRUE)

# set dimensions based on stress plot
DIM <- 3

# ordinate
df.ORD <- metaMDS(comm = t(otu_table(df.filt)), distance = DIST, k = DIM, trymax = 20, autotransform = TRUE)

# check scatter around stress plot
stressplot(df.ORD)

# plot ordinations with 95% confidence (standard deviation) ellipses
plot(df.ORD, main = paste0("NMDS : ", DIST), type = "points", display = "sites", choices = c(1,2))
ordiellipse(ord = df.ORD, groups=sample_data(df.filt)[,"timepoint"][[1]], conf = 0.95, col=c("grey20","grey70"), label = TRUE) # plot 0.95 conf ellipse



###-----------------------------------------------------------------------------------------------------
###----- Run stats with Adonis (permutational multivariate analysis of variance)


# set parameters
PERM <- 9999 # number of permutations
ORDMETHOD <- "PCoA"   # "NMDS","PCoA","CCA"


# calculate distance
set.seed(1234)
df.dist <- phyloseq::distance(df.filt, method = DIST) 


# calculate stats with ADONIS
if(!AGEGROUP %in% c("12mo","24mo")) { # both timepoints
adon.nutrition <- as.data.frame(adonis2(dist(df.dist) ~ center*arm*timepoint*group2, data = as.data.frame(as.matrix(sample_data(df.test))), perm=PERM))
adon.time <- as.data.frame(adonis2(dist(df.dist) ~ center*arm*group2*timepoint, data = as.data.frame(as.matrix(sample_data(df.test))), perm=PERM))
adon.supp <- as.data.frame(adonis2(dist(df.dist) ~ center*timepoint*group2*arm, data = as.data.frame(as.matrix(sample_data(df.test))), perm=PERM))
adon.suppMOD <- as.data.frame(adonis2(dist(df.dist) ~ center*timepoint*group2*arm.mod, data = as.data.frame(as.matrix(sample_data(df.test))), perm=PERM))
adon.loc <- as.data.frame(adonis2(dist(df.dist) ~ timepoint*group2*arm*center, data = as.data.frame(as.matrix(sample_data(df.test))), perm=PERM))

} else if(AGEGROUP %in% c("12mo","24mo")) { # single timepoint
  adon.nutrition <- as.data.frame(adonis2(dist(df.dist) ~ center*arm*group2, data = as.data.frame(as.matrix(sample_data(df.test))), perm=PERM))
  adon.supp <- as.data.frame(adonis2(dist(df.dist) ~ center*group2*arm, data = as.data.frame(as.matrix(sample_data(df.test))), perm=PERM))
  adon.suppMOD <- as.data.frame(adonis2(dist(df.dist) ~ center*group2*arm.mod, data = as.data.frame(as.matrix(sample_data(df.test))), perm=PERM))
  adon.loc <- as.data.frame(adonis2(dist(df.dist) ~ group2*arm*center, data = as.data.frame(as.matrix(sample_data(df.test))), perm=PERM))
}


# prepare results table:

# extract main effects
nuts <- adon.nutrition["group2",]
age <- adon.time["timepoint",]
supp <- adon.supp["arm",]
loc <- adon.loc["center",]

# average interactions
supp.age <- colMeans(rbind(adon.nutrition["arm:timepoint",],adon.time["arm:timepoint",],adon.supp["timepoint:arm",],adon.loc["timepoint:arm",]))
supp.nut <- colMeans(rbind(adon.nutrition["arm:group2",],adon.time["arm:group2",],adon.supp["group2:arm",],adon.loc["group2:arm",]))
age.nut <- colMeans(rbind(adon.nutrition["timepoint:group2",],adon.time["group2:timepoint",],adon.supp["timepoint:group2",],adon.loc["timepoint:group2",]))
loc.arm <- colMeans(rbind(adon.nutrition["center:arm",],adon.time["center:arm",],adon.supp["center:arm",],adon.loc["arm:center",]))
loc.age <- colMeans(rbind(adon.nutrition["center:timepoint",],adon.time["center:timepoint",],adon.supp["center:timepoint",],adon.loc["timepoint:center",]))
loc.nut <- colMeans(rbind(adon.nutrition["center:group2",],adon.time["center:group2",],adon.supp["center:group2",],adon.loc["group2:center",]))

supp.age.nut <- colMeans(rbind(adon.nutrition["arm:timepoint:group2",],adon.time["arm:group2:timepoint",],adon.supp["timepoint:group2:arm",],adon.loc["timepoint:group2:arm",]))
loc.arm.age <- colMeans(rbind(adon.nutrition["center:arm:timepoint",],adon.time["center:arm:timepoint",],adon.supp["center:timepoint:arm",],adon.loc["timepoint:arm:center",]))
loc.age.nut <- colMeans(rbind(adon.nutrition["center:timepoint:group2",],adon.time["center:group2:timepoint",],adon.supp["center:timepoint:group2",],adon.loc["timepoint:group2:center",]))
loc.arm.nut <- colMeans(rbind(adon.nutrition["center:arm:group2",],adon.time["center:arm:group2",],adon.supp["center:group2:arm",],adon.loc["group2:arm:center",]))
loc.age.nut.arm <- colMeans(rbind(adon.nutrition["center:arm:timepoint:group2",],adon.time["center:arm:group2:timepoint",],adon.supp["center:timepoint:group2:arm",],adon.loc["timepoint:group2:arm:center",]))

resid <- colMeans(rbind(adon.nutrition["Residual",],adon.time["Residual",],adon.supp["Residual",],adon.loc["Residual",]))
tot <- colMeans(rbind(adon.nutrition["Total",],adon.time["Total",],adon.supp["Total",],adon.loc["Total",]))

suppMod <- adon.suppMOD[c("arm.mod","Residual","Total"),]

# join into dataframe
adonis.results <- rbind(nuts,age,supp,loc,
                        supp.age,supp.nut,age.nut,loc.arm,loc.age,loc.nut,
                        supp.age.nut,loc.arm.age,loc.age.nut,loc.arm.nut,loc.age.nut.arm,
                        resid,tot,suppMod)
rownames(adonis.results) <- c("NutritionalStatus","Age","Arm","GeographicSetting",
                              "Arm:Age","Arm:NutritionalStatus","Age:NutritionalStatus",
                              "GeoLocation:Arm","GeoLocation:Age","GeoLocation:NutritionalStatus",
                              "Arm:Age:NutritionalStatus","GeoLocation:Arm:Age",
                              "GeoLocation:Age:NutritionalStatus","GeoLocation:Arm:NutritionalStatus",
                              "GeoLocation:Age:NutritionalStatus:Arm",
                              "Residual","Total",
                              "GroupedSupplementation","Residuals.GroupedSupplementation","Total.GroupedSupplementation")
# view results
print(adonis.results)




###----- calculate pairwise adonis (pairwise comparisons for supplementation groups)
source("parwise.adonis.r")
#Martinez Arbizu, P. (2019). pairwiseAdonis: Pairwise multilevel comparison using adonis. R package version 0.3
#https://github.com/pmartinezarbizu/pairwiseAdonis
pairwise.adonis(x = t(otu_table(df.filt)), factors = samplemetadata[rownames(as.matrix(df.dist)),"arm"], 
                sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")



###----- calculate beta-dispersions
VARS <- c("timepoint","group2","arm","center","arm.mod","arm.time","arm.time.group")

bdisp.tbl <- as.data.frame(matrix(ncol = 4, nrow = length(VARS)))
colnames(bdisp.tbl) <- c("var","Bdisp.pval","Fstat","nperm")
bdisp.tbl$nperm <- PERM

for(i in 1:length(VARS)) {
  DF.disp <- betadisper(d = df.dist, group = as.data.frame(as.matrix(sample_data(df.filt)))[rownames(as.matrix(df.dist)),VARS[i]], type = "centroid")
  set.seed(1234)
  DF.disp.p <- permutest(x = DF.disp, pairwise = FALSE, permutations = PERM)
  bdisp.tbl[i,1] <- VARS[i]
  bdisp.tbl[i,2] <- DF.disp.p$tab$`Pr(>F)`[[1]]
  bdisp.tbl[i,3] <- DF.disp.p$tab$F[[1]]
}

print(bdisp.tbl)



# plot dispersions by grouping variable and calculate pairwise differences
DIST <- "bray"

# calculate multivariate dispersions
mod <- betadisper(d = df.dist, as.data.frame(as.matrix(sample_data(df.filt)))[rownames(as.matrix(df.dist)),"arm.time"], type = "centroid")

# view F statistic and overall p val
disp.p <- permutest(x = mod, pairwise = FALSE, permutations = PERM)
disp.p$tab

# view boxplot
p <- boxplot(mod)

# calculate anova
anova(mod)

# get pairwise differences using Tukey's Honest Significant Difference
mod.HSD <- TukeyHSD(mod)
plot(mod.HSD)


# calculate pairwise differences with wilcoxon ranksum
x <- as.data.frame(as.matrix(sample_data(df.filt)))[,c("sample.id","arm.time")]
y <- as.data.frame(mod$distances)
y$sample.id <- rownames(y)
z <- cbind(x,y[rownames(x),])
z$sample.id <- NULL
colnames(z) <- c("var","dist")

stat.mat <- as.data.frame(matrix(ncol = 2, nrow = 6))
colnames(stat.mat) <- c("comparison","pval")
stat.mat$comparison <- c("ctl12-mnp12","ctl12-mnpz12","mnp12-mnpz12",
                         "ctl24-mnp24","ctl24-mnpz24","mnp24-mnpz24")

stat.mat$pval <- c(wilcox.test(x = z[z$var=="control.12mo","dist"], y = z[z$var=="MNP.12mo","dist"])$p.value,
                   wilcox.test(x = z[z$var=="control.12mo","dist"], y = z[z$var=="MNPZ.12mo","dist"])$p.value,
                   wilcox.test(x = z[z$var=="MNP.12mo","dist"], y = z[z$var=="MNPZ.12mo","dist"])$p.value,
                   wilcox.test(x = z[z$var=="control.24mo","dist"], y = z[z$var=="MNP.24mo","dist"])$p.value,
                   wilcox.test(x = z[z$var=="control.24mo","dist"], y = z[z$var=="MNPZ.24mo","dist"])$p.value,
                   wilcox.test(x = z[z$var=="MNP.24mo","dist"], y = z[z$var=="MNPZ.24mo","dist"])$p.value)

stat.mat$padj <- p.adjust(stat.mat$pval, method = "BH")

print(stat.mat)

          

###----- [Partial] Constrained Analysis of Principal Coordinates 
# via capscale, constrained by subject timepoint

ordcap <- ordinate(physeq = df.filt, method = "CAP", 
                   distance = DIST, ~timepoint)
plot_ordination(physeq = df.filt, ordination = ordcap, 
                type = "samples", color="timepoint", axes = c(1,2))

#constrain to take into account repeated /time measures 
capscale(formula = dist.bc ~ subject + arm, data = metadata)





###----- Run partial constrained analysis to view and remove effect of repeated measures

# constrain repeated (timepoint) measures to subject
# see how much variation there is within subjects
DIST <- "bray"
df.dist <- phyloseq::distance(physeq = df.filt, method = DIST)

cap <- capscale(formula = df.dist ~arm + center + group2 + Condition(subject), data=as.data.frame(as.matrix(sample_data(df.filt))))
print(cap)

# view ordination with effects of subject variation constrained
ord.cap <- ordinate(physeq = df.filt, method = "CAP", distance = DIST, ~group2 + arm + Condition(subject))
print(ord.cap)

cap_plot <- plot_ordination(physeq = df.filt, ordination = ord.cap, 
                            color = "timepoint", axes = c(1,2))
print(cap_plot)  


###------ End of document