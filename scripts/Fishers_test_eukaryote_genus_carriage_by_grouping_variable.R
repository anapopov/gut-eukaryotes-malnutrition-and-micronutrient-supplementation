#####--------------  FISHER'S TEST FOR EUKARYOTE GENUS CARRIGE IN GROUPS  --------------


###----- set directories and load data

# load libraries
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(rstatix)

# set working directory
setwd("eukaryote_carriage/")

# load phyloseq object with OTU table, taxonomy and sample metadata
load("microbes.rdata")




###----- set parameters & prepare dataset

# set tax level
TAXLEVEL <- "genus"  #OTU, genus

# set rarefaction depth
MINCOUNT <- 5
MINSAMPLES <- 0.05 #5% samples
MINSAMPLESIZE <- 500
AGEGROUP <- "all" #all, 12mo, 24mo

# select dataset
OBJ <- eukaryotes


# remove read counts for OTUs below minimum threshold
otu_table(OBJ)[otu_table(OBJ) < MINCOUNT] <- 0


# agglomerate if needed
if(TAXLEVEL!="OTU") {
  set.seed(1234)
  OBJ <- tax_glom(OBJ, taxrank = TAXLEVEL)
}


# remove "detection" < minimum threshold of reads
X <- as.data.frame(as.matrix(otu_table(OBJ)))
X[X<MINCOUNT] <- 0
otu_table(OBJ) <- otu_table(X, taxa_are_rows = TRUE)


# filter out low abundance taxa and samples with low read depth 
OBJ <- filter_taxa(OBJ, function(x) sum(x > MINCOUNT) >= (MINSAMPLES*length(x)), TRUE)

OBJ <- prune_samples(sample_sums(OBJ)>=MINSAMPLESIZE, OBJ)


# subset timepoint if needed
if(AGEGROUP %in% c("12mo","24mo")) {
  OBJ <- subset_samples(OBJ, timepoint %in% AGEGROUP)
}




###----- prepare results matrices

# extract grouping information
GROUPING <- c("nutrition.arm","arm","group2","center","timepoint")
X.meta <- as.data.frame(as.matrix(sample_data(OBJ)[,GROUPING]))

# prepare Fisher's test counts matrix for nutritional status
df.fisher.nut <- as.data.frame(matrix(nrow=2,ncol=2))
colnames(df.fisher.nut) <- c("neg","pos")
rownames(df.fisher.nut) <- c("normal","malnourished")

# prepare Fisher's test counts matrix for supplementation
df.fisher.arm <- as.data.frame(matrix(nrow=3,ncol=2))
colnames(df.fisher.arm) <- c("neg","pos")
rownames(df.fisher.arm) <- c("control","MNP","MNP+Zinc")

# prepare Fisher's test counts matrix for supplementation + nutritional status
df.fisher.nut.arm <- as.data.frame(matrix(nrow=6,ncol=2))
colnames(df.fisher.nut.arm) <- c("neg","pos")
rownames(df.fisher.nut.arm) <- c("normal.control","normal.MNP","normal.MNP+Zinc",
                                 "malnourished.control","malnourished.MNP","malnourished.MNP+Zinc")

# prepare Fisher's test counts matrix for location
df.fisher.loc <- as.data.frame(matrix(nrow=2,ncol=2))
colnames(df.fisher.loc) <- c("neg","pos")
rownames(df.fisher.loc) <- c("Rural","Urban")

# prepare Fisher's test counts matrix for age
df.fisher.age <- as.data.frame(matrix(nrow=2,ncol=2))
colnames(df.fisher.age) <- c("neg","pos")
rownames(df.fisher.age) <- c("12mo","24mo")



if(AGEGROUP %in% c("12mo","24mo")) {
  fisher.mat.all <- as.data.frame(matrix(nrow=ntaxa(OBJ),ncol=13))
  colnames(fisher.mat.all) <- c("taxon","p.nutrition","odds.nutrition","CI.nut","p.arm","p.nut.arm",
                                "p.loc","odds.loc","CI.loc","n.nut","n.arm","n.nut.arm","n.loc")
} else {
  fisher.mat.all <- as.data.frame(matrix(nrow=ntaxa(OBJ),ncol=17))
  colnames(fisher.mat.all) <- c("taxon","p.nutrition","odds.nutrition","CI.nut","p.arm","p.nut.arm",
                                "p.loc","odds.loc","CI.loc","p.age","odds.age","CI.age",
                                "n.nut","n.arm","n.nut.arm","n.loc","n.age")
}


# create list for pair-wise fisher results
fisher.arm.pair.list <- list()



###----- loop through organisms, preparing count contigency matrix and testing for differential carriage (Fisher's exact)

for(i in 1:ntaxa(OBJ)) {
  
  print(paste0(i, ": ", tax_table(OBJ)[i,"genus"]))
  
  X <- t(as.data.frame(otu_table(OBJ)[i,]))
  X.comb <- cbind(X.meta, X[rownames(X.meta),])
  colnames(X.comb)[ncol(X.comb)] <- "counts"
  
  # prepare Fisher's test matrix for all results
  if(AGEGROUP %in% c("12mo","24mo")) {
    fishers.tests <- as.data.frame(matrix(nrow = 4, ncol = 4))
    rownames(fishers.tests) <- c("nutritionalstatus","arm","nut.arm","location")
    colnames(fishers.tests) <- c("pval","oddsratio","taxon","n")
  } else {
    fishers.tests <- as.data.frame(matrix(nrow = 5, ncol = 4))
    rownames(fishers.tests) <- c("nutritionalstatus","arm","nut.arm","location","age")
    colnames(fishers.tests) <- c("pval","oddsratio","taxon","n")
  }
  
  
  #--- test by nutritional status
  mat <- df.fisher.nut
  
  # prepare counts matrix for Fisher's exact test
  mat["normal",] <- c(nrow(X.comb[X.comb$group2=="normal"&X.comb$counts==0,]),
                      nrow(X.comb[X.comb$group2=="normal"&X.comb$counts>0,]))
  mat["malnourished",] <- c(nrow(X.comb[X.comb$group2=="malnourished"&X.comb$counts==0,]),
                            nrow(X.comb[X.comb$group2=="malnourished"&X.comb$counts>0,]))
  
  # run test & add results to matrix
  fisher.nut <- fisher.test(mat)
  
  mat$n <- (mat$pos+mat$neg)
  mat$group <- rownames(mat)
  mat$group.n <- paste0(mat$group, " (", mat$n, ")")
  mat$group.n <- as.factor(mat$group.n)
  mat$group.n <- factor(mat$group.n, levels = levels(mat$group.n)[c(2,1)])
 
  fishers.tests["nutritionalstatus","pval"] <-fisher.nut$p.value
  fishers.tests["nutritionalstatus","oddsratio"] <-fisher.nut$estimate
  fishers.tests["nutritionalstatus","CI"] <- paste0(round(fisher.nut$conf.int[1], digits = 2),"-", round(fisher.nut$conf.int[2],digits = 2))
  fishers.tests["nutritionalstatus","taxon"] <- tax_table(OBJ)[i,"genus"]
  fishers.tests["nutritionalstatus","n"] <- paste0("norm. neg ", mat["normal","neg"], ", pos ", mat["normal","pos"],
                                                   "; mal. neg ", mat["malnourished","neg"], ", pos ", mat["malnourished","pos"])
  rm(mat)                                                        
 
   
  
  #---  test by supplementation arm
  mat <- df.fisher.arm
  
  # prepare counts matrix for Fisher's exact test
  mat["control",] <- c(nrow(X.comb[X.comb$arm=="control"&X.comb$counts==0,]),
                       nrow(X.comb[X.comb$arm=="control"&X.comb$counts>0,]))
  mat["MNP",] <- c(nrow(X.comb[X.comb$arm=="MNP"&X.comb$counts==0,]),
                   nrow(X.comb[X.comb$arm=="MNP"&X.comb$counts>0,]))
  mat["MNP+Zinc",] <- c(nrow(X.comb[X.comb$arm=="MNP+Zinc"&X.comb$counts==0,]),
                        nrow(X.comb[X.comb$arm=="MNP+Zinc"&X.comb$counts>0,]))
  
  # run Fisher's exact test
  fisher.arm <- fisher.test(mat, simulate.p.value = TRUE)
  
  # IF SIGNIFICANT, run Fisher's post-hoc pairwise test
  if(fisher.arm$p.value<0.05) {
    fisher.arm.pair <- as.data.frame(pairwise_fisher_test(mat, p.adjust.method = "BH", detailed = TRUE))
    fisher.arm.pair$organism <- as.character(tax_table(OBJ)[i,"genus"])
    fisher.arm.pair$otu <- as.character(rownames(tax_table(OBJ)[i,]))
    fisher.arm.pair.list[[i]] <- fisher.arm.pair
    rm(fisher.arm.pair)
  }
  
  mat$perc <- mat$pos/(mat$pos+mat$neg)
  mat$n <- (mat$pos+mat$neg)
  mat$group <- rownames(mat)
  mat$group.n <- paste0(mat$group, " (", mat$n, ")")
  
  fishers.tests["arm","pval"] <-fisher.arm$p.value
  fishers.tests["arm","oddsratio"] <- "NA.simulated-pval"
  fishers.tests["arm","taxon"] <- tax_table(OBJ)[i,"genus"]
  fishers.tests["arm","n"] <- paste0("ctl. neg ", mat["control","neg"], ", pos ", mat["control","pos"],
                                                   "; MNP. neg ", mat["MNP","neg"], ", pos ", mat["MNP","pos"],
                                     "; MNPZ. neg ", mat["MNP+Zinc","neg"], ", pos ", mat["MNP+Zinc","pos"])
  rm(mat)
  
  
  
  #--- test subgroups (nutritional status + supplementation arm)
  mat <- df.fisher.nut.arm
  
  # prepare counts matrix for Fisher's exact test
  mat["normal.control",] <- c(nrow(X.comb[X.comb$nutrition.arm=="normal.control"&X.comb$counts==0,]),
                              nrow(X.comb[X.comb$nutrition.arm=="normal.control"&X.comb$counts>0,]))
  mat["normal.MNP",] <- c(nrow(X.comb[X.comb$nutrition.arm=="normal.MNP"&X.comb$counts==0,]),
                          nrow(X.comb[X.comb$nutrition.arm=="normal.MNP"&X.comb$counts>0,]))
  mat["normal.MNP+Zinc",] <- c(nrow(X.comb[X.comb$nutrition.arm=="normal.MNP+Zinc"&X.comb$counts==0,]),
                               nrow(X.comb[X.comb$nutrition.arm=="normal.MNP+Zinc"&X.comb$counts>0,]))
  mat["malnourished.control",] <- c(nrow(X.comb[X.comb$nutrition.arm=="malnourished.control"&X.comb$counts==0,]),
                                    nrow(X.comb[X.comb$nutrition.arm=="malnourished.control"&X.comb$counts>0,]))
  mat["malnourished.MNP",] <- c(nrow(X.comb[X.comb$nutrition.arm=="malnourished.MNP"&X.comb$counts==0,]),
                                nrow(X.comb[X.comb$nutrition.arm=="malnourished.MNP"&X.comb$counts>0,]))
  mat["malnourished.MNP+Zinc",] <- c(nrow(X.comb[X.comb$nutrition.arm=="malnourished.MNP+Zinc"&X.comb$counts==0,]),
                                     nrow(X.comb[X.comb$nutrition.arm=="malnourished.MNP+Zinc"&X.comb$counts>0,]))
  
  # run test
  fisher.nutarm <- fisher.test(mat, simulate.p.value = TRUE)
  mat$perc <- mat$pos/(mat$pos+mat$neg)
  mat$n <- (mat$pos+mat$neg)
  mat$group <- rownames(mat)
  mat$group.n <- paste0(mat$group, " (", mat$n, ")")
  mat$group.n <- as.factor(mat$group.n)
  mat$group.n <- factor(x = mat$group.n, levels = levels(mat$group.n)[c(4,5,6,1,2,3)])
  
  fishers.tests["nut.arm","pval"] <- fisher.nutarm$p.value
  fishers.tests["nut.arm","oddsratio"] <- "NA.simulated-pval"
  fishers.tests["nut.arm","taxon"] <- tax_table(OBJ)[i,"genus"]
  fishers.tests["nut.arm","n"] <- paste0("norm.ctl neg ", mat["normal.control","neg"], ", pos ", mat["normal.control","pos"],
                                     "; norm.MNP neg ", mat["normal.MNP","neg"], ", pos ", mat["normal.MNP","pos"],
                                     "; norm.MNPZ neg ", mat["normal.MNP+Zinc","neg"], ", pos ", mat["normal.MNP+Zinc","pos"],
                                     "; mal.ctl neg ", mat["malnourished.control","neg"], ", pos ", mat["malnourished.control","pos"],
                                     "; mal.MNP neg ", mat["malnourished.MNP","neg"], ", pos ", mat["malnourished.MNP","pos"],
                                     "; mal.MNPZ neg ", mat["malnourished.MNP+Zinc","neg"], ", pos ", mat["malnourished.MNP+Zinc","pos"])


  #--- test by location
  mat <- df.fisher.loc
  
  # prepare counts matrix for Fisher's exact test
  mat["Rural",] <- c(nrow(X.comb[X.comb$center=="Rural"&X.comb$counts==0,]),
                      nrow(X.comb[X.comb$center=="Rural"&X.comb$counts>0,]))
  mat["Urban",] <- c(nrow(X.comb[X.comb$center=="Urban"&X.comb$counts==0,]),
                            nrow(X.comb[X.comb$center=="Urban"&X.comb$counts>0,]))
  
  mat <- mat[c(2,1),]
  
  # run test & add results to matrix
  fisher.loc <- fisher.test(mat)
  
  mat$n <- (mat$pos+mat$neg)
  mat$group <- rownames(mat)
  mat$group.n <- paste0(mat$group, " (", mat$n, ")")
  mat$group.n <- as.factor(mat$group.n)
  mat$group.n <- factor(mat$group.n, levels = levels(mat$group.n)[c(2,1)])
  
  fishers.tests["location","pval"] <-fisher.loc$p.value
  fishers.tests["location","oddsratio"] <-fisher.loc$estimate
  fishers.tests["location","CI"] <- paste0(round(fisher.loc$conf.int[1], digits = 2),"-", round(fisher.loc$conf.int[2],digits = 2))
  fishers.tests["location","taxon"] <- tax_table(OBJ)[i,"genus"]
  fishers.tests["location","n"] <- paste0("rural. neg ", mat["Rural","neg"], ", pos ", mat["Rural","pos"],
                                          "; urban. neg ", mat["Urban","neg"], ", pos ", mat["Urban","pos"])
  rm(mat)                                                        

   
  #--- test by age 
  if(!AGEGROUP %in% c("12mo","24mo")) {
    mat <- df.fisher.age
    
    # prepare counts matrix for Fisher's exact test
    mat["12mo",] <- c(nrow(X.comb[X.comb$timepoint=="12mo"&X.comb$counts==0,]),
                       nrow(X.comb[X.comb$timepoint=="12mo"&X.comb$counts>0,]))
    mat["24mo",] <- c(nrow(X.comb[X.comb$timepoint=="24mo"&X.comb$counts==0,]),
                       nrow(X.comb[X.comb$timepoint=="24mo"&X.comb$counts>0,]))
    
    # run test & add results to matrix
    fisher.age <- fisher.test(mat)
    
    mat$n <- (mat$pos+mat$neg)
    mat$group <- rownames(mat)
    mat$group.n <- paste0(mat$group, " (", mat$n, ")")
    mat$group.n <- as.factor(mat$group.n)
    mat$group.n <- factor(mat$group.n, levels = levels(mat$group.n)[c(2,1)])
    
    fishers.tests["age","pval"] <-fisher.age$p.value
    fishers.tests["age","oddsratio"] <- round(fisher.age$estimate, digits = 2)
    fishers.tests["age","CI"] <- paste0(round(fisher.age$conf.int[1], digits = 2),"-", round(fisher.age$conf.int[2],digits = 2))
    fishers.tests["age","taxon"] <- tax_table(OBJ)[i,"genus"]
    fishers.tests["age","n"] <- paste0("12mo. neg ", mat["12mo","neg"], ", pos ", mat["12mo","pos"],
                                            "; 24mo. neg ", mat["24mo","neg"], ", pos ", mat["24mo","pos"])
    rm(mat)                                                        
  }
  
  
  ###--- tabulate all results into single matrix
  fisher.mat.all[i,"taxon"] <- row.names(otu_table(OBJ)[i,])
  fisher.mat.all[i,"p.nutrition"] <- fishers.tests["nutritionalstatus","pval"]
  fisher.mat.all[i,"odds.nutrition"] <- fishers.tests["nutritionalstatus","oddsratio"]
  fisher.mat.all[i,"CI.nut"] <- fishers.tests["nutritionalstatus","CI"]
  fisher.mat.all[i,"p.arm"] <- fishers.tests["arm","pval"]
  fisher.mat.all[i,"p.nut.arm"] <- fishers.tests["nut.arm","pval"]
  fisher.mat.all[i,"p.loc"] <- fishers.tests["location","pval"]
  fisher.mat.all[i,"odds.loc"] <- fishers.tests["location","oddsratio"]
  fisher.mat.all[i,"CI.loc"] <- fishers.tests["location","CI"]
  fisher.mat.all[i,"n.nut"] <- fishers.tests["nutritionalstatus","n"]
  fisher.mat.all[i,"n.arm"] <- fishers.tests["arm","n"]
  fisher.mat.all[i,"n.nut.arm"] <- fishers.tests["nut.arm","n"]
  fisher.mat.all[i,"n.loc"] <- fishers.tests["location","n"]
  
  if(!AGEGROUP %in% c("12mo","24mo")) {
    fisher.mat.all[i,"p.age"] <- fishers.tests["age","pval"]
    fisher.mat.all[i,"odds.age"] <- fishers.tests["age","oddsratio"]
    fisher.mat.all[i,"CI.age"] <- fishers.tests["age","CI"]
    fisher.mat.all[i,"n.age"] <- fishers.tests["age","n"] 
  }
    
  rm(fishers.tests)
}



# adjust p values for multiple testing
fisher.mat.all$p.nut.adj <- p.adjust(fisher.mat.all$p.nutrition, method = "BH")
fisher.mat.all$p.arm.adj <- p.adjust(fisher.mat.all$p.arm, method = "BH")
fisher.mat.all$p.nut.arm.adj <- p.adjust(fisher.mat.all$p.nut.arm, method = "BH")
fisher.mat.all$p.loc.adj <- p.adjust(fisher.mat.all$p.loc, method = "BH")

if(!AGEGROUP %in% c("12mo","24mo")) {
  fisher.mat.all$p.age.adj <- p.adjust(fisher.mat.all$p.age, method = "BH")
}

# add organism genus name
fisher.mat.all$genus <- as.character(tax_table(OBJ)[fisher.mat.all$taxon,"genus"])

# prepare pairwise results matrix
fisher.arm.pair.list <- do.call(rbind, fisher.arm.pair.list)


### End of document