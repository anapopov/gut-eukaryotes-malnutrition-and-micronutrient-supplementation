### Co-occurence of bacterial taxa

#-- steps
# 1. load the data and normalize counts to common scale (min depth)
# 2. (within spieceasi) clr transformation
# 3. inverse covariance estimation along a lambda (sparsity) path
# 4. stability selection using the StARS criterion
# 5. evaluate performance

# lambda.min.ratio = the scaling factor that determines the minimum sparsity/lambda parameter
# rep.num = number of repetitions/subsamples for STARS; 99 or 999
# lambda = soft threshold for setting noisy correlations to zero. Rather than setting this 
#    to some arbitrary number, SPIEC-EASI uses a subsampling scheme (nlambda) to learn the 
#    threshold from data.


###----- load libraries
library(devtools)
install_github("zdk123/SpiecEasi")
install.packages("remotes")
#remotes::install_github("schuyler-smith/phyloschuyler") # helps re-ordering samples
library(phyloseq)
library(ape)
library(SpiecEasi)
library(ggplot2)
library(reshape2)
library(igraph)
library(Matrix)
library(dplyr)
library(tidyverse)



###----- import phyloseq object containing bacterial and eukaryotic OTU tables and taxonomies 
load("microbiome.rdata")



###----- rarefy dataset and filter low abundance taxa
MINSUBJECT <- 0.05 
MINCOUNT <- 5
RAREDEPTH.B <- 25000 # bacterial rarefaction depth
RAREDEPTH.E <- 1000 # eukaryotic rarefaction depth

RAREFY <- "yes"
FILTER <- "no"
TAXRANK <- "OTU"
TAXRANK.EUK <- "genus"


df.B <- bacteria # phyloseq object
df.E <- eukaryotes # phyloseq object


# reorder sample order to match
SAMPLEORDER <- rownames(sample_data(bacteria))

OTU.E <- otu_table(df.E)
META.E <- sample_data(df.E)

otu_table(df.E) <- OTU.E[,SAMPLEORDER]
sample_data(df.E) <- META.E[SAMPLEORDER,]



# rarefy to equal sampling depth and filter
if(RAREFY=="yes") {
  df.B <- rarefy_even_depth(physeq = df.B, sample.size = RAREDEPTH.B, rngseed = 1234)
  df.E <- rarefy_even_depth(physeq = df.E, sample.size = RAREDEPTH.E, rngseed = 1234)
}


# filter low counts
if(FILTER=="yes") {
  df.B.filt <- filter_taxa(df.B, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)
  df.E.filt <- filter_taxa(df.E, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)
} else if(FILTER=="no") {
  df.B.filt <- df.B
  df.E.filt <- df.E
}


# agglomerate to genus
if(TAXRANK=="genus") {
  df.B.filt <- tax_glom(physeq = df.B.filt, taxrank = "genus")
}

if(TAXRANK.EUK=="genus") {
  df.E.filt <- tax_glom(physeq = df.E.filt, taxrank = "genus")
}



###----- subset by age group

# bacteria
df.B.filt.12 <- subset_samples(physeq = df.B.filt, timepoint=="12mo")
df.B.filt.24 <- subset_samples(physeq = df.B.filt, timepoint=="24mo")

df.B.filt.12 <- filter_taxa(df.B.filt.12, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)
df.B.filt.24 <- filter_taxa(df.B.filt.24, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)

# eukaryotes
df.E.filt.12 <- subset_samples(physeq = df.E.filt, timepoint=="12mo")
df.E.filt.24 <- subset_samples(physeq = df.E.filt, timepoint=="24mo")

df.E.filt.12 <- filter_taxa(df.E.filt.12, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)
df.E.filt.24 <- filter_taxa(df.E.filt.24, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)


# filter identical sample sets
SAMPLES.12 <- intersect(rownames(sample_data(df.B.filt.12)), rownames(sample_data(df.E.filt.12)))
SAMPLES.24 <- intersect(rownames(sample_data(df.B.filt.24)), rownames(sample_data(df.E.filt.24)))

df.B.filt.12.p <- prune_samples(samples = SAMPLES.12, x = df.B.filt.12)
df.E.filt.12.p <- prune_samples(samples = SAMPLES.12, x = df.E.filt.12)

df.B.filt.24.p <- prune_samples(samples = SAMPLES.24, x = df.B.filt.24)
df.E.filt.24.p <- prune_samples(samples = SAMPLES.24, x = df.E.filt.24)




###----- choose age group to test, subset samples and filter
AGE <- "12mo" # 12mo, 24mo

if(AGE=="12mo") {
  bac <- df.B.filt.12.p
  euk <- df.E.filt.12.p
} else if(AGE=="24mo") {
  bac <- df.B.filt.24.p
  euk <- df.E.filt.24.p
}




###----- subset by nutritional status
bac.nor <- subset_samples(physeq = bac, group2=="normal")
euk.nor <- subset_samples(physeq = euk, group2=="normal")
bac.mal <- subset_samples(physeq = bac, group2=="malnourished")
euk.mal <- subset_samples(physeq = euk, group2=="malnourished")

# further subdivide by treatment arm
bac.nor.ctl <- subset_samples(physeq = bac.nor, arm=="control")
euk.nor.ctl <- subset_samples(physeq = euk.nor, arm=="control")
bac.nor.mnp <- subset_samples(physeq = bac.nor, arm=="MNP")
euk.nor.mnp <- subset_samples(physeq = euk.nor, arm=="MNP")
bac.nor.mnpz <- subset_samples(physeq = bac.nor, arm=="MNP+Zinc")
euk.nor.mnpz <- subset_samples(physeq = euk.nor, arm=="MNP+Zinc")

bac.mal.ctl <- subset_samples(physeq = bac.mal, arm=="control")
euk.mal.ctl <- subset_samples(physeq = euk.mal, arm=="control")
bac.mal.mnp <- subset_samples(physeq = bac.mal, arm=="MNP")
euk.mal.mnp <- subset_samples(physeq = euk.mal, arm=="MNP")
bac.mal.mnpz <- subset_samples(physeq = bac.mal, arm=="MNP+Zinc")
euk.mal.mnpz <- subset_samples(physeq = euk.mal, arm=="MNP+Zinc")

rm(bac.nor,bac.mal,euk.nor,euk.mal)


# filter out low abundance taxa
bac.nor.ctl <- filter_taxa(bac.nor.ctl, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)
euk.nor.ctl <- filter_taxa(euk.nor.ctl, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)
bac.nor.mnp <- filter_taxa(bac.nor.mnp, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)
euk.nor.mnp <- filter_taxa(euk.nor.mnp, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)
bac.nor.mnpz <- filter_taxa(bac.nor.mnpz, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)
euk.nor.mnpz <- filter_taxa(euk.nor.mnpz, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)

bac.mal.ctl <- filter_taxa(bac.mal.ctl, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)
euk.mal.ctl <- filter_taxa(euk.mal.ctl, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)
bac.mal.mnp <- filter_taxa(bac.mal.mnp, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)
euk.mal.mnp <- filter_taxa(euk.mal.mnp, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)
bac.mal.mnpz <- filter_taxa(bac.mal.mnpz, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)
euk.mal.mnpz <- filter_taxa(euk.mal.mnpz, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)




# filter identical samples
SAMPLES.nor.ctl <- intersect(rownames(sample_data(bac.nor.ctl)), rownames(sample_data(euk.nor.ctl)))
SAMPLES.mal.ctl <- intersect(rownames(sample_data(bac.mal.ctl)), rownames(sample_data(euk.mal.ctl)))
SAMPLES.nor.mnp <- intersect(rownames(sample_data(bac.nor.mnp)), rownames(sample_data(euk.nor.mnp)))
SAMPLES.mal.mnp <- intersect(rownames(sample_data(bac.mal.mnp)), rownames(sample_data(euk.mal.mnp)))
SAMPLES.nor.mnpz <- intersect(rownames(sample_data(bac.nor.mnpz)), rownames(sample_data(euk.nor.mnpz)))
SAMPLES.mal.mnpz <- intersect(rownames(sample_data(bac.mal.mnpz)), rownames(sample_data(euk.mal.mnpz)))


bac.nor.ctl.p <- prune_samples(samples = SAMPLES.nor.ctl, x = bac.nor.ctl)
euk.nor.ctl.p <- prune_samples(samples = SAMPLES.nor.ctl, x = euk.nor.ctl)
bac.nor.mnp.p <- prune_samples(samples = SAMPLES.nor.mnp, x = bac.nor.mnp)
euk.nor.mnp.p <- prune_samples(samples = SAMPLES.nor.mnp, x = euk.nor.mnp)
bac.nor.mnpz.p <- prune_samples(samples = SAMPLES.nor.mnpz, x = bac.nor.mnpz)
euk.nor.mnpz.p <- prune_samples(samples = SAMPLES.nor.mnpz, x = euk.nor.mnpz)

bac.mal.ctl.p <- prune_samples(samples = SAMPLES.mal.ctl, x = bac.mal.ctl)
euk.mal.ctl.p <- prune_samples(samples = SAMPLES.mal.ctl, x = euk.mal.ctl)
bac.mal.mnp.p <- prune_samples(samples = SAMPLES.mal.mnp, x = bac.mal.mnp)
euk.mal.mnp.p <- prune_samples(samples = SAMPLES.mal.mnp, x = euk.mal.mnp)
bac.mal.mnpz.p <- prune_samples(samples = SAMPLES.mal.mnpz, x = bac.mal.mnpz)
euk.mal.mnpz.p <- prune_samples(samples = SAMPLES.mal.mnpz, x = euk.mal.mnpz)



###---------------------------------------------------------------------------------------------------------
###----- test network parameters

# prep result matrix
params.df <- as.data.frame(matrix(ncol=5, nrow=22))
colnames(params.df) <- c("nlambda","lambda.min.ratio","stabilitythreshold","numedges","selectedlambda")
params.df$nlambda[1:5] <- 20
params.df$nlambda[6:10] <- 40
params.df$nlambda[11:15] <- 100
params.df$nlambda[16:20] <- 150
params.df$lambda.min.ratio[c(1,6,11,16)] <- 5e-2
params.df$lambda.min.ratio[c(2,7,12,17)] <- 1e-2
params.df$lambda.min.ratio[c(3,8,13,18)] <- 5e-3
params.df$lambda.min.ratio[c(4,9,14,19)] <- 1e-3
params.df$lambda.min.ratio[c(5,10,15,20)] <- 5e-4
params.df$nlambda[21:22] <- 100
params.df$lambda.min.ratio[21] <- 1e-1
params.df$lambda.min.ratio[22] <- 5e-1

# run tests
for(i in 1:nrow(params.df)) {
  print(paste0("nlambda ", params.df$nlambda[i], "; min lambda ", params.df$lambda.min.ratio[i]))
  
  test.network <- multi.spiec.easi(list(bac.nor.ctl.p, euk.nor.ctl.p), method='mb',
                                   sel.criterion='bstars',
                                   nlambda=params.df$nlambda[i],
                                   lambda.min.ratio=params.df$lambda.min.ratio[i],
                                   pulsar.params = list(thresh = 0.05, rep.num=99))
  
  params.df$stabilitythreshold[i] <- getStability(test.network) # stability threshold
  params.df$numedges[i] <- sum(getRefit(test.network))/2 # number of edges
  params.df$selectedlambda[i] <- getOptInd(test.network) # index of the selected lambda from provided lambda path
  
  print(params.df)
  rm(test.network)
}

# export results
print(params.df)




###---------------------------------------------------------------------------------------------------------
###----- SPIEC-EASI: calculate microbial networks, including cross-domain interactions

# choose parameters - max stability, balanced with time
NLAMBDA <- 100  # 100
MINLAMBDA <- 1e-03  # 1e-02



# nutritional status: normal, control
set.seed(1234)
se.bac.euk.nor.ctl <- multi.spiec.easi(list(bac.nor.ctl.p, euk.nor.ctl.p), method='mb', 
                                   sel.criterion='bstars', nlambda=NLAMBDA,
                                   lambda.min.ratio=MINLAMBDA, pulsar.params = list(thresh = 0.05, rep.num=99))

# nutritional status: malnourished, control
set.seed(1234)
se.bac.euk.mal.ctl <- multi.spiec.easi(list(bac.mal.ctl.p, euk.mal.ctl.p), method='mb', 
                                   sel.criterion='bstars', nlambda=NLAMBDA,
                                   lambda.min.ratio=MINLAMBDA, pulsar.params = list(thresh = 0.05, rep.num=99))



# nutritional status: normal, mnp
set.seed(1234)
se.bac.euk.nor.mnp <- multi.spiec.easi(list(bac.nor.mnp.p, euk.nor.mnp.p), method='mb', 
                                       sel.criterion='bstars', nlambda=NLAMBDA,
                                       lambda.min.ratio=MINLAMBDA, pulsar.params = list(thresh = 0.05, rep.num=99))

# nutritional status: malnourished, mnp
set.seed(1234)
se.bac.euk.mal.mnp <- multi.spiec.easi(list(bac.mal.mnp.p, euk.mal.mnp.p), method='mb', 
                                       sel.criterion='bstars', nlambda=NLAMBDA,
                                       lambda.min.ratio=MINLAMBDA, pulsar.params = list(thresh = 0.05, rep.num=99))



# nutritional status: normal, mnpz
set.seed(1234)
se.bac.euk.nor.mnpz <- multi.spiec.easi(list(bac.nor.mnpz.p, euk.nor.mnpz.p), method='mb', 
                                       sel.criterion='bstars', nlambda=NLAMBDA,
                                       lambda.min.ratio=MINLAMBDA, pulsar.params = list(thresh = 0.05, rep.num=99))

# nutritional status: malnourished, mnpz
set.seed(1234)
se.bac.euk.mal.mnpz <- multi.spiec.easi(list(bac.mal.mnpz.p, euk.mal.mnpz.p), method='mb', 
                                       sel.criterion='bstars', nlambda=NLAMBDA,
                                       lambda.min.ratio=MINLAMBDA, pulsar.params = list(thresh = 0.05, rep.num=99))



###---------------------------------------------------------------------------------------------------------
###----- check network stats
# number of nodes
# number of edges (degrees)
# degrees per node
# degree distribution (distribution of the number of degrees over nodes)


# see sample network stability (target = 0.05)
getStability(se.bac.euk.nor.ctl)
getStability(se.bac.euk.mal.ctl)
getStability(se.bac.euk.nor.mnp)
getStability(se.bac.euk.mal.mnp)
getStability(se.bac.euk.nor.mnpz)
getStability(se.bac.euk.mal.mnpz)



# number of total, positive and negative edges

# getOptBeta: get optimal beta matrix when MB neighborhood selection is run with StARS
# SymBeta: Symmetrize a beta (coefficient) matrix
betaMat.nor.ctl <- as.matrix(symBeta(getOptBeta(se.bac.euk.nor.ctl)))
betaMat.mal.ctl <- as.matrix(symBeta(getOptBeta(se.bac.euk.mal.ctl)))

betaMat.nor.mnp <- as.matrix(symBeta(getOptBeta(se.bac.euk.nor.mnp)))
betaMat.mal.mnp <- as.matrix(symBeta(getOptBeta(se.bac.euk.mal.mnp)))

betaMat.nor.mnpz <- as.matrix(symBeta(getOptBeta(se.bac.euk.nor.mnpz)))
betaMat.mal.mnpz <- as.matrix(symBeta(getOptBeta(se.bac.euk.mal.mnpz)))



total.nor.ctl <- length(betaMat.nor.ctl[betaMat.nor.ctl!=0])/2 
positive.nor.ctl <- length(betaMat.nor.ctl[betaMat.nor.ctl>0])/2 
negative.nor.ctl <- length(betaMat.nor.ctl[betaMat.nor.ctl<0])/2 

total.mal.ctl <- length(betaMat.mal.ctl[betaMat.mal.ctl!=0])/2 
positive.mal.ctl <- length(betaMat.mal.ctl[betaMat.mal.ctl>0])/2 
negative.mal.ctl <- length(betaMat.mal.ctl[betaMat.mal.ctl<0])/2 



total.nor.mnp <- length(betaMat.nor.mnp[betaMat.nor.mnp!=0])/2 
positive.nor.mnp <- length(betaMat.nor.mnp[betaMat.nor.mnp>0])/2 
negative.nor.mnp <- length(betaMat.nor.mnp[betaMat.nor.mnp<0])/2 

total.mal.mnp <- length(betaMat.mal.mnp[betaMat.mal.mnp!=0])/2 
positive.mal.mnp <- length(betaMat.mal.mnp[betaMat.mal.mnp>0])/2 
negative.mal.mnp <- length(betaMat.mal.mnp[betaMat.mal.mnp<0])/2 



total.nor.mnpz <- length(betaMat.nor.mnpz[betaMat.nor.mnpz!=0])/2 
positive.nor.mnpz <- length(betaMat.nor.mnpz[betaMat.nor.mnpz>0])/2 
negative.nor.mnpz <- length(betaMat.nor.mnpz[betaMat.nor.mnpz<0])/2 

total.mal.mnpz <- length(betaMat.mal.mnpz[betaMat.mal.mnpz!=0])/2 
positive.mal.mnpz <- length(betaMat.mal.mnpz[betaMat.mal.mnpz>0])/2 
negative.mal.mnpz <- length(betaMat.mal.mnpz[betaMat.mal.mnpz<0])/2 



# view stats for numbers of associations, positive and negative
print(paste0("total ", total.nor.ctl, "; positive ", positive.nor.ctl, "; negative ", negative.nor.ctl))
print(paste0("total ", total.mal.ctl, "; positive ", positive.mal.ctl, "; negative ", negative.mal.ctl))

print(paste0("total ", total.nor.mnp, "; positive ", positive.nor.mnp, "; negative ", negative.nor.mnp))
print(paste0("total ", total.mal.mnp, "; positive ", positive.mal.mnp, "; negative ", negative.mal.mnp))

print(paste0("total ", total.nor.mnpz, "; positive ", positive.nor.mnpz, "; negative ", negative.nor.mnpz))
print(paste0("total ", total.mal.mnpz, "; positive ", positive.mal.mnpz, "; negative ", negative.mal.mnpz))




###----- prepare stats table for export
stats <- as.data.frame(matrix(ncol=6,nrow=8))
rownames(stats) <- c("edges-total","edges-positive","edges-negative","network-stability","nodes.bac","nodes.euk","edges-per-node","nsamples")
colnames(stats) <- c("normal.ctl","malnourished.ctl","normal.mnp","malnourished.mnp","normal.mnpz","malnourished.mnpz")
stats["edges-total",] <- c(total.nor.ctl, total.mal.ctl, total.nor.mnp, total.mal.mnp, total.nor.mnpz, total.mal.mnpz)
stats["edges-positive",] <- c(positive.nor.ctl, positive.mal.ctl, positive.nor.mnp, positive.mal.mnp, positive.nor.mnpz, positive.mal.mnpz)
stats["edges-negative",] <- c(negative.nor.ctl, negative.mal.ctl, negative.nor.mnp, negative.mal.mnp, negative.nor.mnpz, negative.mal.mnpz)
stats["network-stability",] <- c(getStability(se.bac.euk.nor.ctl), getStability(se.bac.euk.mal.ctl),
                                 getStability(se.bac.euk.nor.mnp), getStability(se.bac.euk.mal.mnp),
                                 getStability(se.bac.euk.nor.mnpz), getStability(se.bac.euk.mal.mnpz))
stats["nodes.bac",] <- c(ntaxa(bac.nor.ctl.p),ntaxa(bac.mal.ctl.p),
                     ntaxa(bac.nor.mnp.p),ntaxa(bac.mal.mnp.p),
                     ntaxa(bac.nor.mnpz.p),ntaxa(bac.mal.mnpz.p))
stats["nodes.euk",] <-  c(ntaxa(euk.nor.ctl.p),ntaxa(euk.mal.ctl.p),
                          ntaxa(euk.nor.mnp.p),ntaxa(euk.mal.mnp.p),
                          ntaxa(euk.nor.mnpz.p),ntaxa(euk.mal.mnpz.p))
stats["edges-per-node",] <- stats["edges-total",]/(stats["nodes.bac",]+stats["nodes.euk",])
stats["nsamples",] <-c(nsamples(bac.nor.ctl.p),nsamples(bac.mal.ctl.p),
                       nsamples(bac.nor.mnp.p),nsamples(bac.mal.mnp.p),
                       nsamples(bac.nor.mnpz.p),nsamples(bac.mal.mnpz.p))

print(stats)
write.csv(x = stats, file = paste0("Network-stats_", TAXRANK, "_age-", AGE, "_nlambda_", NLAMBDA, "_minlambda_", MINLAMBDA, "_", Sys.Date(), "_SUBDIVIDED.csv"), row.names = TRUE)



###---------------------------------------------------------------------------------------------------------
###----- plot degree distribution - connectedness

# getRefit: get the optimal network & related structures when StARS is run.

# plot degree distribution - nutritional status (convert network to igraph format)
se.bac.euk.nor.ctl.ig <- adj2igraph(getRefit(se.bac.euk.nor.ctl), vertex.attr=list(name=c(taxa_names(bac.nor.ctl.p),taxa_names(euk.nor.ctl.p))))
se.bac.euk.mal.ctl.ig <- adj2igraph(getRefit(se.bac.euk.mal.ctl), vertex.attr=list(name=c(taxa_names(bac.mal.ctl.p),taxa_names(euk.mal.ctl.p))))
se.bac.euk.nor.mnp.ig <- adj2igraph(getRefit(se.bac.euk.nor.mnp), vertex.attr=list(name=c(taxa_names(bac.nor.mnp.p),taxa_names(euk.nor.mnp.p))))
se.bac.euk.mal.mnp.ig <- adj2igraph(getRefit(se.bac.euk.mal.mnp), vertex.attr=list(name=c(taxa_names(bac.mal.mnp.p),taxa_names(euk.mal.mnp.p))))
se.bac.euk.nor.mnpz.ig <- adj2igraph(getRefit(se.bac.euk.nor.mnpz), vertex.attr=list(name=c(taxa_names(bac.nor.mnpz.p),taxa_names(euk.nor.mnpz.p))))
se.bac.euk.mal.mnpz.ig <- adj2igraph(getRefit(se.bac.euk.mal.mnpz), vertex.attr=list(name=c(taxa_names(bac.mal.mnpz.p),taxa_names(euk.mal.mnpz.p))))

se.bac.euk.nor.ctl.ig.degree.distribution <- degree.distribution(se.bac.euk.nor.ctl.ig)
se.bac.euk.mal.ctl.ig.degree.distribution <- degree.distribution(se.bac.euk.mal.ctl.ig)
se.bac.euk.nor.mnp.ig.degree.distribution <- degree.distribution(se.bac.euk.nor.mnp.ig)
se.bac.euk.mal.mnp.ig.degree.distribution <- degree.distribution(se.bac.euk.mal.mnp.ig)
se.bac.euk.nor.mnpz.ig.degree.distribution <- degree.distribution(se.bac.euk.nor.mnpz.ig)
se.bac.euk.mal.mnpz.ig.degree.distribution <- degree.distribution(se.bac.euk.mal.mnpz.ig)




# extract degree distributions] values
x1 <- as.data.frame(cbind(degree=0:(length(se.bac.euk.nor.ctl.ig.degree.distribution)-1), freq.nor.ctl=se.bac.euk.nor.ctl.ig.degree.distribution))
x2 <- as.data.frame(cbind(degree=0:(length(se.bac.euk.mal.ctl.ig.degree.distribution)-1), freq.mal.ctl=se.bac.euk.mal.ctl.ig.degree.distribution))

x3 <- as.data.frame(cbind(degree=0:(length(se.bac.euk.nor.mnp.ig.degree.distribution)-1), freq.nor.mnp=se.bac.euk.nor.mnp.ig.degree.distribution))
x4 <- as.data.frame(cbind(degree=0:(length(se.bac.euk.mal.mnp.ig.degree.distribution)-1), freq.mal.mnp=se.bac.euk.mal.mnp.ig.degree.distribution))

x5 <- as.data.frame(cbind(degree=0:(length(se.bac.euk.nor.mnpz.ig.degree.distribution)-1), freq.nor.mnpz=se.bac.euk.nor.mnpz.ig.degree.distribution))
x6 <- as.data.frame(cbind(degree=0:(length(se.bac.euk.mal.mnpz.ig.degree.distribution)-1), freq.mal.mnpz=se.bac.euk.mal.mnpz.ig.degree.distribution))


# join into dataframe
deg.dist <- list(x1, x2, x3, x4, x5, x6) %>% reduce(full_join, by = "degree")
deg.dist.melt <- melt(deg.dist, id.vars = "degree")
deg.dist.melt <- deg.dist.melt[!is.na(deg.dist.melt$value),]


# define line types for plot
LINES <- c("freq.nor.ctl" = "solid", "freq.nor.mnp" = "solid", "freq.nor.mnpz" = "solid", 
           "freq.mal.ctl" = "dashed", "freq.mal.mnp" = "dashed", "freq.mal.mnpz" = "dashed")

# plot on log-log scale (both axes as log scales), to see if scale-free network (if follows Power Law)
p.deg.log <- ggplot(deg.dist.melt[deg.dist.melt$degree!=0,], aes(x=degree, y = value, col=variable, linetype = factor(variable))) + 
  geom_line() + theme_classic() + 
  scale_x_log10() + scale_y_log10() + annotation_logticks() +
  xlab("degree (Log)") + ylab("frequency (Log)") + 
  ggtitle(paste0("Degree Distributions, ", AGE)) +
  scale_color_manual(values = c("black", "black", "seagreen3", "seagreen3", "magenta3","magenta3")) +
  scale_linetype_manual(values = LINES)

print(p.deg.log)

# export
pdf(file = paste0("Degree_distribution_nut-arm_", TAXRANK, "_", AGE, "_nlambda-", NLAMBDA, "_minlambda-", MINLAMBDA, "_LOG-scale_", Sys.Date(), ".pdf"), width = 5, height = 3)
print(p.deg.log)
dev.off()


# plot on regular scale
p.deg <- ggplot(deg.dist.melt[deg.dist.melt$degree!=0,], aes(x=degree, y = value, col=variable, linetype = factor(variable))) + 
  geom_line() + theme_classic() + 
  xlab("degree") + ylab("frequency") + 
  ggtitle(paste0("Degree Distributions, ", AGE)) +
  scale_color_manual(values = c("black", "black", "seagreen3", "seagreen3", "magenta3","magenta3")) +
  scale_linetype_manual(values = LINES)

print(p.deg)

# export
pdf(file = paste0("Degree_distribution_nut-arm_", TAXRANK, "_", AGE, "_nlambda-", NLAMBDA, "_minlambda-", MINLAMBDA, ".pdf"), width = 4, height = 3)
print(p.deg)
dev.off()



# check stats between distributions using Wilcoxon rank sum
wilcox.test(x = deg.dist.melt[deg.dist.melt$variable=="freq.nor.ctl", "value"],
            y = deg.dist.melt[deg.dist.melt$variable=="freq.mal.ctl", "value"], 
            paired =  T)

data <- cbind(norctl=deg.dist.melt[deg.dist.melt$variable=="freq.nor.ctl",c("degree","value")],
      malctl=deg.dist.melt[deg.dist.melt$variable=="freq.mal.ctl","value"])
rownames(data) <- paste0("d",data$norctl.degree)
data$norctl.degree <- NULL

# check stats with chi-squared test
chisq.test(data)





#----- Node (vertex) betweenness and no. node degrees
# node with higher betweenness centrality would have more control over the 
# network, because more information will pass through that node. 

graphs <- c("se.bac.euk.nor.ctl.ig", "se.bac.euk.mal.ctl.ig",
            "se.bac.euk.nor.mnp.ig", "se.bac.euk.mal.mnp.ig",
            "se.bac.euk.nor.mnpz.ig", "se.bac.euk.mal.mnpz.ig")

taxtable <- rbind(as.data.frame(as.matrix(tax_table(bacteria)[,c("phylum","class","genus","species")])),as.data.frame(as.matrix(tax_table(microbes.nocontrol.dichottree)[,c("phylum","class","genus","species")])))
taxtable$phylum <- gsub(pattern = "^.*__", replacement = "", taxtable$phylum)
taxtable$class <- gsub(pattern = "^.*__", replacement = "", taxtable$class)
taxtable$genus <- gsub(pattern = "^.*__", replacement = "", taxtable$genus)
taxtable$species <- gsub(pattern = "^.*__", replacement = "", taxtable$species)

btw.deg.list <- list()

for(i in 1:length(graphs)) {
  ig <- mget(x = graphs[i], envir = globalenv())[[1]]
  
  
  
  # calculate betweenness
  btw <- as.data.frame(betweenness(ig))
  colnames(btw) <- "betweenness"
  btw$otu <- rownames(btw)
  
  
  
  # tabulate no. degrees/node
  deg <- as.data.frame(degree(ig, v = V(ig), normalized = FALSE))
  colnames(deg) <- "degree"
  deg$otu <- rownames(deg)
  deg.norm <- as.data.frame(degree(ig, v = V(ig), normalized = TRUE))
  colnames(deg.norm) <- "degree.normalized"
  deg.norm$otu <- rownames(deg.norm)
  
  btw.deg <- merge(btw, deg, by = "otu")
  btw.deg.norm <- merge(btw.deg, deg.norm, by = "otu")
  rm(btw, deg, deg.norm, btw.deg)
  

  # add taxonomy information
  rownames(btw.deg.norm) <- btw.deg.norm$otu
  btw.deg.norm[,c("phylum","class","genus","species")] <- "NA"
  btw.deg.norm <- btw.deg.norm[order(btw.deg.norm$betweenness, decreasing = TRUE), ]
  
  btw.deg.norm[rownames(btw.deg.norm), "phylum"] <- taxtable[rownames(btw.deg.norm), "phylum"]
  btw.deg.norm[rownames(btw.deg.norm), "class"] <- taxtable[rownames(btw.deg.norm), "class"]
  btw.deg.norm[rownames(btw.deg.norm), "genus"] <- taxtable[rownames(btw.deg.norm), "genus"]
  btw.deg.norm[rownames(btw.deg.norm), "species"] <- taxtable[rownames(btw.deg.norm), "species"]
  
  
  
  # add rankings for betweeness, degree and closeness centralities
  btw.deg.norm$rank.btw <- rank(x = -(btw.deg.norm$betweenness), na.last = TRUE, ties.method = "max")
  btw.deg.norm$rank.deg <- rank(x = -(btw.deg.norm$degree), na.last = TRUE, ties.method = "max")
  #btw.deg.norm$rank.cls <- rank(x = -(btw.deg.norm$closeness), na.last = TRUE, ties.method = "max")
  
  # add ranking based on betweeness AND degree
  btw.deg.norm$rank.btw.deg <- rank(x = (btw.deg.norm$rank.btw + btw.deg.norm$rank.deg), na.last = TRUE, ties.method = "max")
  btw.deg.norm <- btw.deg.norm[order(btw.deg.norm$rank.btw.deg, decreasing = FALSE),]
  
  print(btw.deg.norm)
  
  btw.deg.list[[i]] <- btw.deg.norm
  names(btw.deg.list)[[i]] <- graphs[i]
  
  
  
  # export data
  write.csv(x = btw.deg.norm, file = paste0("Network-degree-betweenness_", graphs[i], "_", TAXRANK, "_", AGE, "_nlambda_", NLAMBDA, "_minlambda_", MINLAMBDA, ".csv"), row.names = TRUE)
  rm(btw.deg.norm)
}



###----- plot betweenness distribution
btw.deg.list

df <- do.call(rbind, btw.deg.list)[,c("degree","betweenness")]
df$type <- gsub(pattern = "\\.ig.*", replacement = "", rownames(df))
df$type <- gsub(pattern = "se.bac.euk.", replacement = "", df$type)
rownames(df) <- NULL


df$type <- factor(df$type, levels = c("nor.ctl", "nor.mnp", "nor.mnpz",
                                      "mal.ctl", "mal.mnp", "mal.mnpz"))

p <- ggplot(df, aes(x = type, y = betweenness, fill = type)) +
  geom_boxplot(lwd=0.4, fatten = 1, outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(), alpha=0.8, size = 0.2) + 
  theme_classic() + scale_y_log10() + annotation_logticks(outside = TRUE) +
  labs(x="", y = "betweenness", title = paste0("Betweenness, ", AGE)) +
  scale_fill_manual(values=c("grey40","seagreen2","magenta2","grey40","seagreen2","magenta2"))

pdf(file = paste0("Betweenness_centrality_distribution_", AGE, "_nlambda_", NLAMBDA, "_minlambda_", MINLAMBDA, ".pdf"), width = 4, height = 2.5)
print(p)
dev.off()



# check pairwise stats using Wilcoxon
stats.btw <- as.data.frame(matrix(ncol = 2, nrow = 9), stringsAsFactors = FALSE)
colnames(stats.btw) <- c("comparison","Wilcoxon.p.val")

stats.btw$comparison <- c("nor.ctl-mal.ctl", "nor.mnp-mal.mnp", "nor.mnpz-mal.mnpz",
                          "nor.ctl-nor.mnp", "nor.ctl-nor.mnpz", "nor.mnp-nor.mnpz",
                          "mal.ctl-mal.mnp", "mal.ctl-mal.mnpz", "mal.mnp-mal.mnpz")

stats.btw[1,"Wilcoxon.p.val"] <- wilcox.test(x = df[df$type=="nor.ctl","betweenness"], y = df[df$type=="mal.ctl","betweenness"], paired = FALSE)$p.value
stats.btw[2,"Wilcoxon.p.val"] <- wilcox.test(x = df[df$type=="nor.mnp","betweenness"], y = df[df$type=="mal.mnp","betweenness"], paired = FALSE)$p.value
stats.btw[3,"Wilcoxon.p.val"] <- wilcox.test(x = df[df$type=="nor.mnpz","betweenness"], y = df[df$type=="mal.mnpz","betweenness"], paired = FALSE)$p.value

stats.btw[4,"Wilcoxon.p.val"] <- wilcox.test(x = df[df$type=="nor.ctl","betweenness"], y = df[df$type=="nor.mnp","betweenness"], paired = FALSE)$p.value
stats.btw[5,"Wilcoxon.p.val"] <- wilcox.test(x = df[df$type=="nor.ctl","betweenness"], y = df[df$type=="nor.mnpz","betweenness"], paired = FALSE)$p.value
stats.btw[6,"Wilcoxon.p.val"] <- wilcox.test(x = df[df$type=="nor.mnp","betweenness"], y = df[df$type=="nor.mnpz","betweenness"], paired = FALSE)$p.value

stats.btw[7,"Wilcoxon.p.val"] <- wilcox.test(x = df[df$type=="mal.ctl","betweenness"], y = df[df$type=="mal.mnp","betweenness"], paired = FALSE)$p.value
stats.btw[8,"Wilcoxon.p.val"] <- wilcox.test(x = df[df$type=="mal.ctl","betweenness"], y = df[df$type=="mal.mnpz","betweenness"], paired = FALSE)$p.value
stats.btw[9,"Wilcoxon.p.val"] <- wilcox.test(x = df[df$type=="mal.mnp","betweenness"], y = df[df$type=="mal.mnpz","betweenness"], paired = FALSE)$p.value

# correct for multiple testing
stats.btw$p.adj <- p.adjust(p = stats.btw$Wilcoxon.p.val, method = "BH")

write.csv(x = stats.btw, file = paste0("Wilcoxon_pairwise_stats_betweenness_nlambda-", NLAMBDA, "_minlambda-", MINLAMBDA, "_", Sys.Date(), ".csv"), row.names = FALSE)




###----- plot degree violin plots
p <- ggplot(df, aes(x = type, y = degree, fill = type)) +
  geom_violin(lwd=0.3, draw_quantiles = c(0.25,0.5,0.75)) + 
  geom_point(position=position_jitterdodge(jitter.width = 1, jitter.height = 0.05), alpha=0.8, size = 0.03) + 
  theme_classic() +
  labs(x="", y = "degree", title = paste0("Degree, ", AGE)) +
  scale_fill_manual(values=c("grey40","seagreen2","magenta2","grey40","seagreen2","magenta2"))

print(p)

pdf(file = paste0("Degree_centrality_violinplots_", AGE, "_nlambda_", NLAMBDA, "_minlambda_", MINLAMBDA, ".pdf"), width = 4, height = 2.5)
print(p)
dev.off()




###----- prepare OTU ranks table by betweenness and degree centrality
btw.deg.list.red <- btw.deg.list

# create new list with subset of columns
btw.deg.list.red <- lapply(btw.deg.list.red, function(x) x%>% select(otu,rank.deg,rank.btw,rank.btw.deg))

# change column names
for(i in names(btw.deg.list.red)) {
  print(i)
  colnames(btw.deg.list.red[[i]])[-1] <- paste0(i, ".", colnames(btw.deg.list.red[[i]])[-1])
}

# join into dataframe
ranks <- btw.deg.list.red %>% reduce(full_join, by = "otu")

# add taxonomy information
rownames(ranks) <- ranks$otu

ranks[,c("phylum","class","genus","species")] <- "NA"
ranks[rownames(ranks), "phylum"] <- taxtable[rownames(ranks), "phylum"]
ranks[rownames(ranks), "class"] <- taxtable[rownames(ranks), "class"]
ranks[rownames(ranks), "genus"] <- taxtable[rownames(ranks), "genus"]
ranks[rownames(ranks), "species"] <- taxtable[rownames(ranks), "species"]

colnames(ranks) <- gsub(pattern = "se.bac.euk.", replacement = "", colnames(ranks))
colnames(ranks) <- gsub(pattern = ".ig", replacement = "", colnames(ranks))

write.csv(x = ranks, file = "OTU_ranking_by_degree_betweenness.csv", row.names = FALSE)



library(pheatmap)
library(RColorBrewer)

RANKTYPE <- c("rank.deg", "rank.btw$", "rank.btw.deg")
maxNA <- 0
CLUSTMETHOD <- "complete" #"ward.D2", "complete"

for(i in 1:length(RANKTYPE)) {
  mat <- ranks[,grep(pattern = RANKTYPE[i], colnames(ranks))]
  
  
  mat <- mat[rowSums(is.na(mat)) <=maxNA,] # remove rows with more than 2 NA
  
  p <- pheatmap(mat = mat,
                cluster_rows = TRUE, cluster_cols = TRUE,
                color = colorRampPalette(rev(brewer.pal(n = 7, name = "BuPu")))(100), #"BuPu", "YlOrRd"
                na_col = "black", #"white", "black"
                clustering_method = CLUSTMETHOD,
                border_color = NA,
                fontsize = 7)
  print(p)
  
  pdf(file = paste0("Heatmap_", RANKTYPE[i], "_maxNA-", maxNA, "_cluster-method-", CLUSTMETHOD, ".pdf"), width = 3, height = 8)
  print(p)
  dev.off()
}



#--------------------------------------------------------------------------------------
#---- graph network

TYPES <- c("nor.ctl","mal.ctl",
           "nor.mnp","mal.mnp",
           "nor.mnpz","mal.mnpz")

graphs <- c("se.bac.euk.nor.ctl.ig", "se.bac.euk.mal.ctl.ig",
            "se.bac.euk.nor.mnp.ig", "se.bac.euk.mal.mnp.ig",
            "se.bac.euk.nor.mnpz.ig", "se.bac.euk.mal.mnpz.ig")

networks <- c("se.bac.euk.nor.ctl", "se.bac.euk.mal.ctl",
            "se.bac.euk.nor.mnp", "se.bac.euk.mal.mnp",
            "se.bac.euk.nor.mnpz", "se.bac.euk.mal.mnpz")

TAXA.bac <- c("bac.nor.ctl.p","bac.mal.ctl.p",
              "bac.nor.mnp.p","bac.mal.mnp.p",
              "bac.nor.mnpz.p","bac.mal.mnpz.p")

TAXA.euk <-  c("euk.nor.ctl.p","euk.mal.ctl.p",
               "euk.nor.mnp.p","euk.mal.mnp.p",
               "euk.nor.mnpz.p","euk.mal.mnpz.p")

# choose dataset
for(i in 1:length(TYPES)) {

  print(TYPES[i])
  
  NETWORK <- mget(x = networks[i], envir = globalenv())[[1]]
  
  BAC <- mget(x = TAXA.bac[i], envir = globalenv())[[1]]
  EUK <- mget(x = TAXA.euk[i], envir = globalenv())[[1]]
  
  GRAPH <- adj2igraph(getRefit(NETWORK), vertex.attr=list(name=c(taxa_names(BAC),taxa_names(EUK))))
  names <- c(taxa_names(BAC),taxa_names(EUK))


  
  mat <- as.data.frame(as.matrix(symBeta(getOptBeta(NETWORK))))
  rownames(mat) <- names
  colnames(mat) <- names

  # prepare edge colors based on positive or negative interaction
  edges <- E(GRAPH)
  edge.colors=c()
  
  for(e.index in 1:length(edges)){
    
    adj.nodes <- ends(GRAPH,edges[e.index])
    xindex <- adj.nodes[1]
    yindex <- adj.nodes[2]
    beta <- mat[xindex,yindex]
    
    if(beta>0){
      edge.colors=append(edge.colors,"grey50")
    } else if(beta<0){
      edge.colors=append(edge.colors,"steelblue2")
    }
  }
  
  E(GRAPH)$color=edge.colors
  
  rm(edges,edge.colors,names)

  # prepare edges for Cytoscape export
  mat.melt <- melt(as.matrix(mat))
  mat.melt <- mat.melt[mat.melt$value!=0,]
  # add taxonomy info
  temp.taxa <- rbind(tax_table(BAC)[,c("phylum","genus")],tax_table(EUK)[,c("phylum","genus")])
  mat.melt$Var1.phylum <- temp.taxa[mat.melt$Var1,"phylum"]
  mat.melt$Var2.phylum <- temp.taxa[mat.melt$Var2,"phylum"]
  mat.melt$Var1.genus <- temp.taxa[mat.melt$Var1,"genus"]
  mat.melt$Var2.genus <- temp.taxa[mat.melt$Var2,"genus"]
  
  write.csv(x = mat.melt, file = paste0("Edge-table_for_Cytoscape_", networks[i], "_", AGE, ".csv"), row.names = FALSE)

  interaction.types <- as.data.frame(matrix(nrow = 3, ncol = 2), stringsAsFactors = FALSE)
  interaction.types$V1 <- c("bacterium-bacterium","bacterium-eukaryote","eukaryote-eukaryote")
  
  X <- mat.melt[grep("^Otu", mat.melt$Var1),]
  Y <- X[grep("^Otu", X$Var2),]
  X <- X[grep("motu|jotu", X$Var2),]
  Z <- mat.melt[grep("motu|jotu", mat.melt$Var1),]
  Z <- Z[grep("motu|jotu", Z$Var2),]
  interaction.types$V2 <- c(nrow(Y)/2, nrow(X), nrow(Z)/2)
  
  write.csv(x = interaction.types, file = paste0("Interaction-type_summary_", networks[i], "_", AGE, ".csv"), row.names = FALSE)
  rm(mat,mat.melt,interaction.types,X,Y,Z)

  
  
  # prepare colour scheme vector based on phyla
  otu.names <- colnames(NETWORK$select$est$data)
  otu.phyla <- c(tax_table(BAC)[,"phylum"],tax_table(EUK)[,"phylum"])
  otu.phyla.num <- length(levels(as.factor(otu.phyla)))
  
  
  # check phylum types represented
  otu.phyla <- gsub(pattern = "p__", replacement = "", otu.phyla)
  otu.phyla.col <- otu.phyla
  
  # set colours manually: bacteria
  otu.phyla.col <- gsub(pattern = "^.*Actinobacteria", replacement = "khaki2", otu.phyla.col)
  otu.phyla.col <- gsub(pattern = "^.*Firmicutes", replacement = "khaki2", otu.phyla.col)
  otu.phyla.col <- gsub(pattern = "^.*Bacteroidetes", replacement = "khaki2", otu.phyla.col)
  otu.phyla.col <- gsub(pattern = "^.*Proteobacteria", replacement = "khaki2", otu.phyla.col)
  otu.phyla.col <- gsub(pattern = "^.*Verrucomicrobia", replacement = "khaki2", otu.phyla.col)
  otu.phyla.col <- gsub(pattern = "^.*Cyanobacteria/Chloroplast", replacement = "khaki2", otu.phyla.col)
  
  # set colours manually: eukaryote
  otu.phyla.col <- gsub(pattern = "Alveolata", replacement = "indianred2", otu.phyla.col)
  otu.phyla.col <- gsub(pattern = "Apicomplexa", replacement = "indianred2", otu.phyla.col)
  otu.phyla.col <- gsub(pattern = "Cercozoa", replacement = "indianred2", otu.phyla.col)
  otu.phyla.col <- gsub(pattern = "Ciliophora", replacement = "indianred2", otu.phyla.col)
  otu.phyla.col <- gsub(pattern = "Parabasalia", replacement = "indianred2", otu.phyla.col)
  otu.phyla.col <- gsub(pattern = "Syndiniales", replacement = "indianred2", otu.phyla.col)
  otu.phyla.col <- gsub(pattern = "Microsporidia", replacement = "indianred2", otu.phyla.col)
  otu.phyla.col <- gsub(pattern = "Amoebozoa|Discosea", replacement = "indianred2", otu.phyla.col)
  otu.phyla.col <- gsub(pattern = "Gracilipodida", replacement = "indianred2", otu.phyla.col)
  otu.phyla.col <- gsub(pattern = "Stramenopiles", replacement = "indianred2", otu.phyla.col)
  
  otu.phyla.col <- gsub(pattern = "Platyhelminthes", replacement = "darkred", otu.phyla.col)
  
  otu.phyla.col <- gsub(pattern = "Ascomycota", replacement = "lightsteelblue4", otu.phyla.col) #"darkolivegreen3"
  otu.phyla.col <- gsub(pattern = "Basidiomycota", replacement = "lightsteelblue4", otu.phyla.col)
  otu.phyla.col <- gsub(pattern = "Mucoromycota", replacement = "lightsteelblue4", otu.phyla.col)
  otu.phyla.col <- gsub(pattern = "Cryptomycota", replacement = "lightsteelblue4", otu.phyla.col)
  
  
  #----- prepare node sizes proportional to betweeness values
  btw <- read.csv(file = paste0("Network-degree-betweenness_", graphs[i], "_", TAXRANK, "_", AGE, "_nlambda_", NLAMBDA, "_minlambda_", MINLAMBDA, ".csv"), row.names = 1, stringsAsFactors = FALSE)
  otu.sizes <- btw[c(rownames(otu_table(BAC)),rownames(otu_table(EUK))),"betweenness"]/130+3 # scale down, and add 0.01 for nodes with no betweenness
  rm(btw)
  
  
  #----- prepare vector of genus names
  otu.genus <- c(tax_table(BAC)[,"genus"],tax_table(EUK)[,"genus"])
  otu.genus <- gsub(pattern = "g__", replacement = "", otu.genus)
  
  
  #----- prepare legend vectors manually
  legend.names <- c("Actinobacteria", "Firmicutes", "Bacteroidetes", "Proteobacteria","Verrucomicrobia", "Cyanobacteria/Chloroplast",
                    "Alveolata","Apicomplexa","Cercozoa","Ciliophora","Parabasalia",
                    "Syndiniales","Microsporidia","Amoebozoa|Discosea","Gracilipodida","Stramenopiles",
                    "Platyhelminthes",
                    "Ascomycota", "Basidiomycota", "Mucoromycota", "Cryptomycota")
  
  legend.cols <- c(rep("khaki2",6), 
                   rep("indianred2",10), "darkred",
                   rep("lightsteelblue4",4))
  
  # extract only those in current network
  subset <- which(legend.names %in% otu.phyla)
  legend.names <- legend.names[subset]
  legend.cols <- legend.cols[subset]
  
  
  #----- plot network with genus names
  
  # change label size
  V(GRAPH)$label.cex = 0.6
  
  # plot with force-based algorithm layout (Fruchterman and Reingold)
  # repelling force calculated only between vertices that are closer to each other than a limit
  GRAPH.coord <- layout.fruchterman.reingold(GRAPH)
  
  pdf(file = paste0("Network_", networks[i], "_", AGE, "_nlambda_", NLAMBDA, "_minlambda_", MINLAMBDA, "_stability_", round(getStability(NETWORK), digits = 3), ".pdf"), width = 10, height = 10)

  set.seed(1234)
  plot.igraph(GRAPH, layout = GRAPH.coord,
              # nodes
              vertex.color=otu.phyla.col,
              vertex.frame.color=otu.phyla.col,
              vertex.size=otu.sizes,
              main=paste0(AGE, ", ", networks[i]),
              vertex.label.color="black", vertex.label.dist=1,
              vertex.label.font=3, vertex.label.degree=-pi/3,
              # edges
              edge.width=2, edge.lty="solid")
  
  legend("bottomleft", legend=legend.names  , col = legend.cols, pch=20, 
         pt.cex = 1.5, cex = 0.8, horiz = FALSE, inset = c(-0.05, 0.02))
  
  dev.off()

  
  
  #----- plot network with OTU names
  
  # change label size
  V(GRAPH)$label.cex = 0.6
  
  # plot with force-based algorithm layout (Fruchterman and Reingold)
  # repelling force calculated only between vertices that are closer to each other than a limit
  GRAPH.coord <- layout.fruchterman.reingold(GRAPH)
  
  pdf(file = paste0("Genus-name_OTU-Network_", networks[i], "_", AGE, "_nlambda_", NLAMBDA, "_minlambda_", MINLAMBDA, "_stability_", round(getStability(NETWORK), digits = 3), ".pdf"), width = 10, height = 10)
  
  set.seed(1234)
  plot.igraph(GRAPH, layout = GRAPH.coord,
              # nodes
              vertex.color=otu.phyla.col,
              vertex.frame.color=otu.phyla.col,
              vertex.size=otu.sizes,
              main=paste0(AGE, ", ", networks[i]),
              vertex.label=otu.genus,
              vertex.label.color="black", vertex.label.dist=1,
              vertex.label.font=3, vertex.label.degree=-pi/3,
              # edges
              edge.width=2, edge.lty="solid")
  
  legend("bottomleft", legend=legend.names  , col = legend.cols, pch=20, 
         pt.cex = 1.5, cex = 0.8, horiz = FALSE, inset = c(-0.05, 0.02))
  
  dev.off()
  
  
  rm(NETWORK, GRAPH, GRAPH.coord, BAC, EUK, legend.names, legend.cols,
     otu.names,otu.phyla,otu.phyla.num,otu.phyla.col,otu.sizes,otu.genus,subset)
}




# plot histogram of edge-weights
elist.mb <- summary(symBeta(getOptBeta(se.bac.euk), mode='maxabs'))
hist(elist.mb[,3], main="", xlab="edge weights")


###----- end of code
