########## Identify variables associated with observed differences in composition using Envfit


###----- load libraries, data and set directory
library(phyloseq)
library(vegan)
library(ggplot2)
library(DESeq2)

setwd("envfit/")


# load phyloseq object with OTU data, taxonomy and sample metadata
load("microbes.rdata")




###----- prepare functions to filter envfit output

# Prep function to filter by minimum correlation (r2): select.envfit(fit, r.select)
select.envfit <- function(fit, r.select){ 
  for (i in 1:length(fit$vectors$r)) { # run loop checking r values in fit$vectors$r
    if (fit$vectors$r[i]<r.select) { 
      fit$vectors$arrows[i,] <- NA # if r below designated threshold, set vector coordinates to NA
      i=i+1 
    }
  } 
  return(fit) #return fit as the result of the function
}

# Prep function to filter by minimum p.val: select.pval.envfit(fit, p.val)
select.pval.envfit <- function(fit, p.val){ 
  for (i in 1:length(fit$vectors$pvals)) {
    
    # filter vectors
    if (fit$vectors$pvals[i]>=p.val) { 
      fit$vectors$arrows[i,] <- NA # if r below designated threshold, set vector coordinates to NA
      i=i+1 
    }
  }
  return(fit) #return fit as the result of the function
}




###----- select dataset, prefilter and transform data

# select microbial dataset
df <- bacteria


# set parameters
MINSUBJECT <- 0.05
MINCOUNT <- 5
DIST <- "bray" # distance
DATATRANSFORM <- "deseq" #rarefy, deseq
RAREDEPTH <- 25000 # 25000 bacteria, 1000 eukaryotes


# pre-process data (transform OTU table)
if(DATATRANSFORM == "rarefy") {
  df.filt <- rarefy_even_depth(df, rngseed = 1234, sample.size = RAREDEPTH)

} else if(DATATRANSFORM == "deseq") {
  
  # first filter out low counts
  df.filt <- filter_taxa(df, function(x) sum(x > MINCOUNT) >= (MINSUBJECT*length(x)), TRUE)

  df.deseq <- phyloseq_to_deseq2(df.filt, ~timepoint)
  dds <- DESeq(df.deseq)
  dds <- estimateSizeFactors(df.deseq) #calculate normalization factors
  NORMCOUNTS <- counts(dds, normalized=TRUE) #extract normalized counts
  
  NORMCOUNTS <- otu_table(NORMCOUNTS, taxa_are_rows = TRUE) # add normalized table to phyloseq object
  otu_table(df.filt) <- NORMCOUNTS
  rm(NORMCOUNTS,df.deseq,dds)

}


# extract pre-processed OTU table
OTU <- t(otu_table(df.filt))



#----- calculate distance for NMDS and ordinate using vegan
OTU.mds <- metaMDS(comm = OTU, distance = DIST, k = 3)

# extract nmds results into dataframe
OTU.nmds <- data.frame(MDS1 = OTU.mds$points[,1], MDS2 = OTU.mds$points[,2])




########################################################################################

#----- ENVFIT - DETERMINE CORRELATED TAXA
# r2 = goodness-of-fit
# output of continuous variables gives the direction cosines = coordinates of vector heads 
# scale vectors by their correlation (square root of r2)  to get coordinates
# permutation to assess significance of fitted vectors

env.species <- envfit(OTU.mds$points, OTU, perm=999)

env.species.df <- as.data.frame(env.species$vectors$arrows*sqrt(env.species$vectors$r))
env.species.df$species <- rownames(env.species.df)


# see scaled vectors
scores(env.species, "vectors")


# filter out low correlations and insignificant p.vals with pre-defined functions
MINR2 <- 0.1
env.species.filt.p <- select.pval.envfit(fit = env.species, p.val = 0.05)
env.species.filt.pr <- select.envfit(fit = env.species.filt.p, r.select = MINR2)



# extract dataframe of coordinates
env.species.filt.df <- as.data.frame(env.species.filt.pr$vectors$arrows*sqrt(env.species.filt.pr$vectors$r))[taxa,]
env.species.filt.df <- env.species.filt.df[!is.na(env.species.filt.df$MDS1),]
# add OTU taxonomy
env.species.filt.df$otu <- rownames(env.species.filt.df)
env.species.filt.df$species <- rownames(env.species.filt.df)
env.species.filt.df$genus <- as.character(as.matrix(tax_table(bacteria)[rownames(env.species.filt.df),"genus"]))
env.species.filt.df$phylum <- as.character(as.matrix(tax_table(bacteria)[rownames(env.species.filt.df),"phylum"]))


# add metadata
OTU.nmds$timepoint <- as.matrix(sample_data(bacteria)[row.names(OTU.nmds), "timepoint"])
OTU.nmds$center <-  as.matrix(sample_data(bacteria)[row.names(OTU.nmds), "center"])
OTU.nmds$arm <-  as.matrix(sample_data(bacteria)[row.names(OTU.nmds), "arm"])



# plot ordination coloured by age, with taxon loading arrows
p <- ggplot(data = OTU.nmds, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = timepoint), size=3) + 
  ggtitle(paste0("NMDS : ", DIST)) + theme_bw() +
  scale_color_manual(values = c("grey70", "grey30")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_segment(data = env.species.filt.df, aes(x=0, xend=MDS1, y=0, yend=MDS2),
               arrow = arrow(length = unit(0.5, "cm")), colour="grey", inherit_aes=FALSE) + 
  geom_text(data = env.species.filt.df, aes(x=MDS1, y=MDS2, label=otu), size=3) +
  coord_fixed()

p <- p + stat_ellipse(data = OTU.nmds, aes(x=MDS1, y=MDS2, fill=timepoint), 
                      geom = "polygon", type="norm", alpha=0.2) +
  scale_fill_manual(values = c("grey70", "grey30"))

print(p)


### End of document