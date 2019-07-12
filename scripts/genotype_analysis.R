# genotype analysis for scarab sp.

library(vcfR)
library(adegenet)
library(ggplot2)
library(wesanderson)
library(tidyverse)
library(reshape2)

# set wd and read in locality data
setwd("/Users/ethanlinck/Dropbox/scarab_migration/")
localities <- read.csv("data/scarab_spp_master.csv")

### d. speciossimum ###

# read in .vcf and sample / pop names; convert to genlight
spec.vcf <- read.vcfR("ipyrad/d_speciossimum_outfiles/d_speciossimum.vcf")

# fix ipyrad vcf missing data issue
spec.vcf@gt <- 
  spec.vcf@gt %>% 
  as_tibble() %>% 
  mutate_all(.funs = function(x) replace(x, which(x == "./.:0:0,0,0,0"| x == "NA"), NA)) %>%
  as.matrix() 

# filter for depth, dropping top 10% and bottom 10% of distribution
dp <- extract.gt(spec.vcf,  element = "DP", as.numeric = TRUE)
class(dp)
dim(dp)
dpf <- melt(dp, varnames = c("Index", "Sample"),
            value.name = "Depth", na.rm = TRUE)
dpf <- dpf[ dpf$Depth > 0, ]
quants <- apply(dp, MARGIN = 2, quantile, probs = c(0.1, 0.9), na.rm = TRUE)
dp2 <- sweep(dp, MARGIN = 2, FUN = "-", quants[1, ])
dp[dp2 < 0] <- NA
dp2 <- sweep(dp, MARGIN = 2, FUN = "-", quants[2, ])
dp[dp2 > 0] <- NA
dp[dp < 4] <- NA
spec.vcf@gt[, -1][ is.na(dp) == TRUE ] <- NA
spec.vcf

## drop individuals with >30% MD
# dp <- extract.gt(spec.vcf,  element = "DP", as.numeric = TRUE)
# spec.miss <- apply(dp, MARGIN = 2, function(x){sum( is.na(x))})
# spec.miss <- spec.miss / nrow(dp)
# spec.vcf@gt <- spec.vcf@gt[, c(TRUE, spec.miss < 0.70)]
# spec.vcf

# drop rows with >20% MD
spec.miss <- apply(dp, MARGIN = 1, function(x){sum(is.na(x))})
spec.miss <- spec.miss / ncol(dp)
spec.vcf <- spec.vcf[spec.miss < 0.3, ]
spec.vcf

# pca 
spec.dna <- vcfR2DNAbin(spec.vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
spec.gen <- DNAbin2genind(spec.dna)
spec.samples <- as.character(read.table("ipyrad/d_speciossimum_outfiles/samples_d_spec.txt")[[1]]) 
spec.pops <- as.factor(as.character(read.table("ipyrad/d_speciossimum_outfiles/populations_d_spec.txt")[[1]])) 
rownames(spec.gen@tab) <- spec.samples
spec.gen@pop <- spec.pops
spec.scaled <- scaleGen(spec.gen,NA.method="mean",scale=F)
spec.pca <- prcomp(spec.scaled,center=F,scale=F)
screeplot(spec.pca)
spec.pc <- data.frame(spec.pca$x[,1:3])
spec.pc$sample <- rownames(spec.pc)
spec.pc$pop <- spec.pops
spec.pc$species <- rep("Deltochilum_speciocissimum",nrow(spec.pc))


### d. tesselatum ###

# read in .vcf and sample / pop names; convert to genlight
tess.vcf <- read.vcfR("ipyrad/d_tesselatum_outfiles/d_tesselatum.vcf")

# fix ipyrad vcf missing data issue
tess.vcf@gt <- 
  tess.vcf@gt %>% 
  as_tibble() %>% 
  mutate_all(.funs = function(x) replace(x, which(x == "./.:0:0,0,0,0"| x == "NA"), NA)) %>%
  as.matrix() 

# filter for depth, dropping top 10% and bottom 10% of distribution
dp <- extract.gt(tess.vcf,  element = "DP", as.numeric = TRUE)
class(dp)
dim(dp)
dpf <- melt(dp, varnames = c("Index", "Sample"),
            value.name = "Depth", na.rm = TRUE)
dpf <- dpf[ dpf$Depth > 0, ]
quants <- apply(dp, MARGIN = 2, quantile, probs = c(0.1, 0.9), na.rm = TRUE)
dp2 <- sweep(dp, MARGIN = 2, FUN = "-", quants[1, ])
dp[dp2 < 0] <- NA
dp2 <- sweep(dp, MARGIN = 2, FUN = "-", quants[2, ])
dp[dp2 > 0] <- NA
dp[dp < 4] <- NA
tess.vcf@gt[, -1][ is.na(dp) == TRUE ] <- NA
tess.vcf

## drop individuals with >30% MD
# dp <- extract.gt(spec.vcf,  element = "DP", as.numeric = TRUE)
# spec.miss <- apply(dp, MARGIN = 2, function(x){sum( is.na(x))})
# spec.miss <- spec.miss / nrow(dp)
# spec.vcf@gt <- spec.vcf@gt[, c(TRUE, spec.miss < 0.70)]
# spec.vcf

# drop rows with >20% MD
tess.miss <- apply(dp, MARGIN = 1, function(x){sum(is.na(x))})
tess.miss <- tess.miss / ncol(dp)
tess.vcf <- tess.vcf[tess.miss < 0.4, ]
tess.vcf

# pca 
tess.dna <- vcfR2DNAbin(tess.vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
tess.gen <- DNAbin2genind(tess.dna)
tess.samples <- as.character(read.table("ipyrad/d_tesselatum_outfiles/samples_d_tess.txt")[[1]]) 
tess.pops <- as.factor(as.character(read.table("ipyrad/d_tesselatum_outfiles/samples_d_tess.txt")[[1]])) 
rownames(tess.gen@tab) <- tess.samples
tess.gen@pop <- tess.pops
tess.scaled <- scaleGen(tess.gen,NA.method="mean",scale=F)
tess.pca <- prcomp(tess.scaled,center=F,scale=F)
screeplot(tess.pca)
tess.pc <- data.frame(tess.pca$x[,1:3])
tess.pc$sample <- rownames(tess.pc)
tess.pc$pop <- tess.pops
tess.pc$species <- rep("Deltochilum_tessellatum",nrow(tess.pc))

pc.master <- rbind(tess.pc,spec.pc)
localities <- localities[,-4]

# merge with sample data
pc.master <- merge(pc.master, localities, by.x = "sample", by.y = "sample_ID")
# names(pc.master)[names(pc.master) == 'species.x'] <- 'species'
# names(pc.master)[names(pc.master) == 'species.y'] <- 'species'

# quick visualization
pal <- wes_palette("Rushmore1", 8, type = "continuous")
ggplot(data=pc.master,aes(x=PC1,y=PC2,col=pc.master$locality)) + 
  geom_point() +
  theme_classic() +
  scale_color_manual(values=pal,name="locality") +
  facet_wrap(~species, scales = "free") +
  theme(strip.background = element_rect(color="white"))

# PC1 by elevation
ggplot(data=pc.master,aes(x=PC1,y=elevation,col=pc.master$locality)) + 
  geom_point(position="jitter") +
  theme_classic() +
  scale_color_manual(values=pal,name="locality") +
  facet_wrap(~species, scales = "free") +
  theme(strip.background = element_rect(color="white"))

# drop bad samples. bad!
# spec.pc.2 <- spec.pc[!spec.pc$pop=="MA2150" & !spec.pc$pop=="2150" & !spec.pc$pop=="YY2150",]





