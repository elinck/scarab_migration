# genotype analysis for scarab sp.

library(vcfR)
library(adegenet)
library(ggplot2)
library(wesanderson)
library(tidyverse)
library(reshape2)
library(patchwork)

# set wd and read in locality data
setwd("/Users/ethanlinck/Dropbox/scarab_migration/")
localities <- read.csv("data/scarab_spp_master.csv")

### d. satanas ###

# read in .vcf and sample / pop names; convert to genlight
satanas.vcf <- read.vcfR("ipyrad/d_satanas_outfiles/d_satanas.vcf")

# fix ipyrad vcf missing data issue
satanas.vcf@gt <- 
  satanas.vcf@gt %>% 
  as_tibble() %>% 
  mutate_all(.funs = function(x) replace(x, which(x == "./.:0:0,0,0,0"| x == "NA"), NA)) %>%
  as.matrix() 

# drop rows with >20% MD
dp <- extract.gt(satanas.vcf,  element = "DP", as.numeric = TRUE)
satanas.miss <- apply(dp, MARGIN = 1, function(x){sum(is.na(x))})
satanas.miss <- satanas.miss / ncol(dp)
satanas.vcf <- satanas.vcf[satanas.miss < 0.4, ]
satanas.vcf

# pca 
satanas.dna <- vcfR2DNAbin(satanas.vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
satanas.gen <- DNAbin2genind(satanas.dna)
satanas.samples <- as.character(read.table("ipyrad/d_satanas_outfiles/samples_d_satanas.txt")[[1]]) 
satanas.pops <- as.factor(as.character(read.table("ipyrad/d_satanas_outfiles/populations_d_satanas.txt")[[1]])) 
rownames(satanas.gen@tab) <- satanas.samples
satanas.gen@pop <- satanas.pops
satanas.scaled <- scaleGen(satanas.gen,NA.method="mean",scale=F)
satanas.pca <- prcomp(satanas.scaled,center=F,scale=F)
screeplot(satanas.pca)
satanas.pc <- data.frame(satanas.pca$x[,1:3])
satanas.pc$sample <- rownames(satanas.pc)
satanas.pc$pop <- satanas.pops
satanas.pc$species <- rep("Dichotomius_satanas",nrow(satanas.pc))

# merge with sample data
satanas.pc <- merge(satanas.pc, localities, by.x = "sample", by.y = "sample_ID")

# quick visualization
pal <- wes_palette("Rushmore1", 8, type = "continuous")

a <- ggplot(data=satanas.pc,aes(x=PC1,y=PC2,col=satanas.pc$pop)) + 
  theme_classic() +
  geom_point()+
  scale_color_manual(values=pal,name="locality")

b <- ggplot(data=satanas.pc,aes(x=elevation,y=PC1,col=satanas.pc$pop)) + 
  theme_classic() +
  geom_point()+
  scale_color_manual(values=pal,name="locality")

a+b




