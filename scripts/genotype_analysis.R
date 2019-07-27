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
satanas.vcf <- read.vcfR("raw_data/d_satanas_filtered.FIL.recode.vcf")

pal <- wes_palette("Darjeeling1", 14, type = "continuous")
scale_color_wes <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(pal, levels(localities$short_locality)), 
    ...
  )
}

# pca 
satanas.dna <- vcfR2DNAbin(satanas.vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
satanas.gen <- DNAbin2genind(satanas.dna)
satanas.pops <-  gsub( "_.*$", "", rownames(satanas.gen@tab))
satanas.gen@pop <- as.factor(satanas.pops)
satanas.scaled <- scaleGen(satanas.gen,NA.method="zero",scale=F)
satanas.pca <- prcomp(satanas.scaled,center=F,scale=F)
screeplot(satanas.pca)
satanas.pc <- data.frame(satanas.pca$x[,1:3])
satanas.pc$sample <- rownames(satanas.pc)
satanas.pc$pop <- satanas.pops
satanas.pc$species <- rep("Dichotomius_satanas",nrow(satanas.pc))

# merge with sample data
satanas.pc <- merge(satanas.pc, localities, by.x = "sample", by.y = "ddocent_ID")

a <- ggplot(data=satanas.pc,aes(x=PC1,y=PC2,col=satanas.pc$short_locality)) + 
  theme_classic() +
  geom_point()+
  scale_color_manual(values=pal,name="locality")

b <- ggplot(data=satanas.pc,aes(x=elevation,y=PC1,col=satanas.pc$short_locality)) + 
  theme_classic() +
  geom_point()+
  scale_color_manual(values=pal,name="locality")

a

b

### d. speciocissimum ###

# read in .vcf and sample / pop names; convert to genlight
spec.vcf <- read.vcfR("raw_data/d_spec_filtered.FIL.recode.vcf")

# pca 
spec.dna <- vcfR2DNAbin(spec.vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
spec.gen <- DNAbin2genind(spec.dna)
spec.pops <-  gsub( "_.*$", "", rownames(spec.gen@tab))
spec.gen@pop <- as.factor(spec.pops)
spec.scaled <- scaleGen(spec.gen,NA.method="mean",scale=F)
spec.pca <- prcomp(spec.scaled,center=F,scale=F)
screeplot(spec.pca)
spec.pc <- data.frame(spec.pca$x[,1:3])
spec.pc$sample <- rownames(spec.pc)
spec.pc$pop <- spec.pops
spec.pc$species <- rep("Deltochilum_speciocissimum",nrow(spec.pc))

# merge with sample data
spec.pc <- merge(spec.pc, localities, by.x = "sample", by.y = "ddocent_ID")

c <- ggplot(data=spec.pc,aes(x=PC1,y=PC2,col=spec.pc$short_locality)) + 
  theme_classic() +
  geom_jitter()+
  scale_color_manual(values=pal,name="locality")

d <- ggplot(data=spec.pc,aes(x=elevation,y=PC1,col=spec.pc$short_locality)) + 
  theme_classic() +
  geom_point()+
  scale_color_manual(values=pal,name="locality")

c

d

### d. tesselatum ###

# read in .vcf and sample / pop names; convert to genlight
tess.vcf <- read.vcfR("raw_data/d_tess_filtered.FIL.recode.vcf")

# pca 
tess.dna <- vcfR2DNAbin(tess.vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
tess.gen <- DNAbin2genind(tess.dna)
tess.pops <-  gsub( "_.*$", "", rownames(tess.gen@tab))
tess.gen@pop <- as.factor(tess.pops)
tess.scaled <- scaleGen(tess.gen,NA.method="mean",scale=F)
tess.pca <- prcomp(tess.scaled,center=F,scale=F)
screeplot(tess.pca)
tess.pc <- data.frame(tess.pca$x[,1:3])
tess.pc$sample <- rownames(tess.pc)
tess.pc$pop <- tess.pops
tess.pc$species <- rep("Deltochilum_tesselatum",nrow(tess.pc))

# merge with sample data
tess.pc <- merge(tess.pc, localities, by.x = "sample", by.y = "ddocent_ID")

e <- ggplot(data=tess.pc,aes(x=PC1,y=PC2,col=tess.pc$pop)) + 
  theme_classic() +
  geom_jitter()+
  scale_color_manual(values=pal,name="locality")

f <- ggplot(data=tess.pc,aes(x=elevation,y=PC1,col=tess.pc$pop)) + 
  theme_classic() +
  geom_point()+
  scale_color_manual(values=pal,name="locality")

### d. podalirius ###

# read in .vcf and sample / pop names; convert to genlight
pod.vcf <- read.vcfR("raw_data/d_pod_filtered.FIL.recode.vcf")

# pca 
pod.dna <- vcfR2DNAbin(pod.vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
pod.gen <- DNAbin2genind(pod.dna)
pod.pops <-  gsub( "_.*$", "", rownames(pod.gen@tab))
pod.gen@pop <- as.factor(pod.pops)
pod.scaled <- scaleGen(pod.gen,NA.method="mean",scale=F)
pod.pca <- prcomp(pod.scaled,center=F,scale=F)
screeplot(pod.pca)
pod.pc <- data.frame(pod.pca$x[,1:3])
pod.pc$sample <- rownames(pod.pc)
pod.pc$pop <- pod.pops
pod.pc$species <- rep("Dichotomius_podalirius",nrow(pod.pc))

# merge with sample data
pod.pc <- merge(pod.pc, localities, by.x = "sample", by.y = "ddocent_ID")

g <- ggplot(data=pod.pc,aes(x=PC1,y=PC2,col=pod.pc$pop)) + 
  theme_classic() +
  geom_jitter()+
  scale_color_manual(values=pal,name="locality")

h <- ggplot(data=pod.pc,aes(x=elevation,y=PC1,col=pod.pc$pop)) + 
  theme_classic() +
  geom_point()+
  scale_color_manual(values=pal,name="locality")

g

h

### e. affin ###

# read in .vcf and sample / pop names; convert to genlight
affin.vcf <- read.vcfR("raw_data/e_affin_filtered.FIL.recode.vcf")

# pca 
affin.dna <- vcfR2DNAbin(affin.vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
affin.gen <- DNAbin2genind(affin.dna)
affin.pops <-  gsub( "_.*$", "", rownames(affin.gen@tab))
affin.gen@pop <- as.factor(affin.pops)
affin.scaled <- scaleGen(affin.gen,NA.method="mean",scale=F)
affin.pca <- prcomp(affin.scaled,center=F,scale=F)
screeplot(affin.pca)
affin.pc <- data.frame(affin.pca$x[,1:3])
affin.pc$sample <- rownames(affin.pc)
affin.pc$pop <- affin.pops
affin.pc$species <- rep("Eurysternus_affin",nrow(affin.pc))

# merge with sample data
affin.pc <- merge(affin.pc, localities, by.x = "sample", by.y = "ddocent_ID")

i <- ggplot(data=affin.pc,aes(x=PC1,y=PC2,col=affin.pc$pop)) + 
  theme_classic() +
  geom_jitter()+
  scale_color_manual(values=pal,name="locality")

j <- ggplot(data=affin.pc,aes(x=elevation,y=PC1,col=affin.pc$pop)) + 
  theme_classic() +
  geom_point()+
  scale_color_manual(values=pal,name="locality")

i

j

# drop outliers
drop <- affin.pc[which(affin.pc$PC1>5),] 
drop <- drop$sample

affin.gen <- affin.gen[indNames(affin.gen)!=drop[1] & indNames(affin.gen)!=drop[2]]
affin.pops <-  gsub( "_.*$", "", rownames(affin.gen@tab))
affin.gen@pop <- as.factor(affin.pops)
affin.scaled <- scaleGen(affin.gen,NA.method="mean",scale=F)
affin.pca <- prcomp(affin.scaled,center=F,scale=F)
screeplot(affin.pca)
affin.pc <- data.frame(affin.pca$x[,1:3])
affin.pc$sample <- rownames(affin.pc)
affin.pc$pop <- affin.pops
affin.pc$species <- rep("Eurysternus_affin",nrow(affin.pc))

# merge with sample data
affin.pc <- merge(affin.pc, localities, by.x = "sample", by.y = "ddocent_ID")

i <- ggplot(data=affin.pc,aes(x=PC1,y=PC2,col=affin.pc$pop)) + 
  theme_classic() +
  geom_jitter()+
  scale_color_manual(values=pal,name="locality")

j <- ggplot(data=affin.pc,aes(x=elevation,y=PC1,col=affin.pc$pop)) + 
  theme_classic() +
  geom_point()+
  scale_color_manual(values=pal,name="locality")

i

j

### plot all with faceting
master.df <- rbind.data.frame(satanas.pc,spec.pc,tess.pc,pod.pc,affin.pc)

k <- ggplot(data=master.df,aes(x=PC1,y=PC2,col=master.df$pop)) + 
  theme_classic() +
  geom_jitter()+
  scale_color_manual(values=pal,name="locality") +
  facet_wrap(~species.x, scales="free")

l <- ggplot(data=master.df,aes(x=elevation,y=PC1,col=master.df$pop)) + 
  theme_classic() +
  geom_jitter()+
  scale_color_manual(values=pal,name="locality") +
  facet_wrap(~species.x, scales="free")

k 

l
