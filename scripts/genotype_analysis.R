# PCA analysis of genotypes for scarab sp.

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

# palette
pal <- wes_palette("Darjeeling1", 14, type = "continuous")
scale_color_wes <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(pal, levels(localities$short_locality)), 
    ...
  )
}

### d. satanas ###

# read in .vcf and sample / pop names; convert to genlight
satanas.vcf <- read.vcfR("raw_data/d_satanas_filtered.FIL.recode.vcf")

# pca 
satanas.dna <- vcfR2DNAbin(satanas.vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
satanas.gen <- DNAbin2genind(satanas.dna)
satanas.pops <-  gsub( "_.*$", "", rownames(satanas.gen@tab))
satanas.gen@pop <- as.factor(satanas.pops)
satanas.scaled <- scaleGen(satanas.gen,NA.method="zero",scale=F)
satanas.pca <- prcomp(satanas.scaled,center=F,scale=F)
summary(satanas.pca) #PC1 6.486%, PC2 4.666%
screeplot(satanas.pca)
satanas.pc <- data.frame(satanas.pca$x[,1:3])
satanas.pc$sample <- rownames(satanas.pc)
satanas.pc$pop <- satanas.pops
satanas.pc$species <- rep("Dichotomius_satanas",nrow(satanas.pc))

# merge with sample data
satanas.pc <- merge(satanas.pc, localities, by.x = "sample", by.y = "ddocent_ID")

# k-means clustering
sat.grp <- find.clusters(satanas.gen, max.n.clust=98)

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
summary(spec.pca) # PC1 78.91%, PC2 1.942%
screeplot(spec.pca)
spec.pc <- data.frame(spec.pca$x[,1:3])
spec.pc$sample <- rownames(spec.pc)
spec.pc$pop <- spec.pops
spec.pc$species <- rep("Deltochilum_speciocissimum",nrow(spec.pc))

# drop outliers!
spec.dna <- vcfR2DNAbin(spec.vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
spec.dna2 <- spec.dna[-c(41:47),]
spec.dna2 <- spec.dna2[-23,]
spec.gen2 <- DNAbin2genind(spec.dna2)
spec.pops2 <-  gsub( "_.*$", "", rownames(spec.gen2@tab))
spec.gen2@pop <- as.factor(spec.pops2)
spec.scaled2 <- scaleGen(spec.gen2,NA.method="mean",scale=F)
spec.pca2 <- prcomp(spec.scaled2,center=F,scale=F)
summary(spec.pca2) # PC1  3.279%, PC2 3.203%
screeplot(spec.pca2)
spec.pc2 <- data.frame(spec.pca2$x[,1:3])
spec.pc2$sample <- as.factor(rownames(spec.pc2))
spec.pc2$pop <- spec.pops2
spec.pc2$species <- rep("Deltochilum_speciocissimum",nrow(spec.pc2))

# merge with sample data
spec.pc2 <- merge(spec.pc2, localities, by.x = "sample", by.y = "ddocent_ID")

# k-means clustering
spec.grp <- find.clusters(spec.gen2, max.n.clust=46)

c <- ggplot(data=spec.pc2,aes(x=PC1,y=PC2,col=spec.pc2$short_locality)) + 
  theme_classic() +
  geom_jitter()+
  scale_color_manual(values=pal,name="locality")

d <- ggplot(data=spec.pc2,aes(x=elevation,y=PC1,col=spec.pc2$short_locality)) + 
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
summary(tess.pca) # PC1 7.606%, PC2 7.138%
screeplot(tess.pca)
tess.pc <- data.frame(tess.pca$x[,1:3])
tess.pc$sample <- rownames(tess.pc)
tess.pc$pop <- tess.pops
tess.pc$species <- rep("Deltochilum_tesselatum",nrow(tess.pc))

# merge with sample data
tess.pc <- merge(tess.pc, localities, by.x = "sample", by.y = "ddocent_ID")

# k-means clustering
tess.grp <- find.clusters(tess.gen, max.n.clust=18)

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
summary(pod.pca) #PC1 23.74%, 8.348%
screeplot(pod.pca)
pod.pc <- data.frame(pod.pca$x[,1:3])
pod.pc$sample <- rownames(pod.pc)
pod.pc$pop <- pod.pops
pod.pc$species <- rep("Dichotomius_podalirius",nrow(pod.pc))

# merge with sample data
pod.pc <- merge(pod.pc, localities, by.x = "sample", by.y = "ddocent_ID")

# k-means clustering
pod.grp <- find.clusters(pod.gen, max.n.clust=25)

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

# k-means clustering
affin.grp <- find.clusters(affin.gen, max.n.clust=23)

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
summary(affin.pca) #PC1 8.911%, #PC2 8.705%
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
master.df <- rbind.data.frame(satanas.pc,spec.pc2,tess.pc,pod.pc,affin.pc)
master.df$species.x <- as.factor(master.df$species.x)
master.df$species.x <- tolower(master.df$species.x)
master.df$species.x <- factor(master.df$species.x,levels=c("eurysternus_affin","dichotomius_podalirius",
                                                 "deltochilum_tesselatum","deltochilum_speciocissimum","dichotomius_satanas"))

k <- ggplot(data=master.df,aes(x=PC1,y=PC2,col=master.df$short_locality)) + 
  theme_bw() +
  geom_jitter()+
  scale_color_wes(name="locality") +
  facet_wrap(~species.x, scales="free") +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank()) 

l <- ggplot(data=master.df,aes(x=elevation,y=PC1,col=master.df$short_locality)) + 
  theme_bw() +
  geom_jitter(size=2.5)+
  scale_color_wes(name="locality") +
  facet_wrap(~species.x, scales="free") +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank()) 

l 

write.csv(master.df, "~/Dropbox/scarab_migration/data/genotypes.df.csv")

# large version for presentation
m <- ggplot(data=master.df,aes(x=elevation,y=PC1,col=master.df$short_locality)) + 
  theme_bw() +
  geom_jitter(size=2.5)+
  scale_color_wes(name="locality") +
  facet_wrap(~species.x, scales="free") +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.text=element_text(size=12),
    strip.text.x=element_text(size=12),
    axis.title.x=element_text(size=12),
    axis.title.y=element_text(size=12),
    legend.title=element_text(size=12),
    legend.text=element_text(size=10)) 

m

# stats for paper
length(unique(satanas.gen@loc.fac)) #372
length(unique(spec.gen@loc.fac)) #27379
length(unique(tess.gen@loc.fac)) #1047
length(unique(pod.gen@loc.fac)) #73
length(unique(affin.gen@loc.fac)) #517

ncol(satanas.gen@tab) #746
ncol(spec.gen@tab) #54890
ncol(tess.gen@tab) #2094
ncol(pod.gen@tab) #147
ncol(affin.gen@tab) #1039

# regressions w/ PC1 and elevation
summary(lm(satanas.pc$PC1 ~ satanas.pc$elevation)) #p-value: 0.6773, Adjusted R-squared:  -0.008498
summary(lm(spec.pc2$PC1 ~ spec.pc2$elevation)) #p-value: 0.1436, Adjusted R-squared:  0.03142
summary(lm(tess.pc$PC1 ~ tess.pc$elevation)) #p-value: 0.7185, Adjusted R-squared:  -0.1395
summary(lm(pod.pc$PC1 ~ pod.pc$elevation)) #p-value: 0.9488, Adjusted R-squared:  -0.04148 
summary(lm(affin.pc$PC1 ~ affin.pc$elevation)) #p-value: 0.3206, Adjusted R-squared:  0.001769 
