## running bedassle

library(adegenet)
library(BEDASSLE)
library(dplyr)
library(ggplot2)
library(broom)
library(data.table)
library(vcfR)
library(fossil)

# set wd and read in locality data
setwd("/Users/ethanlinck/Dropbox/scarab_migration/")
localities <- read.csv("data/scarab_spp_master.csv")
localities$short_locality[which(localities$short_locality=="730")] <- 800
pop.loc <- cbind.data.frame(unique(localities$long),unique(localities$lat),unique(localities$short_locality))
colnames(pop.loc) <- c("long", "lat", "pop")
geodist <- earth.dist(pop.loc[c("long", "lat")], dist = FALSE)

## d satanas

# read in .vcf, turn to genpop
satanas.vcf <- read.vcfR("raw_data/d_satanas_filtered.FIL.recode.vcf", convertNA=TRUE)
satanas.gen <- vcfR2genind(satanas.vcf) 
md <- propTyped(satanas.gen, by=c("both")) # confirm we have expected amount of missing data
table(md)[1]/table(md)[2] # proportion of missing sites
table(satanas.gen@tab)
satanas.pops <-  gsub( "_.*$", "", rownames(satanas.gen@tab)) # add pops
satanas.gen@pop <- as.factor(satanas.pops) # make factor
satanas.b <- genind2genpop(satanas.gen, satanas.pops) # turn into gen pop

# convert to BEDASSLE format
satanas.ac <- as.matrix(satanas.b@tab)
del <- seq(2, ncol(satanas.ac), 2) # sequence of integers to drop non-ref allele
satanas.ac <- satanas.ac[,-del] 

# create equally sized matrix for sample sizes
satanas.n <- matrix(nrow=nrow(satanas.ac), ncol=ncol(satanas.ac))

# name our rows the same thing
rownames(satanas.n) <- rownames(satanas.ac)

# get sample size per population
sample.n <- table(satanas.gen@pop) 

# turn this into a vector
sample.sizes <- as.vector(sample.n)

# populate each row of matrix with sample sizes for pops
for(i in 1:nrow(satanas.n)){
  satanas.n[i,] <- sample.sizes[i]
}
satanas.n <- satanas.n*2 # adjust to account for loss of one allele

# calculate pairwise Fst
satanas.p.fst.all <- calculate.all.pairwise.Fst(satanas.ac, satanas.n)
adj.fst <- satanas.p.fst.all / (1 - satanas.p.fst.all)

# look at global Fst
satanas.p.fst <- calculate.pairwise.Fst(satanas.ac, satanas.n)
drop.pop <- c("800","925","1050","1175")

# drop levels and calc distance
pop.loc.sat <- pop.loc[!pop.loc$pop %in% drop.pop,]
geodist.sat <- earth.dist(pop.loc.sat[c("long", "lat")], dist = FALSE)

# uh what does a negative relationship mean 
plot(geodist.sat,
     satanas.p.fst.all,
     pch=19,
     ylab="pairwise Fst",
     xlab="geographic distance",
     main="isolation by distance")
legend(x="bottomright",pch=19,col=c(1,2))



