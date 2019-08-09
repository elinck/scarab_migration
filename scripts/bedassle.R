## running bedassle

library(adegenet)
library(BEDASSLE)
library(dplyr)
library(ggplot2)
library(broom)
library(data.table)
library(vcfR)
library(fossil)
library(viridis)
library(patchwork)
library(sp)
library(raster)
library(ecodist)
library(radiator)
library(wesanderson)

# set wd and read in locality data
setwd("/Users/ethanlinck/Dropbox/scarab_migration/")
localities <- read.csv("data/scarab_spp_master.csv")
localities$short_locality[which(localities$short_locality=="730")] <- 800
pop.loc <- cbind.data.frame(unique(localities$long),unique(localities$lat),unique(localities$short_locality))
colnames(pop.loc) <- c("long", "lat", "pop")
geodist <- earth.dist(pop.loc[c("long", "lat")], dist = FALSE)

### d satanas

# read in .vcf, turn to genpop
satanas.vcf <- read.vcfR("raw_data/d_satanas.LD.vcf", convertNA=TRUE)
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

# look at global Fst
satanas.p.fst <- calculate.pairwise.Fst(satanas.ac, satanas.n)
drop.pop <- c("800","925","1050","1175")

# drop levels and calc distance
pop.loc.sat <- pop.loc[!pop.loc$pop %in% drop.pop,]
droplevels(pop.loc.sat)
satanas.geo <- earth.dist(pop.loc.sat[c("long", "lat")], dist = FALSE)

# turn to vectors
sat.dist <- as.vector(satanas.geo)
sat.gen <- as.vector(satanas.p.fst.all)

# make data frame with these variables
sat.df <- cbind.data.frame(sat.dist, sat.gen)
sat.df$species <- rep("dichotomius_satanas",nrow(sat.df))
colnames(sat.df) <- c("distance", "fst", "species")

# calculate environmental distance
worldclim <- getData("worldclim",var="bio",res=0.5,lon=localities$long[1],lat=localities$lat[1])
proj <- as.character(worldclim[[2]]@crs)
worldclim <- worldclim[[c(1,12)]]
names(worldclim) <- c("temp","precip")
sp1 <- SpatialPoints(localities[,c('long', 'lat')], proj4string=CRS(proj))
values <- extract(worldclim, sp1)
loc.master <- cbind.data.frame(localities, values)
loc.uniq <- cbind.data.frame(loc.master$locality, loc.master$elevation, loc.master$temp,
                             loc.master$precip) %>% distinct()
colnames(loc.uniq) <- c("population", "elevation", "temp", "precip")
rownames(loc.uniq) <- loc.uniq$population
loc.uniq <- as.data.frame(loc.uniq)
loc.uniq <- loc.uniq[,-which(names(loc.uniq) %in% c("population"))]
env.dist.all <- dist(loc.uniq, diag = TRUE, upper = TRUE)

# drop lower transect
loc.subset <- loc.uniq[-which(rownames(loc.uniq) %in% 
                                       c("CC1-800","CC2-925","CC3-1050","CC4-1175","CC1-730")),]
env.dist.upper <- dist(loc.subset, diag = TRUE, upper = TRUE)

# prep dist matrix
satanas.env <- as.matrix(env.dist.upper)
colnames(satanas.env) <- NULL
rownames(satanas.env) <- NULL

# run MCMC for 10K gens
MCMC(   
  counts = satanas.ac,
  sample_sizes = satanas.n,
  D = satanas.geo,  # geographic distances
  E = satanas.env,  # environmental distances
  k = nrow(satanas.ac), loci = ncol(satanas.ac),  # dimensions of the data
  delta = 0.0001,  # a small, positive, number
  aD_stp = 0.1,   # step sizes for the MCMC
  aE_stp = 0.1,
  a2_stp = 0.025,
  thetas_stp = 0.2,
  mu_stp = 0.35,
  ngen = 10000,        # number of steps (2e6)
  printfreq = 100,  # print progress (10000)
  savefreq = 100,     # save out current state
  samplefreq = 2,     # record current state for posterior (2000)
  prefix = "/Users/ethanlinck/Dropbox/scarab_migration/bedassle/satanas_",   # filename prefix
  continue=FALSE,
  continuing.params=NULL)

# check shit out
show(load("/Users/ethanlinck/Dropbox/scarab_migration/bedassle/satanas_MCMC_output1.Robj"))
layout(t(1:2))
plot(aD[-(1:1000)], xlab="MCMC generation", ylab="value", main="aD")
plot((aD_accept/aD_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="aD acceptance", ylim=c(0,1))
plot(as.vector(aE)[-(1:1000)], xlab="MCMC generation", ylab="value", main="aE")
plot((aE_accept/aE_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="aE acceptance", ylim=c(0,1))
plot(a2[-(1:1000)], xlab="MCMC generation", ylab="value", main="a2")
plot((a2_accept/a2_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="a2 acceptance", ylim=c(0,1))
plot(Prob[-(1:1000)], xlab="MCMC generation", main="log likelihood")
plot((mu_accept/mu_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="mu acceptance", ylim=c(0,1) )
plot((thetas_accept/thetas_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="thetas acceptance", ylim=c(0,1) )
hist((aE/aD)[-(1:1000)],breaks=100,main="posterior of aE/aD ratio")

# convert for ggplotting
sat.ratio <- aE/aD %>% as.data.frame()
sat.bed.df <- cbind.data.frame(sat.ratio, rep("dichotomius_satanas", nrow(sat.ratio)))
colnames(sat.bed.df) <- c("ratio","species")

### d speciocissimum

# read in .vcf, turn to genpop
spec.vcf <- read.vcfR("raw_data/d_spec.LD.vcf", convertNA=TRUE)
spec.gen <- vcfR2genind(spec.vcf) 
md <- propTyped(spec.gen, by=c("both")) # confirm we have expected amount of missing data
table(md)[1]/table(md)[2] # proportion of missing sites
table(spec.gen@tab)
spec.pops <-  gsub( "_.*$", "", rownames(spec.gen@tab)) # add pops
spec.gen@pop <- as.factor(spec.pops) # make factor
spec.b <- genind2genpop(spec.gen, spec.pops) # turn into gen pop

# convert to BEDASSLE format
spec.ac <- as.matrix(spec.b@tab)
spec.ac <- spec.ac[-1,]
del <- seq(2, ncol(spec.ac), 2) # sequence of integers to drop non-ref allele
spec.ac <- spec.ac[,-del] 
dim(spec.ac[,complete.cases(t(spec.ac))])

# create equally sized matrix for sample sizes
spec.n <- matrix(nrow=nrow(spec.ac), ncol=ncol(spec.ac))

# name our rows the same thing
rownames(spec.n) <- rownames(spec.ac)

# get sample size per population
table(spec.gen@pop)
sample.sizes <- as.vector(table(spec.gen@pop))
sample.sizes <- sample.sizes[-1] # drop pop w/ only one individual

# populate each row of matrix with sample sizes for pops
for(i in 1:nrow(spec.n)){
  spec.n[i,] <- sample.sizes[i]
}
spec.n <- spec.n*2 # adjust to account for loss of one allele

# calculate pairwise Fst
spec.p.fst.all <- calculate.all.pairwise.Fst(spec.ac, spec.n)

# look at global Fst
spec.p.fst <- calculate.pairwise.Fst(spec.ac, spec.n)
drop.pop <- c("800","925","1050","1175","1575","Macucaloma","Yanayacu")

# drop levels and calc distance
pop.loc.spec <- pop.loc[!pop.loc$pop %in% drop.pop,]
droplevels(pop.loc.spec)
spec.geo <- earth.dist(pop.loc.spec[c("long", "lat")], dist = FALSE)

# turn to vectors
spec.dist <- as.vector(spec.geo)
spec.gend <- as.vector(spec.p.fst.all)

# make data frame with these variables
spec.df <- cbind.data.frame(spec.dist, spec.gend)
spec.df$species <- rep("deltochilum_speciocissimum",nrow(spec.df))
colnames(spec.df) <- c("distance", "fst", "species")

# prep dist matrix
drop.pop <- c("800","925","1050","1175","HU-1575","Macucaloma","Yanayacu")
spec.env <- as.matrix(env.dist.upper)
spec.env <- spec.env[!rownames(spec.env) %in% drop.pop,]
spec.env <- spec.env[,!colnames(spec.env) %in% drop.pop]
colnames(spec.env) <- NULL
rownames(spec.env) <- NULL

# run MCMC for 100K gens
MCMC(   
  counts = spec.ac,
  sample_sizes = spec.n,
  D = spec.geo,  # geographic distances
  E = spec.env,  # environmental distances
  k = 5, loci = 6899,  # dimensions of the data
  delta = 0.0001,  # a small, positive, number
  aD_stp = 0.1,   # step sizes for the MCMC
  aE_stp = 0.1,
  a2_stp = 0.025,
  thetas_stp = 0.2,
  mu_stp = 0.35,
  ngen = 10000,        # number of steps (2e6)
  printfreq = 100,  # print progress (10000)
  savefreq = 100,     # save out current state
  samplefreq = 2,     # record current state for posterior (2000)
  prefix = "/Users/ethanlinck/Dropbox/scarab_migration/bedassle/spec_",   # filename prefix
  continue=FALSE,
  continuing.params=NULL)

# check shit out
show(load("/Users/ethanlinck/Dropbox/scarab_migration/bedassle/spec_MCMC_output1.Robj"))
layout(t(1:2))
plot(aD[-(1:1000)], xlab="MCMC generation", ylab="value", main="aD")
plot((aD_accept/aD_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="aD acceptance", ylim=c(0,1))
plot(as.vector(aE)[-(1:1000)], xlab="MCMC generation", ylab="value", main="aE")
plot((aE_accept/aE_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="aE acceptance", ylim=c(0,1))
plot(a2[-(1:1000)], xlab="MCMC generation", ylab="value", main="a2")
plot((a2_accept/a2_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="a2 acceptance", ylim=c(0,1))
plot(Prob[-(1:1000)], xlab="MCMC generation", main="log likelihood")
plot((mu_accept/mu_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="mu acceptance", ylim=c(0,1) )
plot((thetas_accept/thetas_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="thetas acceptance", ylim=c(0,1) )
hist((aE/aD)[-(1:1000)],breaks=100,main="posterior of aE/aD ratio")

# convert for ggplotting
spec.ratio <- aE/aD %>% as.data.frame()
spec.bed.df <- cbind.data.frame(spec.ratio, rep("deltochilum_speciocissimum", nrow(spec.ratio)))
colnames(spec.bed.df) <- c("ratio","species")

### d tess

# read in .vcf, turn to genpop
tess.vcf <- read.vcfR("raw_data/d_tess.LD.vcf", convertNA=TRUE)
tess.gen <- vcfR2genind(tess.vcf) 
md <- propTyped(tess.gen, by=c("both")) # confirm we have expected amount of missing data
table(md)[1]/table(md)[2] # proportion of missing sites
table(tess.gen@tab)
tess.pops <-  gsub( "_.*$", "", rownames(tess.gen@tab)) # add pops
tess.gen@pop <- as.factor(tess.pops) # make factor
tess.b <- genind2genpop(tess.gen, tess.pops) # turn into gen pop

# convert to BEDASSLE format
tess.ac <- as.matrix(tess.b@tab)
del <- seq(2, ncol(tess.ac), 2) # sequence of integers to drop non-ref allele
tess.ac <- tess.ac[,-del] 

# create equally sized matrix for sample sizes
tess.n <- matrix(nrow=nrow(tess.ac), ncol=ncol(tess.ac))

# name our rows the same thing
rownames(tess.n) <- rownames(tess.ac)

# get sample size per population
sample.n <- table(tess.gen@pop) 

# turn this into a vector
sample.sizes <- as.vector(sample.n)

# populate each row of matrix with sample sizes for pops
for(i in 1:nrow(tess.n)){
  tess.n[i,] <- sample.sizes[i]
}
tess.n <- tess.n*2 # adjust to account for loss of one allele

# calculate pairwise Fst
tess.p.fst.all <- calculate.all.pairwise.Fst(tess.ac, tess.n)
adj.fst <- tess.p.fst.all / (1 - tess.p.fst.all)

# look at global Fst
tess.p.fst <- calculate.pairwise.Fst(tess.ac, tess.n)
drop.pop <- c("800","925","1050","1175","HU-3","Macucaloma","Yanayacu")

# drop levels and calc distance
pop.loc.tess <- pop.loc[!pop.loc$pop %in% drop.pop,]
droplevels(pop.loc.tess)
tess.geo <- earth.dist(pop.loc.tess[c("long", "lat")], dist = FALSE)

# turn to vectors
tess.dist <- as.vector(tess.geo)
tess.gend <- as.vector(tess.p.fst.all)

# make data frame with these variables
tess.df <- cbind.data.frame(tess.dist, tess.gend)
tess.df$species <- rep("deltochilum_tesselatum",nrow(tess.df))
colnames(tess.df) <- c("distance", "fst", "species")

# prep env matrix
tess.env <- as.matrix(env.dist.upper)
tess.env <- tess.env[!rownames(tess.env) %in% drop.pop,]
tess.env <- tess.env[,!colnames(tess.env) %in% drop.pop]
colnames(tess.env) <- NULL
rownames(tess.env) <- NULL

# run MCMC for 100K gens
MCMC(   
  counts = tess.ac,
  sample_sizes = tess.n,
  D = tess.geo,  # geographic distances
  E = tess.env,  # environmental distances
  k = 5, loci = 3744,  # dimensions of the data
  delta = 0.0001,  # a small, positive, number
  aD_stp = 0.1,   # step sizes for the MCMC
  aE_stp = 0.1,
  a2_stp = 0.025,
  thetas_stp = 0.2,
  mu_stp = 0.35,
  ngen = 10000,        # number of steps (2e6)
  printfreq = 100,  # print progress (10000)
  savefreq = 100,     # save out current state
  samplefreq = 2,     # record current state for posterior (2000)
  prefix = "/Users/ethanlinck/Dropbox/scarab_migration/bedassle/tess_",   # filename prefix
  continue=FALSE,
  continuing.params=NULL)

# check shit out
show(load("/Users/ethanlinck/Dropbox/scarab_migration/bedassle/tess_MCMC_output1.Robj"))
layout(t(1:2))
plot(aD[-(1:1000)], xlab="MCMC generation", ylab="value", main="aD")
plot((aD_accept/aD_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="aD acceptance", ylim=c(0,1))
plot(as.vector(aE)[-(1:1000)], xlab="MCMC generation", ylab="value", main="aE")
plot((aE_accept/aE_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="aE acceptance", ylim=c(0,1))
plot(a2[-(1:1000)], xlab="MCMC generation", ylab="value", main="a2")
plot((a2_accept/a2_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="a2 acceptance", ylim=c(0,1))
plot(Prob[-(1:1000)], xlab="MCMC generation", main="log likelihood")
plot((mu_accept/mu_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="mu acceptance", ylim=c(0,1) )
plot((thetas_accept/thetas_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="thetas acceptance", ylim=c(0,1) )
hist((aE/aD)[-(1:1000)],breaks=100,main="posterior of aE/aD ratio")

# convert for ggplotting
tess.ratio <- aE/aD %>% as.data.frame()
tess.bed.df <- cbind.data.frame(tess.ratio, rep("deltochilum_tesselatum", nrow(tess.ratio)))
colnames(tess.bed.df) <- c("ratio","species")

### d pod

# read in .vcf, turn to genpop
pod.vcf <- read.vcfR("raw_data/d_pod.LD.vcf", convertNA=TRUE)
pod.gen <- vcfR2genind(pod.vcf) 
md <- propTyped(pod.gen, by=c("both")) # confirm we have expected amount of missing data
table(md)[1]/table(md)[2] # proportion of missing sites
table(pod.gen@tab)
pod.pops <-  gsub( "_.*$", "", rownames(pod.gen@tab)) # add pops
pod.gen@pop <- as.factor(pod.pops) # make factor
pod.b <- genind2genpop(pod.gen, pod.pops) # turn into gen pop

# convert to BEDASSLE format
pod.ac <- as.matrix(pod.b@tab)
del <- seq(2, ncol(pod.ac), 2) # sequence of integers to drop non-ref allele
pod.ac <- pod.ac[,-del] 

# create equally sized matrix for sample sizes
pod.n <- matrix(nrow=nrow(pod.ac), ncol=ncol(pod.ac))

# name our rows the same thing
rownames(pod.n) <- rownames(pod.ac)

# get sample size per population
sample.n <- table(pod.gen@pop) 

# turn this into a vector
sample.sizes <- as.vector(sample.n)

# populate each row of matrix with sample sizes for pops
for(i in 1:nrow(pod.n)){
  pod.n[i,] <- sample.sizes[i]
}
pod.n <- pod.n*2 # adjust to account for loss of one allele

# calculate pairwise Fst
pod.p.fst.all <- calculate.all.pairwise.Fst(pod.ac, pod.n)

# look at global Fst
pod.p.fst <- calculate.pairwise.Fst(pod.ac, pod.n)
drop.pop <- c("1175","1575","1700","1825","1950","HU-3","HU-4","Macucaloma","Yanayacu")

# drop levels and calc distance
pop.loc.pod <- pop.loc[!pop.loc$pop %in% drop.pop,]
droplevels(pop.loc.pod)
pod.geo <- earth.dist(pop.loc.pod[c("long", "lat")], dist = FALSE)

# turn to vectors
pod.dist <- as.vector(pod.geo)
pod.gend <- as.vector(pod.p.fst.all)

# make data frame with these variables
pod.df <- cbind.data.frame(pod.dist, pod.gend)
pod.df$species <- rep("dichotomius_podalirius", nrow(pod.df))
colnames(pod.df) <- c("distance", "fst", "species")

# set up env matrix
drop.pop <- c("CC4-1175","CC1-730","HU-1575","HU-1700","HU-1825","HU-1950","HU-3","HU-4","Macucaloma","Yanayacu")
pod.env <- as.matrix(env.dist.all)
pod.env <- pod.env[!rownames(pod.env) %in% drop.pop,]
pod.env <- pod.env[,!colnames(pod.env) %in% drop.pop]
colnames(pod.env) <- NULL
rownames(pod.env) <- NULL

# run MCMC for 100K gens
MCMC(   
  counts = pod.ac,
  sample_sizes = pod.n,
  D = pod.geo,  # geographic distances
  E = pod.env,  # environmental distances
  k = 3, loci = 978,  # dimensions of the data
  delta = 0.0001,  # a small, positive, number
  aD_stp = 0.1,   # step sizes for the MCMC
  aE_stp = 0.1,
  a2_stp = 0.025,
  thetas_stp = 0.2,
  mu_stp = 0.35,
  ngen = 10000,        # number of steps (2e6)
  printfreq = 100,  # print progress (10000)
  savefreq = 100,     # save out current state
  samplefreq = 2,     # record current state for posterior (2000)
  prefix = "/Users/ethanlinck/Dropbox/scarab_migration/bedassle/pod_",   # filename prefix
  continue=FALSE,
  continuing.params=NULL)

# check shit out
show(load("/Users/ethanlinck/Dropbox/scarab_migration/bedassle/tess_MCMC_output1.Robj"))
layout(t(1:2))
plot(aD[-(1:1000)], xlab="MCMC generation", ylab="value", main="aD")
plot((aD_accept/aD_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="aD acceptance", ylim=c(0,1))
plot(as.vector(aE)[-(1:1000)], xlab="MCMC generation", ylab="value", main="aE")
plot((aE_accept/aE_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="aE acceptance", ylim=c(0,1))
plot(a2[-(1:1000)], xlab="MCMC generation", ylab="value", main="a2")
plot((a2_accept/a2_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="a2 acceptance", ylim=c(0,1))
plot(Prob[-(1:1000)], xlab="MCMC generation", main="log likelihood")
plot((mu_accept/mu_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="mu acceptance", ylim=c(0,1) )
plot((thetas_accept/thetas_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="thetas acceptance", ylim=c(0,1) )
hist((aE/aD)[-(1:1000)],breaks=100,main="posterior of aE/aD ratio")

# convert for ggplotting
pod.ratio <- aE/aD %>% as.data.frame()
pod.bed.df <- cbind.data.frame(pod.ratio, rep("dichotomius_podalirius", nrow(pod.ratio)))
colnames(pod.bed.df) <- c("ratio","species")

### e affin

# read in .vcf, turn to genpop
affin.vcf <- read.vcfR("raw_data/e_affin.LD.vcf", convertNA=TRUE)
affin.gen <- vcfR2genind(affin.vcf) 
md <- propTyped(affin.gen, by=c("both")) # confirm we have expected amount of missing data
table(md)[1]/table(md)[2] # proportion of missing sites
table(affin.gen@tab)
affin.pops <-  gsub( "_.*$", "", rownames(affin.gen@tab)) # add pops
affin.gen@pop <- as.factor(affin.pops) # make factor
affin.b <- genind2genpop(affin.gen, affin.pops) # turn into gen pop

# convert to BEDASSLE format
affin.ac <- as.matrix(affin.b@tab)
del <- seq(2, ncol(affin.ac), 2) # sequence of integers to drop non-ref allele
affin.ac <- affin.ac[,-del] 

# create equally sized matrix for sample sizes
affin.n <- matrix(nrow=nrow(affin.ac), ncol=ncol(affin.ac))

# name our rows the same thing
rownames(affin.n) <- rownames(affin.ac)

# get sample size per population
sample.n <- table(affin.gen@pop) 

# turn this into a vector
sample.sizes <- as.vector(sample.n)

# populate each row of matrix with sample sizes for pops
for(i in 1:nrow(affin.n)){
  affin.n[i,] <- sample.sizes[i]
}
affin.n <- affin.n*2 # adjust to account for loss of one allele

# calculate pairwise Fst
affin.p.fst.all <- calculate.all.pairwise.Fst(affin.ac, affin.n)

# look at global Fst
affin.p.fst <- calculate.pairwise.Fst(affin.ac, affin.n)
drop.pop <- c("1575","1700","1825","1950","HU-3","HU-4","Macucaloma","Yanayacu")

# drop levels and calc distance
pop.loc.affin <- pop.loc[!pop.loc$pop %in% drop.pop,]
droplevels(pop.loc.affin)
affin.geo <- earth.dist(pop.loc.affin[c("long", "lat")], dist = FALSE)

# turn to vectors
affin.dist <- as.vector(affin.geo)
affin.gend <- as.vector(affin.p.fst.all)

# make data frame with these variables
affin.df <- cbind.data.frame(affin.dist, affin.gend)
affin.df$species <- rep("eurysternus_affin", nrow(affin.df))
colnames(affin.df) <- c("distance", "fst", "species")

# prep env matrix
drop.pop <- c("CC1-730","HU-1575","HU-1700","HU-1825","HU-1950","HU-3","HU-4","Macucaloma","Yanayacu")
affin.env <- as.matrix(env.dist.all)
affin.env <- affin.env[!rownames(affin.env) %in% drop.pop,]
affin.env <- affin.env[,!colnames(affin.env) %in% drop.pop]
colnames(affin.env) <- NULL
rownames(affin.env) <- NULL

# run MCMC for 100K gens
MCMC(   
  counts = affin.ac,
  sample_sizes = affin.n,
  D = affin.geo,  # geographic distances
  E = affin.env,  # environmental distances
  k = 4, loci = 329,  # dimensions of the data
  delta = 0.0001,  # a small, positive, number
  aD_stp = 0.1,   # step sizes for the MCMC
  aE_stp = 0.1,
  a2_stp = 0.025,
  thetas_stp = 0.2,
  mu_stp = 0.35,
  ngen = 10000,        # number of steps (2e6)
  printfreq = 100,  # print progress (10000)
  savefreq = 100,     # save out current state
  samplefreq = 2,     # record current state for posterior (2000)
  prefix = "/Users/ethanlinck/Dropbox/scarab_migration/bedassle/affin_",   # filename prefix
  continue=FALSE,
  continuing.params=NULL)

# check shit out
show(load("/Users/ethanlinck/Dropbox/scarab_migration/bedassle/affin_MCMC_output1.Robj"))
layout(t(1:2))
plot(aD[-(1:1000)], xlab="MCMC generation", ylab="value", main="aD")
plot((aD_accept/aD_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="aD acceptance", ylim=c(0,1))
plot(as.vector(aE)[-(1:1000)], xlab="MCMC generation", ylab="value", main="aE")
plot((aE_accept/aE_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="aE acceptance", ylim=c(0,1))
plot(a2[-(1:1000)], xlab="MCMC generation", ylab="value", main="a2")
plot((a2_accept/a2_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="a2 acceptance", ylim=c(0,1))
plot(Prob[-(1:1000)], xlab="MCMC generation", main="log likelihood")
plot((mu_accept/mu_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="mu acceptance", ylim=c(0,1) )
plot((thetas_accept/thetas_moves)[-(1:1000)], xlab="MCMC generation", ylab="", main="thetas acceptance", ylim=c(0,1) )
hist((aE/aD)[-(1:1000)],breaks=100,main="posterior of aE/aD ratio")

# convert for ggplotting
affin.ratio <- aE/aD %>% as.data.frame()
affin.bed.df <- cbind.data.frame(affin.ratio, rep("eurysternus_affin", nrow(affin.ratio)))
colnames(affin.bed.df) <- c("ratio","species")

### clean up

# merge FST / distance dfs and write to csv
fst.df <- rbind.data.frame(sat.df, spec.df, tess.df, pod.df, affin.df)
write.csv(fst.df, file="~/Dropbox/scarab_migration/data/fst_dist_species.csv")

# violin plots for aE/aD
load("/Users/ethanlinck/Dropbox/scarab_migration/bedassle/satanas_MCMC_output1.Robj", verbose = FALSE)
sat.ratio <- aE/aD %>% as.data.frame()
sat.bed.df <- cbind.data.frame(sat.ratio, rep("dichotomius_satanas", nrow(sat.ratio)))
colnames(sat.bed.df) <- c("ratio","species")

load("/Users/ethanlinck/Dropbox/scarab_migration/bedassle/spec_MCMC_output1.Robj", verbose = FALSE)
spec.ratio <- aE/aD %>% as.data.frame()
spec.bed.df <- cbind.data.frame(spec.ratio, rep("deltochilum_speciocissimum", nrow(spec.ratio)))
colnames(spec.bed.df) <- c("ratio","species")

load("/Users/ethanlinck/Dropbox/scarab_migration/bedassle/tess_MCMC_output1.Robj", verbose = FALSE)
tess.ratio <- aE/aD %>% as.data.frame()
tess.bed.df <- cbind.data.frame(tess.ratio, rep("deltochilum_tesselatum", nrow(tess.ratio)))
colnames(tess.bed.df) <- c("ratio","species")

load("/Users/ethanlinck/Dropbox/scarab_migration/bedassle/tess_MCMC_output1.Robj", verbose = FALSE)
pod.ratio <- aE/aD %>% as.data.frame()
pod.bed.df <- cbind.data.frame(pod.ratio, rep("dichotomius_podalirius", nrow(pod.ratio)))
colnames(pod.bed.df) <- c("ratio","species")

load("/Users/ethanlinck/Dropbox/scarab_migration/bedassle/affin_MCMC_output1.Robj", verbose = FALSE)
affin.ratio <- aE/aD %>% as.data.frame()
affin.bed.df <- cbind.data.frame(affin.ratio, rep("eurysternus_affin", nrow(affin.ratio)))
colnames(affin.bed.df) <- c("ratio","species")

bed.df <- rbind.data.frame(sat.bed.df, spec.bed.df, tess.bed.df, pod.bed.df, affin.bed.df)
bed.df$species <- factor(bed.df$species,levels=c("eurysternus_affin","dichotomius_podalirius",
                                                 "deltochilum_tesselatum","deltochilum_speciocissimum","dichotomius_satanas"))
ggplot(bed.df, aes(y=ratio, fill=species)) + 
  theme_bw() +
  geom_boxplot(alpha=0.6, binwidth = 0.25) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values=wes_palette(n=5, name="FantasticFox1")) +
  ylim(0,100) +
  ylab("aE/aD ratio") +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank())

