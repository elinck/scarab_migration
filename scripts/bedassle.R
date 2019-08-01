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

# set wd and read in locality data
setwd("/Users/ethanlinck/Dropbox/scarab_migration/")
localities <- read.csv("data/scarab_spp_master.csv")
localities$short_locality[which(localities$short_locality=="730")] <- 800
pop.loc <- cbind.data.frame(unique(localities$long),unique(localities$lat),unique(localities$short_locality))
colnames(pop.loc) <- c("long", "lat", "pop")
geodist <- earth.dist(pop.loc[c("long", "lat")], dist = FALSE)

## d satanas

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
adj.fst <- satanas.p.fst.all / (1 - satanas.p.fst.all)

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

# loop to figure out magnitude
mag <- vector()
for(i in 1:8){
  tmp <- abs(as.vector(satanas.n[i,i]-(satanas.n[,1])))
  mag <- append(mag, tmp)
}

pop <- vector()
for(i in 1:8){
  tmp <- abs(as.vector(satanas.n[i,i]-(satanas.n[,1])))
  mag <- append(mag, tmp)
}


# make data frame with these variables
sat.df <- cbind.data.frame(sat.dist, sat.gen, mag)
colnames(sat.df) <- c("distance", "fst", "mag")

# plot IBD etc
a <- ggplot(data=sat.df, aes(x=distance, y=fst, size = mag)) +
  theme_classic() +
  labs(size = "delta(n)") +
  geom_smooth(method="lm", se=FALSE, color = "black", linetype="dashed") +
  geom_point(fill="#B578EC",color="black",alpha=0.7,pch=21)

b <- ggplot(data=sat.df, aes(x=mag, y=fst, size=distance)) +
  theme_classic() +
  labs(size = "distance") +
  geom_smooth(method="lm", se=FALSE, color = "black", linetype="dashed") +
  geom_point(fill="#B578EC",color="black",alpha=0.7,pch=21) +
  xlab("delta(n)")

# plot locally
a + b

# linear models
mod1 <- lm(fst ~ distance, data=sat.df) 
mod2 <- lm(fst ~ mag, data=sat.df) 
summary(mod1)
summary(mod2)

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

### BEDASSLE for d. satanas

satanas.env <- as.matrix(env.dist.upper)
colnames(satanas.env) <- NULL
rownames(satanas.env) <- NULL

# run MCMC for 100K gens
MCMC(   
  counts = satanas.ac,
  sample_sizes = satanas.n,
  D = satanas.geo,  # geographic distances
  E = satanas.env,  # environmental distances
  k = 8, loci = 1156,  # dimensions of the data
  delta = 0.0001,  # a small, positive, number
  aD_stp = 0.01,   # step sizes for the MCMC
  aE_stp = 0.01,
  a2_stp = 0.025,
  thetas_stp = 0.2,
  mu_stp = 0.35,
  ngen = 100000,        # number of steps (2e6)
  printfreq = 100000,  # print progress (10000)
  savefreq = 1000,     # save out current state
  samplefreq = 5,     # record current state for posterior (2000)
  prefix = "/Users/ethanlinck/Dropbox/scarab_migration/bedassle_test",   # filename prefix
  continue=FALSE,
  continuing.params=NULL)

# check shit out
show(load("/Users/ethanlinck/Dropbox/scarab_migration/bedassle_testMCMC_output1.Robj"))
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
colnames(sat.ratio) <- c("ratio")

# plot w/ ggplot2
ggplot(sat.ratio, aes(x = ratio)) +
  theme_classic() +
  geom_histogram(binwidth = 0.000005) +
  xlim(0, 2e-4) +
  xlab("aE/aD")


