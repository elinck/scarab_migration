library(vcfR)
library(pegas)
library(sGD)
library(cowplot)
source("http://bioconductor.org/biocLite.R")
library(pcadapt)


### read in locality data
localities <- read.csv("data/scarab_spp_master.csv")
localities$population <-  as.factor(gsub( "_.*$", "", localities$ddocent_ID)) # add pops
localities <- localities[!localities$population=="CC1",]
droplevels(localities$population)

# set palette
pal <- wes_palette("Darjeeling1", 14, type = "continuous")
scale_fill_wes <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(wes_palette("Darjeeling1", length(levels(localities$population)), type = "continuous"), levels(localities$population)), 
    ...
  )
}

### estimate theta w/ observed homozygosity

# d. satanas
satanas.vcf <- read.vcfR("raw_data/d_satanas_filtered.FIL.recode.vcf", convertNA=TRUE)
satanas.gen <- vcfR2genind(satanas.vcf) 
satanas.pops <-  gsub( "_.*$", "", rownames(satanas.gen@tab)) # add pops
satanas.gen@pop <- as.factor(satanas.pops) 
sat.S <- as.loci(satanas.gen)

satlist <- list()
pops <- levels(sat.S$population)
for(i in pops){
  tmp <- sat.S[which(sat.S$population==i),]
  sum.tmp <- summary(tmp)
  tmp.theta <- as.vector(sapply(sum.tmp, function(x) theta.h(x$allele)))
  tmp.df <- cbind.data.frame(tmp.theta, rep(i, length(tmp.theta)))
  colnames(tmp.df) <- c("theta", "population")
  satlist[[i]] <- tmp.df
}

sat.df <- do.call(rbind, satlist)
sat.df <- sat.df[sat.df$theta>0,]
sat.df$species <- rep("dichotomius_satanas",nrow(sat.df))
colnames(sat.df) <- c("theta", "population","species")


# d. speciocissimum
spec.vcf <- read.vcfR("raw_data/d_spec_filtered.FIL.recode.vcf", convertNA=TRUE)
spec.gen <- vcfR2genind(spec.vcf) 
spec.pops <-  gsub( "_.*$", "", rownames(spec.gen@tab)) # add pops
spec.gen@pop <- as.factor(spec.pops) 
spec.S <- as.loci(spec.gen)
rownames(spec.S)[23] <- "YY2150_132"
spec.S$population[23] <- "YY2150"
target <- c("1575_41","1700_42","1700_43","1700_44","1700_45","1700_46","1825_47",   
                    "1825_48","1825_49","1825_50","1825_51","1825_86","1950_125","1950_126",  
                    "1950_127","1950_128","1950_129","1950_130","1950_52","1950_53","1950_54",   
                    "1950_55","HU3_107","HU3_108","HU3_109","HU3_110","HU3_111",   
                    "HU3_112","HU3_113","HU3_114","HU4_180","HU4_181","HU4_182","HU4_183",   
                    "HU4_184","HU4_185","HU4_186","HU4_187","HU4_189","MA2150_224","MA2150_226",
                    "MA2150_227","MA2150_230","MA2150_232","MA2150_233","YY2150_235","YY2150_132")
spec.S <- spec.S[match(target, rownames(spec.S)),]
spec.S <- spec.S[-c(39:47),]
levels(spec.S$population) <- droplevels(spec.S$population)
spec.S$population <-  as.factor(gsub( "_.*$", "", rownames(spec.S))) # add pops


speclist <- list()
pops <- levels(spec.S$population)
pops <- pops[-1] # drop pop w/ 1 indv
for(i in pops){
  tmp <- spec.S[which(spec.S$population==i),]
  sum.tmp <- summary(tmp)
  tmp.theta <- as.vector(sapply(sum.tmp, function(x) theta.h(x$allele)))
  tmp.df <- cbind.data.frame(tmp.theta, rep(i, length(tmp.theta)))
  colnames(tmp.df) <- c("theta", "population")
  speclist[[i]] <- tmp.df
}

spec.df <- do.call(rbind, speclist)
spec.df <- spec.df[spec.df$theta>0,]
spec.df$species <- rep("deltochilum_speciocissimum",nrow(spec.df))
colnames(spec.df) <- c("theta", "population","species")

# d. tesselatum

tess.vcf <- read.vcfR("raw_data/d_tess_filtered.FIL.recode.vcf", convertNA=TRUE)
tess.gen <- vcfR2genind(tess.vcf) 
tess.pops <-  gsub( "_.*$", "", rownames(tess.gen@tab)) # add pops
tess.gen@pop <- as.factor(tess.pops) 
tess.S <- as.loci(tess.gen)

tesslist <- list()
pops <- levels(tess.S$population)[-c(1,5)]
for(i in pops){
  tmp <- tess.S[which(tess.S$population==i),]
  sum.tmp <- summary(tmp)
  tmp.theta <- as.vector(sapply(sum.tmp, function(x) theta.h(x$allele)))
  tmp.df <- cbind.data.frame(tmp.theta, rep(i, length(tmp.theta)))
  colnames(tmp.df) <- c("theta", "population")
  tesslist[[i]] <- tmp.df
}

tess.df <- do.call(rbind, tesslist)
tess.df <- tess.df[tess.df$theta>0,]
tess.df$species <- rep("deltochilum_tesselatum",nrow(tess.df))
colnames(tess.df) <- c("theta", "population","species")

# d. podalirius

pod.vcf <- read.vcfR("raw_data/d_pod_filtered.FIL.recode.vcf", convertNA=TRUE)
pod.gen <- vcfR2genind(pod.vcf) 
pod.pops <-  gsub( "_.*$", "", rownames(pod.gen@tab)) # add pops
pod.gen@pop <- as.factor(pod.pops) 
pod.S <- as.loci(pod.gen)

podlist <- list()
pops <- levels(pod.S$population)
for(i in pops){
  tmp <- pod.S[which(pod.S$population==i),]
  sum.tmp <- summary(tmp)
  tmp.theta <- as.vector(sapply(sum.tmp, function(x) theta.h(x$allele)))
  tmp.df <- cbind.data.frame(tmp.theta, rep(i, length(tmp.theta)))
  colnames(tmp.df) <- c("theta", "population")
  podlist[[i]] <- tmp.df
}

pod.df <- do.call(rbind, podlist)
pod.df <- pod.df[pod.df$theta>0,]
pod.df$species <- rep("dichotomius_podalirius",nrow(pod.df))
colnames(pod.df) <- c("theta", "population","species")

# e. affin

affin.vcf <- read.vcfR("raw_data/e_affin_filtered.FIL.recode.vcf", convertNA=TRUE)
affin.gen <- vcfR2genind(affin.vcf) 
affin.pops <-  gsub( "_.*$", "", rownames(affin.gen@tab)) # add pops
affin.gen@pop <- as.factor(affin.pops) 
affin.S <- as.loci(affin.gen)

affinlist <- list()
pops <- levels(affin.S$population)[-2]
for(i in pops){
  tmp <- affin.S[which(affin.S$population==i),]
  sum.tmp <- summary(tmp)
  tmp.theta <- as.vector(sapply(sum.tmp, function(x) theta.h(x$allele)))
  tmp.df <- cbind.data.frame(tmp.theta, rep(i, length(tmp.theta)))
  colnames(tmp.df) <- c("theta", "population")
  affinlist[[i]] <- tmp.df
}

affin.df <- do.call(rbind, affinlist)
affin.df <- affin.df[affin.df$theta>0,]
affin.df$species <- rep("eurysternus_affin",nrow(affin.df))
colnames(affin.df) <- c("theta", "population","species")

master.df <- rbind.data.frame(sat.df, spec.df, tess.df, pod.df, affin.df)

# assign elevation
master.df$elevation <- as.character(master.df$population)
master.df$elevation[master.df$elevation=="HU3"] <- 2025
master.df$elevation[master.df$elevation=="HU4"] <- 1900
master.df$elevation[master.df$elevation=="CC1"] <- 730
master.df$elevation[master.df$elevation=="MA2150"] <- 2125
master.df$elevation[master.df$elevation=="YY2150"] <- 2150
master.df$elevation <- as.numeric(master.df$elevation)

# plot upper and lower transect separately
t.a <- ggplot(master.df[master.df$elevation>1500,], aes(x=elevation, y=theta, fill=population)) +
  theme_classic() +
  geom_boxplot(alpha=0.7) +
  scale_fill_wes() +
  stat_summary(fun.y=median, geom="smooth", aes(group=1), linetype="dashed", color = "black") +
  facet_wrap(~ species, ncol = 1)

t.b <- ggplot(master.df[master.df$elevation<1500,], aes(x=elevation, y=theta, fill=population)) +
  theme_classic() +
  geom_boxplot(alpha=0.7) +
  scale_fill_wes() +
  stat_summary(fun.y=median, geom="smooth", aes(group=1), linetype="dashed", color = "black") +
  facet_wrap(~ species, ncol = 1)

# fst regression
fst.df <- read.csv(file="~/Dropbox/scarab_migration/data/fst_dist_species.csv")[-1]
fst.df <- fst.df[fst.df$distance>0,]
fst.df$fst[fst.df$fst<0] <- 0
fst.df$log_distance <- log2(fst.df$distance)

# d. satanas 
sat.fst <- fst.df[fst.df$species=="dichotomius_satanas",]
sat.fst$log_distance <- log2(sat.fst$distance)
sat.fst$fst.adj <- sat.fst$fst/(1-sat.fst$fst)
sat.mod <- lm(fst.adj ~ log_distance, sat.fst)
summary(sat.mod)
1/as.vector(sat.mod$coefficients[2]) #[1] Inf

# d. spec
spec.fst <- fst.df[fst.df$species=="deltochilum_speciocissimum",]
spec.fst$log_distance <- log2(spec.fst$distance)
spec.fst$fst.adj <- spec.fst$fst/(1-spec.fst$fst)
spec.mod <- lm(fst.adj ~ log_distance, spec.fst)
summary(spec.mod)
1/as.vector(spec.mod$coefficients[2]) #[1] -6557.377

# d. tess
tess.fst <- fst.df[fst.df$species=="deltochilum_tesselatum",]
tess.fst$log_distance <- log2(tess.fst$distance)
tess.fst$fst.adj <- tess.fst$fst/(1-tess.fst$fst)
tess.mod <- lm(fst.adj ~ log_distance, tess.fst)
summary(tess.mod)
1/as.vector(tess.mod$coefficients[2]) # [1] -360.0897

# d. pod
pod.fst <- fst.df[fst.df$species=="dichotomius_podalirius",]
pod.fst$log_distance <- log2(pod.fst$distance)
pod.fst$fst.adj <- pod.fst$fst/(1-pod.fst$fst)
pod.mod <- lm(fst.adj ~ log_distance, pod.fst)
summary(pod.mod)
1/as.vector(pod.mod$coefficients[2]) #[1] -256.9938

# e. affin
affin.fst <- fst.df[fst.df$species=="eurysternus_affin",]
affin.fst$log_distance <- log2(affin.fst$distance)
affin.fst$fst.adj <- affin.fst$fst/(1-affin.fst$fst)
affin.mod <- lm(fst.adj ~ log_distance, affin.fst)
summary(affin.mod)
1/as.vector(pod.mod$coefficients[2]) #[1] -256.9938

upper <- c("dichotomius_satanas", "deltochilum_speciocissimum","deltochilum_tesselatum")
fst.upper <- fst.df[fst.df$species %in% upper,]
fst.lower <- fst.df[!fst.df$species %in% upper,]

f.a <- ggplot(fst.upper, aes(x=log_distance, y=fst)) +
  theme_classic() +
  geom_point() +
  facet_wrap(~ species, ncol=1)

f.b <- ggplot(fst.lower, aes(x=log_distance, y=fst)) +
  theme_classic() +
  geom_point() +
  facet_wrap(~ species, ncol=1)

plot_grid(t.a, f.a, labels = "AUTO", align = 'h', label_size = 12, ncol = 2, rel_widths = c(1,0.75))
plot_grid(t.b, f.b, labels = "AUTO", align = 'h', label_size = 12, ncol = 2, rel_widths = c(1,0.75))

write.csv(master.df, "~/Dropbox/scarab_migration/data/master.theta.df.csv")
master.df <- read.csv("~/Dropbox/scarab_migration/data/master.theta.df.csv")[-1]

master.df$species <- as.factor(master.df$species)
master.df$species <- factor(master.df$species,levels=c("eurysternus_affin","dichotomius_podalirius",
                                                       "deltochilum_tesselatum","deltochilum_speciocissimum","dichotomius_satanas"))
fst.df$species <- as.factor(fst.df$species)
fst.df$species <- factor(fst.df$species,levels=c("eurysternus_affin","dichotomius_podalirius",
                                                 "deltochilum_tesselatum","deltochilum_speciocissimum","dichotomius_satanas"))

# plot
t.c <- ggplot(master.df, aes(x=elevation, y=theta, fill=population)) +
  theme_bw() +
  geom_boxplot(alpha=0.7, width=50) +
  scale_fill_wes() +
  stat_summary(fun.y=median, geom="smooth", aes(group=1), linetype="dashed", color = "black") +
  facet_wrap(~ species, ncol = 1) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank())

f.c <- ggplot(fst.df, aes(x=log_distance, y=fst)) +
  theme_bw() +
  geom_point(pch=1,size=2) +
  facet_wrap(~ species, ncol=1) +
  stat_summary(method=lm, geom="smooth", aes(group=1), linetype="dashed", color = "black") +
  ylab("fst/(1-fst)") +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank()) 

plot_grid(t.c, f.c, labels = "AUTO", align = 'h', label_size = 12, ncol = 2, rel_widths = c(1,0.5))
