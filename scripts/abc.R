library(coala)
library(abc)
library(tidyverse)
library(wesanderson)
library(reshape2)

### d. satanas

# calculate SFS from data
satanas.vcf <- read.vcfR("~/Dropbox/scarab_migration/raw_data/d.satanas.abc.FIL.recode.vcf", convertNA=TRUE)
satanas.gen <- vcfR2genlight(satanas.vcf) 
sat.sfs <- table(glSum(satanas.gen))
sat.sfs <- as.vector(sat.sfs[1:98])

# define two models
sat.null.model <- coal_model(99, 50, 3) +
  feat_mutation(par_prior("theta", runif(1, 1, 6))) +
  sumstat_sfs()

sat.growth.model <- coal_model(99, 50, 3) +
  feat_mutation(par_prior("theta", runif(1, 1, 6))) +
  feat_growth(par_prior("r", runif(1, 0, 1))) +
  sumstat_sfs()

# simulate data
sat.null.sim <- simulate(sat.null.model, nsim = 2000, seed = 69)
sat.growth.sim <- simulate(sat.growth.model, nsim = 2000, seed = 32)

# create params and sumstts for abc
sat.null.param <- create_abc_param(sat.null.sim, sat.null.model)
sat.null.sumstat <- create_abc_sumstat(sat.null.sim, sat.null.model)
sat.growth.param <- create_abc_param(sat.growth.sim, sat.growth.model)
sat.growth.sumstat <- create_abc_sumstat(sat.growth.sim, sat.growth.model)

# look at posterior distributions of params
sat.null.post <- abc(sat.sfs, sat.null.param, sat.null.sumstat, 0.05, method = "rejection")
sat.growth.post <- abc(sat.sfs, sat.growth.param, sat.growth.sumstat, 0.05, method = "rejection")

# prep data for model 
sat.sumstat.merge <- rbind(sat.null.sumstat, sat.growth.sumstat)
sat.index <- c(rep("null",2000),rep("growth",2000))

# check ability to distinguish models
sat.cv.modsel <- cv4postpr(sat.index, sat.sumstat.merge, nval=5, tols=c(.05,.1), method="rejection")
summary(cv.modsel)

# model test
sat.model <- postpr(target=sat.sfs, sumstat=sat.sumstat.merge, index=as.vector(sat.index), tol=.05, method="rejection")
summary(sat.model)

# plotting
sat.null.plot <- as.data.frame(sat.null.post$unadj.values)
colnames(sat.null.plot) <- c("theta")
ggplot(sat.null.plot, aes(x=theta)) +
  theme_classic() +
  geom_density(fill="blue")

sat.growth.plot <- as.data.frame(sat.growth.post$unadj.values)
sat.growth.theta <- cbind.data.frame(sat.growth.plot$theta, rep("theta", nrow(sat.growth.plot)))
colnames(sat.growth.theta) <- c("value","param")
sat.growth.r <- cbind.data.frame(sat.growth.plot$r, rep("r", nrow(sat.growth.plot)))
colnames(sat.growth.r) <- c("value","param")
sat.growth.df <- rbind.data.frame(sat.growth.theta,sat.growth.r)
sat.growth.df$species <- rep("dichotomius_satanas", nrow(sat.growth.df))

# plot both params
ggplot(sat.growth.df, aes(x=value)) +
  theme_classic() +
  geom_density(fill="blue", alpha=0.7) +
  facet_wrap(~ param, scales = "free")

# bayes factor comparison
sat.sum <- summary(sat.model)
sat.bayes <- sat.sum$BayesF
plot(sat.bayes)

### d. spec

# calculate SFS from data
spec.vcf <- read.vcfR("~/Dropbox/scarab_migration/raw_data/d.spec.abc.FIL.recode.vcf", convertNA=TRUE)
spec.gen <- vcfR2genlight(spec.vcf) 
spec.sfs <- table(glSum(spec.gen))
spec.sfs <- as.vector(spec.sfs[2:47])

# define two models
spec.null.model <- coal_model(47, 50, 3) +
  feat_mutation(par_prior("theta", runif(1, 1, 6))) +
  sumstat_sfs()

spec.growth.model <- coal_model(47, 50, 3) +
  feat_mutation(par_prior("theta", runif(1, 1, 6))) +
  feat_growth(par_prior("r", runif(1, 0, 1))) +
  sumstat_sfs()

# simulate data
spec.null.sim <- simulate(spec.null.model, nsim = 2000, seed = 14)
spec.growth.sim <- simulate(spec.growth.model, nsim = 2000, seed = 22)

# create params and sumstts for abc
spec.null.param <- create_abc_param(spec.null.sim, spec.null.model)
spec.null.sumstat <- create_abc_sumstat(spec.null.sim, spec.null.model)
spec.growth.param <- create_abc_param(spec.growth.sim, spec.growth.model)
spec.growth.sumstat <- create_abc_sumstat(spec.growth.sim, spec.growth.model)

# look at posterior distributions of params
spec.null.post <- abc(spec.sfs, spec.null.param, spec.null.sumstat, 0.05, method = "rejection")
spec.growth.post <- abc(spec.sfs, spec.growth.param, spec.growth.sumstat, 0.05, method = "rejection")

# prep data for model 
spec.sumstat.merge <- rbind(spec.null.sumstat, spec.growth.sumstat)
spec.index <- c(rep("null",2000),rep("growth",2000))

# check ability to distinguish models
spec.cv.modsel <- cv4postpr(spec.index, spec.sumstat.merge, nval=5, tols=c(.05,.1), method="rejection")
summary(spec.cv.modsel)

# model test
spec.model <- postpr(target=spec.sfs, sumstat=spec.sumstat.merge, index=as.vector(spec.index), tol=.05, method="rejection")
summary(spec.model)

# plotting
spec.null.plot <- as.data.frame(spec.null.post$unadj.values)
colnames(spec.null.plot) <- c("theta")
ggplot(spec.null.plot, aes(x=theta)) +
  theme_classic() +
  geom_density(fill="blue")

spec.growth.plot <- as.data.frame(spec.growth.post$unadj.values)
spec.growth.theta <- cbind.data.frame(spec.growth.plot$theta, rep("theta", nrow(spec.growth.plot)))
colnames(spec.growth.theta) <- c("value","param")
spec.growth.r <- cbind.data.frame(spec.growth.plot$r, rep("r", nrow(spec.growth.plot)))
colnames(spec.growth.r) <- c("value","param")
spec.growth.df <- rbind.data.frame(spec.growth.theta,spec.growth.r)
spec.growth.df$species <- rep("deltochilum_speciocissimum", nrow(spec.growth.df))

# plot both params
ggplot(spec.growth.df, aes(x=value)) +
  theme_classic() +
  geom_density(fill="blue", alpha=0.7) +
  facet_wrap(~ param, scales = "free")

# bayes factor comparison
spec.sum <- summary(spec.model)
spec.bayes <- spec.sum$BayesF
plot(spec.bayes)

### d. tess

# calculate SFS from data
tess.vcf <- read.vcfR("~/Dropbox/scarab_migration/raw_data/d.tess.abc.FIL.recode.vcf", convertNA=TRUE)
tess.gen <- vcfR2genlight(tess.vcf) 
tess.sfs <- table(glSum(tess.gen))
tess.sfs <- as.vector(tess.sfs[2:19])

# define two models
tess.null.model <- coal_model(19, 50, 3) +
  feat_mutation(par_prior("theta", runif(1, 1, 6))) +
  sumstat_sfs()

tess.growth.model <- coal_model(19, 50, 3) +
  feat_mutation(par_prior("theta", runif(1, 1, 6))) +
  feat_growth(par_prior("r", runif(1, 0, 1))) +
  sumstat_sfs()

# simulate data
tess.null.sim <- simulate(tess.null.model, nsim = 2000, seed = 19)
tess.growth.sim <- simulate(tess.growth.model, nsim = 2000, seed = 32)

# create params and sumstts for abc
tess.null.param <- create_abc_param(tess.null.sim, tess.null.model)
tess.null.sumstat <- create_abc_sumstat(tess.null.sim, tess.null.model)
tess.growth.param <- create_abc_param(tess.growth.sim, tess.growth.model)
tess.growth.sumstat <- create_abc_sumstat(tess.growth.sim, tess.growth.model)

# look at posterior distributions of params
tess.null.post <- abc(tess.sfs, tess.null.param, tess.null.sumstat, 0.05, method = "rejection")
tess.growth.post <- abc(tess.sfs, tess.growth.param, tess.growth.sumstat, 0.05, method = "rejection")

# prep data for model 
tess.sumstat.merge <- rbind(tess.null.sumstat, tess.growth.sumstat)
tess.index <- c(rep("null",2000),rep("growth",2000))

# check ability to distinguish models
tess.cv.modsel <- cv4postpr(tess.index, tess.sumstat.merge, nval=5, tols=c(.05,.1), method="rejection")
summary(tess.cv.modsel)

# model test
tess.model <- postpr(target=tess.sfs, sumstat=tess.sumstat.merge, index=as.vector(tess.index), tol=.05, method="rejection")
summary(tess.model)

# plotting
tess.null.plot <- as.data.frame(tess.null.post$unadj.values)
colnames(tess.null.plot) <- c("theta")
ggplot(tess.null.plot, aes(x=theta)) +
  theme_classic() +
  geom_density(fill="blue")

tess.growth.plot <- as.data.frame(tess.growth.post$unadj.values)
tess.growth.theta <- cbind.data.frame(tess.growth.plot$theta, rep("theta", nrow(tess.growth.plot)))
colnames(tess.growth.theta) <- c("value","param")
tess.growth.r <- cbind.data.frame(tess.growth.plot$r, rep("r", nrow(tess.growth.plot)))
colnames(tess.growth.r) <- c("value","param")
tess.growth.df <- rbind.data.frame(tess.growth.theta,tess.growth.r)
tess.growth.df$species <- rep("deltochilum_tesselatum", nrow(tess.growth.df))

# plot both params
ggplot(tess.growth.df, aes(x=value)) +
  theme_classic() +
  geom_density(fill="blue", alpha=0.7) +
  facet_wrap(~ param, scales = "free")

# bayes factor comparison
tess.sum <- summary(tess.model)
tess.bayes <- tess.sum$BayesF
plot(tess.bayes)

### d. pod

# calculate SFS from data
pod.vcf <- read.vcfR("~/Dropbox/scarab_migration/raw_data/d.pod.abc.FIL.recode.vcf", convertNA=TRUE)
pod.gen <- vcfR2genlight(pod.vcf) 
pod.sfs <- table(glSum(pod.gen))
pod.sfs <- as.vector(pod.sfs[2:26])

# define two models
pod.null.model <- coal_model(26, 50, 3) +
  feat_mutation(par_prior("theta", runif(1, 1, 6))) +
  sumstat_sfs()

pod.growth.model <- coal_model(26, 50, 3) +
  feat_mutation(par_prior("theta", runif(1, 1, 6))) +
  feat_growth(par_prior("r", runif(1, 0, 1))) +
  sumstat_sfs()

# simulate data
pod.null.sim <- simulate(pod.null.model, nsim = 2000, seed = 9)
pod.growth.sim <- simulate(pod.growth.model, nsim = 2000, seed = 51)

# create params and sumstts for abc
pod.null.param <- create_abc_param(pod.null.sim, pod.null.model)
pod.null.sumstat <- create_abc_sumstat(pod.null.sim, pod.null.model)
pod.growth.param <- create_abc_param(pod.growth.sim, pod.growth.model)
pod.growth.sumstat <- create_abc_sumstat(pod.growth.sim, pod.growth.model)

# look at posterior distributions of params
pod.null.post <- abc(pod.sfs, pod.null.param, pod.null.sumstat, 0.05, method = "rejection")
pod.growth.post <- abc(pod.sfs, pod.growth.param, pod.growth.sumstat, 0.05, method = "rejection")

# prep data for model 
pod.sumstat.merge <- rbind(pod.null.sumstat, pod.growth.sumstat)
pod.index <- c(rep("null",2000),rep("growth",2000))

# check ability to distinguish models
pod.cv.modsel <- cv4postpr(pod.index, pod.sumstat.merge, nval=5, tols=c(.05,.1), method="rejection")
summary(pod.cv.modsel)

# model test
pod.model <- postpr(target=pod.sfs, sumstat=pod.sumstat.merge, index=as.vector(pod.index), tol=.05, method="rejection")
summary(pod.model)

# plotting
pod.null.plot <- as.data.frame(pod.null.post$unadj.values)
colnames(pod.null.plot) <- c("theta")
ggplot(pod.null.plot, aes(x=theta)) +
  theme_classic() +
  geom_density(fill="blue")

pod.growth.plot <- as.data.frame(pod.growth.post$unadj.values)
pod.growth.theta <- cbind.data.frame(pod.growth.plot$theta, rep("theta", nrow(pod.growth.plot)))
colnames(pod.growth.theta) <- c("value","param")
pod.growth.r <- cbind.data.frame(pod.growth.plot$r, rep("r", nrow(pod.growth.plot)))
colnames(pod.growth.r) <- c("value","param")
pod.growth.df <- rbind.data.frame(pod.growth.theta,pod.growth.r)
pod.growth.df$species <- rep("dichotomius_podalirius", nrow(pod.growth.df))

# plot both params
ggplot(pod.growth.df, aes(x=value)) +
  theme_classic() +
  geom_density(fill="blue", alpha=0.7) +
  facet_wrap(~ param, scales = "free")

# bayes factor comparison
pod.sum <- summary(pod.model)
pod.bayes <- pod.sum$BayesF
plot(pod.bayes)

### d. affin

# calculate SFS from data
affin.vcf <- read.vcfR("~/Dropbox/scarab_migration/raw_data/e.affin.abc.FIL.recode.vcf", convertNA=TRUE)
affin.gen <- vcfR2genlight(affin.vcf) 
affin.sfs <- table(glSum(affin.gen))
affin.sfs <- as.vector(affin.sfs[2:24])

# define two models
affin.null.model <- coal_model(24, 50, 3) +
  feat_mutation(par_prior("theta", runif(1, 1, 6))) +
  sumstat_sfs()

affin.growth.model <- coal_model(24, 50, 3) +
  feat_mutation(par_prior("theta", runif(1, 1, 6))) +
  feat_growth(par_prior("r", runif(1, 0, 1))) +
  sumstat_sfs()

# simulate data
affin.null.sim <- simulate(affin.null.model, nsim = 5000, seed = 3)
affin.growth.sim <- simulate(affin.growth.model, nsim = 5000, seed = 71)

# create params and sumstts for abc
affin.null.param <- create_abc_param(affin.null.sim, affin.null.model)
affin.null.sumstat <- create_abc_sumstat(affin.null.sim, affin.null.model)
affin.growth.param <- create_abc_param(affin.growth.sim, affin.growth.model)
affin.growth.sumstat <- create_abc_sumstat(affin.growth.sim, affin.growth.model)

# look at posterior distributions of params
affin.null.post <- abc(affin.sfs, affin.null.param, affin.null.sumstat, 0.05, method = "rejection")
affin.growth.post <- abc(affin.sfs, affin.growth.param, affin.growth.sumstat, 0.05, method = "rejection")

# prep data for model 
affin.sumstat.merge <- rbind(affin.null.sumstat, affin.growth.sumstat)
affin.index <- c(rep("null",5000),rep("growth",5000))

# check ability to distinguish models
affin.cv.modsel <- cv4postpr(affin.index, affin.sumstat.merge, nval=, tols=c(.05,.1), method="rejection")
summary(affin.cv.modsel)

# model test
affin.model <- postpr(target=affin.sfs, sumstat=affin.sumstat.merge, index=as.vector(affin.index), tol=.05, method="rejection")
summary(affin.model)

# plotting
affin.null.plot <- as.data.frame(affin.null.post$unadj.values)
colnames(affin.null.plot) <- c("theta")
ggplot(affin.null.plot, aes(x=theta)) +
  theme_classic() +
  geom_density(fill="blue")

affin.growth.plot <- as.data.frame(affin.growth.post$unadj.values)
affin.growth.theta <- cbind.data.frame(affin.growth.plot$theta, rep("theta", nrow(affin.growth.plot)))
colnames(affin.growth.theta) <- c("value","param")
affin.growth.r <- cbind.data.frame(affin.growth.plot$r, rep("r", nrow(affin.growth.plot)))
colnames(affin.growth.r) <- c("value","param")
affin.growth.df <- rbind.data.frame(affin.growth.theta,affin.growth.r)
affin.growth.df$species <- rep("eurysternus_affin", nrow(affin.growth.df))

# plot both params
ggplot(affin.growth.df, aes(x=value)) +
  theme_classic() +
  geom_density(fill="blue", alpha=0.7) +
  facet_wrap(~ param, scales = "free")

# bayes factor comparison
affin.sum <- summary(affin.model)
affin.bayes <- affin.sum$BayesF
plot(affin.bayes)

growth.df <- rbind.data.frame(sat.growth.df, spec.growth.df, tess.growth.df, pod.growth.df, affin.growth.df)
growth.df$species <- factor(growth.df$species,levels=c("eurysternus_affin","dichotomius_podalirius",
                                                                   "deltochilum_tesselatum","deltochilum_speciocissimum","dichotomius_satanas"))
write.csv(growth.df, "~/Dropbox/scarab_migration/data/growth.csv")
r.df <- growth.df[growth.df$param=="r",]
theta.df <- growth.df[growth.df$param=="theta",]

ggplot(r.df, aes(x=value, fill=species)) + 
  theme_classic() +
  geom_density(alpha=0.6) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=9)) +
  scale_fill_manual(values=wes_palette(n=5, name="FantasticFox1")) +
  #xlim(0,20) +
  #ylim(0,2000) +
  xlab("Pr(r)")

ggplot(theta.df, aes(x=value, fill=species)) + 
  theme_classic() +
  geom_density(alpha=0.6) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=9)) +
  scale_fill_manual(values=wes_palette(n=5, name="FantasticFox1")) +
  #xlim(0,20) +
  #ylim(0,2000) +
  xlab("Pr(r)")

ggplot(growth.df, aes(x=value, fill=species)) + 
  theme_bw() +
  geom_density(alpha=0.6) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=9)) +
  scale_fill_manual(values=wes_palette(n=5, name="FantasticFox1")) +
  #xlim(0,20) +
  #ylim(0,2000) +
  facet_wrap(~param, scales="free") +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank()) 
  
melted.sat.bayes <- melt(sat.bayes)
melted.spec.bayes <- melt(spec.bayes)
melted.tess.bayes <- melt(tess.bayes)
melted.pod.bayes <- melt(pod.bayes)
melted.affin.bayes <- melt(affin.bayes)
melted.sat.bayes$species <- rep("dichotomius_satanas", 4)
melted.spec.bayes$species <- rep("deltochilum_speciocissimum", 4)
melted.tess.bayes$species <- rep("deltochilum_tesselatum", 4)
melted.pod.bayes$species <- rep("dichotomius_podalirius", 4)
melted.affin.bayes$species <- rep("eurysternus_affin", 4)

melted.bayes.df <- rbind.data.frame(melted.sat.bayes,melted.spec.bayes,
                                    melted.tess.bayes,melted.pod.bayes,
                                    melted.affin.bayes)
write.csv(melted.bayes.df, "~/Dropbox/scarab_migration/data/melted.bayes.csv")


melted.bayes.df$species <- factor(melted.bayes.df$species,levels=c("eurysternus_affin","dichotomius_podalirius",
                                                 "deltochilum_tesselatum","deltochilum_speciocissimum","dichotomius_satanas"))

ggplot(data = melted.bayes.df, aes(x=Var1, y=Var2, fill=value)) + 
  theme_classic() +
  geom_tile() +
  xlab(element_blank())+
  ylab(element_blank())+
  scale_fill_gradient(name = "Bayes factor",
                      low = "#FFFFFF",
                      high = "#BA4A00") +
  facet_wrap(~species) +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank()) 


