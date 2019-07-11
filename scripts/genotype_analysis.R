# genotype analysis for scarab sp.

library(conStruct)
library(vcfR)
library(adegenet)
library(ggplot2)
library(poppr)
library(ape)
library(igraph)
library(RColorBrewer)
library(pegas)
library(naniar)
library(reshape2)

setwd("/Users/ethanlinck/Dropbox/scarab_migration/")

### analyze d. tesselatum ###

# read in .vcf and sample / pop names; convert to genlight
spec.vcf <- read.vcfR("ipyrad/d_speciossimum_outfiles/d_speciossimum.vcf")

# fix ipyrad vcf missing data issue
spec.vcf@gt <- 
  spec.vcf@gt %>% 
  as.tibble() %>% 
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
spec.vcf <- spec.vcf[spec.miss < 0.15, ]
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

# quick visualization
ggplot(data=spec.pc,aes(x=PC1,y=PC2,col=spec.pc$pop)) + 
  geom_text(aes(label=samples))

# drop bad samples. bad!
spec.pc.2 <- spec.pc[!spec.pc$pop=="MA2150" & !spec.pc$pop=="2150" & !spec.pc$pop=="YY2150",]

# quick visualization -- still panmixia!
ggplot(data=spec.pc.2,aes(x=PC1,y=PC2,col=spec.pc.2$pop)) + 
  geom_text(aes(label=rownames(spec.pc.2)))



