# recreate plots in linck and sheldon in prep.

library(ggplot2)
library(ggmap)
library(wesanderson)
library(gridExtra)
library(cowplot)
library(patchwork)
library(reshape2)

# figure 1

# read in locality data
localities <- read.csv("/Users/ethanlinck/Dropbox/scarab_migration/data/scarab_spp_master.csv")
levels(localities$short_locality) # see levels 
localities$transect <- ifelse(grepl("CC",localities$locality),'colonso','pipeline') # add levels for transects
localities$population <-  as.factor(gsub( "_.*$", "", localities$ddocent_ID)) # add pops

# set palette
scale_fill_wes <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(wes_palette("BottleRocket2", 13, type = "continuous"), levels(localities$short_locality)), 
    ...
  )
}

scale_color_wes <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(wes_palette("BottleRocket2", 13, type = "continuous"), levels(localities$short_locality)), 
    ...
  )
}

# pick the center point for our maps, and an appropriate zoom level
pipeline <- get_map(location=c(lon=-77.85, lat=-0.615), zoom = 13, color = "bw")
colonso <- get_map(location=c(lon=-77.89, lat=-0.937), zoom = 14, color = "bw")

# subset locality data by transect
df.pipeline <- localities[localities$transect=='pipeline',]
df.colonso <- localities[localities$transect=='colonso',]

p1 <- ggmap(pipeline) +
  theme_bw() +
  geom_jitter(data = df.pipeline,
              aes(x = long, y = lat, color = short_locality, size=elevation),
              pch=1)+
  scale_color_wes(name="Locality") +
  ggtitle("Pipeline") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(size = "Elevation (m)") +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank()) +
  guides(fill = guide_legend(order=1),
         size = guide_legend(order=2))

p2 <- ggmap(colonso) +
  theme_bw() +
  geom_jitter(data = df.colonso,
             aes(x = long, y = lat, color = short_locality, size=elevation),
             pch=1)+
  scale_color_wes(name="Locality") +
  ggtitle("Colonso Chalupas") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(size = "Elevation (m)") +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank()) +
  guides(fill = guide_legend(order=1),
         size = guide_legend(order=2))


pdf("~/Dropbox/scarab_migration/manuscript/fig1.pdf",width=12,height=5)
p1 + p2
dev.off()

# figure 2

pc.df <- read.csv("~/Dropbox/scarab_migration/data/genotypes.df.csv")[-1]
pc.df$species.x <- factor(pc.df$species.x,levels=c("eurysternus_affin","dichotomius_podalirius",
                                                           "deltochilum_tesselatum","deltochilum_speciocissimum","dichotomius_satanas"))
levels(pc.df$species.x) <- c("Eurysternus affin. flocossus", "Dichotomius podalirius", "Deltochilum tessellatum", "Deltochilum speciosissimum", "Dichotomius satanas")

p3 <- ggplot(data=pc.df,aes(x=elevation,y=PC1,color=pc.df$short_locality)) + 
  theme_bw() +
  geom_jitter(pch=1,size=2.5)+
  scale_color_wes(name="locality") +
  facet_wrap(~species.x, scales="free", ncol=1) +
  xlab("Elevation (m)") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "italic"),
    legend.position="none",
    panel.grid = element_blank()) 

# read coancestry data
melted.coanc.df <- read.csv(file="~/Dropbox/scarab_migration/data/melted.coancestry.df.csv")[-1]
melted.coanc.df$species <- factor(melted.coanc.df$species,levels=c("eurysternus_affin","dichotomius_podalirius",
                                                   "deltochilum_tesselatum","deltochilum_speciocissimum","dichotomius_satanas"))
levels(melted.coanc.df$species) <- c("Eurysternus affin. flocossus", "Dichotomius podalirius", "Deltochilum tessellatum", "Deltochilum speciosissimum", "Dichotomius satanas")

p4 <- ggplot(data=melted.coanc.df,aes(x=Var1, y=Var2, fill=relative_coancestry)) + 
  theme_bw() +
  geom_tile() +
  xlab("Sample")+
  ylab("Sample")+
  scale_fill_gradient(low = "white", high = "#BA4A00",name="coancestry") +
  facet_wrap(~species, scales="free", ncol=1) +
  theme(
    strip.background = element_blank(),
    legend.position="none",
    panel.grid = element_blank(),
    strip.text = element_text(face = "italic"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()) 

# fst regression
fst.df <- read.csv(file="~/Dropbox/scarab_migration/data/fst_dist_species.csv")[-1]
fst.df <- fst.df[fst.df$distance>0,]
fst.df$fst[fst.df$fst<0] <- 0
fst.df$log_distance <- log2(fst.df$distance)
fst.df$species <- as.factor(fst.df$species)
fst.df$species <- factor(fst.df$species,levels=c("eurysternus_affin","dichotomius_podalirius",
                                                 "deltochilum_tesselatum","deltochilum_speciocissimum","dichotomius_satanas"))
levels(fst.df$species) <- c("Eurysternus affin. flocossus", "Dichotomius podalirius", "Deltochilum tessellatum", "Deltochilum speciosissimum", "Dichotomius satanas")

p5 <- ggplot(fst.df, aes(x=log_distance, y=fst)) +
  theme_bw() +
  geom_point(pch=1,size=2.5) +
  facet_wrap(~ species, ncol=1,scales ="free") +
  geom_smooth(method="lm",linetype="dashed",color="black",se = FALSE) +
  labs(y=expression(F[ST]/(1-F[ST])), x="ln(Distance)") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "italic"),
    panel.grid = element_blank())

pdf("~/Dropbox/scarab_migration/manuscript/fig2.pdf",width=8,height=10)
plot_grid(p3,p4,p5,ncol=3,labels = "AUTO")
dev.off()

# figure 3
bed.df <- read.csv("~/Dropbox/scarab_migration/data/bedassle_species.csv")[-1]
bed.df$species <- factor(bed.df$species,levels=c("eurysternus_affin","dichotomius_podalirius",
                                                 "deltochilum_tesselatum","deltochilum_speciocissimum","dichotomius_satanas"))
levels(bed.df$species) <- c("Eurysternus affin. flocossus", "Dichotomius podalirius", "Deltochilum tessellatum", "Deltochilum speciosissimum", "Dichotomius satanas")

p6 <- ggplot(bed.df, aes(y=ratio, x=species, fill=species)) + 
  theme_bw() +
  geom_boxplot(alpha=0.7) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values=wes_palette(n=5, name="Cavalcanti1")) +
  ylim(0,0.075) +
  labs(fill= "Species") +
  ylab("aE/aD ratio") +
  theme(
    strip.background = element_blank(),
    axis.title.x=element_blank(),
    legend.text = element_text(face = "italic"),
    panel.grid = element_blank())

pdf("~/Dropbox/scarab_migration/manuscript/fig3.pdf",width=8,height=5)
p6
dev.off()

# figure 4

sat.bayes <- read.csv("~/Dropbox/scarab_migration/data/sat.bayes.csv")
spec.bayes <- read.csv("~/Dropbox/scarab_migration/data/spec.bayes.csv")
tess.bayes <- read.csv("~/Dropbox/scarab_migration/data/tess.bayes.csv")
pod.bayes <- read.csv("~/Dropbox/scarab_migration/data/pod.bayes.csv")
affin.bayes <- read.csv("~/Dropbox/scarab_migration/data/affin.bayes.csv")

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

colnames(melted.bayes.df) <- c("Var1","Var2","value","species")
write.csv(melted.bayes.df, "~/Dropbox/scarab_migration/data/melted.bayes.csv")
melted.bayes.df <- read.csv("~/Dropbox/scarab_migration/data/melted.bayes.csv")[-1]
melted.bayes.df$species <- factor(melted.bayes.df$species,levels=c("eurysternus_affin","dichotomius_podalirius",
                                                                   "deltochilum_tesselatum","deltochilum_speciocissimum","dichotomius_satanas"))
levels(melted.bayes.df$species) <- c("Eurysternus affin. flocossus", "Dichotomius podalirius", "Deltochilum tessellatum", "Deltochilum speciosissimum", "Dichotomius satanas")

p7 <- ggplot(data = melted.bayes.df, aes(x=Var1, y=Var2, fill=value)) + 
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
    strip.text = element_text(face = "italic"),
    panel.grid = element_blank()) 

sat.growth <- read.csv("~/Dropbox/scarab_migration/data/sat.growth.csv")[-1]
spec.growth <- read.csv("~/Dropbox/scarab_migration/data/spec.growth.csv")[-1]
tess.growth <- read.csv("~/Dropbox/scarab_migration/data/tess.growth.csv")[-1]
pod.growth <- read.csv("~/Dropbox/scarab_migration/data/pod.growth.csv")[-1]
affin.growth <- read.csv("~/Dropbox/scarab_migration/data/affin.growth.csv")[-1]

growth.df <- rbind.data.frame(sat.growth, spec.growth, tess.growth, pod.growth, affin.growth)
growth.df$species <- as.factor(growth.df$species)
growth.df$species <- factor(growth.df$species,levels=c("eurysternus_affin","dichotomius_podalirius",
                                                      "deltochilum_tesselatum","deltochilum_speciocissimum","dichotomius_satanas"))
levels(growth.df$species) <- c("Eurysternus affin. flocossus", "Dichotomius podalirius", "Deltochilum tessellatum", "Deltochilum speciosissimum", "Dichotomius satanas")

r.df <- growth.df[growth.df$param=="r",]
theta.df <- growth.df[growth.df$param=="theta",]

p8 <- ggplot(r.df, aes(x=value, fill=species)) + 
  theme_classic() +
  geom_density(alpha=0.6) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=9),
        legend.text = element_text(face = "italic")) +
  labs(fill= "Species") +
  scale_fill_manual(values=wes_palette(n=5, name="Cavalcanti1")) +
  #theme(legend.position = "bottom") +
  #xlim(0,20) +
  #ylim(0,2000) +
  xlab("Pr(r)") +
  ylab("Density")

pdf("~/Dropbox/scarab_migration/manuscript/fig4.pdf",width=13,height=5)
plot_grid(p7,p8,ncol=2,labels = "AUTO", rel_widths = c(1.25,1))
dev.off()

# figure 5

master.df <- read.csv("~/Dropbox/scarab_migration/data/master.theta.df.csv")
master.df$species <- as.factor(master.df$species)
master.df$species <- factor(master.df$species,levels=c("eurysternus_affin","dichotomius_podalirius",
                                                       "deltochilum_tesselatum","deltochilum_speciocissimum","dichotomius_satanas"))
levels(master.df$species) <- c("Eurysternus affin. flocossus", "Dichotomius podalirius", "Deltochilum tessellatum", "Deltochilum speciosissimum", "Dichotomius satanas")

master.df$population <- as.character(master.df$population)
master.df$population[which(master.df$population=="MA2150")] <- "Macucaloma"
master.df$population[which(master.df$population=="YY2150")] <- "Yanayacu"
master.df$population[which(master.df$population=="HU3")] <- "HU-3"
master.df$population[which(master.df$population=="HU4")] <- "HU-4"
master.df$population[which(master.df$population=="CC1")] <- "730"
colnames(master.df)[3] <- "short_locality"

p9 <- ggplot(master.df, aes(x=elevation, y=theta, fill=short_locality)) +
  theme_bw() +
  geom_boxplot(alpha=0.7, width=50) +
  scale_fill_wes(name="Locality") +
  stat_summary(fun.y=median, geom="smooth", aes(group=1), linetype="dashed", color = "black") +
  facet_wrap(~ species, ncol = 1) +
  labs(y=expression(theta),x="Elevation (m)") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "italic"),
    panel.grid = element_blank())

pdf("~/Dropbox/scarab_migration/manuscript/fig5.pdf",width=6,height=7)
p9
dev.off()

