### plots for publication

library(ggplot2)
library(viridis)
library(reshape2)
library(wesanderson)
library(cowplot)

### fineRADstructure 

# d. satanas

# load and format correlation matrix
satanas <- as.matrix(read.table("~/Dropbox/scarab_migration/raw_data/d_satanas_painter_chunks.out"))
colnames(satanas) <- NULL
names <- as.vector(satanas[,1])[-1]
satanas <- satanas[-1,]
satanas <- satanas[,-1]
colnames(satanas) <- names
rownames(satanas) <- names
class(satanas) <- "numeric"

# melt df
melted_satanas<- melt(satanas)
melted_satanas$pop.x <- gsub( "_.*$", "", melted_satanas$Var1) %>% as.factor()
melted_satanas$pop.y <- gsub( "_.*$", "", melted_satanas$Var2) %>% as.factor()

# plot correlation
pop.cols <- as.character(pal)
names(pop.cols) <- levels(melted_satanas$pop.x)
a <- ggplot(data = melted_satanas, aes(x=Var1, y=Var2, fill=value)) + 
  theme_classic() +
  geom_tile() +
  xlab(element_blank())+
  ylab(element_blank())+
  scale_fill_gradient(name = "coancestry",
                      low = "#FFFFFF",
                      high = "#BA4A00") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

# d. speciocissimum

# load and format correlation matrix
spec <- as.matrix(read.table("~/Dropbox/scarab_migration/raw_data/d_spec_painter_chunks.out"))
colnames(spec) <- NULL
names <- as.vector(spec[,1])[-1]
spec <- spec[-1,]
spec <- spec[,-1]
colnames(spec) <- names
rownames(spec) <- names
class(spec) <- "numeric"

# melt df
melted_spec<- melt(spec)
melted_spec$pop.x <- gsub( "_.*$", "", melted_spec$Var1) %>% as.factor()
melted_spec$pop.y <- gsub( "_.*$", "", melted_spec$Var2) %>% as.factor()
#melted_spec<-melted_spec[!(melted_spec$Var1=="2150_132"),]
#melted_spec<-melted_spec[!(melted_spec$Var2=="2150_132"),]

# plot correlation
b <- ggplot(data = melted_spec, aes(x=Var1, y=Var2, fill=value)) + 
  theme_classic() +
  geom_tile() +
  xlab(element_blank())+
  ylab(element_blank())+
  scale_fill_gradient(name = "coancestry",
                      low = "#FFFFFF",
                      high = "#BA4A00") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

# d tesselatum

# load and format correlation matrix
tess <- as.matrix(read.table("~/Dropbox/scarab_migration/raw_data/d_tess_painter_chunks.out"))
colnames(tess) <- NULL
names <- as.vector(tess[,1])[-1]
tess <- tess[-1,]
tess <- tess[,-1]
colnames(tess) <- names
rownames(tess) <- names
class(tess) <- "numeric"

# melt df
melted_tess<- melt(tess)
melted_tess$pop.x <- gsub( "_.*$", "", melted_tess$Var1) %>% as.factor()
melted_tess$pop.y <- gsub( "_.*$", "", melted_tess$Var2) %>% as.factor()

# plot correlation
c <- ggplot(data = melted_tess, aes(x=Var1, y=Var2, fill=value)) + 
  theme_classic() +
  geom_tile() +
  xlab(element_blank())+
  ylab(element_blank())+
  scale_fill_gradient(name = "coancestry",
                      low = "#FFFFFF",
                      high = "#BA4A00") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

# d. podalirius

# load and format correlation matrix
pod<- as.matrix(read.table("~/Dropbox/scarab_migration/raw_data/d_pod_painter_chunks.out"))
colnames(pod) <- NULL
names <- as.vector(pod[,1])[-1]
pod <- pod[-1,]
pod <- pod[,-1]
colnames(pod) <- names
rownames(pod) <- names
class(pod) <- "numeric"

# melt df
melted_pod<- melt(pod)
melted_pod$pop.x <- gsub( "_.*$", "", melted_pod$Var1) %>% as.factor()
melted_pod$pop.y <- gsub( "_.*$", "", melted_pod$Var2) %>% as.factor()

# plot correlation
d <- ggplot(data = melted_pod, aes(x=Var1, y=Var2, fill=value)) + 
  theme_classic() +
  geom_tile() +
  xlab(element_blank())+
  ylab(element_blank())+
  scale_fill_gradient(name = "coancestry",
                      low = "#FFFFFF",
                      high = "#BA4A00") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

# load and format correlation matrix
affin <- as.matrix(read.table("~/Dropbox/scarab_migration/raw_data/e_affin_painter_chunks.out"))
colnames(affin) <- NULL
names <- as.vector(affin[,1])[-1]
affin <- affin[-1,]
affin <- affin[,-1]
colnames(affin) <- names
rownames(affin) <- names
class(affin) <- "numeric"

# melt df
melted_affin<- melt(affin)
melted_affin$pop.x <- gsub( "_.*$", "", melted_affin$Var1) %>% as.factor()
melted_affin$pop.y <- gsub( "_.*$", "", melted_affin$Var2) %>% as.factor()

# plot correlation
e <- ggplot(data = melted_affin, aes(x=Var1, y=Var2, fill=value)) + 
  theme_classic() +
  geom_tile() +
  xlab(element_blank())+
  ylab(element_blank())+
  scale_fill_gradient(name = "coancestry",
                      low = "#FFFFFF",
                      high = "#BA4A00") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

plot_grid(a, b, c, d, e, labels = "AUTO", align = 'h', label_size = 12, nrow = 2, label_x = 0.5)

melted_satanas$relative_coancestry <- melted_satanas$value/max(melted_satanas$value)
melted_spec$relative_coancestry <- melted_spec$value/max(melted_spec$value)
melted_tess$relative_coancestry <- melted_tess$value/max(melted_tess$value)
melted_pod$relative_coancestry <- melted_pod$value/max(melted_pod$value)
melted_affin$relative_coancestry <- melted_affin$value/max(melted_affin$value)

melted_satanas$species <- rep("dichotomius_satanas", nrow(melted_satanas))
melted_spec$species <- rep("deltochilum_speciocissimum", nrow(melted_spec))
melted_tess$species <- rep("deltochilum_tesselatum", nrow(melted_tess))
melted_pod$species <- rep("dichotomius_podalirius", nrow(melted_pod))
melted_affin$species <- rep("eurysternus_affin", nrow(melted_affin))

melted.coanc.df <- rbind.data.frame(melted_satanas,melted_spec,melted_tess,melted_pod,melted_affin)
melted.coanc.df$species <- factor(melted.coanc.df$species,levels=c("eurysternus_affin","dichotomius_podalirius",
                                                           "deltochilum_tesselatum","deltochilum_speciocissimum","dichotomius_satanas"))

l <- ggplot(data=melted.coanc.df,aes(x=Var1, y=Var2, fill=relative_coancestry)) + 
  theme_bw() +
  geom_tile() +
  xlab(element_blank())+
  ylab(element_blank())+
  scale_fill_gradient(name = "coancestry",
                      low = "#FFFFFF",
                      high = "#BA4A00") +
  facet_wrap(~species, scales="free") +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()) 






