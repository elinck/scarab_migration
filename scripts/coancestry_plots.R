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



