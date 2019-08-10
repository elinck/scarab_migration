### plot sampling localities

library(ggmap)
library(wesanderson)
library(gridExtra)

# read in locality data
localities <- read.csv("/Users/ethanlinck/Dropbox/scarab_migration/data/scarab_spp_master.csv")
levels(localities$locality) # see levels 
localities$transect <- ifelse(grepl("CC",localities$locality),'colonso','pipeline') # add levels for transects
localities$population <-  as.factor(gsub( "_.*$", "", localities$ddocent_ID)) # add pops

# set palette
scale_fill_wes <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(wes_palette("Darjeeling1", 13, type = "continuous"), levels(localities$population)), 
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
  geom_point(data = df.pipeline,
             aes(x = long, y = lat,size=elevation, fill = population),
             pch=21)+
  scale_fill_wes(name="locality") +
  ggtitle("pipeline") +
  xlab("longitude") +
  ylab("latitude") +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank()) +
  guides(fill = guide_legend(order=1),
         size = guide_legend(order=2))

p2 <- ggmap(colonso, alpha=0.2) +
  theme_bw() +
  geom_point(data = df.colonso,
             aes(x = long, y = lat, fill = population, size=elevation),
             pch=21)+
  scale_fill_wes(name="locality") +
  ggtitle("colonso de chalupas") +
  xlab("longitude") +
  ylab("latitude") +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank()) +
  guides(fill = guide_legend(order=1),
         size = guide_legend(order=2))

grid.arrange(p1,p2,ncol=2)
p1 + p2






