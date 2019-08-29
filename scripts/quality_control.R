# basic sequencing and assembly quality control 

install.packages("ggplot2")
library(ggplot2)

setwd("~/Dropbox/scarab_migration/ipyrad/")
file <- "s1_demultiplex_stats.txt"
nsamples <- 239
raw <- readLines("s1_demultiplex_stats.txt",warn=FALSE)
a <- grep("sample_name",raw) # where the data begins
a <- a[1] 
reads <- read.table(file,skip=(a),nrow=(nsamples))
colnames(reads) <- c("sample_ID", "read_count")

reads <- read.table("reads.txt")[-1,]
colnames(reads) <- c("sample_ID", "read_count")
reads$sample_ID <- gsub("EU","Eurysternus_affin",reads$sample_ID)
reads$sample_ID <- gsub("DI1","Dichotomious_satanas",reads$sample_ID)
reads$sample_ID <- gsub("DI5","Dichotomious_satanas",reads$sample_ID)
reads$sample_ID <- gsub("DI7","Dichotomious_podalirius",reads$sample_ID)
reads$sample_ID <- gsub("DE2","Deltochilum_tesselatum",reads$sample_ID)
reads$sample_ID <- gsub("DE1","Deltochilum_speciosissimum",reads$sample_ID)
reads$sample_ID <- gsub("DE9","Deltochilum_speciosissimum",reads$sample_ID)
reads$sample_ID <- gsub("239","Dichotomious_protectus",reads$sample_ID)
reads$sample_ID <- gsub("238","Dichotomious_protectus",reads$sample_ID)
reads$sample_ID <- gsub("219","Deltochilum_amazonicum",reads$sample_ID)
reads$sample_ID <- gsub("220","Coprophanaeus_telamon",reads$sample_ID)
reads$sample_ID <- gsub("217","Phanaeus_meleagris",reads$sample_ID)
reads$sample_ID <- gsub("218","Phanaeus_meleagris",reads$sample_ID)
reads$sample_ID <- gsub("223","Oxysternon_silenus",reads$sample_ID)
reads$sample_ID <- gsub("222","Oxysternon_silenus",reads$sample_ID)
reads$sample_ID <- gsub("221","Oxysternon_conspicullatum",reads$sample_ID)
reads$sample_ID <- gsub('-.*', '', reads$sample_ID) # drop numbers
colnames(reads) <- c("Species","Read_count")
reads$Read_count <- as.numeric(as.character(reads$Read_count))

pdf("~/Dropbox/scarab_migration/figures/readcount.pdf",7,7)
ggplot(reads, aes(x=Species, y=Read_count)) +
  theme_classic() +
  geom_boxplot()+
  geom_jitter(height = 0, width = 0.1) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=9)) +
  ylab("Number of reads")
dev.off()

# stats for paper
levels(as.factor(reads$Species))
drop <- c("Phanaeus_meleagris","Deltochilum_amazonicum","Coprophanaeus_telamon",
          "Oxysternon_silenus","Dichotomious_protectus","Oxysternon_conspicullatum")
reads.subset <- reads[!reads$Species %in% drop,]
levels(as.factor(reads.subset$Species))

# average coverage
m1 <- reads.subset[reads.subset$Species=="Deltochilum_speciosissimum",]$Read_count %>% mean() #884498.7
m2 <- reads.subset[reads.subset$Species=="Deltochilum_tesselatum",]$Read_count %>% mean() #890779.1
m3 <- reads.subset[reads.subset$Species=="Dichotomious_podalirius",]$Read_count %>% mean() #905850.6
m4 <- reads.subset[reads.subset$Species=="Dichotomious_satanas",]$Read_count %>% mean() #917563.1
m5 <- reads.subset[reads.subset$Species=="Eurysternus_affin",]$Read_count %>% mean() #824833.5
mdf <- as.vector(c(m1,m2,m3,m4,m5))
mean(mdf) #884705
sd(mdf) #35875.92
mean(reads.subset$Read_count)
sd(reads.subset$Read_count)



