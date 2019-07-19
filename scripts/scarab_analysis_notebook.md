Scarab ddRADseq data analysis
================

## Demultiplexing and basic quality control

We’ve just recieved A LOT (9.39 GB) of data: 100 bp single-end read
sequencing data from 239 samples of 11 species of dung beetle
(Scarabaeinae), multiplexed on a single Illumina lane, to be precise.
Before we dig in to the all the interesting questions we want to ask
with it, we should do a basic check that the sequencing reaction was
successful, and that we don’t need to start from scratch (lolsob).
First, we need to demultiplex the original file,
`raw_data/UTN-DUNG_S167_L007_R1_001.fastq.gz1`, using the list of
barcodes we appended to each sample during library preparation
(`ipyrad/barcodes`). We can look at it and see it matches our
expectations for .fastq data, which I won’t display here. (Enter all
shell commands in the appropriate directory via terminal.)

``` bash
gunzip -c ~/Dropbox/scarab_migration/raw_data/UTN-DUNG_S167_L007_R1_001.fastq.gz | head -n 12
```

We then switch to using program `ipyrad`, [which has extensive
documentation you can view here](https://ipyrad.readthedocs.io/). We
first generate a params file, which we will use for the remainder of
data assembly.

``` bash
ipyrad -n scarab
```

We modify the file `ipyrad/params-scarab.txt` to match our data and fill
in the appropriate paths to our data and barcode file. We then run the
first two steps of the ipyrad pipeline, which demultiplex the single
large gunzipped fastq file and quality trims and filters reads.

``` bash
ipyrad -p params-scarab.txt -s 12
```

When this completes, you’ll see the relevant directory fill with fastq
files, and a text file labeled `s1_demultiplex_stats.txt` generated. We
want to extract data on the total number of reads per library, and group
libraries by species. First, we’ll read in the (unformatted) stats file
and select only the information we want:

``` r
setwd("~/Dropbox/scarab_migration/ipyrad/")
file <- "s1_demultiplex_stats.txt"
nsamples <- 239
raw <- readLines("s1_demultiplex_stats.txt",warn=FALSE)
a <- grep("sample_name",raw) # where the data begins
a <- a[1] 
reads <- read.table(file,skip=(a),nrow=(nsamples))
colnames(reads) <- c("sample_ID", "read_count")
```

We know have a simple two column data frame to play around with. We will
do some tedious substitution of our field-based morphospecies assignment
shorthand with the real species IDs:

``` r
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
```

Kind of a mess, but that’s biology\! Now we’ll plot it.

![](scarab_analysis_notebook_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

So far so good: pretty even sequencing effort across the five species
we’re interested in with the most data.

## Data assembly

Now that we’re comfortable we sequenced *something*, we’re going to try
and assemble the remainder of our data. We’re only going to analyze the
five species with reasonable sampling: *Deltochilum speciocissimum*,
*Deltochilum tesselatum*, *Dichotomius podalirius*, *Dichotomius
satanas*, and *Eurysternus affin*. Because I’m not confident in how
`ipyrad` is assembling these data right now, I’m going to keep this
brief as it may change, but essentially we want to make branches of our
original pipeline for each species and then run these independently.
We’ll need to move our files to unique folders first as well

``` bash

# make directories
mkdir d_speciocissimum
mkdir d_tesselatum
mkdir d_podalirirus
mkdir d_satanas
mkdir e_affin

# move files
mv DE1* d_speciocissimum/
mv DE9* d_speciocissimum/ # two morphospecies lumped into one species
mv DE2* d_tesselatum/
mv DI7* d_podalirius/ 
mv DI1* d_satanas/
mv DI5* d_satanas/
mv EU* e_affin/ 

# make ipyrad branches
ipyrad -p params-scarab.txt -b d_speciocissimum
ipyrad -p params-scarab.txt -b d_tesselatum
ipyrad -p params-scarab.txt -b d_podalirirus
ipyrad -p params-scarab.txt -b d_satanas
ipyrad -p params-scarab.txt -b e_affin

# run independently; be sure to edit params files for paths
ipyrad -p params-d_speciocissimum.txt -s 34567
ipyrad -p params-d_tesselatum.txt -s 34567
ipyrad -p params-d_podalirirus.txt -s 34567
ipyrad -p params-d_satanas.txt -s 34567
ipyrad -p params-e_affin.txt -s 34567
```

The output we’re interested in are .vcf files; we’ll be working with
these for the rest of the analysis.

## Additional filtering and testing for genetic differentiation

Our driving question is whether dispersal (and gene flow) is reduced
across elevational gradients relative to within an elevational band. We
are going to try and answer this question in a number of ways. An
important first step is running principle component analysis to ask
whether there is population genetic structure, and whether it relates to
elevation. For now, let’s only look at *D. satanas*, as we have the best
sampling. We’ll read in locality data, our .vcf file, and then use a
custom chunk of code to replace `ipyrad's` weird .vcf syntax with
something `R` can use.

``` r
# load libraries
library(vcfR)
library(adegenet)
library(ggplot2)
library(wesanderson)
library(tidyverse)
library(reshape2)
library(poppr)

# set wd and read in locality data
setwd("/Users/ethanlinck/Dropbox/scarab_migration/")
localities <- read.csv("data/scarab_spp_master.csv")

# read in .vcf and sample / pop names; convert to genlight
satanas.vcf <- read.vcfR("ipyrad/d_satanas_outfiles/d_satanas.vcf", verbose=FALSE)

# fix ipyrad vcf missing data issue
satanas.vcf@gt <- 
  satanas.vcf@gt %>% 
  as_tibble() %>% 
  mutate_all(.funs = function(x) replace(x, which(x == "./.:0:0,0,0,0"| x == "NA"), NA)) %>%
  as.matrix() 
```

Now, let’s check the .vcf file:

``` r
satanas.vcf
```

    ## ***** Object of Class vcfR *****
    ## 100 samples
    ## 8364 CHROMs
    ## 61,912 variants
    ## Object size: 52.8 Mb
    ## 49.18 percent missing data
    ## *****        *****         *****

49.18% missing data is typical for ddRADseq data, so that’s reassuring.
Let’s drop loci missing data for lots of individuals (\>40%).

``` r
# drop rows with >40% MD
dp <- vcfR::extract.gt(satanas.vcf,  element = "DP", as.numeric = TRUE)
satanas.miss <- apply(dp, MARGIN = 1, function(x){sum(is.na(x))})
satanas.miss <- satanas.miss / ncol(dp)
satanas.vcf <- satanas.vcf[satanas.miss < 0.4, ]
satanas.vcf
```

    ## ***** Object of Class vcfR *****
    ## 100 samples
    ## 1500 CHROMs
    ## 10,092 variants
    ## Object size: 9 Mb
    ## 29.71 percent missing data
    ## *****        *****         *****

10,000 SNPs on 1500 chromosomes is pretty reasonable. Let’s check the
average depth of coverage for a handful of individuals:

``` r
# randomly sample ten indiduals to check depth, arrange data frame
dp2 <- dp[,sample(ncol(dp),size=10,replace=FALSE)]
dpf <- melt(dp2, varnames = c("Index", "Sample"),
            value.name = "Depth", na.rm = TRUE)
dpf <- dpf[ dpf$Depth > 0, ]
```

Now, let’s plot it:

![](scarab_analysis_notebook_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Variable, but we can see our basic filter worked (no calls \<3 reads).
Let’s consider it okay for now and run a PCA. We’ll have to change the
file format to something `adegenet` can read, then do some manipulation
and assign new individual and population labels from our directory:

``` r
satanas.dna <- vcfR2DNAbin(satanas.vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
```

    ## After extracting indels, 10092 variants remain.

``` r
satanas.gen <- DNAbin2genind(satanas.dna)
satanas.samples <- as.character(read.table("/Users/ethanlinck/Dropbox/scarab_migration/ipyrad/d_satanas_outfiles/samples_d_satanas.txt")[[1]]) 
satanas.pops <- as.factor(as.character(read.table("/Users/ethanlinck/Dropbox/scarab_migration/ipyrad/d_satanas_outfiles/populations_d_satanas.txt")[[1]])) 
rownames(satanas.gen@tab) <- satanas.samples
satanas.gen@pop <- satanas.pops
satanas.scaled <- scaleGen(satanas.gen,NA.method="mean",scale=F)
satanas.pca <- prcomp(satanas.scaled,center=F,scale=F)
satanas.pc <- data.frame(satanas.pca$x[,1:3])
satanas.pc$sample <- rownames(satanas.pc)
satanas.pc$pop <- satanas.pops
satanas.pc$species <- rep("Dichotomius_satanas",nrow(satanas.pc))
satanas.pc <- merge(satanas.pc, localities, by.x = "sample", by.y = "sample_ID")
```

We’ll now plot the first two PCs against each other, and plot PC1
against elevation:

![](scarab_analysis_notebook_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Hmm, looks a lot like panmixia\! I’ll add the other species next week,
after I’ve checked our assembly.