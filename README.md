# Analysis of scarab beetle genotypes for Linck & Sheldon *in prep*
    
### Scripts  
  
(in `scripts/` subdirectory)  
  
`abc.R:` run approximate Bayesian computation for demographic inference  
  
`bedassle.R:` run BEDASSLE to estimate relative impact of environment on genetic differentiation   
  
`genotype_analysis.R:` handle `.vcf` files and prepare for plotting; run PCA etc.     
  
`quality_control.R:` visualize sequencing results (readcount by spp.)   
  
`scarab_analysis_notebook.Rmd:` Rmarkdown digital notebook file for all assembly and analyses  
  
`scarab_analysis_notebook.md:` knit digital notebook file, with images  
  
`snp_filtering.sh`: SNP filtering parameters  
  
`stats.R:` calculate summary statistics and mixed models  
  
(in `data/` subdirectory)  
  
`fst_dist_species.csv:`  pairwise FST by geographic distance data  
  
`genotypes.df.csv:` genotype PCA and sample data  
  
`growth.csv:` posterior distribution of parameters from ABC "growth" model  
  
`master.theta.df:` etsimtates of theta across elevation  
  
`melted.bayes.csv:` Bayes factor data from ABC  
  
`melted.coancestry.df.csv:` fineRADstructure results  
  
`scarab_spp_master.csv:` full sample and locality data  
  
(in `dDocent/` subdirectory)  
  
`*_config.file:` species-specific dDocent assembly parameters  
  
`*_rename.txt:` name manipulations to meet program reqs  
  
(in `bedassle/` subdirectory) 
  
`*.Robj:` BEDASSLE output for each species  

  



