# look at data
gunzip -c ~/Dropbox/scarab_migration/raw_data/UTN-DUNG_S167_L007_R1_001.fastq.gz | head -n 12

# make ipyrad params file
ipyrad -n scarab

# demultiplex; filter
ipyrad -p params-scarab.txt -s 12

# make branch for dichotomius satanas
ipyrad -p params-scarab.txt -b d_satanas

# make branch for dichotomius satanas
ipyrad -p params-scarab.txt -b d_speciossimum

# move samples to new subfolder
mkdir scarab_fastqs/d_satanas
mv scarab_fastqs/DI1* scarab_fastqs/d_satanas
mv scarab_fastqs/DI5* scarab_fastqs/d_satanas

# run assembly for dichotomius satanas
ipyrad -p params-d_satanas.txt -s 34567
