# look at data
gunzip -c ~/Dropbox/scarab_migration/raw_data/UTN-DUNG_S167_L007_R1_001.fastq.gz | head -n 12

# make ipyrad params file
ipyrad -n scarab

#demultiplex
ipyrad -p params-scarab.txt -s 1
