#!/bin/bash

cd d_satanas
echo "now filtering D. satanas"
mv Final.recode.vcf d_satanas.raw.vcf
vcftools --vcf d_satanas.raw.vcf --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out d.satanas.raw.g5mac3
vcftools --vcf d.satanas.raw.g5mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out d.satanas.raw.g5mac3dp3
vcftools --vcf d.satanas.raw.g5mac3dp3.recode.vcf --missing-indv
mawk '$5 > 0.3' out.imiss | cut -f1 > lowDP.indv
vcftools --vcf d.satanas.raw.g5mac3dp3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out d.satanas.raw.g5mac3dplm # drop one individual
vcftools --vcf d.satanas.raw.g5mac3dplm.recode.vcf --max-missing 0.75 --maf 0.05 --recode --recode-INFO-all --out d.satanas.DP3g95maf05 --min-meanDP 5
bash dDocent_filters d.satanas.DP3g95maf05.recode.vcf d_satanas_filtered

cd ../d_speciocissimum
echo "now filtering D. speciocissimum"
mv Final.recode.vcf d_spec.raw.vcf
vcftools --vcf d_spec.raw.vcf --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out d.spec.raw.g5mac3
vcftools --vcf d.spec.raw.g5mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out d.spec.raw.g5mac3dp3
vcftools --vcf d.spec.raw.g5mac3dp3.recode.vcf --missing-indv
mawk '$5 > 0.3' out.imiss | cut -f1 > lowDP.indv
vcftools --vcf d.spec.raw.g5mac3dp3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out d.spec.raw.g5mac3dplm # drop two individuals
vcftools --vcf d.spec.raw.g5mac3dplm.recode.vcf --max-missing 0.75 --maf 0.05 --recode --recode-INFO-all --out d.spec.DP3g95maf05 --min-meanDP 5
bash dDocent_filters d.spec.DP3g95maf05.recode.vcf d_spec_filtered

cd ../d_tesselatum
echo "now filtering D. tesselatum"
mv Final.recode.vcf d_tess.raw.vcf
vcftools --vcf d_tess.raw.vcf --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out d.tess.raw.g5mac3
vcftools --vcf d.tess.raw.g5mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out d.tess.raw.g5mac3dp3
vcftools --vcf d.tess.raw.g5mac3dp3.recode.vcf --missing-indv
mawk '$5 > 0.3' out.imiss | cut -f1 > lowDP.indv
vcftools --vcf d.tess.raw.g5mac3dp3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out d.tess.raw.g5mac3dplm # drop one individual
vcftools --vcf d.tess.raw.g5mac3dplm.recode.vcf --max-missing 0.75 --maf 0.05 --recode --recode-INFO-all --out d.tess.DP3g95maf05 --min-meanDP 5
bash dDocent_filters d.tess.DP3g95maf05.recode.vcf d_tess_filtered

cd ../e_affin
echo "now filtering E. affin"
mv Final.recode.vcf e_affin.raw.vcf
vcftools --vcf e_affin.raw.vcf --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out e.affin.raw.g5mac3
vcftools --vcf e.affin.raw.g5mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out e.affin.raw.g5mac3dp3
vcftools --vcf e.affin.raw.g5mac3dp3.recode.vcf --missing-indv
mawk '$5 > 0.3' out.imiss | cut -f1 > lowDP.indv
vcftools --vcf e.affin.raw.g5mac3dp3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out e.affin.raw.g5mac3dplm # droped nine individual
vcftools --vcf e.affin.raw.g5mac3dplm.recode.vcf --max-missing 0.75 --maf 0.05 --recode --recode-INFO-all --out e.affin.DP3g95maf05 --min-meanDP 5
bash dDocent_filters e.affin.DP3g95maf05.recode.vcf e_affin_filtered

cd ../d_podalirius
echo "now filtering D. podalirius"
mv Final.recode.vcf d_pod.raw.vcf
vcftools --vcf d_pod.raw.vcf --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out d.pod.raw.g5mac3
vcftools --vcf d.pod.raw.g5mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out d.pod.raw.g5mac3dp3
vcftools --vcf d.pod.raw.g5mac3dp3.recode.vcf --missing-indv
mawk '$5 > 0.3' out.imiss | cut -f1 > lowDP.indv
vcftools --vcf d.pod.raw.g5mac3dp3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out d.pod.raw.g5mac3dplm # drop no individuals
vcftools --vcf d.pod.raw.g5mac3dplm.recode.vcf --max-missing 0.75 --maf 0.05 --recode --recode-INFO-all --out d.pod.DP3g95maf05 --min-meanDP 5
bash dDocent_filters d.pod.DP3g95maf05.recode.vcf d_pod_filtered
