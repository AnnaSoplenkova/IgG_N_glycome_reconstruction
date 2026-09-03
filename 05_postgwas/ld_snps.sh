# creates a list of snps that are in high linkage disequilibrium with replicated snps
plink --bfile /home/common/DataStorage/1kg_503_EUR_82M_SNP_ref/1kg_503_eur_82M_ref --chr 1,4,7,10,19 --out /mnt/polyomica/projects/GPD/PostGWAS_v2/SNPs_in_LD_with_5_SNPs.txt --show-tags /mnt/polyomica/projects/GPD/PostGWAS_v2/5_SNPs.txt --tag-r2 0.8 --threads 5
