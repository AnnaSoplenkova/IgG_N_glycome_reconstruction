#!/usr/bin/python3
import re



with open('SNPs_in_LD_with_3_repl_SNPs.txt.tags', 'r') as required_rsid:
    list_rsid = [line.strip() for line in required_rsid]
  
    
with open('/home/common/DataStorage/1kg_503_EUR_82M_SNP_ref/1kg_503_eur_82M_ref.bim', mode='r') as in_file:
    for line in in_file:
        list_line = re.findall(r'\w+', line)
        if list_line[1] in list_rsid:
            with open('required_rsid_3.bim', mode='a') as out_file:
                new_line = list_line[1] + ',' + list_line[0] + ',' + list_line[3] + ',' + list_line[4] + ',' + list_line[5]
                out_file.write(f'{new_line}\n')
