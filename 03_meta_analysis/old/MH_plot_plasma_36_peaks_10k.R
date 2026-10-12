# MH plot plasma 10k

library(data.table)
library(qqman)

setwd("/mnt/polyomica/projects/GPD/peak_prediction")
coefs = read.xlsx("OLS/coeffs.xlsx")

IgG_coef = names(coefs)[startsWith(names(coefs), "Estimate")]

names = paste0("/mnt/polyomica/projects/plasma_v2/11_EUR_META/01_GWAS_MAP/01_upload_data/02_unification_out/pgp", 1:36, "/pgp", 1:36, "_done.csv")

dataFiles <- lapply(names, fread)
print("downloading files now finished")

dataFiles[[1]]

snp_info = subset(dataFiles[[1]], select = c(rs_id,chr,bp))

pgp = function(i) {
  dataFiles[[i]]$gwas_id = rep(paste0("pgp", i), nrow(dataFiles[[i]]) )
  return(dataFiles[[i]])
}

dataFiles = lapply(1:36, pgp)

p.select <- function(mydata){
  pvals <- subset(mydata, select = c( rs_id, p))
  gwas.id <- unique(mydata$gwas_id)
  colnames(pvals)[2] <- paste("p", gwas.id, sep ="_")
  return(pvals)
}

list.pvals  <- lapply(dataFiles, p.select)
tmp <- Reduce(left_join, list.pvals)

tmp_1 <- left_join(snp_info, tmp)

pvals = names(tmp_1)[startsWith(names(tmp_1), "p_")]
tmp_1$p_min <- do.call(pmin, tmp_1[, ..pvals])

fwrite(tmp_1, "/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/Plasma_36peaks_pvals_10k.csv")

plasma_gwas = data.frame(BP = tmp_1$bp, CHR = tmp_1$chr, P = tmp_1$p_min, SNP = tmp_1$rs_id)
plasma_gwas = plasma_gwas[!is.na(plasma_gwas$P),]
plasma_gwas_lite = plasma_gwas[plasma_gwas$P < 0.001,]
nrow(plasma_gwas_lite)


pdf("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MH_plasma_10k.pdf", width = 16, height = 9)
manhattan(plasma_gwas_lite, chr="CHR", bp="BP", snp="SNP", p="P", col = c("#606060FF", "#B1009B"),
          suggestiveline = F, genomewideline = -log10(5e-08/28))
dev.off() 

svg("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MH_plasma_10k.svg", width = 16, height = 9)
manhattan(plasma_gwas_lite, chr="CHR", bp="BP", snp="SNP", p="P", col = c("#606060FF", "#B1009B"),
          suggestiveline = F, genomewideline = -log10(5e-08/28))
dev.off() 
