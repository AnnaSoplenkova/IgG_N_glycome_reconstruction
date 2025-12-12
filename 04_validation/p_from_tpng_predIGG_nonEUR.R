setwd("/mnt/polyomica/projects/GPD/result_ma/02_IgG_new_models_with_non_EUR")


IgG_gwas_pred = data.frame(BP = predicted_IgG$bp, CHR = predicted_IgG$chr, P = predicted_IgG$min_p, SNP = predicted_IgG$rs_id)
IgG_gwas_pred_lite = IgG_gwas_pred[IgG_gwas_pred$P < 0.00005,]
# list of snps for 16 IGPs with small sample size n=170 CHR1
rs_vector = c("rs79574235", "rs144090205", "rs77999053")
IgG_gwas_pred_lite = IgG_gwas_pred_lite[!(IgG_gwas_pred_lite$SNP %in% rs_vector),]

# list of snps for 16 IGPs with small sample size n=350 CHR9
rs_vector = c("rs79574235", "rs144090205", "rs77999053")
IgG_gwas_pred_lite = IgG_gwas_pred_lite[!(IgG_gwas_pred_lite$SNP %in% rs_vector),]

# list of snps for 16 IGPs with small sample size n=170 CHR1
rs_vector = c("rs116467783","rs9931970", "rs9932247")
IgG_gwas_pred_lite = IgG_gwas_pred_lite[!(IgG_gwas_pred_lite$SNP %in% rs_vector),]


snp = IgG_gwas_pred_lite$SNP
#IgG_gwas_pred_lite = IgG_gwas_pred
nrow(IgG_gwas_pred_lite)
#head(IgG_gwas_pred_lite)



real_IgG = fread("/mnt/polyomica/projects/GPD/result_ma/02_IgG_new_models_with_non_EUR/min_p_values_across_36_PGPs.csv")
#names(real_IgG)

#head(IgG_gwas_real)
IgG_gwas_real = data.frame(BP = real_IgG$bp, CHR = real_IgG$chr, P = real_IgG$min_p, SNP = real_IgG$rs_id)
IgG_gwas_real = IgG_gwas_real[!is.na(IgG_gwas_real$BP),]
IgG_gwas_real = IgG_gwas_real[!is.na(IgG_gwas_real$P),]
IgG_gwas_real_lite = IgG_gwas_real[IgG_gwas_real$P < 0.00005,]


# list of snps for 16 IGPs with small sample size n=170
rs_vector = c("rs79574235", "rs144090205", "rs77999053")
IgG_gwas_real_lite = IgG_gwas_real_lite[!(IgG_gwas_real_lite$SNP %in% rs_vector),]
# list of snps for 16 IGPs with small sample size n=350 CHR9
rs_vector = c("rs79574235", "rs144090205", "rs77999053")
IgG_gwas_real_lite = IgG_gwas_real_lite[!(IgG_gwas_real_lite$SNP %in% rs_vector),]


#IgG_gwas_real_lite = IgG_gwas_real
snp = c(snp, IgG_gwas_real_lite$SNP)
nrow(IgG_gwas_real_lite)


IgG_gwas_real_lite = IgG_gwas_real[IgG_gwas_real$SNP %in% snp,]
IgG_gwas_pred_lite = IgG_gwas_pred[IgG_gwas_pred$SNP %in% snp,]


library(openxlsx)
loci = read.xlsx("p_across_tpng_igg.xlsx")
names(loci) = c("gene", "chr", "bp")
head(IgG_gwas_real_lite)

get_min_p <- function(chr, bp_min, bp_max, data) {
  subset_data <- data[data$CHR == chr & data$BP >= bp_min & data$BP <= bp_max, ]
  min(subset_data$P)
}

loci$P_igg = NA
loci$P_plasma = NA

p_func = function(i, data) {
  get_min_p(chr = loci$chr[i], bp_min = loci$bp[i]-250000, bp_max = loci$bp[i] + 250000, data = data)
}

loci$P_plasma = sapply(1:nrow(loci), p_func, data =  IgG_gwas_real_lite)
loci$P_igg = sapply(1:nrow(loci), p_func, data =  IgG_gwas_pred_lite)

loci$P_plasma[loci$gene == "HLA"] = get_min_p(chr = 6, bp_min = 28000000, bp_max = 32500000, data = IgG_gwas_real_lite)
loci$P_igg[loci$gene == "HLA"] = get_min_p(chr = 6, bp_min = 28000000, bp_max = 32500000, data = IgG_gwas_pred_lite)

loci$P_plasma[loci$gene == "FUT8"] = get_min_p(chr = 14, bp_min = 65472891, bp_max = 66284991, data = IgG_gwas_real_lite)
loci$P_igg[loci$gene == "FUT8"] = get_min_p(chr = 14, bp_min = 65472891, bp_max = 66284991, data = IgG_gwas_pred_lite)


write.xlsx(loci,"p_across_tpng_igg_full.xlsx")


