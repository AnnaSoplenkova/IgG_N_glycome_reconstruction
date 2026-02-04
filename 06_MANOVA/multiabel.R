library(MultiABEL)
load('/mnt/polyomica/projects/MV_TPNG/MultiABEL/indep.snps.RData')
data <- load.summary(
  file.path("/mnt/polyomica/projects/GPD/result_ma/02_IgG_new_models_with_non_EUR/unification_out",
           paste0("igp", c(1:15, 20)),
           paste0("igp", c(1:15, 20), "_done.csv")), columnNames = c("rs_id", "ea", "ra", "eaf", "beta", "se", "n"),
  indep.snps = indep.snps
)
result <- MultiSummary(data, type = "outbred")
fwrite(result$scan, file = "/mnt/polyomica/projects/GPD/01_manova/result_scan_16", sep = '\t')
fwrite(result$coef, file = "/mnt/polyomica/projects/GPD/01_manova/result_coef_16", sep = '\t')
#add chr/pos/ea/ra
info_df = fread("/mnt/polyomica/projects/GPD/result_ma/02_IgG_new_models_with_non_EUR/unification_out/igp1/igp1_done.csv")
names(info_df)

info_cols = c("rs_id", "chr", "bp", "ea", "ra")
SNP_info = as.data.frame(info_df)[, info_cols]
scan <- merge(result$scan, SNP_info, by.x = 'marker', by.y = 'rs_id', all.x = T)

scan <- scan[, -c(7:54)]
scan <- scan[, c(1, 7, 8, 9, 10, 2, 3, 4, 5, 6)]
fwrite(scan, file = "/mnt/polyomica/projects/GPD/01_manova/result_scan_16_coord")
