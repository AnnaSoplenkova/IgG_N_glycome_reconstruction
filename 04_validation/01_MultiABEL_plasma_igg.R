# скрипт multiABEL
library(MultiABEL)

load('/mnt/polyomica/projects/MV_TPNG/MultiABEL/indep.snps.RData')

names_plasma = paste0("/mnt/polyomica/projects/plasma_v2/11_EUR_META/01_GWAS_MAP/01_upload_data/02_unification_out/pgp", 1:36, "/pgp", 1:36, "_done.csv")

good_igg = c(2:5,7:12,15,20)
names_igg = rep(NA, 12)
for (i in 1:12) {
  names_igg[i] = paste0("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/unification_out/igp",good_igg[i],"/igp",good_igg[i],"_done.csv")
}


names_all = c(names_plasma, names_igg)

#load list of independent snps
load('/mnt/polyomica/projects/MV_TPNG/MultiABEL/indep.snps.RData')

stats <- load.summary(files = names_all, indep.snps = indep.snps, columnNames = c("rs_id", "ea", "ra", "eaf", "beta", "se", "n"))


save(stats, file = "/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MultiABEL/stats_load_summary.rda")



library(corrplot)
corrplot(stats$cor.pheno)
names(as.matrix(stats$cor.pheno))
corr = stats$cor.pheno
names(corr[,1])

fwrite(as.matrix(stats$cor.pheno), "MultiABEL/cor.pheno.csv")
corDT = fread("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MultiABEL/cor.pheno.csv")
names(corDT) = c(paste0("plasma_", 1:36), paste0("igg_", good_igg))
corMT = as.matrix(corDT)
rownames(corMT) = names(corDT)
corrplot(abs(as.matrix(corMT)) )
getwd()
