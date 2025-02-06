

#file.copy("/mnt/polyomica/projects/GPD/result_ma/PGP93_FUC_C.csv", "/mnt/polyomica/projects/GPD/check_lincomb_lambda")

realPGP.name = "/mnt/polyomica/projects/plasma_v2/11_EUR_META/01_GWAS_MAP/01_upload_data/02_unification_out/pgp93/pgp93_done.csv"
lincombPGP.name = "/mnt/polyomica/projects/GPD/check_lincomb_lambda/PGP93_FUC_C.csv"


dt.real = fread(realPGP.name)
dt.linkomb = fread(lincombPGP.name)


library(GWASTools)

df("/mnt/polyomica/projects/GPD/check_lincomb_lambda/PGP93_FUC_C_plot.real.pdf")
qqPlot(dt.real$p)
dev.off()

pdf("/mnt/polyomica/projects/GPD/check_lincomb_lambda/PGP93_FUC_C_plot.lincomb.pdf")
qqPlot(dt.linkomb$Pval)
dev.off()

library(gap)
gc.lambda(dt.real$p)
gc.lambda(dt.linkomb$Pval)
compare_GWASs(list(realPGP.name,lincombPGP.name), "/mnt/polyomica/projects/GPD/result_ma/")

