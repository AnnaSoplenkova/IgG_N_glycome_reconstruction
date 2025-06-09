#код, который выдает минимальное р-значение для локусов из IgG/TF GWASes в реальном IgG на 23 пиках и в TPNG на 36 пиках

library(dplyr)

setwd("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/IgG_real/clean")
names = list.files(pattern = "*done.csv")
n = length(names)


dt = fread(names[1])
df = data.frame(snp = dt$SNP, p_1 = dt$p)

for(i in 2:n) {
  dt = fread(names[i])
  df_1 = data.frame(snp = dt$SNP, p = dt$p)
  df = left_join(df, df_1, by = c("snp"))
}

head(df)
mmin = function(x) {
  min(x, na.rm = T)
}

p_real_igg = apply(df[,-1], 1, mmin)
res_real_igg = data.frame(snp = df[,1], p_real_igg)
snp_info = data.frame(snp = dt$rs_id, chr = dt$chr, bp = dt$bp)

res = left_join(snp_info,res_real_igg)

head(res)

min(res$p_real_igg[res$chr == 2 & res$bp > 26113120-250000 & res$bp < 26113120+250000], na.rm = T)
min(res$p_real_igg[res$chr == 2 & res$bp > 135014078-250000 & res$bp < 135014078+250000], na.rm = T)

min(res$p_real_igg[res$chr == 3 & res$bp > 133466457-250000 & res$bp < 133466457+250000], na.rm = T)
min(res$p_real_igg[res$chr == 3 & res$bp > 142901537-250000 & res$bp < 142901537+250000], na.rm = T)
min(res$p_real_igg[res$chr == 3 & res$bp > 196203979-250000 & res$bp < 196203979+250000], na.rm = T)

11	114381448
min(res$p_real_igg[res$chr == 11 & res$bp > 114381448-250000 & res$bp < 114381448+250000], na.rm = T)
min(res$p_real_igg[res$chr == 11 & res$bp > 126232337-250000 & res$bp < 126232337+250000], na.rm = T)
min(res$p_real_igg[res$chr == 11 & res$bp > 134264227-250000 & res$bp < 134264227+250000], na.rm = T)

min(res$p_real_igg[res$chr == 12 & res$bp > 121424574-250000 & res$bp < 121424574+250000], na.rm = T)

min(res$p_real_igg[res$chr == 17 & res$bp > 16848750-250000 & res$bp < 16848750+250000], na.rm = T)





###################

igp = c(2,3,4,5,7,8,9,10,11,12,15,20)
i = 2
name = paste0("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/unification_out/igp",igp[i],"/igp",igp[i],"_done.csv")
dt = fread(name)
df = data.frame(snp = dt$rs_id, p_1 = dt$p)



for (i in 2:length(igp)) {
  name = paste0("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/unification_out/igp",igp[i],"/igp",igp[i],"_done.csv")
  dt = fread(name)
  df_1 = data.frame(snp = dt$rs_id, p = dt$p)
  df = left_join(df, df_1, by = c("snp"))
}

head(df)
ncol(df)

p_real_igg = apply(df[,-1], 1, mmin)
res_real_igg = data.frame(snp = df[,1], p_real_igg)
snp_info = data.frame(snp = dt$rs_id, chr = dt$chr, bp = dt$bp)

res = left_join(snp_info,res_real_igg)

min(res$p_real_igg[res$chr == 2 & res$bp > 135014078-250000 & res$bp < 135014078+250000], na.rm = T)

min(res$p_real_igg[res$chr == 3 & res$bp > 133466457-250000 & res$bp < 133466457+250000], na.rm = T)
min(res$p_real_igg[res$chr == 3 & res$bp > 142901537-250000 & res$bp < 142901537+250000], na.rm = T)
min(res$p_real_igg[res$chr == 3 & res$bp > 196203979-250000 & res$bp < 196203979+250000], na.rm = T)

min(res$p_real_igg[res$chr == 5 & res$bp > 95273410-250000 & res$bp < 95273410+250000], na.rm = T)
min(res$p_real_igg[res$chr == 5 & res$bp > 131725329-250000 & res$bp < 131725329+250000], na.rm = T)

min(res$p_real_igg[res$chr == 11 & res$bp > 65551957-250000 & res$bp < 65551957+250000], na.rm = T)
min(res$p_real_igg[res$chr == 11 & res$bp > 114381448-250000 & res$bp < 114381448+250000], na.rm = T)
min(res$p_real_igg[res$chr == 11 & res$bp > 126232337-250000 & res$bp < 126232337+250000], na.rm = T)
min(res$p_real_igg[res$chr == 11 & res$bp > 134264227-250000 & res$bp < 134264227+250000], na.rm = T)

min(res$p_real_igg[res$chr == 12 & res$bp > 121424574-250000 & res$bp < 121424574+250000], na.rm = T)

min(res$p_real_igg[res$chr == 14 & res$bp > 105998544-250000 & res$bp < 105998544+250000], na.rm = T)

min(res$p_real_igg[res$chr == 17 & res$bp > 16848750-250000 & res$bp < 16848750+250000], na.rm = T)
min(res$p_real_igg[res$chr == 17 & res$bp > 45835278-250000 & res$bp < 45835278+250000], na.rm = T)

min(res$p_real_igg[res$chr == 19 & res$bp > 5835677-250000 & res$bp < 5835677+250000], na.rm = T)



names = paste0("/mnt/polyomica/projects/plasma_v2/11_EUR_META/01_GWAS_MAP/01_upload_data/02_unification_out/pgp", 1:36, "/pgp", 1:36, "_done.csv")
dt = fread(names[1])
df = data.frame(snp = dt$rs_id, p_1 = dt$p)
n = length(names)

for(i in 2:n) {
  dt = fread(names[i])
  df_1 = data.frame(snp = dt$rs_id, p = dt$p)
  df = left_join(df, df_1, by = c("snp"))
}

head(df)
ncol(df)

p_real_igg = apply(df[,-1], 1, mmin)
res_real_igg = data.frame(snp = df[,1], p_real_igg)
snp_info = data.frame(snp = dt$rs_id, chr = dt$chr, bp = dt$bp)

res = left_join(snp_info,res_real_igg)
head(res)


min(res$p_real_igg[res$chr == 2 & res$bp > 26113120-250000 & res$bp < 26113120+250000], na.rm = T)

min(res$p_real_igg[res$chr == 3 & res$bp > 133466457-250000 & res$bp < 133466457+250000], na.rm = T)
min(res$p_real_igg[res$chr == 3 & res$bp > 196203979-250000 & res$bp < 196203979+250000], na.rm = T)

min(res$p_real_igg[res$chr == 5 & res$bp > 95273410-250000 & res$bp < 95273410+250000], na.rm = T)
min(res$p_real_igg[res$chr == 5 & res$bp > 131725329-250000 & res$bp < 131725329+250000], na.rm = T)

min(res$p_real_igg[res$chr == 6 & res$bp > 139625074-250000 & res$bp < 139625074+250000], na.rm = T)

min(res$p_real_igg[res$chr == 11 & res$bp > 65551957-250000 & res$bp < 65551957+250000], na.rm = T)
min(res$p_real_igg[res$chr == 11 & res$bp > 114381448-250000 & res$bp < 114381448+250000], na.rm = T)

min(res$p_real_igg[res$chr == 14 & res$bp > 105998544-250000 & res$bp < 105998544+250000], na.rm = T)

min(res$p_real_igg[res$chr == 17 & res$bp > 38069949-250000 & res$bp < 38069949+250000], na.rm = T)
min(res$p_real_igg[res$chr == 17 & res$bp > 44797919-250000 & res$bp < 44797919+250000], na.rm = T)
min(res$p_real_igg[res$chr == 17 & res$bp > 45835278-250000 & res$bp < 45835278+250000], na.rm = T)
min(res$p_real_igg[res$chr == 17 & res$bp > 79192430-250000 & res$bp < 79192430+250000], na.rm = T)

min(res$p_real_igg[res$chr == 19 & res$bp > 19267990-250000 & res$bp < 19267990+250000], na.rm = T)

min(res$p_real_igg[res$chr == 20 & res$bp > 17829280-250000 & res$bp < 17829280+250000], na.rm = T)

min(res$p_real_igg[res$chr == 21 & res$bp > 36565278-250000 & res$bp < 36565278+250000], na.rm = T)

