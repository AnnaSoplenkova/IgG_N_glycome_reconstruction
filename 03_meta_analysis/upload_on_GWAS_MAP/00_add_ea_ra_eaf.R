# скрипт чтобы добавить eaf, ea, ra

setwd("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS")

names = paste0("10k_IGP_",1:21,".csv")

list = list()
for (i in 1:21) {
  list[[i]] = fread(names[i])
}

extra = fread("/mnt/polyomica/projects/plasma_v2/11_EUR_META/01_GWAS_MAP/PGP1/pgp1.csv")
extra = extra[,1:6]


for (i in 1:21) {
  list[[i]] = left_join(list[[i]], extra)
}

for (i in 1:21) {
  fwrite(list[[i]], paste0("10k_IGP_",i,".csv"))
}
