#код перетаскивает в GPD коллекции реального ИгГ (не все пики, потому что код повторный). Остальные кода, перетаскивающие коллекции, утеряны.

#пики в анализе:
#Наши	Луция
#IGP1	IGP1
#IGP2	IGP2
#IGP3	IGP3
#IGP4	IGP4+IGP5
#IGP5	IGP6
#IGP6	IGP7
#IGP7	IGP8
#IGP8	IGP9
#IGP9	IGP10
#IGP10	IGP11
#IGP11	IGP12+IGP13
#IGP12	IGP14
#IGP13	IGP15
#IGP14	IGP16
#IGP15	IGP17
#IGP20	IGP22



peaks_lucia = c(1,7,8,15,16,17)
peaks_gpd =   c(1,6,7,13,14,15)


#Korsula
for (i in 1:6) {
  df = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Korcula/Korcula.IGP",peaks_lucia[i],".10122014.txt"))
  fwrite(df,paste0("/mnt/polyomica/projects/GPD/Real_IgG_collections/Korsula/IGP",peaks_gpd[i],".csv"))
}
#Vis
for (i in 1:6) {
  df = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Vis/Vis.IGP",peaks_lucia[i],".16122014.txt"))
  fwrite(df,paste0("/mnt/polyomica/projects/GPD/Real_IgG_collections/Vis/IGP",peaks_gpd[i],".csv"))
}
#TwinsUK
for (i in 1:6) {
  df = fread(paste0("/home/common/projects/IgGGlycans_GWAS/TWINSUK/TWINSUK.IGP",peaks_lucia[i],".20141208.txt"))
  fwrite(df,paste0("/mnt/polyomica/projects/GPD/Real_IgG_collections/TwinsUK/IGP",peaks_gpd[i],".csv"))
}
#ORCADES
for (i in 1:6) {
  df = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Orkney2000/Orkney2000.IGP",peaks_lucia[i],".17122014.txt"))
  fwrite(df,paste0("/mnt/polyomica/projects/GPD/Real_IgG_collections/ORCADES/IGP",peaks_gpd[i],".csv"))
}
#mv_igg_repl_3k
for (i in 1:6) {
  df = fread(paste0("/mnt/polyomica/projects/MV_igg_replication/3k_meta_23traits/meta_results/3K_meta_IGP",peaks_lucia[i],".csv"))
  fwrite(df,paste0("/mnt/polyomica/projects/GPD/Real_IgG_collections/mv_igg_repl_3k/IGP",peaks_gpd[i],".csv"))
}


