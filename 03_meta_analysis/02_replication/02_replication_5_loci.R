#replication

setwd("/mnt/polyomica/projects/GPD/meta_real_IgG_without_TwinsUK")

#rs67411091	1	161593962	IGP_7
#rs747559	4	103414175	IGP_3	
#rs6796	7	6502367	IGP_12
#rs7918084	10	94429467	IGP_14
#rs6510615	19	1649363	IGP_8



dt7 = fread("IGP7/igp7.csv")
dt3 = fread("IGP3/igp3.csv")
dt12 = fread("IGP12/igp12.csv")
dt14 = fread("IGP14/igp14.csv")
dt8 = fread("IGP8/igp8.csv")

dt9 = fread("IGP9/igp9.csv")

rep1 = dt7[dt7$rs_id %in% "rs67411091",]
rep_new = dt3[dt3$rs_id %in% "rs747559",]
rep2 = dt12[dt12$rs_id %in% "rs6796",]
rep3 = dt14[dt14$rs_id %in% "rs7918084",]
rep4 = dt8[dt8$rs_id %in% "rs6510615",]



res = rbind(rbind(rbind(rbind(rep1, rep_new), rep2),rep3), rep4)
write.xlsx(res, "/mnt/polyomica/projects/GPD/result_ma/02_IgG_new_models_with_non_EUR/replication_with_nonEUR.xlsx")
