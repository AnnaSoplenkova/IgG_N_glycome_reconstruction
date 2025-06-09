#код для подготовки к импутации

"/mnt/polyomica/projects/GPD/Real_IgG_collections/Vis/02_unification_out"
library(data.table)

for (i in c(2,3,4,5,7:12,15,20)) {
  
  dt = fread(
    paste0("/mnt/polyomica/projects/GPD/Real_IgG_collections/Korsula/02_unification_out/igp",i,"/igp",i,"_done.csv") )
  
  names(dt)[c(3,5:8,11:15)] = c("rsid", "chr", "bp", "effect_allele", 
                                "other_allele", "beta", "se", "p", "n", "z")
  dt1 = dt[,c(3,5:8,11:15)]
  dt1
  
  fwrite(dt1, paste0("/mnt/polyomica/projects/GPD/imputation_output/imput_data/Korsula/igp",i,".csv"))
  
}



dt = fread("/mnt/polyomica/projects/GPD/Real_IgG_collections/ORCADES/02_unification_out/igp3/igp3_done.csv")


names(dt)[c(3,5:8,11:15)] = c("rsid", "chr", "bp", "effect_allele", 
                              "other_allele", "beta", "se", "p", "n", "z")
dt1 = dt[,c(3,5:8,11:15)]
dt1

fwrite(dt1, "/mnt/polyomica/projects/GPD/imputation_output/imput_data/ORCADES/igp3.csv")


dt = fread("/mnt/polyomica/projects/GPD/imputation_output/replication/ORCADES/igp2_chr1.csv")
table(dt$r2.pred > 0.7)
table(dt$source)
max(dt$N.imp)


setwd("/mnt/polyomica/projects/GPD/imputation_output/replication/ORCADES/igp2")
names = list.files("/mnt/polyomica/projects/GPD/imputation_output/replication/ORCADES/igp2")

list_dt = lapply(names, fread)

dt = do.call(rbind, list_dt)

#9507306 - всего
#216471 - вообще случайные (r2 = 0)
#4511084 - больше чем 0.8
#2001020 - больше чем 0.9
#1994149 - исходно было снипов
table(dt$r2.pred > 0.9)
table(dt$source)
