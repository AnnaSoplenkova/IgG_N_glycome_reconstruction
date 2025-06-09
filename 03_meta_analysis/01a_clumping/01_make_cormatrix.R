#This script makes correlation matrix of 17 synthetic IGG peaks with help of MultiABEL package
#copy files to another directory

peaks = c(1:15,20)



for (i in peaks) {
  path <- paste0("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/unification_out/igp",i,"/igp",i,"_done.csv")
  file.copy(path, "/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MultiABEL/files/")
}

fread(path)

#run load.summary() and save as .RData

library(MultiABEL)
setwd("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MultiABEL/files/")
load('/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MultiABEL/indep.snps.RData')
sst <- load.summary(files = list.files(), cor.pheno = NULL, indep.snps = indep.snps, est.var = FALSE, columnNames = c("rs_id", "ea", "ra", "eaf", "beta", "se", "n"), fixedN = NULL)
save(file='/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MultiABEL/16_synth_IgG_sst.Rdata', list=c("sst"))
