#This script makes correlation matrix of 143 GWAMA with help of MultiABEL package
#copy files to another directory

flist <- c("H5N2","H6N2","H5N3","H3N4F1","H7N2","H3N3F1E1","H4N3E1","H4N4F1","H5N4","H3N5F1","H4N5","H8N2","H5N3E1","H5N4F1","H4N5F1","H5N5","H9N2","H6N3E1","H4N4F1E1","H5N4E1","H5N5F1","H5N4F1L1","H5N4F1E1","H4N5F1E1","H5N5E1","H5N4L2","H5N4E1L1","H5N4E2","H5N5F1E1","H6N5E1","H5N4F1L2","H5N4F1E1L1","H5N4F1E2","H6N5F1E1","H5N5E2","H5N5F1E1L1","H6N5E1L1","H5N5F1E2","H6N5E2","H6N5F1E1L1","H6N5F1E2","H6N5E1L2","H6N5E2L1","H6N5E3","H6N5F1E1L2","H6N5F1E2L1","H6N5F2E2L1","H7N6E1L2","H7N6E2L1","H7N6F1E1L2","H7N6F1E2L1","H7N6E1L3","H7N6E2L2","H7N6E3L1","H7N6F1E1L3","H7N6F1E2L2","H7N6F2E2L2","A1F","A2B","A2E","A2E0F","A2EF","A2F","A2F0B","A2F0E","A2F0G","A2F0GE","A2F0GL","A2F0GS","A2F0L","A2F0S","A2F0S0B","A2F0S0G","A2F0SB","A2FB","A2FE","A2FG","A2FGE","A2FGL","A2FGS","A2FL","A2FS","A2FS0B","A2FS0G","A2FSB","A2FSG","A2G","A2GE","A2GL","A2GS","A2L","A2L0F","A2LF","A2S","A2S0B","A2S0F","A2S0G","A2SB","A2SF","A2SG","A3E","A3F","A3F0E","A3F0L","A3F0S","A3Fa","A3FE","A3FL","A3FS","A3L","A3L0F","A3LF","A3S","A4E","A4F","A4F0E","A4F0GE","A4F0GL","A4F0L","A4F0S","A4Fa","A4FE","A4FGE","A4FGL","A4FGS","A4FL","A4FS","A4L","A4S","CA1","CA2","CA3","CA4","CB","CF","CFa","CG","CS","MHy","MM","TA2FS0","THy","TM")

for (i in flist) {
  path <- paste0("/mnt/polyomica/projects/ms41loci/03_GWAMA/", i, "/", i)
  file.copy(path, "/mnt/polyomica/projects/ms41loci/05_MultiABEL/")
}

#run load.summary() and save as .RData

library(MultiABEL)
setwd("/mnt/polyomica/projects/ms41loci/05_MultiABEL/")
load('/mnt/polyomica/projects/ms41loci/indep.snps.RData')
sst <- load.summary(files = list.files(), cor.pheno = NULL, indep.snps = indep.snps, est.var = FALSE, columnNames = c("rs_id", "ea", "ra", "eaf_ma", "beta_ma", "se_ma", "n_ma"), fixedN = NULL)
save(file='/mnt/polyomica/projects/ms41loci/143_tpng_sst.Rdata', list=c("sst"))
