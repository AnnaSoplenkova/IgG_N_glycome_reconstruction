library(data.table)

list.files("/home/common/projects/IgGGlycans_GWAS/Orkney2000")
fread("/home/common/projects/IgGGlycans_GWAS/Orkney2000/Orkney2000.IGP1.17122014.txt")


list.files("/home/common/projects/IgGGlycans_GWAS/TWINSUK")
fread("/home/common/projects/IgGGlycans_GWAS/TWINSUK/TWINSUK.IGP1.20141208.txt")


list.files("/home/common/projects/IgGGlycans_GWAS/Vis")
fread("/home/common/projects/IgGGlycans_GWAS/Vis/Vis.IGP1.16122014.txt")


list.files("/home/common/projects/IgGGlycans_GWAS/Korcula")
fread("/home/common/projects/IgGGlycans_GWAS/Korcula/Korcula.IGP1.10122014.txt")


list.files("/mnt/polyomica/projects/MV_igg_replication/3k_meta_23traits/meta_results")
fread("/mnt/polyomica/projects/MV_igg_replication/3k_meta_23traits/meta_results/3K_meta_IGP1.csv")


IGPs_real = c("IGP2","IGP3","IGP6","IGP8","IGP9","IGP10","IGP11","IGP14","IGP17","IGP22")
IGPs_our = c("IGP2","IGP3","IGP5","IGP7","IGP8","IGP9","IGP10","IGP12","IGP15","IGP20")


for (i in 1:length(IGPs_real)) {
  x = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Orkney2000/Orkney2000.",IGPs_real[i],".17122014.txt"))
  fwrite(x, paste0("/mnt/polyomica/projects/GPD/Real IgG collections/ORCADES/", IGPs_our[i],".csv"))
  
  x = fread(paste0("/home/common/projects/IgGGlycans_GWAS/TWINSUK/TWINSUK.",IGPs_real[i],".20141208.txt"))
  fwrite(x, paste0("/mnt/polyomica/projects/GPD/Real IgG collections/TwinsUK/", IGPs_our[i],".csv"))

  x = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Vis/Vis.",IGPs_real[i],".16122014.txt"))
  fwrite(x, paste0("/mnt/polyomica/projects/GPD/Real IgG collections/Vis/", IGPs_our[i],".csv"))
  
  x = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Korcula/Korcula.",IGPs_real[i],".10122014.txt"))
  fwrite(x, paste0("/mnt/polyomica/projects/GPD/Real IgG collections/Korsula/", IGPs_our[i],".csv"))
  
  x = fread(paste0("/mnt/polyomica/projects/MV_igg_replication/3k_meta_23traits/meta_results/3K_meta_",IGPs_real[i],".csv"))
  fwrite(x, paste0("/mnt/polyomica/projects/GPD/Real IgG collections/mv_igg_repl_3k/", IGPs_our[i],".csv"))
  
}



GWAS_linear_combination=function(a,beta,se1,vary1=1,covm,N){
  vary2=sum(covm*(a%o%a))
  b=(beta%*%a)
  
  vary1_varg_n=se1^2+(beta[,1]^2)/N 
  
  se2=vary1_varg_n*(vary2/vary1)
  
  se2=se2-b^2/N
  
  se=sqrt(se2)
  
  b=b/sqrt(vary2)
  se=se/sqrt(vary2)
  
  out=cbind(b=b,se=se)
  colnames(out)=c("b","se")
  out=as.data.frame(out)
  return(out)
}



IGP4	IGP4+IGP5

cor = fread("/mnt/polyomica/projects/GPD/Real IgG collections/cor_matrix_real_IgG.txt")
cor = as.data.frame(cor)
rownames(cor) = colnames(cor)

cor_45 = matrix(c(1,0.3033895,0.3033895,1), nrow = 2)
cor_1213 = matrix(c(1,0.3965152,0.3965152,1), nrow = 2)

betas.select <- function(mydata, name){
  betas <- subset(mydata, select = c( rs_id, beta))
  gwas.id <- name
  colnames(betas)[2] <- paste("beta", gwas.id, sep ="_")
  return(betas)
}


main_func = function(x,y, names, cor) {
  print(names(x))
  print(names(y))
  dataFiles = list(x,y)
  nrow(x)
  
  snp_info = subset(x, select = c(rs_id,chr,bp))
  nrow(snp_info)
  
  list.betas  <- list(betas.select(x, names[1]), betas.select(y, names[2]))
  tmp <- Reduce(left_join, list.betas)
  nrow(tmp)
  
  se_gwas1 <- subset(x, select = c(rs_id, se))
  tmp_1 <- left_join(tmp, se_gwas1)
  colnames(tmp_1)[ncol(tmp_1)] <- "se_gwas1"
  nrow(tmp_1)
  
  ## add a column with minimal sample size across all GWASes for each snp 
  size.list <- lapply(dataFiles, subset, select = c(rs_id, n))
  names(size.list) <- seq(1:length(size.list))
  size.list_2 <- lapply(seq_along(size.list), function(i) {colnames(size.list[[i]])[2] <- paste("n", names(size.list)[i]) 
  size.list[[i]]})
  size.list.m <-Reduce(left_join, size.list_2)
  #table(betas.se.N$N_min, useNA = "always")
  size.list.m$N_min <- do.call(pmin, size.list.m[,2:ncol(size.list.m)])
  nrow(size.list)
  
  betas.se.N <- left_join(tmp_1, size.list.m[, c("rs_id", "N_min")])
  beta.matrix <- as.matrix(betas.se.N[,2:3])
  
  trait_lin_comb <- GWAS_linear_combination(a=c(1,1), beta = beta.matrix, se1 = betas.se.N$se_gwas1, covm=cor , N = betas.se.N$N_min)
  
  res = mutate(trait_lin_comb, Z=b/se,Pval=pchisq(Z^2,df=1,low=F))
  res$rs_id <- betas.se.N$rs_id
  
  res = left_join(res, snp_info)
  res$n = betas.se.N$N_min
  return(res)
}







###### 4 5 пики ########

x = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Orkney2000/Orkney2000.IGP4.17122014.txt"))
names(x)[1:3] = c("rs_id","chr", "bp")
y = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Orkney2000/Orkney2000.IGP5.17122014.txt"))
names(y)[1:3] = c("rs_id","chr", "bp")

res = main_func(x,y,c("IGP4","IGP5"), cor = cor_45)
eaf_dt = subset(x, select = c(rs_id,chr,bp, effect_allele,other_allele,EAF))

res = left_join(res,eaf_dt)
names(res) = c("beta","se","Z", "p", "name","chromosome","position", "n", "effect_allele","other_allele","EAF")
fwrite(res, paste0("/mnt/polyomica/projects/GPD/Real IgG collections/ORCADES/IGP4.csv"))



x = fread(paste0("/home/common/projects/IgGGlycans_GWAS/TWINSUK/TWINSUK.IGP4.20141208.txt"))
names(x)[1:3] = c("rs_id","chr", "bp")
y = fread(paste0("/home/common/projects/IgGGlycans_GWAS/TWINSUK/TWINSUK.IGP5.20141208.txt"))
names(y)[1:3] = c("rs_id","chr", "bp")

res = main_func(x,y,c("IGP4","IGP5"), cor = cor_45)
eaf_dt = subset(x, select = c(rs_id,chr,bp, effect_allele,other_allele,EAF))

res = left_join(res,eaf_dt)
names(res) = c("beta","se","Z", "p", "name","chromosome","position","n", "effect_allele","other_allele","EAF")
fwrite(res, paste0("/mnt/polyomica/projects/GPD/Real IgG collections/TwinsUK/IGP4.csv"))



x = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Vis/Vis.IGP4.16122014.txt"))
names(x)[1:3] = c("rs_id","chr", "bp")
y = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Vis/Vis.IGP5.16122014.txt"))
names(y)[1:3] = c("rs_id","chr", "bp")

res = main_func(x,y,c("IGP4","IGP5"), cor = cor_45)
eaf_dt = subset(x, select = c(rs_id,chr,bp, effect_allele,other_allele,EAF))

res = left_join(res,eaf_dt)
names(res) = c("beta","se","Z", "p", "name","chromosome","position", "n", "effect_allele","other_allele","EAF")
fwrite(res, paste0("/mnt/polyomica/projects/GPD/Real IgG collections/Vis/IGP4.csv"))



x = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Korcula/Korcula.IGP4.10122014.txt"))
names(x)[1:3] = c("rs_id","chr", "bp")
y = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Korcula/Korcula.IGP5.10122014.txt"))
names(y)[1:3] = c("rs_id","chr", "bp")

res = main_func(x,y,c("IGP4","IGP5"), cor = cor_45)
eaf_dt = subset(x, select = c(rs_id,chr,bp, effect_allele,other_allele,EAF))

res = left_join(res,eaf_dt)
names(res) = c("beta","se","Z", "p", "name","chromosome","position","n", "effect_allele","other_allele","EAF")
fwrite(res, paste0("/mnt/polyomica/projects/GPD/Real IgG collections/Korsula/IGP4.csv"))



x = fread(paste0("/mnt/polyomica/projects/MV_igg_replication/3k_meta_23traits/meta_results/3K_meta_IGP4.csv"))
names(x)[c(8:9,12) ] = c("beta","se","n")
y = fread(paste0("/mnt/polyomica/projects/MV_igg_replication/3k_meta_23traits/meta_results/3K_meta_IGP5.csv"))
names(y)[c(8:9,12)] = c("beta","se", "n")

#table(is.na(x$rs_id))
#table(is.na(res$b))

#head(res[is.na(res$b),])
#x[x$rs_id =="rs181866923",]
#y[y$rs_id =="rs181866923",]

res = main_func(x,y,c("IGP4","IGP5"), cor = cor_45)
eaf_dt = subset(x, select = c(rs_id,chr,bp, ea,ra,eaf_ma))


res = left_join(res,eaf_dt)
names(res) = c("beta_ma","se_ma","z_ma", "p_ma", "rs_id","chr","bp", "n_ma","ea","ra","eaf_ma")
fwrite(res, paste0("/mnt/polyomica/projects/GPD/Real IgG collections/mv_igg_repl_3k/IGP4.csv"))







######### 12 13 пики ################

#IGP11	IGP12+IGP13

x = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Orkney2000/Orkney2000.IGP12.17122014.txt"))
names(x)[1:3] = c("rs_id","chr", "bp")
y = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Orkney2000/Orkney2000.IGP13.17122014.txt"))
names(y)[1:3] = c("rs_id","chr", "bp")

res = main_func(x,y,c("IGP12","IGP13"), cor = cor_1213)
eaf_dt = subset(x, select = c(rs_id,chr,bp, effect_allele,other_allele,EAF))

res = left_join(res,eaf_dt)
names(res) = c("beta","se","Z", "p", "name","chromosome","position", "n", "effect_allele","other_allele","EAF")
fwrite(res, paste0("/mnt/polyomica/projects/GPD/Real IgG collections/ORCADES/IGP11.csv"))



x = fread(paste0("/home/common/projects/IgGGlycans_GWAS/TWINSUK/TWINSUK.IGP12.20141208.txt"))
names(x)[1:3] = c("rs_id","chr", "bp")
y = fread(paste0("/home/common/projects/IgGGlycans_GWAS/TWINSUK/TWINSUK.IGP13.20141208.txt"))
names(y)[1:3] = c("rs_id","chr", "bp")

res = main_func(x,y,c("IGP12","IGP13"), cor = cor_1213)
eaf_dt = subset(x, select = c(rs_id,chr,bp, effect_allele,other_allele,EAF))

res = left_join(res,eaf_dt)
names(res) = c("beta","se","Z", "p", "name","chromosome","position","n", "effect_allele","other_allele","EAF")
fwrite(res, paste0("/mnt/polyomica/projects/GPD/Real IgG collections/TwinsUK/IGP11.csv"))



x = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Vis/Vis.IGP12.16122014.txt"))
names(x)[1:3] = c("rs_id","chr", "bp")
y = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Vis/Vis.IGP13.16122014.txt"))
names(y)[1:3] = c("rs_id","chr", "bp")

res = main_func(x,y,c("IGP12","IGP13"), cor = cor_1213)
eaf_dt = subset(x, select = c(rs_id,chr,bp, effect_allele,other_allele,EAF))

res = left_join(res,eaf_dt)
names(res) = c("beta","se","Z", "p", "name","chromosome","position", "n", "effect_allele","other_allele","EAF")
fwrite(res, paste0("/mnt/polyomica/projects/GPD/Real IgG collections/Vis/IGP11.csv"))



x = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Korcula/Korcula.IGP12.10122014.txt"))
names(x)[1:3] = c("rs_id","chr", "bp")
y = fread(paste0("/home/common/projects/IgGGlycans_GWAS/Korcula/Korcula.IGP13.10122014.txt"))
names(y)[1:3] = c("rs_id","chr", "bp")

res = main_func(x,y,c("IGP12","IGP13"), cor = cor_1213)
eaf_dt = subset(x, select = c(rs_id,chr,bp, effect_allele,other_allele,EAF))

res = left_join(res,eaf_dt)
names(res) = c("beta","se","Z", "p", "name","chromosome","position","n", "effect_allele","other_allele","EAF")
fwrite(res, paste0("/mnt/polyomica/projects/GPD/Real IgG collections/Korsula/IGP11.csv"))



x = fread(paste0("/mnt/polyomica/projects/MV_igg_replication/3k_meta_23traits/meta_results/3K_meta_IGP12.csv"))
names(x)[c(8:9,12) ] = c("beta","se","n")
y = fread(paste0("/mnt/polyomica/projects/MV_igg_replication/3k_meta_23traits/meta_results/3K_meta_IGP13.csv"))
names(y)[c(8:9,12)] = c("beta","se", "n")

res = main_func(x,y,c("IGP12","IGP13"), cor = cor_1213)
eaf_dt = subset(x, select = c(rs_id,chr,bp, ea,ra,eaf_ma))

res = left_join(res,eaf_dt)
names(res) = c("beta_ma","se_ma","z_ma", "p_ma", "rs_id","chr","bp", "n_ma","ea","ra","eaf_ma")
fwrite(res, paste0("/mnt/polyomica/projects/GPD/Real IgG collections/mv_igg_repl_3k/IGP11.csv"))


#fread("/mnt/polyomica/projects/GPD/plasma_EUR_ma_withoutTwinsUK/PGP99/pgp99.csv")

