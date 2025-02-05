
#Make PGP93_FUC_C trait as linear combination of summary statistics
#PGP93	FUC-C	The percentage of core fucosylated structures in total plasma glycans	PGP1+PGP2+PGP4+PGP5+PGP6+PGP10+PGP11+PGP13+PGP15+PGP16+PGP20+PGP21+PGP28+PGP31+PGP32

library(data.table)

pgps = c(1,2,4,5,6,10,11,13,15,16,20,21,28,31,32)
names = paste0("/mnt/polyomica/projects/plasma_v2/11_EUR_META/01_GWAS_MAP/01_upload_data/02_unification_out/pgp", pgps, "/pgp", pgps, "_done.csv")

setwd("/mnt/polyomica/projects/GPD/peak_prediction")

cor = fread("data/cor_38_38.txt")
cor.matrix.filtered = as.matrix(cor[3:38,4:39])

cor.matrix.small = cor.matrix.filtered[pgps,pgps]

list.files = lapply(names,fread)




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

dataFiles = list.files


pgp = function(i) {
  dataFiles[[i]]$gwas_id = rep(paste0("pgp", pgps[i]), nrow(dataFiles[[i]]))
  return(dataFiles[[i]])
}

dataFiles = lapply(1:length(pgps), pgp)
dataFiles[[5]]

## combine all snp effects into one dataframe
betas.select <- function(mydata){
  betas <- subset(mydata, select = c( rs_id, beta))
  gwas.id <- unique(mydata$gwas_id)
  colnames(betas)[2] <- paste("beta", gwas.id, sep ="_")
  return(betas)
}

list.betas  <- lapply(dataFiles, betas.select)
tmp <- Reduce(left_join, list.betas)

se_gwas1 <- subset(dataFiles[[1]], select = c(rs_id, se))
tmp_1 <- left_join(tmp, se_gwas1)
colnames(tmp_1)[ncol(tmp_1)] <- "se_gwas1"


## add a column with minimal sample size across all GWASes for each snp 
size.list <- lapply(dataFiles, subset, select = c(rs_id, n))
names(size.list) <- seq(1:length(size.list))
size.list_2 <- lapply(seq_along(size.list), function(i) {colnames(size.list[[i]])[2] <- paste("n", names(size.list)[i]) 
size.list[[i]]})
size.list.m <-Reduce(left_join, size.list_2)
size.list.m$N_min <- do.call(pmin, size.list.m[,2:ncol(size.list.m)])

betas.se.N <- left_join(tmp_1, size.list.m[, c("rs_id", "N_min")])



beta.matrix <- as.matrix(betas.se.N[,2:(length(pgps)+1)])

snp_info = subset(dataFiles[[1]], select = c(rs_id,chr,bp))


  a1 = rep(1, length(pgps))
  trait_lin_comb <- GWAS_linear_combination(a=a1, beta = beta.matrix, se1 = betas.se.N$se_gwas1, covm=cor.matrix.small , N = betas.se.N$N_min)
  
  
  res = mutate(trait_lin_comb, Z=b/se,Pval=pchisq(Z^2,df=1,low=F))
  res$rs_id <- betas.se.N$rs_id
  res = left_join(res, snp_info) 
  res$trait <- "PGP93_FUC-C"
  res = res[,c(5:8, 1:4)]


rm(list.betas)
rm(dataFiles)
rm(size.list)
rm(size.list_2)

res1 = as.data.table(left_join(res, as.data.frame(betas.se.N)[,c(1, ncol(betas.se.N))] ))

file_name = paste0("/mnt/polyomica/projects/GPD/result_ma/PGP93_FUC_C.csv")
fwrite(res1, file_name)
