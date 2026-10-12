#бенчмарк: создаю сумстаты по GWIS

library(openxlsx)
library(data.table)
library(dplyr)
library(tidyr)


#PGP92 = PGP24+PGP30+PGP32+PGP36
coefs  = c(1,1,1,1)
IDs = c(24,30,32,36)
cor_manova = fread("/mnt/home/timoschuk/cormatrix_new.csv")

cor_manova = as.matrix(cor_manova)
cor_manova.need = cor_manova[IDs,IDs]



GWAS_linear_combination=function(a, beta, ses, covm){
  b <- beta %*% a
  v <- rowSums(ses * (ses %*% (covm * outer(a, a))))
  result = data.frame(b = b, se = sqrt(v))
  
  return(result)
}




setwd("/mnt/polyomica/projects/GPD/peak_prediction")

IgG_coef = names(coefs)

# все когорты с не европейцами
names = paste0("/mnt/polyomica/projects/plasma_v2/14_all_META/01_GWAS_MAP/pgp", IDs, ".csv")


a = fread(names[1])
a
eaf_dt = subset(a, select = c(rs_id,chr,bp, ea,ra,eaf_ma))

dataFiles <- lapply(names, fread)
print("downloading files now finished")

dataFiles[[1]]

snp_info = subset(dataFiles[[1]], select = c(rs_id,chr,bp))


pgp = function(i) {
  dataFiles[[i]]$gwas_id = rep(paste0("pgp", IDs[i]), nrow(dataFiles[[i]]))
  return(dataFiles[[i]])
}

dataFiles = lapply(1:length(IDs), pgp)
dataFiles[[3]]



ses.select <- function(mydata){
  ses <- subset(mydata, select = c( rs_id, se_ma))
  gwas.id <- unique(mydata$gwas_id)
  colnames(ses)[2] <- paste("se", gwas.id, sep ="_")
  return(ses)
}

list.ses  <- lapply(dataFiles, ses.select)
tmp <- Reduce(left_join, list.ses)


## combine all snp effects into one dataframe
betas.select <- function(mydata){
  betas <- subset(mydata, select = c( rs_id, beta_ma))
  gwas.id <- unique(mydata$gwas_id)
  colnames(betas)[2] <- paste("beta", gwas.id, sep ="_")
  return(betas)
}

list.betas  <- lapply(dataFiles, betas.select)
tmp_beta <- Reduce(left_join, list.betas)



tmp_1 <- left_join(tmp_beta, tmp)


## add a column with minimal sample size across all GWASes for each snp 
size.list <- lapply(dataFiles, subset, select = c(rs_id, n_ma))
names(size.list) <- seq(1:length(size.list))
size.list_2 <- lapply(seq_along(size.list), function(i) {colnames(size.list[[i]])[2] <- paste("n", names(size.list)[i]) 
size.list[[i]]})
size.list.m <-Reduce(left_join, size.list_2)
size.list.m$N_min <- do.call(pmin, size.list.m[,2:ncol(size.list.m)])

betas.se.N <- left_join(tmp_1, size.list.m[, c("rs_id", "N_min")])

names(betas.se.N)

beta.matrix <- as.matrix(betas.se.N[,2:5])
se.matrix <- as.matrix(betas.se.N[,6:9])





rm(list.betas)
rm(dataFiles)
rm(size.list)
rm(size.list_2)

trait_lin_comb <- GWAS_linear_combination(a=coefs, beta = beta.matrix, ses = se.matrix, covm=cor_manova.need)
res = mutate(trait_lin_comb, Z=b/se,Pval=pchisq(Z^2,df=1,low=F))
res$rs_id <- betas.se.N$rs_id
res = left_join(res, snp_info) 
res$trait <- "PGP107"

names(betas.se.N)

res_full = as.data.table(left_join(res[,c(5:8, 1:4)], betas.se.N[,c(1, 10)]))
PGP92_final = as.data.table(left_join(res_full, eaf_dt))




file_name = paste0("/mnt/polyomica/projects/GPD/benchmark/PGP92_GWIS.csv")
fwrite(PGP92_final, file_name)

