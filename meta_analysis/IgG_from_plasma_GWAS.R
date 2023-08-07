# generates IgG summary statistics from a linear combination of plasma summary statistics
# main function from https://github.com/Defrag1236/MV_IGG/blob/master/3K_replication/MANOVA_6snps/04-replicating_linear_combinations_coefficeints.R

library(openxlsx)
library(data.table)
library(dplyr)
library(tidyr)

setwd("/mnt/polyomica/projects/GPD/peak_prediction")
coefs = read.xlsx("OLS/coeffs.xlsx")

IgG_coef = names(coefs)[startsWith(names(coefs), "Estimate")]
names = paste0("/mnt/polyomica/projects/plasma_v2/04_meta_analisys/01_meta_analisys/meta_analysis_filtered/pgp", 1:36,"_done.csv")

dataFiles <- lapply(names, fread)
print("downloading files now finished")


snp_info = subset(dataFiles[[1]], select = c(rs_id,chr,bp))


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

cor = fread("data/cor_38_38.txt")
cor.matrix.filtered = as.matrix(cor[3:38,4:39])


beta.matrix <- as.matrix(betas.se.N[,2:(36+1)])



IgG_peak = function (i) {
  a1 = coefs[,names(coefs) == IgG_coef[i]][-1]
  trait_lin_comb <- GWAS_linear_combination(a=a1, beta = beta.matrix, se1 = betas.se.N$se_gwas1, covm=cor.matrix.filtered , N = betas.se.N$N_min)
  

  res = mutate(trait_lin_comb, Z=b/se,Pval=pchisq(Z^2,df=1,low=F))
  res$rs_id <- betas.se.N$rs_id
  res = left_join(res, snp_info) 
  res$trait <- paste0("IGP_", i)
  print(paste0("IGP_", i))
  return(res[,c(5:8, 1:4)])
}

rm(list.betas)
rm(dataFiles)
rm(size.list)
rm(size.list_2)
IgG_peak_list = lapply(1:22, IgG_peak)


make_res_list = function(i) {
  IgG_peak_list[[i]] = as.data.table(left_join(IgG_peak_list[[i]], betas.se.N[,c(1, 39)]))
}

res_list = lapply(1:22, make_res_list)

save_df = function(i) {
  file_name = paste0("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/IGP_",i,".csv")
  fwrite(res_list[[i]], file_name)
}

lapply(1:22,save_df)
