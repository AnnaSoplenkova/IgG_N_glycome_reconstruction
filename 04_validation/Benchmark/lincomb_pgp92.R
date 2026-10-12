# Бенчмарк: сравниваем работу GWIS и Tsepilov
# создаем наш реконструированный гвас

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




setwd("/mnt/polyomica/projects/GPD/peak_prediction")

IgG_coef = names(coefs)

# все когорты с не европейцами
names = paste0("/mnt/polyomica/projects/plasma_v2/14_all_META/01_GWAS_MAP/pgp", IDs, ".csv")


a = fread(names[1])
eaf_dt = subset(a, select = c(rs_id,chr,bp, ea,ra,eaf_ma))

dataFiles <- lapply(names, fread)
print("downloading files now finished")

dataFiles[[1]]

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

pgp = function(i) {
  dataFiles[[i]]$gwas_id = rep(paste0("pgp", IDs[i]), nrow(dataFiles[[i]]))
  return(dataFiles[[i]])
}

dataFiles = lapply(1:length(IDs), pgp)
dataFiles[[3]]

## combine all snp effects into one dataframe
betas.select <- function(mydata){
  betas <- subset(mydata, select = c( rs_id, beta_ma))
  gwas.id <- unique(mydata$gwas_id)
  colnames(betas)[2] <- paste("beta", gwas.id, sep ="_")
  return(betas)
}

list.betas  <- lapply(dataFiles, betas.select)
tmp <- Reduce(left_join, list.betas)

se_gwas1 <- subset(dataFiles[[1]], select = c(rs_id, se_ma))
tmp_1 <- left_join(tmp, se_gwas1)
colnames(tmp_1)[ncol(tmp_1)] <- "se_gwas1"


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
head(betas.se.N)


rm(list.betas)
rm(dataFiles)
rm(size.list)
rm(size.list_2)

print("aaoaoa")

trait_lin_comb <- GWAS_linear_combination(a=coefs, beta = beta.matrix, se1 = betas.se.N$se_gwas1, covm=cor_manova.need , N = betas.se.N$N_min)
res = mutate(trait_lin_comb, Z=b/se,Pval=pchisq(Z^2,df=1,low=F))
res$rs_id <- betas.se.N$rs_id
res = left_join(res, snp_info) 
res$trait <- "PGP107"

names(betas.se.N)

res_full = as.data.table(left_join(res[,c(5:8, 1:4)], betas.se.N[,c(1, 7)]))
PGP92_final = as.data.table(left_join(res_full, eaf_dt))




file_name = paste0("/mnt/polyomica/projects/GPD/benchmark/PGP92_linear_combination.csv")
fwrite(PGP92_final, file_name)

