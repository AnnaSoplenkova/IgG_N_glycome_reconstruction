library(data.table)
library(qqman)

setwd("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS")

x = fread("/mnt/polyomica/projects/MV_TPNG_v2/Data_discovery/multi_traits_unfilt/B/meta_disc11")
x = x[,1:3]

temp = list.files(path = "/mnt/polyomica/projects/prs_glycomics_igg_project/gwas/gwas_11K/", pattern="*.csv")

read_files = function(i) {
  file_name = paste0("/mnt/polyomica/projects/prs_glycomics_igg_project/gwas/gwas_11K/",i)
  res = fread(file_name)
  return(res)
}
df_list = lapply(temp, read_files)
df_list = lapply(df_list, function(y) {left_join(y,x, by = c("SNP" = "rs_id"))})
table(is.na(df_list[[1]]$bp))


## add a column with minimal sample size across all GWASes for each snp 
size.list <- lapply(df_list, subset, select = c(SNP, chr, bp, p))
names(size.list) <- seq(1:length(size.list))
size.list_2 <- lapply(seq_along(size.list), function(i) {colnames(size.list[[i]])[4] <- paste("pval", names(size.list)[i])
size.list[[i]]})
size.list.m <-Reduce(left_join, size.list_2)

pvals = names(size.list.m)[startsWith(names(size.list.m), "pval")]

size.list.m$p_min <- do.call(pmin, size.list.m[, ..pvals])
fwrite(size.list.m, "/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/IgG_24peaks_pvals_real.csv")


IgG_gwas = data.frame(BP = size.list.m$bp, CHR = size.list.m$chr, P = size.list.m$p_min, SNP = size.list.m$SNP)
IgG_gwas = IgG_gwas[!is.na(IgG_gwas$BP),]
IgG_gwas = IgG_gwas[!is.na(IgG_gwas$P),]
IgG_gwas_lite = IgG_gwas[IgG_gwas$P < 0.001,]

head(IgG_gwas_lite)
nrow(IgG_gwas)



pdf("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MH_real_igg.pdf", width = 16, height = 9)
manhattan(IgG_gwas_lite[,c(2,1,4,3)], chr="CHR", bp="BP", snp="SNP", p="P", col = c("#5a613e", "#27d1a6"),
          suggestiveline = F, genomewideline = -log10(5e-08/21), ylim=c(0,70))
dev.off() 

svg("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MH_real_igg.svg", width = 16, height = 9)
manhattan(IgG_gwas_lite, chr="CHR", bp="BP", snp="SNP", p="P", col = c("#5a613e", "#27d1a6"),
          suggestiveline = F, genomewideline = -log10(5e-08/21), , ylim=c(0,70))
dev.off() 
