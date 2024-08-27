# MH plot for predicted IgG on MA plasma 10k EUR

library(data.table)
library(qqman)

#renv::install("MASS@7.3-60")

setwd("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS")


#порог - пирсон 80%
peaks = c(2,3,4,5,7,8,9,10,11,12,15,20)
n_peaks = length(peaks)

read_files = function(i) {
  #file_name = paste0("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/clean/10k_IGP_",i,"_clean_done.csv")
  file_name = paste0("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/unification_out/igp",i,"/igp",i,"_done.csv")
  res = fread(file_name)
  return(res)
}

df_list = lapply(peaks, read_files)

for (i in 1:length(peaks)) {
  df_list[[i]]$trait = paste0("IGP_", peaks[i])
}

df_list[[1]]

## add a column with minimal sample size across all GWASes for each snp 
size.list <- lapply(df_list, subset, select = c(rs_id, chr, bp,trait, p))
names(size.list) <- seq(1:length(size.list))
size.list_2 <- lapply(seq_along(size.list), function(i) {colnames(size.list[[i]])[4] <- paste("trait", names(size.list)[i])
colnames(size.list[[i]])[5] <- paste("pval", names(size.list)[i])
size.list[[i]]})
size.list.m <-Reduce(left_join, size.list_2)

pvals = names(size.list.m)[startsWith(names(size.list.m), "pval")]
traits_0 = names(size.list.m)[startsWith(names(size.list.m), "trait")]
traits = unlist(size.list.m[, ..traits_0][1])

size.list.m$p_min <- do.call(pmin, size.list.m[, ..pvals])
size.list.m[is.na(size.list.m$p_min),]
size.list.m = size.list.m[!is.na(size.list.m$p_min),]
size.list.m$main_trait<-traits[unlist(apply(size.list.m[, ..pvals],1,which.min))]

fwrite(size.list.m, "/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/IgG_8peaks_pvals_predicted_on_10k_pearson_0.8.csv")

IgG_gwas = data.frame(BP = size.list.m$bp, CHR = size.list.m$chr, P = size.list.m$p_min, SNP = size.list.m$rs_id)
IgG_gwas_lite = IgG_gwas[IgG_gwas$P < 0.001,]
nrow(IgG_gwas_lite)

pdf("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MH_10k_pearson_0.8.pdf", width = 16, height = 9)
manhattan(IgG_gwas_lite, chr="CHR", bp="BP", snp="SNP", p="P", col = c("#606060FF", "#D6ED17FF"),
          suggestiveline = F, genomewideline = -log10(5e-08/n_peaks))
dev.off() 

svg("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MH_10k_pearson_0.8.svg", width = 16, height = 9)
manhattan(IgG_gwas_lite, chr="CHR", bp="BP", snp="SNP", p="P", col = c("#606060FF", "#D6ED17FF"),
          suggestiveline = F, genomewideline = -log10(5e-08/n_peaks))
dev.off() 
