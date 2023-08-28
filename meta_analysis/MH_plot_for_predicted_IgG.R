# MH plot for predicted IgG on MA plasma 7k EUR

library(data.table)
library(qqman)

setwd("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS")


#порог - пирсон 88%
peaks = c(3,6,9,10,11,12,13,16)

read_files = function(i) {
  file_name = paste0("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/IGP_",i,".csv")
  res = fread(file_name)
  return(res)
}

df_list = lapply(peaks, read_files)

## add a column with minimal sample size across all GWASes for each snp 
size.list <- lapply(df_list, subset, select = c(rs_id, chr, bp,trait, Pval))
names(size.list) <- seq(1:length(size.list))
size.list_2 <- lapply(seq_along(size.list), function(i) {colnames(size.list[[i]])[4] <- paste("trait", names(size.list)[i])
colnames(size.list[[i]])[5] <- paste("pval", names(size.list)[i])
size.list[[i]]})
size.list.m <-Reduce(left_join, size.list_2)

pvals = names(size.list.m)[startsWith(names(size.list.m), "pval")]
traits_0 = names(size.list.m)[startsWith(names(size.list.m), "trait")]
traits = unlist(size.list.m[, ..traits_0][1])

size.list.m$p_min <- do.call(pmin, size.list.m[, ..pvals])
size.list.m$main_trait<-traits[apply(size.list.m[, ..pvals],1,which.min)]

IgG_gwas = data.frame(BP = size.list.m$bp, CHR = size.list.m$chr, P = size.list.m$p_min, SNP = size.list.m$rs_id)
IgG_gwas_lite = IgG_gwas[IgG_gwas$P < 0.001,]
nrow(IgG_gwas_lite)

pdf("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MH.pdf", width = 16, height = 9)
manhattan(IgG_gwas_lite, chr="CHR", bp="BP", snp="SNP", p="P", col = c("#606060FF", "#D6ED17FF"),
          suggestiveline = F, genomewideline = -log10(5e-08/8))
dev.off() 

svg("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MH.svg", width = 16, height = 9)
manhattan(IgG_gwas_lite, chr="CHR", bp="BP", snp="SNP", p="P", col = c("#606060FF", "#D6ED17FF"),
          suggestiveline = F, genomewideline = -log10(5e-08/8))
dev.off() 
