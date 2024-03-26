#Miami plot

library(data.table)
library(qqman)

predicted_IgG = fread("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/IgG_8peaks_pvals_predicted_on_10k_pearson_0.8.csv")

IgG_gwas_pred = data.frame(BP = predicted_IgG$bp, CHR = predicted_IgG$chr, P = predicted_IgG$p_min, SNP = predicted_IgG$rs_id)
IgG_gwas_pred_lite = IgG_gwas_pred[IgG_gwas_pred$P < 0.001,]
nrow(IgG_gwas_pred_lite)


real_IgG = fread("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/Plasma_36peaks_pvals_10k.csv")
#names(real_IgG)

#head(IgG_gwas_real)
IgG_gwas_real = data.frame(BP = real_IgG$bp, CHR = real_IgG$chr, P = real_IgG$p_min, SNP = real_IgG$rs_id)
IgG_gwas_real = IgG_gwas_real[!is.na(IgG_gwas_real$BP),]
IgG_gwas_real = IgG_gwas_real[!is.na(IgG_gwas_real$P),]
IgG_gwas_real_lite = IgG_gwas_real[IgG_gwas_real$P < 0.001,]
nrow(IgG_gwas_real_lite)

pdf("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/Maiami_plot_plasma_IgG_predicted.pdf", width = 16, height = 9)
par(mfrow=c(2,1))
par(mar=c(0,5,3,3))
manhattan(IgG_gwas_pred_lite,ylim=c(0,70),cex=2.2,cex.lab=2.5,font.lab=2,font.axis=2,cex.axis=1.6,las=2,font=4,
          chr="CHR", bp="BP", snp="SNP", p="P", col = c("#606060FF", "#D6ED17FF"),
          suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08/8))
par(mar=c(5,5,3,3))
manhattan(IgG_gwas_real_lite,ylim=c(70,0),cex=2.2,cex.lab=2.5,font.lab=2,font.axis=2,cex.axis=1.6,
          las=2,font=4,xlab="",xaxt="n",
          chr="CHR", bp="BP", snp="SNP", p="P", col = c("#5a613e", "#B1009B"),
          suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08/21))
dev.off()


svg("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/Maiami_plot_plasma_IgG__predicted.svg", width = 16, height = 9)
par(mfrow=c(2,1))
par(mar=c(0,5,3,3))
manhattan(IgG_gwas_pred_lite,ylim=c(0,70),cex=2.2,cex.lab=2.5,font.lab=2,font.axis=2,cex.axis=1.6,las=2,font=4,
          chr="CHR", bp="BP", snp="SNP", p="P", col = c("#606060FF", "#D6ED17FF"),
          suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08/8))
par(mar=c(5,5,3,3))
manhattan(IgG_gwas_real_lite,ylim=c(70,0),cex=2.2,cex.lab=2.5,font.lab=2,font.axis=2,cex.axis=1.6,
          las=2,font=4,xlab="",xaxt="n",
          chr="CHR", bp="BP", snp="SNP", p="P", col = c("#5a613e", "#B1009B"),
          suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08/21))
dev.off()
