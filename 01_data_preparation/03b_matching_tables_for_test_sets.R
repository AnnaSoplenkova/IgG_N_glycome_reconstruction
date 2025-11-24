#matching_tables_for_test_sets

igg<-data.table::fread('/mnt/polyomica/projects/GPD/test_set/QC_results/TestIgG_merge45_2025_QC/TestIgG_merge45_2025_QC_qced_glyco.csv')
plasma<-data.table::fread('/mnt/polyomica/projects/GPD/test_set/QC_results/TestPlasma_QC/TestPlasma_QC_qced_glyco.csv')


#intersection by ids
length(intersect(igg$Sample, plasma$Sample))

#taking ids of igg presented in plasma
ind<-match(igg$Sample,plasma$Sample)
igg<-igg[!is.na(ind),]

res_table_igg<-igg[,-2]
colnames(res_table_igg)[1]<-'ID'
data.table::fwrite(res_table_igg, '/mnt/polyomica/projects/GPD/test_set/IgG_test_set_merge45_2025.csv')

colnames(plasma)[1]<-'ID'
ind2<-match(plasma$ID,res_table_igg$ID)
plasma<-plasma[!is.na(ind2),]
res_table_plasma<-plasma[,-2]

#matching again
ind3<-match(res_table_igg$ID,res_table_plasma$ID)
res_table_plasma<-res_table_plasma[ind3,]

data.table::fwrite(res_table_plasma, '/mnt/polyomica/projects/GPD/test_set/plasma_test_set_36peaks_2023.csv')

colnames(res_table_plasma)<-paste0('Plasma_', colnames(res_table_plasma))
colnames(res_table_igg)<-paste0('IgG_', colnames(res_table_igg))
cor1<-cor(cbind(res_table_plasma[,4:(42-3)],res_table_igg[,4:24]))
library(corrplot)
png('/mnt/polyomica/projects/GPD/test_set/corr_heatmap_CEDAR.png')
corrplot(cor1,tl.cex=0.5)
dev.off()
