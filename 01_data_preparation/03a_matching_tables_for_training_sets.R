#создает IgG и plasma файлы с одинаковыми айдишниками и одинаковыми по структуре (ID, пол/позраст, пики)

info<-readxl::read_excel('/home/common/projects/IgG_prediction_from_Plasma/data/Twins_Info.xlsx')
igg<-data.table::fread('/mnt/polyomica/projects/GPD/training_set/QC_results/TrainIgG_merge45_2025_QC/TrainIgG_merge45_2025_QC_qced_glyco.csv')
plasma<-data.table::fread('/mnt/polyomica/projects/GPD/training_set/QC_results/TrainPlasma_QC/TrainPlasma_QC_qced_glyco.csv')
plasma_ids<-data.table::fread('/home/common/projects/IgG_prediction_from_Plasma/data/TwinsLongitudinalPlasma_newIDs.txt')
ind<-match(plasma$Sample, plasma_ids$random_ID)
plasma$true_id<-plasma_ids$true_ID[ind]


#Plasma

ind2<-match(plasma$true_id, info$LAB_NO)
plasma$study_no<-info$STUDY_NO[ind2]
plasma$Sex<-info$Sex[ind2]
plasma$Age<-info$Age[ind2]

#Remove samples are not desribed in the info file
plasma<-plasma[!is.na(plasma$study_no),]

#Take the data for the youngest age of each person
sub_plasma<-lapply(unique(plasma$study_no), function(n){
  min_age<-min(plasma$Age[plasma$study_no==n])
  return(plasma[plasma$Age==min_age & plasma$study_no==n,])
})
res<-data.table::rbindlist(sub_plasma)
#ID, пол/возраст, пики (раньше было 39 пиков, теперь 36, поэтому добавилось "-3")
res_table<-res[,c(44-3,46-3,47-3,5:(43-3))]
colnames(res_table)[1]<-'ID'
res_table$ID<-as.numeric(res_table$ID)



#IgG
#intersection by ids
length(intersect(igg$Sample, res_table$ID))

#taking ids of igg presented in res_table (in plasma)
ind<-match(igg$Sample,res_table$ID)
igg<-igg[!is.na(ind),]

#Add age and sex info 
ind2<-match(igg$Sample, info$LAB_NO)
igg$study_no<-info$STUDY_NO[ind2]
igg$Sex<-info$Sex[ind2]
igg$Age<-info$Age[ind2]

ncol(igg)
res_table_igg<-igg[,c(1,28,29,6:26)]
colnames(res_table_igg)[1]<-'ID'
res_table_igg$ID<-as.numeric(res_table_igg$ID)
data.table::fwrite(res_table_igg, '/mnt/polyomica/projects/GPD/training_set/IgG_training_set_merge45_2025.csv')
ind3<-match(res_table$ID,res_table_igg$ID)
#removing that are not matched
res_table<-res_table[!is.na(ind3),]
#matching again
ind3<-match(res_table_igg$ID,res_table$ID)
res_table<-res_table[ind3,]
data.table::fwrite(res_table, '/mnt/polyomica/projects/GPD/training_set/plasma_training_set_36peaks_2023.csv')
colnames(res_table)<-paste0('Plasma_', colnames(res_table))
colnames(res_table_igg)<-paste0('IgG_', colnames(res_table_igg))
cor1<-cor(cbind(res_table[,4:(42-3)],res_table_igg[,4:24]))
library(corrplot)
png('/mnt/polyomica/projects/GPD/training_set/corr_heatmap.png')
corrplot(cor1,tl.cex=0.5)
dev.off()
