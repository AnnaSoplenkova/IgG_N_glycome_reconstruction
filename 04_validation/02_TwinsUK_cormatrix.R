#corrplot


#код рисует хитмап для тестовой выборки (плазма и ИгГ)

library(corrplot)
library(openxlsx)
library(dplyr)

setwd("/mnt/polyomica/projects/GPD/training_set/")


plasma_train = read.csv("/mnt/polyomica/projects/GPD/peak_prediction/data/plasma_training_set_117peaks_2025_sex_age.csv")
IgG_train = read.csv("/mnt/polyomica/projects/GPD/peak_prediction/data/IgG_training_set_merge45_2025_sex_age.csv")
#plasma_train = read.csv("plasma_training_set.csv")
#IgG_train = read.csv("IgG_training_set.csv")

plasma_train = plasma_train[,1:39]
head(plasma_train)

names(IgG_train)[-(1:3)] = paste0("IGP",1:21)
names(plasma_train)[-(1:3)] = paste0("PGP", 1:36)

train_set = left_join(plasma_train, IgG_train)

nrow(train_set)
head(train_set)

corrected_to_gender_age <- train_set[, -(1:3)]

table(is.na(corrected_to_gender_age$IgG_7))

names = c(paste0("PGP", 1:36), paste0("IGP", 1:21))
corrected_to_gender_age = corrected_to_gender_age[,names(corrected_to_gender_age) %in% names]

cor_matrix = as.matrix(cor(corrected_to_gender_age, use = "complete.obs"))

fwrite(cor_matrix, '/mnt/polyomica/projects/GPD/result_ma/02_IgG_new_models_with_non_EUR/corrplot_plasma_IgG_TwinsUK.csv')


pdf('/mnt/polyomica/projects/GPD/result_ma/02_IgG_new_models_with_non_EUR/corrplot_plasma_IgG_TwinsUK.pdf')
corrplot(cor_matrix, tl.cex = 0.7)
dev.off()

pdf('/mnt/polyomica/projects/GPD/result_ma/02_IgG_new_models_with_non_EUR/corrplot_abs_plasma_IgG_TwinsUK.pdf')
corrplot(abs(cor_matrix))
dev.off()
