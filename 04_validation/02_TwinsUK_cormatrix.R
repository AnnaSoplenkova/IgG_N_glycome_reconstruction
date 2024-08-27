код рисует хитмап для тестовой выборки (плазма и ИгГ)

library(corrplot)
library(openxlsx)
library(dplyr)

setwd("/mnt/polyomica/projects/GPD/training_set/")

plasma_train = read.csv("/mnt/polyomica/projects/GPD/peak_prediction/data/plasma_training_set_36peaks_2023_sex_age_merge45.csv")
IgG_train = read.csv("/mnt/polyomica/projects/GPD/peak_prediction/data/IgG_training_set_2023_sex_age_merge45.csv")
#plasma_train = read.csv("plasma_training_set.csv")
#IgG_train = read.csv("IgG_training_set.csv")


names(IgG_train)[-(1:3)] = paste0("IgG_",1:21)
names(plasma_train)[-(1:3)] = paste0("plasma_", 1:36)

train_set = left_join(plasma_train, IgG_train)

nrow(train_set)
head(train_set)

corrected_to_gender_age <- train_set[, -(1:3)]

cor_matrix = as.matrix(cor(corrected_to_gender_age))

fwrite(cor_matrix, '/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MultiABEL/corrplot_plasma_IgG_TwinsUK.csv')


pdf('/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MultiABEL/corrplot_plasma_IgG_TwinsUK.pdf')
corrplot(cor_matrix)
dev.off()

pdf('/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MultiABEL/corrplot_abs_plasma_IgG_TwinsUK.pdf')
corrplot(abs(cor_matrix))
dev.off()
