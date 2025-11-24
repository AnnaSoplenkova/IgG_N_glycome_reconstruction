#Тренировочный и тестовый набор правились каждый на свои параметры (среднее и сд)
library(data.table)
library(dplyr)


setwd("/mnt/polyomica/projects/GPD/peak_prediction")

#file.copy("/mnt/polyomica/projects/GPD/test_set/IgG_test_set_merge45_2025.csv", "/mnt/polyomica/projects/GPD/peak_prediction/data")
#file.copy("/mnt/polyomica/projects/GPD/training_set/IgG_training_set_merge45_2025.csv", "/mnt/polyomica/projects/GPD/peak_prediction/data")


# беру сеты
plasma_train = as.data.frame(fread("data/plasma_training_set_117pgp_2025.csv"))
IgG_train = as.data.frame(fread("data/IgG_training_set_merge45_2025.csv"))
plasma_test = as.data.frame(fread("data/plasma_test_set_117pgp_2025.csv"))
IgG_test = as.data.frame(fread("data/IgG_test_set_merge45_2025.csv"))


# correction to sex and age
for(i in 4:ncol(plasma_train)){
  plasma_train[, i] <- residuals(lm(plasma_train[,i] ~ plasma_train$Sex + plasma_train$Age, 
                                    na.action = na.exclude))
  plasma_test[, i] <- residuals(lm(plasma_test[,i] ~ plasma_test$Sex + plasma_test$Age, 
                                   na.action = na.exclude))
}

for(i in 4:ncol(IgG_train)){
  IgG_train[, i] <- residuals(lm(IgG_train[,i] ~ IgG_train$Sex + IgG_train$Age, 
                                 na.action = na.exclude))
  IgG_test[, i] <- residuals(lm(IgG_test[,i] ~ IgG_test$Sex + IgG_test$Age, 
                                na.action = na.exclude))
}



# унификация

mean_plasma = unlist(lapply(4:120, function(i) {mean(plasma_train[, i])}))
sd_plasma = unlist(lapply(4:120, function(i) {sd(plasma_train[, i])}))
mean_IgG = unlist(lapply(4:24, function(i) {mean(IgG_train[, i])}))
sd_IgG = unlist(lapply(4:24, function(i) {sd(IgG_train[, i])}))

plasma_train[, 4:120] = t(apply(plasma_train[, 4:120], 1, function(x) {(x - mean_plasma)/sd_plasma}))
IgG_train[, 4:24] = t(apply(IgG_train[, 4:24], 1, function(x) {(x - mean_IgG)/sd_IgG}))


mean_plasma = unlist(lapply(4:120, function(i) {mean(plasma_test[, i])}))
sd_plasma = unlist(lapply(4:120, function(i) {sd(plasma_test[, i])}))
mean_IgG = unlist(lapply(4:24, function(i) {mean(IgG_test[, i])}))
sd_IgG = unlist(lapply(4:24, function(i) {sd(IgG_test[, i])}))

plasma_test[, 4:120] = t(apply(plasma_test[, 4:120], 1, function(x) {(x - mean_plasma)/sd_plasma}))
IgG_test[, 4:24] = t(apply(IgG_test[, 4:24], 1, function(x) {(x - mean_IgG)/sd_IgG}))



plasma_train = as.data.table(plasma_train)
IgG_train = as.data.table(IgG_train)
plasma_test = as.data.table(plasma_test)
IgG_test = as.data.table(IgG_test)



fwrite(plasma_train, "data/plasma_training_set_117peaks_2025_sex_age.csv")
fwrite(IgG_train, "data/IgG_training_set_merge45_2025_sex_age.csv")
fwrite(plasma_test, "data/plasma_test_set_117peaks_2025_sex_age.csv")
fwrite(IgG_test, "data/IgG_test_set_merge45_2025_sex_age.csv")
