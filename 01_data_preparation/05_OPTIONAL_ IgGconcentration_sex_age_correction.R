library(data.table)
library(dplyr)


setwd("/mnt/polyomica/projects/GPD/training_set")

# беру сеты
plasma_train = as.data.frame(fread("plasma_training_set_36peaks_merge45_2023.csv"))
IgG_train = as.data.frame(fread("IgG_training_set_merge45_2023.csv"))

setwd("/mnt/polyomica/projects/GPD/test_set")
plasma_test = as.data.frame(fread("plasma_test_set_36peaks_merge45_2023.csv"))
IgG_test = as.data.frame(fread("IgG_test_set_merge45_2023.csv"))

df = plasma_train
df[,-(1:3)] = exp(df[,-(1:3)])

denominator = df$GP1 + df$GP2 + df$GP3 + df$GP4 + df$GP5 + df$GP6 + df$GP7 + df$GP8 + df$GP9 + df$GP10 +
  df$GP11 + df$GP12 + df$GP13 + df$GP14 + df$GP15 + df$GP16 + df$GP17 + df$GP18 + df$GP19 + df$GP20 +
  df$GP21 + df$GP22 + df$GP23 + df$GP24 + df$GP25 + df$GP26 + df$GP27 + df$GP28 + df$GP29 + df$GP30 +
  df$GP31 + df$GP32 + df$GP33 + df$GP34 + df$GP35 + df$GP36
nominator = df$GP1 + df$GP4 + df$GP5 + df$GP6 + 1/2*df$GP16
df$GP37 = nominator/denominator
df[,-(1:3)] = log(df[,-(1:3)])
plasma_train = df

df = plasma_test
df[,-(1:3)] = exp(df[,-(1:3)])
denominator = df$GP1 + df$GP2 + df$GP3 + df$GP4 + df$GP5 + df$GP6 + df$GP7 + df$GP8 + df$GP9 + df$GP10 +
  df$GP11 + df$GP12 + df$GP13 + df$GP14 + df$GP15 + df$GP16 + df$GP17 + df$GP18 + df$GP19 + df$GP20 +
  df$GP21 + df$GP22 + df$GP23 + df$GP24 + df$GP25 + df$GP26 + df$GP27 + df$GP28 + df$GP29 + df$GP30 +
  df$GP31 + df$GP32 + df$GP33 + df$GP34 + df$GP35 + df$GP36
nominator = df$GP1 + df$GP4 + df$GP5 + df$GP6 + 1/2*df$GP16
df$GP37 = nominator/denominator
df[,-(1:3)] = log(df[,-(1:3)])
plasma_test = df


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

mean_plasma = unlist(lapply(4:40, function(i) {mean(plasma_train[, i])}))
sd_plasma = unlist(lapply(4:40, function(i) {sd(plasma_train[, i])}))
mean_IgG = unlist(lapply(4:24, function(i) {mean(IgG_train[, i])}))
sd_IgG = unlist(lapply(4:24, function(i) {sd(IgG_train[, i])}))


plasma_train[, 4:40] = t(apply(plasma_train[, 4:40], 1, function(x) {(x - mean_plasma)/sd_plasma}))
plasma_test[, 4:40] = t(apply(plasma_test[, 4:40], 1, function(x) {(x - mean_plasma)/sd_plasma}))
IgG_train[, 4:24] = t(apply(IgG_train[, 4:24], 1, function(x) {(x - mean_IgG)/sd_IgG}))
IgG_test[, 4:24] = t(apply(IgG_test[, 4:24], 1, function(x) {(x - mean_IgG)/sd_IgG}))



plasma_train = as.data.table(plasma_train)
IgG_train = as.data.table(IgG_train)
plasma_test = as.data.table(plasma_test)
IgG_test = as.data.table(IgG_test)

setwd("/mnt/polyomica/projects/GPD/peak_prediction")

fwrite(plasma_train, "data/plasma_training_set_36peaks_2023_sex_age_IgGconcentration.csv")
fwrite(IgG_train, "data/IgG_training_set_2023_sex_age_IgGconcentration.csv")
fwrite(plasma_test, "data/plasma_test_set_36peaks_2023_sex_age_IgGconcentration.csv")
fwrite(IgG_test, "data/IgG_test_set_2023_sex_age_IgGconcentration.csv")
fwrite(plasma_test, "data/plasma_test_set_36peaks_2023_sex_age_IgGconcentration.csv")
fwrite(IgG_test, "data/IgG_test_set_2023_sex_age_IgGconcentration.csv")
