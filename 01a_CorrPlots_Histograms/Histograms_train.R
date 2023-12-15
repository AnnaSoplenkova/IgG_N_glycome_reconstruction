#код рисует гистограммы c плотностью для тренировочной выборки (плазма и ИгГ)
#с поправкой на пол и возраст и без

library(openxlsx)
library(dplyr)
library(ggplot2)
library(gridExtra)

setwd("/mnt/polyomica/projects/GPD/training_set/")

plasma_train = read.csv("plasma_training_set_ta.csv")
IgG_train = read.csv("IgG_training_set_ta.csv")

#plasma_train = read.csv("plasma_training_set.csv")
#IgG_train = read.csv("IgG_training_set.csv")

names(IgG_train)[-(1:3)] = paste0(names(IgG_train),"_IgG")[-(1:3)]
names(plasma_train)[-(1:3)] = paste0(names(plasma_train),"_plasma")[-(1:3)]

train_set = left_join(plasma_train, IgG_train)
names(train_set)

table(train_set$ID == 10002)

corrected_to_gender_age <- train_set[, -(2:3)]
train_set$Sex = ifelse(train_set$Sex == "F", 0, 1)
for(i in 2:ncol(corrected_to_gender_age)){
  corrected_to_gender_age[, i] <- residuals(lm(train_set[,i+2] ~ train_set$Sex + train_set$Age, 
                                               na.action = na.exclude))
}

#рисую гистограммы для чистых данных
gplot_train = function(i) {
  a = ggplot(train_set, aes(x = .data[[names(train_set)[i] ]])) + 
    geom_histogram(aes(y = ..density..), colour = "gray32", fill = "gray50", alpha = 0.25) +
    geom_density(lwd = 0.6, colour = 4, fill = 4, alpha = 0.2) +
    theme_classic()
  
  return(a)
}

plot_list = lapply(c(4:ncol(train_set)), gplot_train)

pdf('hist_training_set.pdf')
marrangeGrob(grobs = plot_list, ncol = 2, nrow = 2)
dev.off()


#рисую гистограммы для поправленных на пол и возраст данных
head(corrected_to_gender_age)
gplot_train_corrected_sex_age = function(i) {
  a = ggplot(corrected_to_gender_age, aes(x = .data[[names(corrected_to_gender_age)[i] ]])) + 
    geom_histogram(aes(y = ..density..), colour = "gray32", fill = "gray50", alpha = 0.25) +
    geom_density(lwd = 0.6, colour = 4, fill = 4, alpha = 0.2) +
    theme_classic()
  
  return(a)
}

plot_list = lapply(c(2:ncol(corrected_to_gender_age)), gplot_train_corrected_sex_age)

pdf('hist_training_set_corrected_sex_age.pdf')
marrangeGrob(grobs = plot_list, ncol = 2, nrow = 2)
dev.off()

