#код рисует хитмап для тренировочной выборки (плазма и ИгГ)

library(corrplot)
library(openxlsx)
library(dplyr)

setwd("/mnt/polyomica/projects/GPD/training_set/")

plasma_train = read.csv("plasma_training_set_ta.csv")
IgG_train = read.csv("IgG_training_set_ta.csv")

#plasma_train = read.csv("plasma_training_set.csv")
#IgG_train = read.csv("IgG_training_set.csv")

head(plasma_train)
head(IgG_train)

names = read.xlsx("Description_GP.xlsx", sheet = 1)
names(plasma_train)[4:ncol(plasma_train)] = names$trait

names = read.xlsx("Description_GP.xlsx", sheet = 2)
names(IgG_train)[4:ncol(IgG_train)] = names$trait

nrow(IgG_train)
nrow(plasma_train)

names(IgG_train)[-(1:3)] = paste0(names(IgG_train),"_IgG")[-(1:3)]
names(plasma_train)[-(1:3)] = paste0(names(plasma_train),"_plasma")[-(1:3)]



train_set = left_join(plasma_train, IgG_train)

row(train_set)
head(train_set)

cor_matrix = as.matrix(cor(train_set[,-(1:3)]))

head(train_set_matrix)


pdf('corrplot_plasma_IgG_training.pdf')
 corrplot(cor_matrix, method = "square", is.corr = T, order='hclust', hclust.method='ward.D2',
          tl.cex = 0.7, title = "Plasma_IgG_training_set", mar=c(0,0,1,0), , tl.col = "black")
 
dev.off()

pdf('corrplot_abs_plasma_IgG_training.pdf')
corrplot(abs(cor_matrix), method = "square", is.corr = T, order='hclust', hclust.method='ward.D2',
         tl.cex = 0.7, title = "Plasma_IgG_training_set", mar=c(0,0,1,0), , tl.col = "black")

dev.off()
