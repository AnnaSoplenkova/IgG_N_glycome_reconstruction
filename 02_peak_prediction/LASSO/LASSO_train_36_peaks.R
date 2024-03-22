#LASSO + plot actual/predict on training set with 36 plasma peaks

library(glmnet)
library(data.table)
library(dplyr)
library(Metrics)
library(ggplot2)
library(gridExtra)
library(tibble)

setwd("/mnt/polyomica/projects/GPD/peak_prediction")

# беру тренировочные сеты
#plasma_train = fread("data/plasma_training_set_36peaks_2023.csv")
#IgG_train = fread("data/IgG_training_set_2023.csv")
plasma_train = fread("data/plasma_training_set_36peaks_2023_sex_age_merge45.csv")
IgG_train = fread("data/IgG_training_set_2023_sex_age_merge45.csv")
#plasma_train = fread("data/plasma_training_set_36peaks_2023_sex_age.csv")
#IgG_train = fread("data/IgG_training_set_2023_sex_age.csv")

n_igg = 21

# merge 2 datasets in one
names(IgG_train)[-(1:3)] = paste0(names(IgG_train),"_IgG")[-(1:3)]
names(plasma_train)[-(1:3)] = paste0(names(plasma_train),"_plasma")[-(1:3)]
train_set = left_join(plasma_train, IgG_train)

IgG_names = names(train_set)[endsWith(names(train_set), "IgG")]
plasma_names = names(train_set)[endsWith(names(train_set), "plasma")]

###############################
##### functions for LASSO #####
###############################

BestLambda = function(i, train) {
  plasma = as.matrix(train[, ..plasma_names])
  IgG = as.matrix(train[, ..IgG_names])[,i]
  
  set.seed(111)
  #perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(x = plasma, y = IgG, alpha = 1)
  #find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.1se
  
  return(best_lambda)
}


LASSO = function(i, train, test, lambda) {
  plasma = as.matrix(train[, ..plasma_names])
  IgG = as.matrix(train[, ..IgG_names])[,i]
  best_lambda = lambda[i]
  
  test_plasma = as.matrix(test[, ..plasma_names])
  
  
  best_model <- glmnet(x = plasma, y = IgG, alpha = 1, lambda = best_lambda)
  coef(best_model)
  x = as.matrix(test[, ..IgG_names])[,i]
  y = predict(best_model, newx = test_plasma)
  
  return(c(rmse = rmse(x,y), rel_err = mean(abs(x-y)/x),
           pearson = cor(x,y), spearman = cor(x,y,method = "spearman")) )
} 


#draw LASSO plot (actual/predict)
LASSO_plot = function(i, train, test, lambda) {
  
  plasma = as.matrix(train[, ..plasma_names])
  IgG = as.matrix(train[, ..IgG_names])[,i]
  best_lambda = lambda[i]
  
  test_plasma = as.matrix(test[, ..plasma_names])
  
  best_model <- glmnet(x = plasma, y = IgG, alpha = 1, lambda = best_lambda)
  coef(best_model)
  x = as.matrix(test[, ..IgG_names])[,i]
  y = predict(best_model, newx = test_plasma)
  
  
  df1 = data.frame(x,y)
  df1 = df1[order(x),]
  n = length(x)
  df = data.frame(element_of_sample = c(1:n, 1:n),
                  type = c(rep("actual", n), rep("predict", n)),
                  value = c(df1[,1],df1[,2]))
  plot = ggplot(df, aes(x = element_of_sample, y = value, color = type)) +
    geom_point() + ggtitle(IgG_names[i])
  return(plot)
}


LASSO_coefs = function(i, train, lambda) {
  plasma = as.matrix(train[, ..plasma_names])
  IgG = as.matrix(train[, ..IgG_names])[,i]
  best_lambda = lambda[i]
  
  best_model <- glmnet(x = plasma, y = IgG, alpha = 1, lambda = best_lambda)
  coef = coef(best_model)
  coef = rbind(lambda = best_lambda, coef)
  a = rownames(coef)
  coef = data.frame(coeff = coef[,1])
  rownames(coef) = a
  
  names(coef) = paste0("coef_", IgG_names[i])
  return(coef)
}

####################################
##### кросс-валидация 5 частей #####
####################################

set.seed(111)
N = floor(nrow(train_set)*.2)

dt = sort(sample(nrow(train_set), N))
df1<-train_set[dt,]
rest<-train_set[-dt,]

dt = sort(sample(nrow(rest), N))
df2<-rest[dt,]
rest<-rest[-dt,]

dt = sort(sample(nrow(rest), N))
df3<-rest[dt,]
rest<-rest[-dt,]

dt = sort(sample(nrow(rest), N))
df4<-rest[dt,]
df5<-rest[-dt,]

df = list(df1,df2,df3,df4,df5)


#######################################
##### запускаю с кросс-валидацией #####
#######################################

list_coefs = list()
for (i in 1:5) {
  test = df[[i]]
  train = do.call(rbind, df[-i])
  
  lambda = lapply(1:n_igg, BestLambda, train = train)
  lambda = unlist(lambda)
  
  coeffs = lapply(1:n_igg, LASSO_coefs, train = train, lambda = lambda)
  coeffs = do.call(cbind,coeffs)
  
  res = sapply(1:n_igg, LASSO, train = train, test = test, lambda = lambda)
  res = as.data.frame(t(res))
  
  if (i == 1) {
    RES = res
    names(RES) = paste0(names(RES),"_set",i)
    plot_list = lapply(1:n_igg, LASSO_plot, train = train, test = test, lambda = lambda) }
  if (i %in% 2:5) {
    names(res) = paste0(names(res),"_set",i)
    RES = cbind(RES, res)
  }
  
  list_coefs[[i]] = coeffs
}

Mean_coefs = bind_rows(list_coefs[[1]] %>% rownames_to_column(), list_coefs[[2]] %>% rownames_to_column(), 
                       list_coefs[[3]] %>% rownames_to_column(), list_coefs[[4]] %>% rownames_to_column(), 
                       list_coefs[[5]] %>% rownames_to_column()) %>%  group_by(rowname) %>% 
  summarise_all(mean)

target = rownames(list_coefs[[1]])
Mean_coefs = as.data.frame(Mean_coefs)
Mean_coefs = Mean_coefs[match(target, Mean_coefs$rowname),
                        !(startsWith(names(Mean_coefs), "Pr"))]

list_coefs[[6]] = Mean_coefs
names(list_coefs) = c(1:5, "Mean_coefs")

#####################################################
##### вывожу средние от статистик по 5 запускам #####
#####################################################

RES$rmse = apply(RES[,startsWith(names(RES), "rmse")],1,mean)
RES$rel_err = apply(RES[,startsWith(names(RES), "rel_err")],1,function(x) {mean(abs(x))})
RES$pearson = apply(RES[,startsWith(names(RES), "pearson")],1,mean)
RES$spearman = apply(RES[,startsWith(names(RES), "spearman")],1,mean)
rownames(RES) = IgG_names

write.xlsx(RES, "LASSO/stats.xlsx", rowNames = T)

wb <- createWorkbook()
sheetnames <- as.character(names(list_coefs))
sheets <- lapply(sheetnames, addWorksheet, wb = wb)
for (i in 1:length(sheetnames)) {
  writeData(wb, sheet = sheetnames[i], x = list_coefs[[i]], rowNames = T)
}
saveWorkbook(wb, 'LASSO/coeffs.xlsx', overwrite = T)

pdf('LASSO/plot_actual_predict_mq.pdf')
marrangeGrob(grobs = plot_list, ncol = 1, nrow = 2)
dev.off()
