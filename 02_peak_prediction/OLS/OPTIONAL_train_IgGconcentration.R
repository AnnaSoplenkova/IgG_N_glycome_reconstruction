#OLS + plot actual/predict on training set with 36 plasma peaks + IgGconcentration (37 peak)

library(data.table)
library(dplyr)
library(Metrics)
library(ggplot2)
library(gridExtra)
library(tibble)
library(openxlsx)

setwd("/mnt/polyomica/projects/GPD/peak_prediction")

# беру тренировочные сеты
#plasma_train = fread("data/plasma_training_set_36peaks_2023.csv")
#IgG_train = fread("data/IgG_training_set_2023.csv")
#plasma_train = fread("data/plasma_training_set_36peaks_2023_sex_age_merge45.csv")
#IgG_train = fread("data/IgG_training_set_2023_sex_age_merge45.csv")
plasma_train = fread("data/plasma_training_set_36peaks_2023_sex_age_IgGconcentration_merge45.csv")
IgG_train = fread("data/IgG_training_set_2023_sex_age_IgGconcentration_merge45.csv")

# merge 2 datasets in one
names(IgG_train)[-(1:3)] = paste0(names(IgG_train),"_IgG")[-(1:3)]
names(plasma_train)[-(1:3)] = paste0(names(plasma_train),"_plasma")[-(1:3)]
train_set = left_join(plasma_train, IgG_train)

#OLS
IgG_names = names(train_set)[endsWith(names(train_set), "IgG")]
plasma_names = names(train_set)[endsWith(names(train_set), "plasma")]


# function for OLS
OLS = function(name, train, test) {
  f <- as.formula(paste(name, "~", paste(plasma_names, collapse="+")))
  mod = lm(f, data = train)
  x = unlist(test[,..name])
  y = predict(mod,newdata=test)
  
  return(c(rmse = rmse(x,y), rel_err = mean(abs(x-y)/x),
           pearson = cor(x,y), spearman = cor(x,y,method = "spearman")) )
} 

#draw OLS plot (actual/predict)
OLS_plot = function(name, train, test) {
  f <- as.formula(paste(name, "~", paste(names(train)[c(1:39)], collapse="+")))
  mod = lm(f, data = train)
  x = unlist(test[,..name])
  y = predict(mod,newdata=test)
  
  df1 = data.frame(x,y)
  df1 = df1[order(x),]
  n = length(x)
  df = data.frame(element_of_sample = c(1:n, 1:n),
                  type = c(rep("actual", n), rep("predict", n)),
                  value = c(df1[,1],df1[,2]))
  plot = ggplot(df, aes(x = element_of_sample, y = value, color = type)) +
    geom_point() + ggtitle(name)
  return(plot)
}


OLS_coefs = function(name, train) {
  f <- as.formula(paste(name, "~", paste(plasma_names, collapse="+")))
  mod = lm(f, data = train)
  coef = as.data.frame(summary(mod)$coef)
  names(coef) = paste0(names(coef),"_", name)
  return(coef[,-(2:3)])
}


#кросс-валидация 5 частей
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

#запускаю с кросс-валидацией
list_coefs = list()
for (i in 1:5) {
  test = df[[i]]
  train = do.call(rbind, df[-i])
  coeffs = lapply(IgG_names, OLS_coefs, train = train)
  coeffs = do.call(cbind,coeffs)
  
  res = sapply(IgG_names, OLS, train = train, test = test)
  res = as.data.frame(t(res))
  
  if (i == 1) {
    RES = res
    names(RES) = paste0(names(RES),"_set",i)
    plot_list = lapply(IgG_names, OLS_plot, train = train, test = test) }
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

#вывожу средние от статистик по 5 запускам
RES$rmse = apply(RES[,startsWith(names(RES), "rmse")],1,mean)
RES$rel_err = apply(RES[,startsWith(names(RES), "rel_err")],1,function(x) {mean(abs(x))})
RES$pearson = apply(RES[,startsWith(names(RES), "pearson")],1,mean)
RES$spearman = apply(RES[,startsWith(names(RES), "spearman")],1,mean)

write.xlsx(RES, "OLS/stats_IgGconcentration.xlsx", rowNames = T)


write.xlsx(list_coefs, "OLS/coeffs_IgGconcentration.xlsx", rowNames = T)
length(list_coefs)

pdf('OLS/plot_actual_predict_mq_IgGconcentration.pdf')
marrangeGrob(grobs = plot_list, ncol = 1, nrow = 2)
dev.off()
