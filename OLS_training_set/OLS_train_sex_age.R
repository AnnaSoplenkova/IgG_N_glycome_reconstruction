#OLS with sex and age on training set

library(openxlsx)
library(dplyr)
library(Metrics)

setwd("/mnt/polyomica/projects/GPD/training_set/")

# беру тренировочные сеты в зависимости от нормализации
plasma_train = read.csv("plasma_training_set_ta.csv")
IgG_train = read.csv("IgG_training_set_ta.csv")
#plasma_train = read.csv("plasma_training_set.csv")
#IgG_train = read.csv("IgG_training_set.csv")

# merge 2 datasets in one
names(IgG_train)[-(1:3)] = paste0(names(IgG_train),"_IgG")[-(1:3)]
names(plasma_train)[-(1:3)] = paste0(names(plasma_train),"_plasma")[-(1:3)]
train_set = left_join(plasma_train, IgG_train)


# correction to sex and age
corrected_to_gender_age <- train_set[, -(1:3)]
train_set$Sex = ifelse(train_set$Sex == "F", 0, 1)
for(i in 1:ncol(corrected_to_gender_age)){
  corrected_to_gender_age[, i] <- residuals(lm(train_set[,i+3] ~ train_set$Sex + train_set$Age, 
                                               na.action = na.exclude))
}

names(corrected_to_gender_age)

#create test and train samples (20% and 80%)
dt = sort(sample(nrow(corrected_to_gender_age), nrow(corrected_to_gender_age)*.8))
train<-corrected_to_gender_age[dt,]
test<-corrected_to_gender_age[-dt,]
train_set = corrected_to_gender_age

# info about first IgG/Plasma column and N columns 
Start_plasma = 4
Start_IgG = 43
N_plasma_GP = 39
N_IgG_GP = 22

#OLS
IgG_names = names(train)[1:22+39]

# function for OLS
OLS = function(name) {
  f1 <- as.formula(paste(name, "~", paste(names(train)[c(1:39)], collapse="+")))
  mod = lm(f1, data = train)
  x = test[,names(test) == name]
  y = predict(mod,newdata=test)
  
  return(c(rmse = rmse(x,y), rel_err = mean(abs(x-y)/x),
           R2 = cor(x,y)^2, spearman = cor(x,y,method = "spearman")) )
} 

#draw OLS plot (actual/predict)
OLS_plot = function(name) {
  f1 <- as.formula(paste(name, "~", paste(names(train)[c(1:39)], collapse="+")))
  mod = lm(f1, data = train)
  x = test[,names(test) == name]
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

#get coefs from OLS
OLS_coefs = function(name) {
  f1 <- as.formula(paste(name, "~", paste(names(train)[c(1:39)], collapse="+")))
  mod = lm(f1, data = train)
  coef = as.data.frame(summary(mod)$coef)
  names(coef) = paste0(names(coef),"_", name)
  return(coef[,-(2:3)])
}

coeffs = lapply(IgG_names, OLS_coefs)
coeffs = do.call(cbind,coeffs)
res = sapply(IgG_names, OLS)
res = t(res)


#кросс-валидация 5 частей
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
  coeffs = lapply(IgG_names, OLS_coefs)
  coeffs = do.call(cbind,coeffs)
  res = sapply(IgG_names, OLS)
  res = as.data.frame(t(res))

  
  if (i == 1) {
    RES = res
    names(RES) = paste0(names(RES),"_set",i) }
  if (i %in% 2:5) {
    names(res) = paste0(names(res),"_set",i)
    RES = cbind(RES, res)
    
    plot_list = lapply(IgG_names, OLS_plot)
  }
  
  list_coefs[[i]] = coeffs
}


#вывожу средние от статистик по 5 запускам
names(RES)
rmse = c(1,5,9,13,17)
RES$rmse = apply(RES[,rmse],1,mean)
RES$rel_err = apply(RES[,rmse+1],1,function(x) {mean(abs(x))})
RES$R2 = apply(RES[,rmse+2],1,mean)
RES$spearman = apply(RES[,rmse+3],1,mean)

write.xlsx(RES, "cross_validation/ta/stats_all_sets_sex_age.xlsx",overwrite = T, rowNames = T)
write.xlsx(list_coefs, "cross_validation/ta/coeffs_all_sets_sex_age.xlsx",overwrite = T, rowNames = T)


pdf('plot_actual_predict_sex_age_ta.pdf')
marrangeGrob(grobs = plot_list, ncol = 1, nrow = 2)
dev.off()
