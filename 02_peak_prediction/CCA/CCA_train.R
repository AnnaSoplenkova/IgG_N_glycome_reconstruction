#CCA + plot actual/predict on training set with 36 plasma peaks

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

# merge 2 datasets in one
names(IgG_train)[-(1:3)] = paste0(names(IgG_train),"_IgG")[-(1:3)]
names(plasma_train)[-(1:3)] = paste0(names(plasma_train),"_plasma")[-(1:3)]
train_set = left_join(plasma_train, IgG_train)

IgG_names = names(train_set)[endsWith(names(train_set), "IgG")]
plasma_names = names(train_set)[endsWith(names(train_set), "plasma")]






cca.predict = function(X, Y, X.new, k.x, k.y) {
  
  nXp = scale(t(X))
  nYq = scale(t(Y))
  
  nXp.new = sweep(t(X.new), 2, attr(nXp, 'scaled:center'), '-')
  nXp.new = sweep(nXp.new, 2, attr(nXp, 'scaled:scale'), '/')
  
  n = nrow(nXp)
  p = ncol(nXp)
  q = ncol(nYq)
  
  cor.eigen = function(X, scale = F) {
    # efficiently compute non-zero eigenvectors of the correlation matrix for X
    
    if(scale) { X = scale(X) }
    
    if(nrow(X) < ncol(X)) {
      X = X / sqrt(n-1)
      e = eigen(X %*% t(X))
      e$values[e$values<0]=0
      list(vectors = sweep(t(X) %*% e$vectors, 2, sqrt(e$values), FUN = "/"), 
           values = e$values)
    } else {
      eigen(t(X)%*%X/(n-1))
    }
  }
  
  pRp.e = cor.eigen(nXp)
  qRq.e = cor.eigen(nYq)
  
  pLp = diag(pRp.e$values)
  pEp = pRp.e$vectors
  qMq = diag(qRq.e$values)
  qFq = qRq.e$vectors
  
  nUpp = scale(nXp %*% pEp[,1:k.x])
  nVqq = scale(nYq %*% qFq[,1:k.y])
  
  nUpp.new = sweep(nXp.new %*% pEp[,1:k.x], 2, attr(nUpp, 'scaled:center'), '-')
  nUpp.new = sweep(nUpp.new, 2, attr(nUpp, 'scaled:scale'), '/')
  
  pRp = t(nUpp) %*% nUpp / (n-1)
  pRq = t(nUpp) %*% nVqq / (n-1)
  qRq = t(nVqq) %*% nVqq / (n-1)
  
  # Cook eq. 19
  B.e = eigen(solve(qRq) %*% t(pRq) %*% solve(pRp) %*% pRq)
  B = B.e$vectors
  B.e$values[B.e$values<0]=0
  qLambdaq = diag(B.e$values)
  
  # Glahn (1968) eq. 13
  A = solve(pRp) %*% pRq %*% B %*% diag(1/sqrt(B.e$values))
  
  V.hat = nUpp.new %*% A %*% qLambdaq %*% solve(B)
  
  # unscale and backtransform V.hat
  V.hat = sweep(V.hat, 2, attr(nVqq, 'scaled:scale'), "*")
  V.hat = sweep(V.hat, 2, attr(nVqq, 'scaled:center'), "+")
  Y.hat = V.hat %*% t(qFq[,1:k.y])
  
  t(Y.hat)
}

my.cca.predict = function(X, Y, X.new, k.x, k.y) {
  res = cca.predict(t(X), t(Y), t(X.new), k.x, k.y)
  return(t(res))
}


main_function = function(train, test) {
  
  plasma = as.matrix(train[, ..plasma_names])
  IgG = as.matrix(train[, ..IgG_names])
  
  test_plasma = as.matrix(test[, ..plasma_names])
  test_IgG = as.matrix(test[, ..IgG_names])
  
  pred_IgG = my.cca.predict(X = plasma, Y = IgG, X.new = test_plasma, 
                            k.x = ncol(plasma), k.y = ncol(IgG))
  
  
  metrics = function(i) {
    x = pred_IgG[,i]
    y = test_IgG[,i]
    return(c(rmse = rmse(x,y), rel_err = mean(abs(x-y)/x),
             pearson = cor(x,y), spearman = cor(x,y,method = "spearman")) )
  }
  
  N = ncol(pred_IgG)
  res_list = lapply(1:N, metrics)
  res = do.call(rbind, res_list)
  
  rownames(res) = IgG_names
  return(res)
  
}



#####################

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



for (i in 1:5) {
  test = df[[i]]
  train = do.call(rbind, df[-i])
  
  res = as.data.frame(main_function(train, test))
  
  if (i == 1) {
    RES = res
    names(RES) = paste0(names(RES),"_set",i) }
  if (i %in% 2:5) {
    names(res) = paste0(names(res),"_set",i)
    RES = cbind(RES, res)
  }
  
}

RES$rmse = apply(RES[,startsWith(names(RES), "rmse")],1,mean)
RES$rel_err = apply(RES[,startsWith(names(RES), "rel_err")],1,function(x) {mean(abs(x))})
RES$pearson = apply(RES[,startsWith(names(RES), "pearson")],1,mean)
RES$spearman = apply(RES[,startsWith(names(RES), "spearman")],1,mean)


write.xlsx(RES, "CCA/stats.xlsx", rowNames = T)

