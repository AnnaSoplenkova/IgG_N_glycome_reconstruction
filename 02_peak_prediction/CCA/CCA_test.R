#CCA + plot actual/predict on test set with 36 plasma peaks

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
#plasma_train = fread("data/plasma_training_set_36peaks_2023_sex_age.csv")
#IgG_train = fread("data/IgG_training_set_2023_sex_age.csv")
plasma_train = fread("data/plasma_training_set_36peaks_2023_sex_age_merge45.csv")
IgG_train = fread("data/IgG_training_set_2023_sex_age_merge45.csv")


# беру тестовые сеты
#plasma_test = fread("data/plasma_test_set_36peaks_2023.csv")
#IgG_test = fread("data/IgG_test_set_2023.csv")
#plasma_test = fread("data/plasma_test_set_36peaks_2023_sex_age.csv")
#IgG_test = fread("data/IgG_test_set_2023_sex_age.csv")
plasma_test = fread("data/plasma_test_set_36peaks_2023_sex_age_merge45.csv")
IgG_test = fread("data/IgG_test_set_2023_sex_age_merge45.csv")

# merge 2 datasets in one
names(IgG_train)[-(1:3)] = paste0(names(IgG_train),"_IgG")[-(1:3)]
names(plasma_train)[-(1:3)] = paste0(names(plasma_train),"_plasma")[-(1:3)]
train_set = left_join(plasma_train, IgG_train)


# merge 2 datasets in one
names(IgG_test)[-(1:3)] = paste0(names(IgG_test),"_IgG")[-(1:3)]
names(plasma_test)[-(1:3)] = paste0(names(plasma_test),"_plasma")[-(1:3)]
test_set = left_join(plasma_test, IgG_test)

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


res = as.data.frame(main_function(train_set, test_set))

write.xlsx(res, "CCA/stats_test.xlsx", rowNames = T)
