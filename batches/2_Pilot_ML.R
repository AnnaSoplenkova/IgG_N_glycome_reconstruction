## Imports and settings
library(ggplot2)
library(patchwork)
library(Cairo)
library(ggsci) #!
library(forestplot)

library(caret)
library(glmnet)
library(lars)

library(huxtable)

library(factoextra)


WORKDIR <- ''

setwd(WORKDIR)

### Definitions
## General purpose
zip <- function(...){
  mapply(list, ..., SIMPLIFY = F)
}


enumerate <- function(...){
  zip(ix = seq_along(..1), ...)
}


sum.gglist <- function(gglist, na.rm = TRUE){
  Reduce('+', gglist)
}


## CV routines

ModelCVRunner <- function(data, fold_col, method, tuneGrid = NULL, indep_vars, dep_vars){
  
  folds <- lapply(sort(unique(data[,fold_col])), function(i){
    match(rownames(data[(data[,fold_col] != i) & (data[,fold_col] > 0),]), 
          rownames(data))
  })
  names(folds) <- paste0('Fold', sort(unique(data[,fold_col])))
  fitControl <- trainControl(method = "cv", index = folds, number = length(unique(data[,fold_col])))
  
  ans <- lapply(dep_vars, function(y){
    train(data[,indep_vars], data[, y], method = method, trControl = fitControl,
          tuneGrid = tuneGrid)
  })
  names(ans) <- dep_vars
  class(ans) <- c('carettrainlist', 'list')
  return(ans)
}


EnetCVRunner <- function(data, fold_col, method, alpha = NULL, indep_vars, dep_vars){
  
  stopifnot(!is.null(alpha), max(alpha) > 1, min(alpha) < 0)
  
  X <- as.matrix(data[indep_vars])
  y <- data[, dep_vars[1]]
  ## collecting best lambdas for each alphas
  lapply(dep_vars, function(var_name){
    y <- data[,var_name]
    cvs <- lapply(alpha, function(a) cv.glmnet(X, y, alpha = a, 
                                               type.measure = "mse", 
                                               foldid = data[fold_col,]))
    lambdas <- do.call(rbind, lapply(cvs, function(cv){
      df <- cv$lambda
      df$alpha <- a
      return(df[, c(ncol(df)), ])
    }))
    for (a in alpha){
      cv <- 
        lambdas <- cv$lambda
      up <- cv$up
    }
  })
  for (a in alpha){
    cv <- cv.glmnet(X, )
  }
  cv <- cv.glmnet(X, y, alpha = 1, type.measure = c("mse", "mae"),
                  foldid = data[, fold_col])
  
}


lambda1sd <- function(all.mods){
  
  do.call(rbind, lapply(unique(all.mods$Trait), function(trait){
    dft <- all.mods[all.mods$Trait == trait,]
    best2.lambdas <- lapply(unique(dft$alpha), function(a){
      dfa <- dft[dft$alpha == a,]
      lam.rmse.sd <- sd(dfa$RMSE)
      lam.min <- max(dfa[dfa$RMSE == min(dfa$RMSE), 'lambda']) 
      lam.1sd <- max(dfa[dfa$RMSE < (min(dfa$RMSE) + lam.rmse.sd), 'lambda'])
      dfa <- dfa[dfa$lambda %in% c(lam.min, lam.1sd),]
      dfa$RMSE_sd_lambda <- lam.rmse.sd
      # Explicity
      return(dfa[, c(1:4, ncol(dfa), 5:(ncol(dfa)- 1))])
    })
    # Explicity
    return(best2.lambdas[[which.min(sapply(best2.lambdas, function(x) x[2, 'RMSE']))]])
  }))
}


CVReport.carettrainlist <- function(cv_list, print.best = NULL, 
                                    metric_cols = c(RMSE = 'RMSE', Rsquared ='Rsq', MAE = 'MAE', 
                                                    RMSESD ='RMSE_SD', RsquaredSD = 'Rsq_SD', MAESD = 'MAE_SD')){
  
  method <- cv_list[[1]]$method
  tune_idx <- sort(match(names(cv_list[[1]]$bestTune), names(cv_list[[1]]$results)))
  all.params <- do.call(rbind, lapply(names(cv_list), function(var_name){
    tab <- cv_list[[var_name]]$results
    names(tab)[-tune_idx] <- metric_cols[names(tab)[-tune_idx]] 
    tab$Trait <- var_name
    tab$Method <- method
    return(tab[, c(ncol(tab)-1, ncol(tab), 1:(ncol(tab)-2))])
  }))
  best.params <- do.call(rbind, lapply(names(cv_list), function(var_name){
    bt <- cv_list[[var_name]]$bestTune
    tab <- cv_list[[var_name]]$results
    tab <- tab[apply(tab[c(tune_idx)], 1, function(row){
      all(row == bt[1,])
    }), ]
    names(tab)[-tune_idx] <- metric_cols[names(tab)[-tune_idx]] 
    tab$Trait <- var_name
    tab$Method <- method
    return(tab[, c(ncol(tab)-1, ncol(tab), 1:(ncol(tab)-2))])
  }))
  
  if (!is.null(print.best)){
    if (print.best == TRUE){
      print(best.params)
    } else if (typeof(print.best) == "character"){
      print(best.params)
    }
  }
  
  return(list(all.params, best.params))
}


CVReport <- function(cv_list, print.best = NULL){
  UseMethod("CVReport")
}


## Formatters

TestModel <- function(params, model, Xtrain, Xtest, Ytrain, Ytest){
  res <- lapply(colnames(Ytrain), function(trait){
    df <- params[[trait]]
    fit <- model(Xtrain, Ytrain[trait], df)
    preds <- tryCatch(predict(fit, as.matrix(Xtest)), 
                      error = function(cond){
                        message(cond)
                        message("\nProccessing data as dataframe instead")
                        return(predict(fit, Xtest))
                      })
    # preds <- predict(fit, as.matrix(Xtest))
    df <- c(list(Trait = trait), df, 
            list(RMSE = sqrt(mean((preds - Ytest[,trait])^2)), 
                 Rsq = 1 - sum((preds - Ytest[,trait])^2) / sum((Ytest[,trait] - mean(Ytest[,trait]))^2),
                 MAE = mean(abs(preds - Ytest[,trait])), pearson = cor(preds, Ytest[,trait]),
                 spearman = cor(preds, Ytest[,trait], method = 'spearman')))
    
    return(list(data.frame(df), fit))
  })
  
  names(res) <- colnames(Ytrain)
  return(res)
}


params_fromtab <- function(tab, cols){
  ans <-lapply(rownames(tab), function(n){
    as.list(tab[n, cols])
  })
  names(ans) <- tab$Trait
  return(ans)
}


OutTable <- function(df, file = NULL, caption, number, num_hyp, method = 'pdf'){
  require('magrittr')
  df <- as_hux(df)
  df <- df %>%  set_background_color(evens, everywhere, "grey95") %>%
    insert_row(fill = "", after = 0) %>% set_contents(1, c(3, 3 + num_hyp), 
                                                      c('Hyper parameters', 'Quality Metrics')) %>%
    merge_cells(1, 3:(2+num_hyp)) %>% merge_cells(1, (3+num_hyp):ncol(df)) %>%
    set_align(1, everywhere, "center") %>% set_tb_padding(1, everywhere, 0) %>% 
    set_bold(1, everywhere) %>%
    set_caption(paste0("Table ", number, ". ", caption)) %>% theme_article()
  
  if (!is.null(file)){
    if (method == 'pdf')
      quick_pdf(df, file = file, open = F)
    else quick_pptx(df, file = file, open = F)
  } else 
    return(df)
}



## Evaluation on the out-test
TestOutSet <- function(X, Y, fits){
  do.call(rbind, lapply(names(Y), function(trait){
    preds <- tryCatch(predict(fits[[trait]], X),
                      error = function(cond){
                        warning(cond)
                        warning("\nUsing `as.matrix` on X")
                        return(predict(fits[[trait]], as.matrix(X)))
                      })  
    df <- data.frame(Trait = trait, RMSE = sqrt(mean((preds - Y[,trait])^2)),
                     Rsq = 1 - sum((preds - Y[,trait])^2) / sum((Y[,trait] - mean(Y[,trait]))^2),
                     MAE = mean(abs(preds - Y[,trait])), pearson = cor(preds, Y[,trait]),
                     spearman = cor(preds, Y[,trait], method = 'spearman'))
    return(df)
  }))
}

predMat <-  function(X, Y, fits){
  preds <- do.call(cbind, lapply(names(Y), function(trait){
    preds <- tryCatch(predict(fits[[trait]], X),
                      error = function(cond){
                        warning(cond)
                        warning("\nUsing `as.matrix` on X")
                        return(predict(fits[[trait]], as.matrix(X)))
                      })
  }))
  colnames(preds) <- colnames(Y)
  return(preds)
}


### Analysis
## Data Import and initial processing
# TwinsUK
IgG_train <- read.csv('IgG_training_set_ta.csv')
Pl_train <- read.csv('plasma_training_set_ta.csv')
names(IgG_train) <- paste0('i', names(IgG_train))
names(Pl_train)[4:length(Pl_train)] <- paste0('p', names(Pl_train)[4:length(Pl_train)]) 
all(Pl_train$ID == IgG_train$iID)# Little check
# [1] TRUE
twins <- cbind(Pl_train, IgG_train[, 4:length(IgG_train)])

# Dummy variable for splits: Sex-Age_quantile
quis <- quantile(twins$Age)
twins$Age_Sex <- ifelse(twins$Age < quis[2], 25, 
                        ifelse(twins$Age < quis[3], 50,
                               ifelse(twins$Age < quis[4], 75, 100)))

twins$Age_Sex <- paste(twins$Sex, twins$Age_Sex, sep = '_')

set.seed(1337)
trainIndex <- createDataPartition(twins$Age_Sex, p = 0.8, list = F)
twins$Fold <- -1
twins$Fold[trainIndex[,1]] <- 12345 
interIndex <- createDataPartition(twins[twins$Fold == -1, 'Age_Sex'], p = 0.5, list = F)
twins[twins$Fold == -1,][interIndex, 'Fold'] <- 0

# Now we have these partions: test - 10%, intermediate-test - 10%, train-cv - 80%
# Inner test and intermediate-test were merged

# Assigning folds
twins[twins$Fold > 0, 'Fold'] <- createFolds(twins[twins$Fold > 0, 'Age_Sex'], 5, list = F)

# Adding sex-age interaction
twins$Sex <- ifelse(twins$Sex == 'M', 1, 0)
twins[, 'Age:Sex'] <- twins$Sex * twins$Age
rownames(twins) <- twins$ID; twins <- twins[,-1]


# Cedar
Pl_test <- read.csv('plasma_test_set_ta.csv')
IgG_test <- read.csv('IgG_test_set_ta.csv')

names(IgG_test) <- paste0('i', names(IgG_test))
names(Pl_test)[4:length(Pl_test)] <- paste0('p', names(Pl_test)[4:length(Pl_test)]) 
all(Pl_test$ID == IgG_test$iID)
# [1] TRUE
cedar <- cbind(Pl_test, IgG_test[, 4:length(IgG_test)])

cedar$Sex <- ifelse(cedar$Sex == 'M', 1, 0)
cedar[, 'Age:Sex'] <- cedar$Sex * cedar$Age
rownames(cedar) <- cedar$ID; cedar <- cedar[,-1]


## CV-selection
indep_vars <- c('Sex', 'Age', 'Age:Sex', names(twins)[startsWith(names(twins), 'p')])
dep_vars <- c(names(twins)[startsWith(names(twins), 'i')])

# ordinary least squares
linreg <- ModelCVRunner(twins[twins$Fold > 0,], 'Fold', method = 'lm', tuneGrid = NULL,
                        indep_vars = indep_vars, dep_vars = dep_vars)

reg <- CVReport(linreg)
reg[[2]] <- tibble::add_column(reg[[2]], Intercept = T, .after = 2)

# Ridge
ridge <- ModelCVRunner(twins[twins$Fold > 0,], 'Fold', method = 'glmnet', 
                       tuneGrid = expand.grid(alpha = seq(0, 0.05, length.out = 10),
                                              lambda = exp(seq(-6, 2, length.out = 1000))),
                       indep_vars = indep_vars, dep_vars = dep_vars)
rig <- CVReport(ridge)
rig[[2]] <- lambda1sd(rig[[1]])

# Lasso
lasso <- ModelCVRunner(twins[twins$Fold > 0,], 'Fold', method = 'glmnet', 
                       tuneGrid = expand.grid(alpha = seq(0.95, 1, length.out = 10),
                                              lambda = exp(seq(-6, 2, length.out = 1000))),
                       indep_vars = indep_vars, dep_vars = dep_vars)
las <- CVReport(lasso)
las[[2]] <- lambda1sd(las[[1]])

# Enet
enet <- ModelCVRunner(twins[twins$Fold > 0,], 'Fold', method = 'glmnet', 
                      tuneGrid = expand.grid(alpha = seq(0.5, 0.95, length.out = 46),
                                             lambda = exp(seq(-6, 2, length.out = 1000))),
                      indep_vars = indep_vars, dep_vars = dep_vars)
en <- CVReport(enet)
en[[2]] <- lambda1sd(en[[1]])

# Least angle

lars1 <- ModelCVRunner(twins[twins$Fold > 0,], 'Fold', method = 'lars2', 
                       tuneGrid = expand.grid(step = seq(length(indep_vars), 1, by = -1)),
                       indep_vars = indep_vars, dep_vars = dep_vars)

lars_ <- CVReport(lars1)


## Choosing the champion on the inner test

# Overused expressions
Xtrain <- twins[twins$Fold > 0, indep_vars]
Xtest <- twins[twins$Fold < 1, indep_vars]
Ytrain <- twins[twins$Fold > 0, dep_vars]
Ytest <- twins[twins$Fold < 1, dep_vars]

# Inner Test Performance
lreg <- TestModel(params_fromtab(reg[[2]], 'Intercept'), 
                  function(x, y, params){
                    yn <- names(y)
                    lm(formula(paste(yn, '~ .')), cbind(x, y))
                    }, Xtrain, Xtest, Ytrain, Ytest)

rid_1sd <- TestModel(params_fromtab(rig[[2]][c(F,T),], c('alpha', 'lambda')), 
                     function(x, y, params){
                       glmnet(as.matrix(x), y[,1], alpha = params[['alpha']], 
                              lambda = params[['lambda']], family = "gaussian", 
                              nlambda = 0)
                       }, Xtrain, Xtest, Ytrain, Ytest)

rid_min <- TestModel(params_fromtab(rig[[2]][c(T,F),], c('alpha', 'lambda')), 
                     function(x, y, params){
                       glmnet(as.matrix(x), y[,1], alpha = params[['alpha']], 
                              lambda = params[['lambda']], family = "gaussian", 
                              nlambda = 0)
                       }, Xtrain, Xtest, Ytrain, Ytest)

las_1sd <- TestModel(params_fromtab(las[[2]][c(F,T),], c('alpha', 'lambda')), 
                     function(x, y, params){
                       glmnet(as.matrix(x), y[,1], alpha = params[['alpha']], 
                              lambda = params[['lambda']], family = "gaussian",
                              nlambda = 0)
                       }, Xtrain, Xtest, Ytrain, Ytest)

las_min <- TestModel(params_fromtab(las[[2]][c(T,F),], c('alpha', 'lambda')),
                     function(x, y, params){
                       glmnet(as.matrix(x), y[,1], alpha = params[['alpha']],
                              lambda = params[['lambda']], family = "gaussian",
                              nlambda = 0)}, Xtrain, Xtest, Ytrain, Ytest)

en_1sd <- TestModel(params_fromtab(en[[2]][c(F,T),], c('alpha', 'lambda')), 
                    function(x, y, params){
                      glmnet(as.matrix(x), y[,1], alpha = params[['alpha']], 
                             lambda = params[['lambda']], family = "gaussian", 
                             nlambda = 0)
                      }, Xtrain, Xtest, Ytrain, Ytest)

en_min <- TestModel(params_fromtab(las[[2]][c(T,F),], c('alpha', 'lambda')), 
                    function(x, y, params){
                      glmnet(as.matrix(x), y[,1], alpha = params[['alpha']],
                             lambda = params[['lambda']], family = "gaussian",
                             nlambda = 0)}, Xtrain, Xtest, Ytrain, Ytest)

all <- list(lm = lreg, ridge.1sd = rid_1sd, ridge.min = rid_min, lasso.1sd = las_1sd, 
            lasso.min = las_min, enet.1sd = en_1sd, enet.min = en_min) # Redefinition of builtins

best_of <- sapply(1:length(all[[1]]), function(i){
  which.min(sapply(all, function(x) x[[i]][[1]][1, 'RMSE']))
})

best.mods <- lapply(1:length(best_of), function(i){
  df <- all[[best_of[i]]][[i]][[1]][, c('Trait', 'RMSE', 'Rsq', 'MAE', 'pearson', 'spearman')]
  df$Method <- names(all)[best_of[i]]
  fit <- all[[best_of[i]]][[i]][[2]]
  return(list(df[, c(1, ncol(df), 2:(ncol(df)-1))], fit))
})

best.mods.tab <- do.call(rbind, lapply(best.mods, function(x) x[[1]]))
best.mods.fit <- lapply(best.mods, function(x) x[[2]])
names(best.mods.fit) <- best.mods.tab$Trait


## Final evaluation

fin_perf <- TestOutSet(cedar[,indep_vars], cedar[,dep_vars], best.mods.fit)
preds <- predMat(cedar[,indep_vars], cedar[,dep_vars], best.mods.fit)
fdf <- data.frame(cedar = unlist(cedar[, dep_vars]), preds = preds[1:length(preds)],
                  sex = rep(cedar$Sex, length(dep_vars)))

### Tables and plots

## Scatterplot GTrue vs Predicted with sexwise colorfill
scs <- ggplot(fdf, aes(x = cedar, y = preds, fill = as.factor(sex))) + 
  geom_point(shape = 21, size = 2.7, alpha = 0.65) + 
  geom_abline(intercept = 0, slope = 1, color = "#ffeda0", linetype = "dashed", size = 1.6) +
  scale_fill_manual(labels = c("F", "M"), values = c("tomato","blue")) +
  theme_bw() + labs(title = "iGP Predicted vs. Ground Truth", 
                    subtitle = paste("corr", round(cor(fdf$preds, fdf$cedar), 4), sep = " = "),
                    x = "CEDAR", y = "Predicted", fill = "Sex") +
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12))


## Scatterplot with peakwise clusters
# All clusters
res.km <- kmeans(fdf[, -3], 6, nstart = 25)
res.km2 <- res.km
res.km2$cluster <- factor(rep(dep_vars, each = nrow(cedar)),
                          levels = dep_vars)


# All peaks in one plot
all_clust <- fviz_cluster(res.km2, data = fdf[, -3], # It is faster and cleaner
                                                     # To just use factoextra's
                                                     # function to draw poligons
                                                     # then calculating and drawin them manually
                                                     # The latter is more preferable, though,
                                                     # Since the picture is prettier that way
                                                     # Here, for the sake of simplicity
                                                     # And readability we use shorter way
                          geom = 'point', eliipse.type = "convex", stand = F, ggtheme = theme_bw()) + 
  labs(title = "iGP Predicted vs. Ground Truth", subtitle = paste("Each iGP highlighted"),
       x = "CEDAR", y = "Predicted", colour = "iGP", shape = "iGP", fill = "iGP") +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", size = 1.6) +
  geom_abline(intercept = sd(rot[,2])/2, slope = 1, color = "#005a32", linetype = "dotted", size = 1.3) +
  geom_abline(intercept = -sd(rot[,2])/2, slope = 1, color = "#005a32", linetype = "dotted", size = 1.3) +
  geom_abline(intercept = sd(rot[,2]), slope = 1, color = "#b2182b", linetype = "dotted", size = 1.3) +
  geom_abline(intercept = -sd(rot[,2]), slope = 1, color = "#b2182b", linetype = "dotted", size = 1.3) +
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 12), 
        plot.title = element_text(size = 14, face = "bold"))


# Rotating coordinates by pi/4 counterclockwise
rot <- as.data.frame(matrix(
  apply(fdf[,1:2], 1, function(vec) {
    (vec %*% matrix(c(cos(pi/4), -sin(pi/4), sin(pi/4), cos(pi/4)), 
                    nrow=2, byrow = T))[1,]
  }), ncol = 2, byrow = T
  ))

rot$Trait <- res.km2$cluster
out_perc <- tapply(rot$V2, rot$Trait, function(group){
  sum(abs(group) > sd(rot[,2]))/length(group)
})

# We need these to draw asteriks in the boxplot
annotations <- data.frame(V2 = -3, Trait = factor(levels(rot$Trait), levels(rot$Trait)),
                          lab = c('', '', '', '**', '**', '', '', '', '', '', '', 
                                  '', '', '', '*', '', '*', '', '**', '**', '', '*'))

# Rotated boxplots themselves

boxrot <- ggplot(rot, aes(x = -V2, y = Trait, colour = Trait)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(x = -V2, y = Trait, colour = Trait), width = 0.2, height = 0.2, alpha = 0.15) +
  geom_vline(xintercept = sd(rot[,2])/2, color = "#005a32", linetype = "dotted", size = 1.3) +
  geom_vline(xintercept = -sd(rot[,2])/2, color = "#005a32", linetype = "dotted", size = 1.3) +
  geom_vline(xintercept = sd(rot[,2]), color = "#b2182b", linetype = "dotted", size = 1.3) +
  geom_vline(xintercept = -sd(rot[,2]), color = "#b2182b", linetype = "dotted", size = 1.3) +
  geom_text(data = annotations, aes(x = V2, y = Trait, label = lab), color = '#b2182b',
            size = 8, fontface = "bold") +
  theme_bw() + labs(title = "Deviation from the (1,1) diagonal", y = 'Trait',
                    x = 'Deviation') + 
  scale_x_continuous(breaks = c(-3, -2, -1, -sd(rot[,2]), -sd(rot[,2])/2, 0, sd(rot[,2])/2,
                                sd(rot[,2]), 1, 2), 
                     labels = c('3', '2', '1', '1SD', '0.5SD', "0", "-0.5SD", '-1SD', '-1', '-2')) +
  theme(legend.position = "none", axis.text.x = element_text(size = 10), 
        axis.title.x = element_text(size = 12), axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 12), plot.title = element_text(size = 14, face = "bold"))

# Scatter for good peaks
res.kmgood <- res.km2
res.kmgood$cluster <- res.kmgood$cluster[!(res.kmgood$cluster %in% c('iGP22', 'iGP20', 'iGP19',
                                                                     'iGP17', 'iGP15', 'iGP5', 'iGP4'))]
res.kmgood$cluster <- droplevels(res.kmgood$cluster)
good_clust <- fviz_cluster(res.kmgood, data = fdf[!(res.km2$cluster %in% c(
  'iGP22', 'iGP20', 'iGP19', 'iGP17', 'iGP15', 'iGP5', 'iGP4')), -3],
  palette = pal_material('green', n = 20, rev = T)(17), stand = F,
  geom = 'point', eliipse.type = "convex", ggtheme = theme_bw()) + 
  labs(title = "iGP Predicted vs. Ground Truth", subtitle = paste('"Good" iGPs highlighted'),
       x = "CEDAR", y = "Predicted", colour = "iGP", shape = "iGP", fill = "iGP") +
  xlim(layer_scales(all_clust)$x$range$range[1], layer_scales(all_clust)$x$range$range[2]) + 
  ylim(layer_scales(all_clust)$y$range$range[1], layer_scales(all_clust)$y$range$range[2]) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", size = 1.6) +
  geom_abline(intercept = sd(rot[,2])/2, slope = 1, color = "#005a32", linetype = "dotted", size = 1.3) +
  geom_abline(intercept = -sd(rot[,2])/2, slope = 1, color = "#005a32", linetype = "dotted", size = 1.3) +
  geom_abline(intercept = sd(rot[,2]), slope = 1, color = "#b2182b", linetype = "dotted", size = 1.3) +
  geom_abline(intercept = -sd(rot[,2]), slope = 1, color = "#b2182b", linetype = "dotted", size = 1.3) +
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 12), 
        plot.title = element_text(size = 14, face = "bold"))

# Scatter for med peaks
res.kmmed <- res.km2
res.kmmed$cluster <- res.kmmed$cluster[res.kmmed$cluster %in% c('iGP22', 'iGP17', 'iGP15')]
res.kmmed$cluster <- droplevels(res.kmmed$cluster)
med_clust <- fviz_cluster(res.kmmed, data = fdf[res.km2$cluster %in% c('iGP22', 'iGP17', 'iGP15'), -3],
                          palette = pal_material("yellow", n = 7, rev = T)(7)[c(1,2,5)],
                          geom = 'point', eliipse.type = "convex",stand = F,  ggtheme = theme_bw()) + 
  labs(title = "iGP Predicted vs. Ground Truth", subtitle = paste('"Concerning" iGPs highlighted'),
       x = "CEDAR", y = "Predicted", colour = "iGP", shape = "iGP", fill = "iGP") +
  xlim(layer_scales(all_clust)$x$range$range[1], layer_scales(all_clust)$x$range$range[2]) + 
  ylim(layer_scales(all_clust)$y$range$range[1], layer_scales(all_clust)$y$range$range[2]) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", size = 1.6) +
  geom_abline(intercept = sd(rot[,2])/2, slope = 1, color = "#005a32", linetype = "dotted", size = 1.3) +
  geom_abline(intercept = -sd(rot[,2])/2, slope = 1, color = "#005a32", linetype = "dotted", size = 1.3) +
  geom_abline(intercept = sd(rot[,2]), slope = 1, color = "#b2182b", linetype = "dotted", size = 1.3) +
  geom_abline(intercept = -sd(rot[,2]), slope = 1, color = "#b2182b", linetype = "dotted", size = 1.3) +
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 12), 
        plot.title = element_text(size = 14, face = "bold"))

# Scatter for bad peaks
res.kmbad <- res.km2
res.kmbad$cluster <- res.kmbad$cluster[res.kmbad$cluster %in% c('iGP20', 'iGP19', 'iGP5', 'iGP4')]
res.kmbad$cluster <- droplevels(res.kmbad$cluster)

bad_clust <- fviz_cluster(res.kmbad, data = fdf[res.km2$cluster %in% c('iGP20', 'iGP19','iGP5', 'iGP4'), -3],
                          palette = pal_material("pink", n =9, rev = T)(9)[c(1,3,6,8)],
                          geom = 'point', eliipse.type = "convex", stand = F, ggtheme = theme_bw()) + 
  labs(title = "iGP Predicted vs. Ground Truth", subtitle = paste('"Bad" iGPs highlighted'),
       x = "CEDAR", y = "Predicted", colour = "iGP", shape = "iGP", fill = "iGP") +
  xlim(layer_scales(all_clust)$x$range$range[1], layer_scales(all_clust)$x$range$range[2]) + 
  ylim(layer_scales(all_clust)$y$range$range[1], layer_scales(all_clust)$y$range$range[2]) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", size = 1.6) +
  geom_abline(intercept = sd(rot[,2])/2, slope = 1, color = "#005a32", linetype = "dotted", size = 1.3) +
  geom_abline(intercept = -sd(rot[,2])/2, slope = 1, color = "#005a32", linetype = "dotted", size = 1.3) +
  geom_abline(intercept = sd(rot[,2]), slope = 1, color = "#b2182b", linetype = "dotted", size = 1.3) +
  geom_abline(intercept = -sd(rot[,2]), slope = 1, color = "#b2182b", linetype = "dotted", size = 1.3) +
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 12), 
        plot.title = element_text(size = 14, face = "bold"))

## Forest plot 
# t-test for difference for between-cohort means
ts <- do.call(rbind, lapply(dep_vars, function(v){
  tt <- t.test(twins[,v], cedar[,v], conf.level = 1 - 0.002272727) # Bonferroni-adjusted
  data.frame(mean = tt$estimate[1] - tt$estimate[2], lower = tt$conf.int[1], 
             upper = tt$conf.int[2], label = v)
  }))

bads <- c('iGP20', 'iGP19', 'iGP5', 'iGP4')
meds <- c('iGP22', 'iGP17', 'iGP15')
ts$Qual <- ifelse(ts$label %in% bads, "#870D4EFF", 
                  ifelse(ts$label %in% meds, "#F47F17FF", "#3A923EFF"))

# Drawing to a file straightaway
Cairo(1800, 1200, file = 'forest.png', dpi = 200)
ts %>% forestplot(mean = 'mean', lower = 'lower', upper = 'upper',
                  labeltext = 'label', clip = c(-0.5, 0.5), 
                  xlab = 'Difference of Means [log(ta)]', xticks = seq(-0.5, 0.5, 0.1),
                  boxsize = 0.5, shapes_gp = fpShapesGp(box = ts$Qual %>%
                                                          lapply(function(x) gpar(fill = x, col = "#555555")),
                                                        lines = ts$Qual %>%
                                                          lapply(function(x) gpar(col  = x)))) %>%
  fp_set_zebra_style("#F5F9F9")
dev.off()

## Boxplots for peaks in two cohorts
boxes <- lapply(dep_vars, function(v){
  df <- data.frame(value = c(twins[,v], cedar[,v]), lab = c(rep('Twins', nrow(twins)),
                                                            rep('CEDAR', nrow(cedar))))
  p <- ggplot(df, aes(y = value, x = lab, color = lab)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(y = value, x = lab, color = lab), width = 0.1, height = 0.01, alpha = 0.05) +
    theme_bw() + labs(title = v, y = 'log(ta)', x = '') + 
    scale_color_manual(values = c("#00468BFF", "#ED0000FF")) +
    theme(legend.position = "none", axis.text.x = element_text(size = 10), 
          axis.title.x = element_blank(), axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12), 
          plot.title = element_text(size = 12, face = "bold"))
  return(p)
})


class(boxes) <- 'gglist'
names(boxes) <- dep_vars
bboxes_toviz <- boxes[c(rev(bads), rev(meds))]
class(bboxes_toviz) <- 'gglist'
sum(bboxes_toviz) + plot_layout(ncol = 4)

gboxes_toviz <- boxes[!names(boxes) %in% c(rev(bads), rev(meds))]
class(gboxes_toviz) <- 'gglist'
sum(gboxes_toviz) + plot_layout(ncol = 4)

## Scatterplot GTruth vs Predicted colored by peak quality

fdf$Status <- ifelse(res.km2$cluster %in% bads, 'Bad', ifelse(res.km2$cluster %in% meds, 'Medium', 'Good'))
statussc <- ggplot(fdf, aes(x = cedar, y = preds, fill = as.factor(Status))) + 
  geom_point(shape = 21, size = 2.7, alpha = 0.65) + 
  geom_abline(intercept = 0, slope = 1, color = "#ffeda0", linetype = "dashed", size = 1.6) +
  # geom_smooth(aes(cedar, preds), method = "loess", se = F, fullrange = T) 
  scale_fill_manual(labels = c("Bad", "Good", 'Medium'), 
                    values = c("#870D4EFF","#3A923EFF", "#F47F17FF")) +
  theme_bw() + labs(title = "iGP Predicted vs. Ground Truth", 
                    subtitle = paste("corr", round(cor(fdf$preds, fdf$cedar), 4), sep = " = "),
                    x = "CEDAR", y = "Predicted", fill = "Quality") +
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12))

## Saving plots

ggsave(statussc, file = 'status_sc.png', width = 1600, height = 1200, units ='px', dpi = 180)

ggsave(sum(gboxes_toviz) + plot_layout(ncol = 4), file = 'goodboxes.png',
       width = 2000, height = 1600, units ='px', dpi = 180)

ggsave(sum(bboxes_toviz) + plot_layout(ncol = 2, byrow = F), file = 'badboxes.png',
       width = 1000, height = 1600, units ='px', dpi = 180)

ggsave('sex_scatter.png', scs, width = 1600, height = 1200,
       units ='px', dpi = 180)

ggsave('alliGP_scatter.png', all_clust, width = 13, height = 7.5,
       units ='in', dpi = 220)

ggsave('boxrot.png', boxrot, width = 13, height = 7.5,
       units ='in', dpi = 220)

ggsave('good_scatter.png', good_clust, width = 1600, height = 1200,
       units ='px', dpi = 180)

ggsave('med_scatter.png', med_clust, width = 1600, height = 1200,
       units ='px', dpi = 180)

ggsave('bad_scatter.png', bad_clust, width = 1600, height = 1200,
       units ='px', dpi = 180)

## Saving formatted tables

rig2 <- rig[[2]][c(T, F),]
rig2[,c(4, 6:11)] <- round(rig2[,c(4, 6:11)], 4)
rig2[,c(4, 6:11)] <- matrix(paste(as.matrix(rig2[,c(4,6:11)]), round(as.matrix(rig[[2]][c(F,T), c(4,6:11)]),4), sep = '/'),
                            ncol = length(c(4, 6:11)))
las2 <- las[[2]][c(T, F),]
las2[,c(4, 6:11)] <- round(las2[,c(4, 6:11)], 4)
las2[,c(4, 6:11)] <- matrix(paste(as.matrix(las2[,c(4,6:11)]), round(as.matrix(las[[2]][c(F,T), c(4,6:11)]),4), sep = '/'),
                            ncol = length(c(4, 6:11)))
enet2 <- en[[2]][c(T, F),]
enet2[,c(4, 6:11)] <- round(enet2[,c(4, 6:11)], 4)
enet2[,c(4, 6:11)] <- matrix(paste(as.matrix(enet2[,c(4,6:11)]), round(as.matrix(en[[2]][c(F,T), c(4,6:11)]),4), sep = '/'),
                             ncol = length(c(4, 6:11)))

quick_xlsx(OutTable(reg[[2]], caption = 'CV linear regression', number = 1, num_hyp = 1), 
           OutTable(rig2, caption = 'CV Ridge regression, [min lambda/1sd lambda]', number = 2, num_hyp = 3),
           OutTable(las2, caption = 'CV Lasso regression, [min lambda/1sd lambda]', number = 3, num_hyp = 3),
           OutTable(enet2, caption = 'CV Elastic Net [min lambda/1sd lambda]', number = 4, num_hyp = 3),
           OutTable(lars_[[2]], caption = 'CV Least Angle Regression', number = 5, num_hyp = 1),
           theme_article(set_caption(set_background_color(as_hux(best.mods.tab), evens, 
                                                          everywhere, "grey95"), 
                                     "Table 6. Best models evaluation (in-test)")), 
           theme_article(set_caption(set_background_color(as_hux(fin_perf), evens, 
                                                          everywhere, "grey95"), 
                                     "Table 7. Best models evaluation (out-test)")),
           file = 'tables.xlsx', open = F)



