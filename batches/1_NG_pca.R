## Imports and settings
library(factoextra)

library(ggplot2)
library(patchwork)
library(Cairo)

WORKDIR <- ''

setwd(WORKDIR)

## Function definitions
zip <- function(...){
  mapply(list, ..., SIMPLIFY = F)
}

enumerate <- function(...){
  zip(ix = seq_along(..1), ...)
}

sum.gglist <- function(gglist, na.rm = TRUE){
  Reduce('+', gglist)
}


## Data Import and initial processing
IgG_train <- read.csv('IgG_training_set_ta.csv')
IgG_train$st <- 'Train'
IgG_test <- read.csv('IgG_test_set_ta.csv')
IgG_test$st <- 'Test'
Ig <- rbind(IgG_train, IgG_test)
rownames(Ig) <- Ig[,1]; Ig <- Ig[,-1]

Pl_train <- read.csv('plasma_training_set_ta.csv')
Pl_train$st <- 'Train'
Pl_test <- read.csv('plasma_test_set_ta.csv')
Pl_test$st <- 'Test'
Pl <- rbind(Pl_train, Pl_test)
rownames(Pl) <- Pl[,1]; Pl <- Pl[,-1]

## PCA
# Plasma
pl_pca <- princomp(Pl[, -c(1:2, ncol(Pl))], cor = T, scores = T)
eig <- get_eigenvalue(pl_pca)$variance.percent
pl_plc_ind <- as.data.frame(get_pca_ind(pl_pca)$coord)
pl_plc_ind$status <- Pl$st

plots <- list()
class(plots) <- 'gglist'
for (pair in enumerate(combn(6, 2, simplify = F))){
  i <- pair[[1]]
  pair <- pair[[2]]
  df <- pl_plc_ind[, c(pair, ncol(pl_plc_ind))]
  names(df) <- c('PC1', 'PC2', 'Status')
  p <- ggplot(data = df, aes(x = PC1, y = PC2, colour = Status)) +
    stat_density2d(size = 1.2, alpha = 0.75) + theme_bw() + 
    theme(legend.key.height = unit(1, 'cm'),  legend.key.width = unit(1.25, 'cm'),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12)) +
    ggtitle(paste0('PC', pair[1], ' vs. PC', pair[2])) +
    xlab(paste0('PC', pair[1], ': ', round(eig[pair[1]],1), '% var explained')) + 
    ylab(paste0('PC', pair[2], ': ', round(eig[pair[2]],1), '% var explained'))
  plots[[i]] <- p
}

pp <- sum(plots) + fviz_eig(pl_pca, choice = "variance", ncp = 15, 
                            ggtheme = theme_bw()) +
  plot_layout(guides = 'collect',
              design = "
              ABC
              DEF
              GHI
              JKL
              MNO
              PP#
              " )

CairoPDF(width = 16, height = 25, 'Plasma_train_test.pdf')
pp
dev.off()

# IgG
Ig_pca <- princomp(Ig[, -c(1:2, ncol(Ig))], cor = T, scores = T)
eig <- get_eigenvalue(Ig_pca)$variance.percent
Ig_pca_ind <- as.data.frame(get_pca_ind(Ig_pca)$coord)
Ig_pca_ind$status <- Ig$st

plots <- list()
class(plots) <- 'gglist'
for (pair in enumerate(combn(6, 2, simplify = F))){
  i <- pair[[1]]
  pair <- pair[[2]]
  df <- Ig_pca_ind[, c(pair, ncol(Ig_pca_ind))]
  names(df) <- c('PC1', 'PC2', 'Status')
  p <- ggplot(data = df, aes(x = PC1, y = PC2, colour = Status)) +
    stat_density2d(size = 1.2, alpha = 0.75) + theme_bw() + 
    theme(legend.key.height = unit(1, 'cm'),  legend.key.width = unit(1.25, 'cm'),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12)) +
    ggtitle(paste0('PC', pair[1], ' vs. PC', pair[2])) +
    xlab(paste0('PC', pair[1], ': ', round(eig[pair[1]],1), '% var explained')) + 
    ylab(paste0('PC', pair[2], ': ', round(eig[pair[2]],1), '% var explained'))
  plots[[i]] <- p
}

pp <- sum(plots) + fviz_eig(Ig_pca, choice = "variance", ncp = 15, 
                            ggtheme = theme_bw()) +
  plot_layout(guides = 'collect',
              design = "
              ABC
              DEF
              GHI
              JKL
              MNO
              PP#
              " )

CairoPDF(width = 16, height = 25, 'IgG_train_test.pdf')
pp
dev.off()

# All
data <- cbind(Pl[, -c(1:2, ncol(Pl))], Ig[,-(1:2)])

data_pca <- princomp(data[, -ncol(data)], cor = T, scores = T)
eig <- get_eigenvalue(data_pca)$variance.percent
data_pca_ind <- as.data.frame(get_pca_ind(data_pca)$coord)
data_pca_ind$status <- data$st

plots <- list()
class(plots) <- 'gglist'
for (pair in enumerate(combn(6, 2, simplify = F))){
  i <- pair[[1]]
  pair <- pair[[2]]
  df <- data_pca_ind[, c(pair, ncol(data_pca_ind))]
  names(df) <- c('PC1', 'PC2', 'Status')
  p <- ggplot(data = df, aes(x = PC1, y = PC2, colour = Status)) +
    stat_density2d(size = 1.2, alpha = 0.75) + theme_bw() + 
    theme(legend.key.height = unit(1, 'cm'),  legend.key.width = unit(1.25, 'cm'),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12)) +
    ggtitle(paste0('PC', pair[1], ' vs. PC', pair[2])) +
    xlab(paste0('PC', pair[1], ': ', round(eig[pair[1]],1), '% var explained')) + 
    ylab(paste0('PC', pair[2], ': ', round(eig[pair[2]],1), '% var explained'))
  plots[[i]] <- p
}

pp <- sum(plots) + fviz_eig(Ig_pca, choice = "variance", ncp = 15, 
                            ggtheme = theme_bw()) +
  plot_layout(guides = 'collect',
              design = "
              ABC
              DEF
              GHI
              JKL
              MNO
              PP#
              " )

CairoPDF(width = 16, height = 25, 'All_train_test.pdf')
pp
dev.off()
