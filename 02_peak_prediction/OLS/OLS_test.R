#Plot actual/predict + stats on test set

library(data.table)
library(dplyr)
library(Metrics)
library(ggplot2)
library(gridExtra)
library(tibble)

setwd("/mnt/polyomica/projects/GPD/peak_prediction")

# беру тестовые сеты
#plasma_test = fread("data/plasma_test_set_36peaks_2023.csv")
#IgG_test = fread("data/IgG_test_set_2023.csv")
plasma_test = fread("data/plasma_test_set_36peaks_2023_sex_age_merge45.csv")
IgG_test = fread("data/IgG_test_set_2023_sex_age_merge45.csv")
#plasma_test = fread("data/plasma_test_set_36peaks_2023_sex_age.csv")
#IgG_test = fread("data/IgG_test_set_2023_sex_age.csv")

# беру коэффициенты
coeffs = openxlsx::read.xlsx("OLS/coeffs.xlsx", sheet = 6)[-1]
rownames(coeffs) = coeffs$rowname
coeffs = coeffs[-1]



# merge 2 datasets in one
names(IgG_test)[-(1:3)] = paste0(names(IgG_test),"_IgG")[-(1:3)]
names(plasma_test)[-(1:3)] = paste0(names(plasma_test),"_plasma")[-(1:3)]




test_set = left_join(plasma_test, IgG_test)

#OLS
IgG_names = names(test_set)[endsWith(names(test_set), "IgG")]
plasma_names = names(test_set)[endsWith(names(test_set), "plasma")]


#пытаюсь сама че то написать
predict_IgG = data.frame(rep(NA, nrow(test_set)))

for (i in 1:length(IgG_names)) {
predictions <- as.matrix(as.matrix(cbind(1,test_set[,..plasma_names]) )) %*% t(coeffs)[i,]
predict_IgG[,i] = predictions
}

names(predict_IgG) = IgG_names

i=1
test_set[,plasma_names[i]]

OLS_plot = function(name, x, y) {
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

# function for stats OLS
OLS_stats = function(name, x, y) {
  return(c(rmse = rmse(x,y), rel_err = mean(abs(x-y)/x),
           pearson = cor(x,y), spearman = cor(x,y,method = "spearman")) )
} 




plots = lapply(1:length(IgG_names), function(i) {OLS_plot(IgG_names[i],
                                                          x = as.data.frame(test_set[,..IgG_names])[,i],
                                                          y = predict_IgG[,i])})

stats_list = lapply(1:length(IgG_names), function(i) {OLS_stats(IgG_names[i],
                                                          x = as.data.frame(test_set[,..IgG_names])[,i],
                                                          y = predict_IgG[,i])})

stats = do.call(rbind, stats_list)
rownames(stats) = IgG_names
openxlsx::write.xlsx(as.data.frame(stats), "OLS/stats_test.xlsx", rowNames = T)

pdf('OLS/plot_actual_predict_mq_test.pdf')
marrangeGrob(grobs = plots, ncol = 1, nrow = 2)
dev.off()
