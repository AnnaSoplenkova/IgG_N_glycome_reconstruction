# stats on test set ELASTIC_NET

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
coeffs = openxlsx::read.xlsx("ELASTIC_NET/coeffs.xlsx", sheet = 6)[-1]
rownames(coeffs) = coeffs$rowname
coeffs = coeffs[-1]



# merge 2 datasets in one
names(IgG_test)[-(1:3)] = paste0(names(IgG_test),"_IgG")[-(1:3)]
names(plasma_test)[-(1:3)] = paste0(names(plasma_test),"_plasma")[-(1:3)]




test_set = left_join(plasma_test, IgG_test)

#ELASTIC_NET
IgG_names = names(test_set)[endsWith(names(test_set), "IgG")]
plasma_names = names(test_set)[endsWith(names(test_set), "plasma")]


#пытаюсь сама че то написать
predict_IgG = data.frame(rep(NA, nrow(test_set)))



for (i in 1:length(IgG_names)) {
  predictions <- as.matrix(as.matrix(cbind(1,test_set[,..plasma_names]) )) %*% t(coeffs)[i,][-1]
  predict_IgG[,i] = predictions
}

names(predict_IgG) = IgG_names

i=1
test_set[,plasma_names[i]]


# function for stats 
stats = function(name, x, y) {
  return(c(rmse = rmse(x,y), rel_err = mean(abs(x-y)/x),
           pearson = cor(x,y), spearman = cor(x,y,method = "spearman")) )
} 



stats_list = lapply(1:length(IgG_names), function(i) {stats(IgG_names[i],
                                                            x = as.data.frame(test_set[,..IgG_names])[,i],
                                                            y = predict_IgG[,i])})

stats = do.call(rbind, stats_list)
rownames(stats) = IgG_names
openxlsx::write.xlsx(as.data.frame(stats), "ELASTIC_NET/stats_test.xlsx", rowNames = T)
