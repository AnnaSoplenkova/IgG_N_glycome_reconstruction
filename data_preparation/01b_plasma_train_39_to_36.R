#переводит пики плазмы в тренировочной выборке из 39 пиков в 36

library(data.table)
library(tidyr)
library(dplyr)

mydata <- fread('/mnt/polyomica/projects/GPD/training_set/TwinsLongPlasma.csv')
file_final = "/mnt/polyomica/projects/GPD/training_set/TwinsPlasma_36.csv"

head(mydata)
my_data_harm <- mydata

#перевожу в 36 пиков
my_data_harm$GP14 <- my_data_harm$GP14 + my_data_harm$GP15
my_data_harm <- my_data_harm %>% select(-GP15)
my_data_harm$GP20 <- my_data_harm$GP20 + my_data_harm$GP21
my_data_harm <- my_data_harm %>% select(-GP21)
my_data_harm$GP24 <- my_data_harm$GP24 + my_data_harm$GP25
my_data_harm <- my_data_harm %>% select(-GP25)
head(my_data_harm)

sum(grepl('GP',colnames(my_data_harm)))

colnames(my_data_harm)[grepl('GP',colnames(my_data_harm))] <- paste0('GP',1:36)
head(my_data_harm)

print(cor(my_data_harm[,'GP14'],mydata[,'GP14']+mydata[,'GP15']))
print(cor(my_data_harm[,'GP19'],mydata[,'GP20']+mydata[,'GP21']))
print(cor(my_data_harm[,'GP22'],mydata[,'GP24']+mydata[,'GP25']))


cor_list <- NULL
for(i in 1:13){
  GP <- paste0('GP',i)
  cor_list <- c(cor_list,cor(my_data_harm[,..GP],mydata[,..GP]))
}

for(i in 15:18){
  GPx <- paste0('GP',i)
  GPy <- paste0('GP',i+1)
  cor_list <- c(cor_list,cor(my_data_harm[,..GPx],mydata[,..GPy]))
}

for(i in 20:21){
  GPx <- paste0('GP',i)
  GPy <- paste0('GP',i+2)
  cor_list <- c(cor_list,cor(my_data_harm[,..GPx],mydata[,..GPy]))
}

for(i in 23:36){
  GPx <- paste0('GP',i)
  GPy <- paste0('GP',i+3)
  cor_list <- c(cor_list,cor(my_data_harm[,..GPx],mydata[,..GPy]))
}
table(cor_list)

fwrite(my_data_harm, file= file_final)
