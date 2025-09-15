# Copyright 2019 Sodbo Sharapov, Elgaeva Elizaveta, Lucija Klaric, Maja Pučić Baković, Yurii Aulchenko, Gordan Lauc
library(data.table)

plasma_derived_traits <- function(my_data){
  
  my_data <- my_data[,grep(colnames(my_data), pattern='GP'), with = FALSE]
  
  colnames(my_data) <- c('PGP1', 'PGP2', 'PGP3', 'PGP4', 'PGP5', 'PGP6', 'PGP7', 'PGP8', 'PGP9', 'PGP10',
                         'PGP11', 'PGP12', 'PGP13', 'PGP14', 'PGP15', 'PGP16', 'PGP17', 'PGP18', 'PGP19', 'PGP20',
                         'PGP21', 'PGP22', 'PGP23', 'PGP24', 'PGP25', 'PGP26', 'PGP27', 'PGP28', 'PGP29', 'PGP30',
                         'PGP31', 'PGP32', 'PGP33', 'PGP34', 'PGP35', 'PGP36')
  
  my_data$PGP37 <- with(my_data, (PGP13 + PGP15 + PGP20 + PGP28 + PGP31 + PGP32) / (PGP4 + PGP5 + PGP10 + PGP13 + PGP15 + PGP20 + PGP28 + PGP31 + PGP32))
  
  my_data$PGP38 <- with(my_data, (PGP16 + PGP21) / (PGP6 + PGP11 + PGP16 + PGP21))
  
  my_data$PGP39 <- with(my_data, (PGP13 + PGP15 + PGP20 + PGP28 + PGP31 + PGP32) / (PGP1 + PGP4 + PGP5 + PGP10 + PGP13 + PGP15 + PGP20 + PGP28 + PGP31 + PGP32))
  
  my_data$PGP40 <- with(my_data, (PGP16 + PGP21) / (PGP2 + PGP6 + PGP11 + PGP16 + PGP21))
  
  my_data$PGP41 <- with(my_data, PGP13 / (PGP4 + PGP5 + PGP13))
  
  my_data$PGP42 <- with(my_data, PGP15 / (PGP10 + PGP15 + PGP20))
  
  my_data$PGP43 <- with(my_data, PGP20 / (PGP10 + PGP15 + PGP20))
  
  my_data$PGP44 <- with(my_data, PGP16 / (PGP11 + PGP16 + PGP21))
  
  my_data$PGP45 <- with(my_data, PGP21 / (PGP11 + PGP16 + PGP21))
  
  my_data$PGP46 <- with(my_data, (PGP13 + PGP15 + PGP16) / (PGP20 + PGP21 + PGP24))
  
  my_data$PGP47 <- with(my_data, (PGP13 + PGP15) / (PGP20 + PGP24))
  
  my_data$PGP48 <- with(my_data, PGP16 / PGP21)
  
  my_data$PGP49 <- with(my_data, (PGP13 + PGP15 + PGP16) / (PGP28 + PGP31 + PGP32))
  
  my_data$PGP50 <- with(my_data, (PGP13 + PGP15) / (PGP28 + PGP31 + PGP32))
  
  my_data$PGP51 <- with(my_data, (PGP20 + PGP21) / (PGP28 + PGP31 + PGP32))
  
  my_data$PGP52 <- with(my_data, PGP20 / (PGP28 + PGP31 + PGP32))
  
  my_data$PGP53 <- with(my_data, (PGP16 + PGP21) / (PGP13 + PGP15 + PGP20 + PGP28 + PGP31 + PGP32))
  
  my_data$PGP54 <- with(my_data, PGP16 / (PGP13 + PGP15))
  
  my_data$PGP55 <- with(my_data, PGP16 / (PGP13 + PGP15 + PGP16))
  
  my_data$PGP56 <- with(my_data, PGP21 / (PGP20 + PGP24))
  
  my_data$PGP57 <- with(my_data, PGP21 / (PGP20 + PGP21 + PGP24))
  
  my_data$PGP58 <- with(my_data, PGP1 / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP59 <- with(my_data, PGP2 / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP60 <- with(my_data, PGP3 / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP61 <- with(my_data, PGP4 / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP62 <- with(my_data, PGP5 / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP63 <- with(my_data, PGP6 / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP64 <- with(my_data, PGP7 / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP65 <- with(my_data, PGP8 / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP66 <- with(my_data, PGP9 / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP67 <- with(my_data, PGP10 / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP68 <- with(my_data, PGP11 / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP69 <- with(my_data, PGP18 / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP70 <- with(my_data, (PGP1 + PGP2) / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP71 <- with(my_data, (PGP3 + PGP4 + PGP5 + PGP6) / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP72 <- with(my_data, (PGP8 + PGP9 + PGP10 + PGP11) / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP73 <- with(my_data, (PGP1 + PGP2 + PGP4 + PGP5 + PGP6 + PGP10 + PGP11) / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP74 <- with(my_data, (PGP4 + PGP5 + PGP6) / (PGP3 + PGP4 + PGP5 + PGP6))
  
  my_data$PGP75 <- with(my_data, (PGP10 + PGP11) / (PGP8 + PGP9 + PGP10 + PGP11))
  
  my_data$PGP76 <- with(my_data, (PGP1 + PGP4 + PGP5 + PGP10) / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP77 <- with(my_data, PGP1/(PGP1 + PGP2) )
  
  my_data$PGP78 <- with(my_data, (PGP4 + PGP5) / (PGP3 + PGP4 + PGP5 + PGP6))
  
  my_data$PGP79 <- with(my_data, PGP10 / (PGP8 + PGP9 + PGP10 + PGP11))
  
  my_data$PGP80 <- with(my_data, (PGP2 + PGP6 + PGP11) / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP81 <- with(my_data, PGP2/(PGP1 + PGP2) )
  
  my_data$PGP82 <- with(my_data, PGP6 / (PGP3 + PGP4 + PGP5 + PGP6))
  
  my_data$PGP83 <- with(my_data, PGP11 / (PGP8 + PGP9 + PGP10 + PGP11))
  
  my_data$PGP84 <- with(my_data, (PGP2 + PGP6 + PGP11) / (PGP1 + PGP4 + PGP5 + PGP10))
  
  my_data$PGP85 <- with(my_data, (PGP2 + PGP6 + PGP11) / (PGP1 + PGP2 + PGP4 + PGP5 + PGP6 + PGP10 + PGP11))
  
  my_data$PGP86 <- with(my_data, (PGP1 + PGP4 + PGP5 + PGP10) / (PGP2 + PGP3 + PGP6 + PGP9 + PGP11))
  
  my_data$PGP87 <- with(my_data, (PGP3 + PGP9) / (PGP1 + PGP2 + PGP4 + PGP5 + PGP6 + PGP10 + PGP11))
  
  my_data$PGP88 <- with(my_data, PGP11 / PGP10)
  
  my_data$PGP89 <- with(my_data, PGP11 / (PGP10 + PGP11))
  
  my_data$PGP90 <- with(my_data, PGP10 / (PGP9 + PGP11))
  
  my_data$PGP91 <- with(my_data, PGP9 / (PGP10 + PGP11))
  
  my_data$PGP92 <- with(my_data, PGP24 + PGP30 + PGP32 + PGP36)
  
  my_data$PGP93 <- with(my_data, PGP1 + PGP2 + PGP4 + PGP5 + PGP6 + PGP10 + PGP11 + PGP13 + PGP15 + PGP16 + PGP20 + PGP21 + PGP28 + PGP31 + PGP32)
  
  my_data$PGP94 <- with(my_data, PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18)
  
  my_data$PGP95 <- with(my_data, 0.5 * PGP12 + PGP13 + PGP14 + PGP15 + PGP16)
  
  my_data$PGP96 <- with(my_data, PGP17 + PGP19 + PGP20 + PGP21 + PGP22 + PGP23 + PGP24)
  
  my_data$PGP97 <- with(my_data, PGP25 + PGP26 + PGP27 + PGP28 + PGP29 + PGP30 + PGP31 + PGP32 + PGP33)
  
  my_data$PGP98 <- with(my_data, PGP34 + PGP35 + PGP36)
  
  my_data$PGP99 <- with(my_data, PGP1 + PGP2)
  
  my_data$PGP100 <- with(my_data, PGP3 + PGP4 + PGP5 + PGP6 + PGP13)
  
  my_data$PGP101 <- with(my_data, PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP14 + PGP15 + PGP16 + PGP17 + PGP19 + PGP20 + PGP21)
  
  my_data$PGP102 <- with(my_data, PGP22 + PGP23 + PGP24 + PGP25 + PGP26 + PGP27 + PGP28 + PGP29 + PGP30 + PGP31 + PGP32)
  
  my_data$PGP103 <- with(my_data, PGP33 + PGP34 + PGP35 + PGP36)
  
  my_data$PGP104 <- with(my_data, PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP13 + PGP14 + PGP15 + PGP16 + PGP17 + PGP19 + PGP20 + PGP21)
  
  my_data$PGP105 <- with(my_data, PGP22 + PGP23 + PGP24 + PGP25 + PGP26 + PGP27 + PGP28 + PGP29 + PGP30 + PGP31 + PGP32)
  
  my_data$PGP106 <- with(my_data, PGP33 + PGP34 + PGP35 + PGP36)
  
  my_data$PGP107 <- with(my_data, PGP7 + 0.5 * PGP12 + PGP18)
  
  my_data$PGP108 <- with(my_data, PGP2 + PGP3 + PGP6 + PGP9 + PGP11 + PGP16 + PGP21)
  
  my_data$PGP109 <- with(my_data, (PGP22 + PGP23 + PGP24) / (PGP25 + PGP26 + PGP27 + PGP28 + PGP29 + PGP30 + PGP31 + PGP32))
  
  my_data$PGP110 <- with(my_data, PGP33 / (PGP34 + PGP35 + PGP36))
  
  my_data$PGP111 <- with(my_data, (PGP28 + PGP31 + PGP32) / (PGP22 + PGP23 + PGP24 + PGP25 + PGP26 + PGP27 + PGP28 + PGP29 + PGP30 + PGP31 + PGP32))
  
  my_data$PGP112 <- with(my_data, (PGP24 + PGP30 + PGP32) / (PGP22 + PGP23 + PGP24 + PGP25 + PGP26 + PGP27 + PGP28 + PGP29 + PGP30 + PGP31 + PGP32))
  
  my_data$PGP113 <- with(my_data, PGP36 / (PGP33 + PGP34 + PGP35 + PGP36))
  
  my_data$PGP114 <- with(my_data, (0.5 * PGP12) / (PGP1 + PGP2 + PGP3 + PGP4 + PGP5 + PGP6 + PGP7 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP18))
  
  my_data$PGP115 <- with(my_data, (0.5 * PGP12 + PGP13 + PGP14 + PGP15 + PGP16 + PGP17 + PGP19 + PGP20 + PGP21 + PGP22 + PGP23 + PGP24 + PGP25 + PGP26 + PGP27 + PGP28 + PGP29 + PGP30 + PGP31 + PGP32 + PGP33 + PGP34 + PGP35 + PGP36))
  
  my_data$PGP116 <- with(my_data, (PGP3 + PGP4 + PGP5 + PGP6 + PGP8 + PGP9 + PGP10 + PGP11 + 0.5 * PGP12 + PGP13 + PGP14 + PGP15 + PGP16 + PGP17 + PGP19 + PGP20 + PGP21 + PGP22 + PGP23 + PGP24 + PGP25 + PGP26 + PGP27 + PGP28 + PGP29 + PGP30 + PGP31 + PGP32 + PGP33 + PGP34 + PGP35 + PGP36))
  
  my_data$PGP117 <- with(my_data, (PGP3 + PGP4 + PGP5 + PGP6 + PGP8 + PGP9 + PGP10 + PGP11))
  
  my_data <- my_data[, -c(1:36)]
  
  return(my_data)
}



plasma_train = fread("/mnt/polyomica/projects/GPD/peak_prediction/data/plasma_training_set_36peaks_2023.csv")
plasma_test = fread("/mnt/polyomica/projects/GPD/peak_prediction/data/plasma_test_set_36peaks_2023.csv")

grep(colnames(plasma_train), pattern='GP')
plasma_train[,4:39] = exp(plasma_train[,4:39])
grep(colnames(plasma_test), pattern='GP')
plasma_test[,4:39] = exp(plasma_test[,4:39])

derived_train = plasma_derived_traits(plasma_train)
derived_test = plasma_derived_traits(plasma_test)



plasma_train[,4:39] = log(plasma_train[,4:39])
plasma_test[,4:39] = log(plasma_test[,4:39])

derived_train = log(derived_train)
derived_test = log(derived_test)


names(plasma_train)[4:39] = c('PGP1', 'PGP2', 'PGP3', 'PGP4', 'PGP5', 'PGP6', 'PGP7', 'PGP8', 'PGP9', 'PGP10',
                                'PGP11', 'PGP12', 'PGP13', 'PGP14', 'PGP15', 'PGP16', 'PGP17', 'PGP18', 'PGP19', 'PGP20',
                                'PGP21', 'PGP22', 'PGP23', 'PGP24', 'PGP25', 'PGP26', 'PGP27', 'PGP28', 'PGP29', 'PGP30',
                                'PGP31', 'PGP32', 'PGP33', 'PGP34', 'PGP35', 'PGP36')
names(plasma_test)[4:39]  = c('PGP1', 'PGP2', 'PGP3', 'PGP4', 'PGP5', 'PGP6', 'PGP7', 'PGP8', 'PGP9', 'PGP10',
                       'PGP11', 'PGP12', 'PGP13', 'PGP14', 'PGP15', 'PGP16', 'PGP17', 'PGP18', 'PGP19', 'PGP20',
                       'PGP21', 'PGP22', 'PGP23', 'PGP24', 'PGP25', 'PGP26', 'PGP27', 'PGP28', 'PGP29', 'PGP30',
                       'PGP31', 'PGP32', 'PGP33', 'PGP34', 'PGP35', 'PGP36')


res_train = cbind(plasma_train, derived_train)
res_test = cbind(plasma_test, derived_test)


fwrite(res_train, "/mnt/polyomica/projects/GPD/peak_prediction/data/plasma_training_set_117pgp_2025.csv")
fwrite(res_test, "/mnt/polyomica/projects/GPD/peak_prediction/data/plasma_test_set_117pgp_2025.csv")

