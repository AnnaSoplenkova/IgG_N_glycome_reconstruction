library(MultiABEL)
load('/mnt/polyomica/projects/MV_TPNG/MultiABEL/indep.snps.RData')
data <- load.summary(
  file.path("/mnt/polyomica/projects/GPD/result_ma/02_IgG_new_models_with_non_EUR/unification_out",
           paste0("igp", c(1:15, 20)),
           paste0("igp", c(1:15, 20), "_done.csv")), columnNames = c("rs_id", "ea", "ra", "eaf", "beta", "se", "n"),
  indep.snps = indep.snps
)
