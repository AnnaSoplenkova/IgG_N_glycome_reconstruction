#лямбда геномного контроля

library(gap)
library(openxlsx)

igp = c(2,3,4,5,7,8,9,10,11,12,15,20)
lambda_discovery = rep(NA, length(igp))
lambda_replication = rep(NA, length(igp))

for (i in 1:length(igp)) {
  name = paste0("/mnt/polyomica/projects/GPD/discovery_meta/unification_out/igp",igp[i],"/igp",igp[i],"_done.csv")
  dt = fread(name)
  lambda_discovery[i] = gc.lambda(dt$p)
}


for (i in 1:length(igp)) {
  name = paste0("/mnt/polyomica/projects/GPD/replication_meta/unification_out/igp",igp[i],"/igp",igp[i],"_done.csv")
  dt = fread(name)
  lambda_replication[i] = gc.lambda(dt$p)
}

res = data.frame(igp, lambda_discovery, lambda_replication)
write.xlsx(res, "/mnt/polyomica/projects/GPD/lambda_gc.xlsx")
