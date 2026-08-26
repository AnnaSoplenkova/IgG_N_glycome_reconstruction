library(MultiABEL)
library(data.table)

# 1. Загружаем список независимых SNP
load('/mnt/polyomica/projects/MV_TPNG/MultiABEL/indep.snps.RData')

# 2. Пути к 36 файлам (только pgp1–pgp36)
base_path <- "/mnt/polyomica/projects/plasma_v2/14_all_META/01_GWAS_MAP"
files <- paste0(base_path, "/pgp", 1:36, ".csv")

# 3. Проверка наличия всех файлов (необязательно, но полезно)
stopifnot(all(file.exists(files)))

# 4. Загрузка сводной статистики
#    Если имена колонок в ваших CSV другие — замените вектор columnNames
stats <- load.summary(
  files        = files,
  indep.snps   = indep.snps,
  columnNames  = c("rs_id", "ea", "ra", "eaf_ma", "beta_ma", "se_ma", "n_ma")
)

# 5. Проверка размерности (должно быть SNP × 36)
print(dim(stats$beta))

# 6. Сохранение объекта для дальнейшего использования
save(stats, file = "/mnt/polyomica/projects/MV_TPNG/DATA_discovery_7541/disc_multi_traits/N1/N1_stats.RData")

library(data.table)
cormatrix_new <- as.matrix(stats$cor.pheno)

fwrite(cormatrix_new,"/mnt/home/timoschuk/cormatrix_new.csv")
