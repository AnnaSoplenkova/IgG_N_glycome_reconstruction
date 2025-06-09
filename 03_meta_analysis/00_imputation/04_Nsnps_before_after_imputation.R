#число снипов в каждой когорте до/после импутации


setwd("/mnt/polyomica/projects/GPD/imputation_output/imput_data/")

dt = fread("ORCADES/igp2.csv")
nrow(dt)


dt = fread("/mnt/polyomica/projects/GPD/synthetic_IgG_TPNG_without_TwinsUK/unification_out/igp2/igp2_done.csv")

dt = fread("/mnt/polyomica/projects/GPD/imputation_output/final/ORCADES/02_unification_out/igp2/igp2_done.csv")
nrow(dt)

table(dt$n > 6000)
max(dt$n)


# Функция для построения графика зависимости числа маркеров от объема выборки
plot_markers_vs_sample_size <- function(sample_sizes) {
  # sample_sizes - вектор объемов выборки для каждого маркера
  
  # Упорядочиваем объемы выборки по убыванию
  sorted_sizes <- sort(sample_sizes, decreasing = TRUE)
  
  # Создаем данные для графика:
  # x - объем выборки (по убыванию)
  # y - количество маркеров с объемом выборки >= x
  x <- sorted_sizes
  y <- 1:length(sorted_sizes)
  
  # Строим график
  plot(x, y, 
       type = "l", 
       #log = "x",  # логарифмическая шкала по X (часто полезно для GWAS)
       xlab = "Sample size", 
       ylab = "Number of snps with sample size >= X",
       main = "After imputation",
       col = "blue", 
       lwd = 2)
  
  # Добавляем сетку для лучшей читаемости
  grid()
  
  # Возвращаем данные графика неявно
  invisible(list(x = x, y = y))
}


dt = fread("/mnt/polyomica/projects/GPD/discovery_meta/IGP5/igp5.csv")
sample_sizes = dt$n_ma
# Строим график
jpeg("markers_vs_sample_size.jpeg")
plot_markers_vs_sample_size(sample_sizes)
dev.off()
