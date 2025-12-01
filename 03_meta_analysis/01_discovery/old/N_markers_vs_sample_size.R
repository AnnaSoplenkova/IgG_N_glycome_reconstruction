# Скрипт для построения графика зависимости числа маркеров от объема выборки

dt = fread("/mnt/polyomica/projects/GPD/discovery_meta/IGP5/igp5.csv")
sample_sizes = dt$n_ma


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
       log = "x",  # логарифмическая шкала по X (часто полезно для GWAS)
       xlab = "Объем выборки (человек)", 
       ylab = "Число маркеров с выборкой >= X",
       main = "Зависимость числа маркеров от объема выборки",
       col = "blue", 
       lwd = 2)
  
  # Добавляем сетку для лучшей читаемости
  grid()
  
  # Возвращаем данные графика неявно
  invisible(list(x = x, y = y))
}


# Строим график
jpeg("markers_vs_sample_size.jpeg")
plot_markers_vs_sample_size(sample_sizes)
dev.off()

setwd("/mnt/polyomica/projects/GPD/discovery_meta/")
