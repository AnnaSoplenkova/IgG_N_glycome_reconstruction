
setwd("/mnt/polyomica/projects/GPD/result_ma/02_IgG_new_models_with_non_EUR")

library(MethComp)
library(data.table)
library(ggplot2)
library(deming)
#install.packages("MethComp")
#install.packages("deming")



# Создаем список для хранения графиков
plots_list <- list()

peaks <- c(1:15, 20)

for (i_idx in seq_along(peaks)) {
  i <- peaks[i_idx]
  
  file_name = paste0("/mnt/polyomica/projects/GPD/result_ma/02_IgG_new_models_with_non_EUR/10k_IGP_",i,".csv")
  dt_pred = fread(file_name)
  file_name = paste0("/mnt/polyomica/projects/GPD/meta_real_IgG_without_TwinsUK/", "IGP",i,"/igp",i,".csv")
  dt_real = fread(file_name)
  
  # Объединение данных
  # Фильтрация снипов
  common_snps <- intersect(dt_pred$rs_id, dt_real$rs_id)
  
  # Фильтрация по размеру выборки
  dt_pred_filtered <- dt_pred[rs_id %in% common_snps & N_min > 1500]
  dt_real_filtered <- dt_real[rs_id %in% common_snps & n_ma > 1500]
  
  # Находим снипы, которые есть в ОБОИХ исследованиях после фильтрации по N
  #common_snps_with_N <- intersect(dt_pred_filtered$rs_id, dt_real_filtered$rs_id)
  
  # Фильтруем данные, оставляя только снипы с N>5000 в обоих исследованиях
  dt_pred_filtered <- dt_pred_filtered[rs_id %in% common_snps_with_N]
  dt_real_filtered <- dt_real_filtered[rs_id %in% common_snps_with_N]
  
  # Фильтрация по p-value < 1e-8 хотя бы в одном исследовании
  # Определяем значимые снипы
  sig_snps_pred <- dt_pred_filtered[Pval < 5e-8, rs_id]
  sig_snps_real <- dt_real_filtered[p_ma < 5e-8, rs_id]
  sig_snps <- union(sig_snps_pred, sig_snps_real)  # объединение - хотя бы в одном
  
  # Фильтруем данные по значимым снипам
  dt_pred_sig <- dt_pred_filtered[rs_id %in% sig_snps]
  dt_real_sig <- dt_real_filtered[rs_id %in% sig_snps]
  
  common_snps_final <- intersect(dt_pred_sig$rs_id, dt_real_sig$rs_id)
  dt_merged <- merge(
    dt_pred_sig[rs_id %in% common_snps_final, .(rs_id, b_pred = b, se_pred = se, p_pred = Pval)],
    dt_real_sig[rs_id %in% common_snps_final, .(rs_id, b_real = beta_ma, se_real = se_ma, p_real = p_ma)],
    by = "rs_id"
  )
  
  dt_plot <- dt_merged
  
  # 1. Рассчитываем регрессию Деминга
  # Выполняем регрессию Деминга
  deming_fit <- deming(
    b_pred ~ b_real, 
    data = dt_merged,
    xstd = se_real,  # стандартные ошибки для X
    ystd = se_pred   # стандартные ошибки для Y
  )
  
  # Извлекаем коэффициенты
  intercept_deming <- deming_fit$coefficients[1]
  slope_deming <- deming_fit$coefficients[2]
  
  beta_text <- paste0("Deming beta = ", round(slope_deming, 3))
  
  # 2. Создаем график - ВАЖНО: правильный порядок слоев
  p <- ggplot(dt_plot, aes(x = b_real, y = b_pred)) +
    # 1. Сначала рисуем УСЫ (погрешности) - они будут внизу
    geom_errorbarh(
      aes(xmin = b_real - se_real, xmax = b_real + se_real),
      height = 0,
      color = "gray40", 
      alpha = 0.4,
      linewidth = 0.3
    ) +
    geom_errorbar(
      aes(ymin = b_pred - se_pred, ymax = b_pred + se_pred),
      width = 0,
      color = "gray40", 
      alpha = 0.4,
      linewidth = 0.3
    ) +
    # 2. Затем точки - они поверх усов
    geom_point(
      alpha = 0.6, 
      size = 0.8, 
      shape = 16,
      color = "black"
    ) +
    # 3. Линия регрессии Деминга - поверх точек
    geom_abline(
      slope = slope_deming, 
      intercept = intercept_deming, 
      color = "blue",
      linewidth = 1.0,  # чуть толще для лучшей видимости
      alpha = 0.8
    ) +
    # 4. Прямая y=x для сравнения - также поверх
    geom_abline(
      slope = 1, 
      intercept = 0, 
      color = "red", 
      linetype = "dashed",
      linewidth = 0.7,
      alpha = 0.8
    ) +
    # 5. Аннотации - самый верхний слой
    annotate(
      "text", 
      x = min(dt_plot$b_real, na.rm = TRUE), 
      y = max(dt_plot$b_pred, na.rm = TRUE),
      label = beta_text, 
      hjust = 0, 
      vjust = 1, 
      color = "blue", 
      fontface = "bold", 
      size = 3
    ) +
    labs(
      x = "Beta (real)", 
      y = "Beta (predicted)",
      title = paste0(
        "IGP", i, 
        " (N snps = ", format(nrow(dt_merged), big.mark = ","), ")"
      )
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8)
    )
  
  # 3. Если слишком много точек, упрощаем визуализацию
  if (nrow(dt_plot) > 1000) {
    p <- p +
      # Обновляем усы (уже нарисованы, но можно переопределить)
      geom_errorbarh(
        aes(xmin = b_real - se_real, xmax = b_real + se_real),
        height = 0,
        color = "gray60", 
        alpha = 0.3,  # еще прозрачнее
        linewidth = 0.15
      ) +
      geom_errorbar(
        aes(ymin = b_pred - se_pred, ymax = b_pred + se_pred),
        width = 0,
        color = "gray60", 
        alpha = 0.3,
        linewidth = 0.15
      ) +
      # Обновляем точки
      geom_point(
        alpha = 0.3,  # более прозрачные
        size = 0.5,   # меньше
        shape = 16,
        color = "black"
      ) +
      # Линии делаем немного ярче, чтобы были видны
      geom_abline(
        slope = slope_deming, 
        intercept = intercept_deming, 
        color = "blue",
        linewidth = 1.2,  # толще
        alpha = 0.9       # менее прозрачная
      ) +
      geom_abline(
        slope = 1, 
        intercept = 0, 
        color = "red", 
        linetype = "dashed",
        linewidth = 0.8,
        alpha = 0.9
      )
  }
  
  print(p)
  
  # Добавляем график в список
  plots_list[[i_idx]] <- p
}


png("all_peaks_comparison.png", 
    width = 3508, height = 2480, res = 300) # A4 альбомный: высота и ширина поменяны местами

# Создаем сетку 4x4 (4 ряда, 4 столбца)
library(gridExtra)
# Разбиваем на 4 ряда по 4 графика
grid.arrange(grobs = plots_list, nrow = 4, ncol = 4)

dev.off()

