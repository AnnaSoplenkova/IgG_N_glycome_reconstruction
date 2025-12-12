
setwd("/mnt/polyomica/projects/GPD/result_ma/02_IgG_new_models_with_non_EUR")


library(data.table)
library(ggplot2)

# Создаем список для хранения графиков
plots_list <- list()

peaks <- c(1:15, 20)

# Цикл по всем значениям i
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
  
  # Фильтрация по p-value < 1e-8 хотя бы в одном исследовании
  # Определяем значимые снипы
  sig_snps_pred <- dt_pred_filtered[Pval < 5e-8, rs_id]
  sig_snps_real <- dt_real_filtered[p_ma < 5e-8, rs_id]
  sig_snps <- union(sig_snps_pred, sig_snps_real)  # объединение - хотя бы в одном
  
  # Фильтруем данные по значимым снипам
  dt_pred_sig <- dt_pred_filtered[rs_id %in% sig_snps]
  dt_real_sig <- dt_real_filtered[rs_id %in% sig_snps]
  
  # Объединение данных
  common_snps_final <- intersect(dt_pred_sig$rs_id, dt_real_sig$rs_id)
  dt_merged <- merge(dt_pred_sig[rs_id %in% common_snps_final, .(rs_id, b_pred = b, p_pred = Pval)],
                     dt_real_sig[rs_id %in% common_snps_final, .(rs_id, b_real = beta_ma, p_real = p_ma)],
                     by = "rs_id")
  
 
    dt_plot <- dt_merged

  
  # Рассчитываем регрессию отдельно для получения коэффициентов
  lm_fit <- lm(b_pred ~ b_real, data = dt_merged)
  beta_coef <- coef(lm_fit)[["b_real"]]
  beta_text <- paste0("beta_regr = ", round(beta_coef, 3))
  
  # Создаем график
  p <- ggplot(dt_plot, aes(x = b_real, y = b_pred)) +
    geom_point(alpha = 0.3, size = 0.5, shape = 16) +  # оптимизированные точки
    geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue", 
                linewidth = 0.5) +  # линия регрессии с ДИ
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", 
                linewidth = 0.5) +  # прямая y=x
    # Добавляем аннотацию с бета-коэффициентом
    annotate("text", x = min(dt_plot$b_real), y = max(dt_plot$b_pred),
             label = beta_text, hjust = 0, vjust = 1, 
             color = "blue", fontface = "bold", size = 2) +
    labs(x = "Beta (real)", 
         y = "Beta (predicted)",
         title = paste0("IGP", i, " (N snps = ", 
                        format(nrow(dt_merged), big.mark = ","), 
                        ")")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 8),
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 5))
  
  # Добавляем график в список
  plots_list[[i_idx]] <- p
}

# Сохраняем все графики на одну PDF-страницу A4
pdf("all_peaks_comparison.pdf", width = 11.7, height = 8.3)  # A4 landscape

# Создаем сетку 4x4 (4 ряда, 4 столбца)
library(gridExtra)
# Разбиваем на 4 ряда по 4 графика
grid.arrange(grobs = plots_list, nrow = 4, ncol = 4)

dev.off()

# Альтернативный вариант на 2 страницы, если нужен (раскомментировать при необходимости):
# pdf("all_peaks_comparison_2pages.pdf", width = 11.7, height = 8.3)
# # Первая страница: 8 графиков
# grid.arrange(grobs = plots_list[1:8], nrow = 4, ncol = 2)
# # Вторая страница: остальные 8 графиков
# grid.arrange(grobs = plots_list[9:16], nrow = 4, ncol = 2)
# dev.off()
