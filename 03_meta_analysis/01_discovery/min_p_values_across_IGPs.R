#min p-value for 13 IGPs for MH plot

setwd("/mnt/polyomica/projects/GPD/result_ma/01_IgG_new_models")

library(data.table)
library(dplyr)

# Вектор IGP индексов
igp_indices <- c(2:5, 7:13, 15, 20)

# Создаем пустой список для хранения данных
all_data <- list()

# Читаем и обрабатываем каждый файл
for (i in igp_indices) {
  filename <- sprintf("10k_IGP_%d.csv", i)
  
  if (file.exists(filename)) {
    cat("Обрабатываю файл:", filename, "\n")
    
    # Читаем файл
    data <- fread(filename)
    
    # Выбираем нужные колонки и добавляем IGP индекс
    data_subset <- data %>%
      select(rs_id, chr, bp, Pval) %>%
      mutate(IGP = paste0("IGP_", i)) %>%
      rename(p_value = Pval)
    
    # Добавляем в список
    all_data[[length(all_data) + 1]] <- data_subset
    
  } else {
    warning(paste("Файл", filename, "не найден"))
  }
}

# Объединяем все данные
combined_data <- rbindlist(all_data)

# Находим минимальное p-значение для каждого rs_id
min_p_values <- combined_data %>%
  group_by(rs_id) %>%
  arrange(p_value) %>%
  slice(1) %>%
  ungroup() %>%
  select(rs_id, min_p = p_value, min_igp = IGP)

# Теперь добавляем информацию о chr и bp (берем из первого вхождения)
chr_bp_info <- combined_data %>%
  select(rs_id, chr, bp) %>%
  distinct(rs_id, .keep_all = TRUE)

# Создаем финальную таблицу
final_table <- min_p_values %>%
  left_join(chr_bp_info, by = "rs_id") %>%
  select(rs_id, chr, bp, min_p, min_igp)

# Сортируем по chr и bp для удобства
final_table <- final_table %>%
  arrange(chr, bp)


# Сохраняем результат
output_file <- "min_p_values_across_13_IGPs.csv"
fwrite(final_table, output_file)

