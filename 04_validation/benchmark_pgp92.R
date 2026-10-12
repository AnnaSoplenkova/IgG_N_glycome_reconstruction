library(data.table)
library(qqman)
library(dplyr)

setwd("/mnt/polyomica/projects/GPD/benchmark/")

gwas_lincomb = fread("/mnt/polyomica/projects/GPD/benchmark/PGP92_linear_combination.csv")
gwas_gwis = fread("/mnt/polyomica/projects/GPD/benchmark/PGP92_GWIS.csv")
gwas_real = fread("/mnt/polyomica/projects/plasma_v2/14_all_META/01_GWAS_MAP/pgp92.csv")

cor(gwas_lincomb$b, gwas_real$beta_ma, use = "complete.obs")

df = inner_join(gwas_lincomb[, .(rs_id, b)], gwas_real[, .(rs_id, beta_ma)], by = "rs_id")
cor(df$b, df$beta_ma, use = "complete.obs")

df = inner_join(gwas_gwis[, .(rs_id, b)], df, by = "rs_id")
head(df)


names(gwas_real)


make_manhattan <- function(dt, pval_col, out_file, title) {
  dt <- as.data.frame(dt)
  dt <- dt[dt[[pval_col]] < 1e-3 & !is.na(dt[[pval_col]]), ]  # фильтр по порогу
  
    jpeg(out_file, width = 12, height = 5, units = "in", res = 150, type = "cairo")
  manhattan(dt,
            chr = "chr",
            bp = "bp",
            p = pval_col,
            snp = "rs_id",
            main = title,
            suggestiveline = FALSE,
            genomewideline = FALSE)
  dev.off()
}

unique( gwas_lincomb$chr)
table(is.na(gwas_lincomb$bp))

gwas_lincomb$chr <- as.integer(gwas_lincomb$chr)
table(is.na(gwas_lincomb$Pval), useNA = "always")
gwas_lincomb$bp  <- as.numeric(gwas_lincomb$bp)

# Строим и сохраняем все 3 графика
make_manhattan(gwas_lincomb, "Pval", "manhattan_lincomb_pgp92.jpg", "Linear Combination")
make_manhattan(gwas_gwis,    "Pval", "manhattan_gwis_pgp92.jpg",    "GWIS")
make_manhattan(gwas_real,    "p_ma", "manhattan_real_pgp92.jpg",    "Real GWAS")




lambda_lincomb <- median(qchisq(gwas_lincomb$Pval, df = 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, df = 1)
lambda_gwis <- median(qchisq(gwas_gwis$Pval, df = 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, df = 1)
lambda_real <- median(qchisq(gwas_real$p_ma, df = 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, df = 1)

lambda_lincomb
lambda_gwis
lambda_real

names(gwas_lincomb)
names(gwas_real)

cor(gwas_lincomb$b, gwas_real$beta_ma,use = "na.or.complete")
cor(gwas_gwis$b, gwas_real$beta_ma,use = "na.or.complete")


library(data.table)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(writexl)

# ===== Объединение датасетов =====
df_lincomb <- inner_join(
  gwas_lincomb[, .(rs_id, b, se, Z)],
  gwas_real[, .(rs_id, beta_ma, se_ma, z_ma)],
  by = "rs_id"
)

df_gwis <- inner_join(
  gwas_gwis[, .(rs_id, b, se, Z)],
  gwas_real[, .(rs_id, beta_ma, se_ma, z_ma)],
  by = "rs_id"
)

# ===== Функция построения scatter-плота =====
make_scatter <- function(df, x_col, y_col, x_lab, y_lab, title, color) {
  fit <- lm(as.formula(paste(y_col, "~", x_col)), data = df)
  co <- coef(fit)
  cor_val <- cor(df[[x_col]], df[[y_col]], use = "complete.obs")
  
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(color = color, alpha = 0.3, size = 0.8) +
    geom_abline(intercept = co[1], slope = co[2],
                color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "grey40", linewidth = 0.3) +
    labs(title = title, x = x_lab, y = y_lab) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      axis.title = element_text(size = 8),
      axis.text  = element_text(size = 7)
    ) +
    annotate("text",
             x = min(df[[x_col]], na.rm = TRUE) + 0.05 * diff(range(df[[x_col]], na.rm = TRUE)),
             y = max(df[[y_col]], na.rm = TRUE) - 0.05 * diff(range(df[[y_col]], na.rm = TRUE)),
             label = sprintf("Cor = %.3f\ny = %.3f + %.3f·x", cor_val, co[1], co[2]),
             size = 2.8, hjust = 0, vjust = 1)
}

# ===== 1) Beta-Beta plots =====
p_beta_lincomb <- make_scatter(df_lincomb, "b", "beta_ma",
                               "Beta (LinComb)", "Beta (Real)",
                               "Beta: LinComb vs Real", "#E41A1C")
p_beta_gwis    <- make_scatter(df_gwis, "b", "beta_ma",
                               "Beta (GWIS)", "Beta (Real)",
                               "Beta: GWIS vs Real", "#377EB8")

jpeg("beta_beta_comparison.jpg", width = 10, height = 5, units = "in", res = 200, type = "cairo")
grid.arrange(p_beta_lincomb, p_beta_gwis, ncol = 2)
dev.off()

# ===== 2) SE-SE plots =====
p_se_lincomb <- make_scatter(df_lincomb, "se", "se_ma",
                             "SE (LinComb)", "SE (Real)",
                             "SE: LinComb vs Real", "#E41A1C")
p_se_gwis    <- make_scatter(df_gwis, "se", "se_ma",
                             "SE (GWIS)", "SE (Real)",
                             "SE: GWIS vs Real", "#377EB8")

jpeg("se_se_comparison.jpg", width = 10, height = 5, units = "in", res = 200, type = "cairo")
grid.arrange(p_se_lincomb, p_se_gwis, ncol = 2)
dev.off()

# ===== 3) Z-Z plots =====
p_z_lincomb <- make_scatter(df_lincomb, "Z", "z_ma",
                            "Z (LinComb)", "Z (Real)",
                            "Z-score: LinComb vs Real", "#E41A1C")
p_z_gwis    <- make_scatter(df_gwis, "Z", "z_ma",
                            "Z (GWIS)", "Z (Real)",
                            "Z-score: GWIS vs Real", "#377EB8")

jpeg("z_z_comparison.jpg", width = 10, height = 5, units = "in", res = 200, type = "cairo")
grid.arrange(p_z_lincomb, p_z_gwis, ncol = 2)
dev.off()

# ===== 4) Лямбды геномного контроля → Excel =====
lambda <- data.table(
  method = c("Linear Combination", "GWIS", "Real GWAS"),
  lambda = c(lambda_lincomb, lambda_gwis, lambda_real)
)

write_xlsx(lambda, "genomic_control_lambda_pgp92.xlsx")
cat("Lambda saved:\n")
print(lambda)










