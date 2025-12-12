
#Делаю красивые майами плоты 16 пиков ИгГ, 36 пиков плазсы (не европейцы)


library(data.table)
library(ggplot2)
library(ggrepel)

#загружаю дата фреймы
setwd("/mnt/polyomica/projects/GPD/result_ma/02_IgG_new_models_with_non_EUR")

predicted_IgG = fread("min_p_values_across_16_IGPs.csv")



IgG_gwas_pred = data.frame(BP = predicted_IgG$bp, CHR = predicted_IgG$chr, P = predicted_IgG$min_p, SNP = predicted_IgG$rs_id)
IgG_gwas_pred_lite = IgG_gwas_pred[IgG_gwas_pred$P < 0.00005,]
# list of snps for 16 IGPs with small sample size n=170 CHR1
rs_vector = c("rs79574235", "rs144090205", "rs77999053")
IgG_gwas_pred_lite = IgG_gwas_pred_lite[!(IgG_gwas_pred_lite$SNP %in% rs_vector),]

# list of snps for 16 IGPs with small sample size n=350 CHR9
rs_vector = c("rs79574235", "rs144090205", "rs77999053")
IgG_gwas_pred_lite = IgG_gwas_pred_lite[!(IgG_gwas_pred_lite$SNP %in% rs_vector),]

# list of snps for 16 IGPs with small sample size n=170 CHR1
rs_vector = c("rs116467783","rs9931970", "rs9932247")
IgG_gwas_pred_lite = IgG_gwas_pred_lite[!(IgG_gwas_pred_lite$SNP %in% rs_vector),]


snp = IgG_gwas_pred_lite$SNP
#IgG_gwas_pred_lite = IgG_gwas_pred
nrow(IgG_gwas_pred_lite)
#head(IgG_gwas_pred_lite)



real_IgG = fread("/mnt/polyomica/projects/GPD/result_ma/02_IgG_new_models_with_non_EUR/min_p_values_across_36_PGPs.csv")
#names(real_IgG)

#head(IgG_gwas_real)
IgG_gwas_real = data.frame(BP = real_IgG$bp, CHR = real_IgG$chr, P = real_IgG$min_p, SNP = real_IgG$rs_id)
IgG_gwas_real = IgG_gwas_real[!is.na(IgG_gwas_real$BP),]
IgG_gwas_real = IgG_gwas_real[!is.na(IgG_gwas_real$P),]
IgG_gwas_real_lite = IgG_gwas_real[IgG_gwas_real$P < 0.00005,]


# list of snps for 16 IGPs with small sample size n=170
rs_vector = c("rs79574235", "rs144090205", "rs77999053")
IgG_gwas_real_lite = IgG_gwas_real_lite[!(IgG_gwas_real_lite$SNP %in% rs_vector),]
# list of snps for 16 IGPs with small sample size n=350 CHR9
rs_vector = c("rs79574235", "rs144090205", "rs77999053")
IgG_gwas_real_lite = IgG_gwas_real_lite[!(IgG_gwas_real_lite$SNP %in% rs_vector),]


#IgG_gwas_real_lite = IgG_gwas_real
snp = c(snp, IgG_gwas_real_lite$SNP)
nrow(IgG_gwas_real_lite)


IgG_gwas_real_lite = IgG_gwas_real[IgG_gwas_real$SNP %in% snp,]
IgG_gwas_pred_lite = IgG_gwas_pred[IgG_gwas_pred$SNP %in% snp,]

as.data.table(IgG_gwas_pred_lite[IgG_gwas_pred_lite$CHR == 16 & IgG_gwas_pred_lite$P < 5e-8 & IgG_gwas_pred_lite$BP < 78772743 + 300000,])

as.data.table(IgG_gwas_real_lite[IgG_gwas_real_lite$CHR == 5 & IgG_gwas_real_lite$P < 5e-8/28 ,])

#выбираю ОНП интереса 


IgG_snp = function(chr, bp, n = 250000) {
  start = (bp-n)
  end = (bp+n)
  cond = (IgG_gwas_pred_lite$CHR == chr) & (IgG_gwas_pred_lite$BP %in% (start:end))
  snp = IgG_gwas_pred_lite$SNP[cond]
  return(snp)
}

list_igg = list()



list_igg[[1]] = IgG_snp(chr = 1, bp = 25316924, n = 1000000)
list_igg[[2]] = IgG_snp(chr = 2,	bp = 26148521, n = 350000)
list_igg[[3]] = IgG_snp(chr = 3,	bp = 186725887)
list_igg[[4]] = IgG_snp(chr = 6,	bp = 30000000, n = 2500000)
list_igg[[5]] = IgG_snp(chr = 6,	bp = 139629524)
list_igg[[6]] = IgG_snp(chr = 6,	bp = 143169723)
#list_igg[[7]] = IgG_snp(chr = 7,	bp = 6502367)

list_igg[[7]] = NA
list_igg[[8]] = IgG_snp(chr = 7,	bp = 50353144)
list_igg[[9]] = IgG_snp(chr = 7,	bp = 150942349)
list_igg[[10]] = IgG_snp(chr = 9,	bp = 33128021)

list_igg[[11]] = IgG_snp(chr = 14,	bp = 106000231)
list_igg[[12]] = IgG_snp(chr = 17,	bp = 37912377)
list_igg[[13]] = IgG_snp(chr = 17,	bp = 43824360)
list_igg[[14]] = IgG_snp(chr = 17,	bp = 79217478)
list_igg[[15]] = IgG_snp(chr = 19,	bp = 19277296)

list_igg[[16]] = IgG_snp(chr = 20,	bp = 17832658)
list_igg[[17]] = IgG_snp(chr = 21,	bp = 36588180)
list_igg[[18]] = IgG_snp(chr = 22,	bp = 24170996)
list_igg[[19]] = IgG_snp(chr = 22,	bp = 39831986, n = 350000)

list_igg[[20]] = IgG_snp(chr = 3,	bp = 142902719)
list_igg[[21]] = IgG_snp(chr = 17,	bp = 16854480)

vect_igg = unlist(list_igg)

TF_snp = function(chr, bp, n = 250000) {
  start = (bp-n)
  end = (bp+n)
  cond = (IgG_gwas_real_lite$CHR == chr) & (IgG_gwas_real_lite$BP %in% (start:end))
  snp = IgG_gwas_real_lite$SNP[cond]
  return(snp)
}

list_tf = list()

list_tf[[1]] = TF_snp(chr = 2,	bp = 135014685)
list_tf[[2]] = TF_snp(chr = 11,	bp = 126234820)
list_tf[[3]] = TF_snp(chr = 11,	bp = 134264529)
list_tf[[4]] = TF_snp(chr = 12,	bp = 121420263)

vect_tf = unlist(list_tf)


# List of SNPs to highlight are in the snpsOfInterest object
names(IgG_gwas_pred_lite)



other_snp = function(chr, bp, n = 250000) {
  start = (bp-n)
  end = (bp+n)
  cond = (IgG_gwas_real_lite$CHR == chr) & (IgG_gwas_real_lite$BP %in% (start:end))
  snp = IgG_gwas_real_lite$SNP[cond]
  return(snp)
}

list_other = list()

list_other[[1]] = other_snp(chr = 1, bp = 161666633)
list_other[[2]] = other_snp(chr = 4, bp =	103411099)
list_other[[3]] = other_snp(chr = 10,	bp = 94429467)
list_other[[4]] = other_snp(chr = 14,	bp = 65855232)

list_other[[5]] = other_snp(chr = 1, bp =	93344692)
list_other[[6]] = other_snp(chr = 1, bp =	161590277)
list_other[[7]] = other_snp(chr = 2, bp =	27598097)
list_other[[8]] = other_snp(chr = 4, bp =	103422504)
list_other[[9]] = other_snp(chr = 7, bp =	73013269)
list_other[[10]] = other_snp(chr = 8, bp = 126500031)
list_other[[11]] = other_snp(chr = 11, bp =	61602003)
list_other[[12]] = other_snp(chr = 11, bp =	123344959)

list_other[[13]] = other_snp(chr = 12, bp =	31997635)
list_other[[14]] = other_snp(chr = 13, bp =	41611271)
list_other[[15]] = other_snp(chr = 14, bp =	65879354)
list_other[[16]] = other_snp(chr = 15, bp =	60932837)
list_other[[17]] = other_snp(chr = 16, bp =	72106024)
list_other[[18]] = other_snp(chr = 19, bp =	5848000)
list_other[[19]] = other_snp(chr = 19, bp =	45412079)

vect_other = unlist(list_other)







# Prepare the dataset
don_plasma <- IgG_gwas_real_lite %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len) )-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(IgG_gwas_real_lite, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  # Add highlight information
  mutate( is_highlight_igg=ifelse(SNP %in% vect_igg, "yes", "no")) %>%
  mutate( is_highlight_tf=ifelse(SNP %in% vect_tf, "yes", "no")) %>%
  mutate( is_highlight_other=ifelse(SNP %in% vect_other, "yes", "no"))




# Prepare X axis для plasma
axisdf <- don_plasma %>% group_by(CHR) %>% summarize(center = (max(BPcum) + min(BPcum)) / 2)

# Make the plot для plasma - с увеличенными шрифтами и точками
ggplasma = ggplot(don_plasma, aes(x = BPcum, y = -log10(P))) +
  
  # Show all points - УВЕЛИЧИВАЕМ РАЗМЕР ТОЧЕК
  geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 2.0) +
  scale_color_manual(values = rep(c("grey60", "gray73"), 22)) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center,
                     expand = ggplot2::expansion(mult = 0.01),
                     position = "top") +
  labs(x = "", y = "Total blood plasma N-glycome") +
  
  # Add highlighted points - УВЕЛИЧИВАЕМ РАЗМЕР ВЫДЕЛЕННЫХ ТОЧЕК
  geom_point(data = subset(don_plasma, is_highlight_igg == "yes"), 
             color = "limegreen", size = 2.0) +
  geom_point(data = subset(don_plasma, is_highlight_tf == "yes"), 
             color = "magenta3", size = 2.0) +
  
  # genome significance line - ЕЩЕ ТОЛЩЕ
  geom_hline(yintercept = -log10(5e-8/28), color = "red3", 
             linetype = "dashed", size = 1.2, alpha = 0.9) +
  
  # reverse MH plot
  scale_y_reverse(limits = c(70, 0), expand = ggplot2::expansion(mult = c(0, 0.02))) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = ggplot2::margin(t = 0, b = 10, l = 10, r = 10),
    
    # ЗНАЧИТЕЛЬНО УВЕЛИЧИВАЕМ ШРИФТЫ
    axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 10)),
    axis.text.y = element_text(size = 13, face = "bold")
  )


# Prepare the dataset
don_igg <- IgG_gwas_pred_lite %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len) )-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(IgG_gwas_pred_lite, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  # Add highlight information
  mutate( is_highlight_igg=ifelse(SNP %in% vect_igg, "yes", "no")) %>%
  mutate( is_highlight_tf=ifelse(SNP %in% vect_tf, "yes", "no"))



# Prepare X axis для IgG
axisdf <- don_igg %>% group_by(CHR) %>% summarize(center = (max(BPcum) + min(BPcum)) / 2)

# Make the plot для IgG - с увеличенными шрифтами и точками
ggigg = ggplot(don_igg, aes(x = BPcum, y = -log10(P))) +
  
  # Show all points - УВЕЛИЧИВАЕМ РАЗМЕР ТОЧЕК
  geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 2.0) +
  scale_color_manual(values = rep(c("grey60", "gray73"), 22)) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center,
                     expand = ggplot2::expansion(mult = 0.01),
                     guide = ggplot2::guide_axis(check.overlap = TRUE)) +
  scale_y_continuous(limits = c(0, 70),
                     expand = expansion(mult = c(0.02, 0))) +
  labs(x = "", y = "Reconstructed IgG N-glycome") +    
  
  # Add highlighted points - УВЕЛИЧИВАЕМ РАЗМЕР ВЫДЕЛЕННЫХ ТОЧЕК
  geom_point(data = subset(don_igg, is_highlight_igg == "yes"), 
             color = "limegreen", size = 2.0) +
  geom_point(data = subset(don_igg, is_highlight_tf == "yes"), 
             color = "magenta3", size = 2.0) +
  
  # genome significance line - ЕЩЕ ТОЛЩЕ
  geom_hline(yintercept = -log10(5e-8/16), color = "red3", 
             linetype = "dashed", size = 1.2, alpha = 0.9) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = ggplot2::margin(t = 10, b = 0, l = 10, r = 10),
    
    # ЗНАЧИТЕЛЬНО УВЕЛИЧИВАЕМ ШРИФТЫ
    axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 13, face = "bold")
  )

# Создаем композицию графиков
gridExtra::grid.arrange(ggigg, ggplasma, nrow = 2)

# Сохраняем в PDF с ОПТИМИЗАЦИЕЙ ДЛЯ УМЕНЬШЕНИЯ РАЗМЕРА
pdf("miami_plasma.pdf", width = 16, height = 9,
    useDingbats = FALSE,  # уменьшает размер файла
    compress = TRUE)      # включает сжатие

gridExtra::grid.arrange(ggigg, ggplasma, nrow = 2)
dev.off()
