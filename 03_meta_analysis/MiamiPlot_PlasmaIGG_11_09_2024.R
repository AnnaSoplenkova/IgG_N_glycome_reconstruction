#Делаю красивые майами плоты

library(data.table)
library(ggplot2)
library(ggrepel)

#загружаю дата фреймы

predicted_IgG = fread("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/IgG_8peaks_pvals_predicted_on_10k_pearson_0.8.csv")

IgG_gwas_pred = data.frame(BP = predicted_IgG$bp, CHR = predicted_IgG$chr, P = predicted_IgG$p_min, SNP = predicted_IgG$rs_id)
IgG_gwas_pred_lite = IgG_gwas_pred[IgG_gwas_pred$P < 0.0005,]
snp = IgG_gwas_pred_lite$SNP
#IgG_gwas_pred_lite = IgG_gwas_pred
nrow(IgG_gwas_pred_lite)


real_IgG = fread("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/Plasma_36peaks_pvals_10k.csv")
#names(real_IgG)

#head(IgG_gwas_real)
IgG_gwas_real = data.frame(BP = real_IgG$bp, CHR = real_IgG$chr, P = real_IgG$p_min, SNP = real_IgG$rs_id)
IgG_gwas_real = IgG_gwas_real[!is.na(IgG_gwas_real$BP),]
IgG_gwas_real = IgG_gwas_real[!is.na(IgG_gwas_real$P),]
IgG_gwas_real_lite = IgG_gwas_real[IgG_gwas_real$P < 0.0005,]
#IgG_gwas_real_lite = IgG_gwas_real
snp = c(snp, IgG_gwas_real_lite$SNP)
nrow(IgG_gwas_real_lite)

IgG_gwas_real_lite = IgG_gwas_real[IgG_gwas_real$SNP %in% snp,]
IgG_gwas_pred_lite = IgG_gwas_pred[IgG_gwas_pred$SNP %in% snp,]

#выбираю ОНП интереса 


IgG_snp = function(chr, bp, n = 250000) {
  start = (bp-n)
  end = (bp+n)
  cond = (IgG_gwas_pred_lite$CHR == chr) & (IgG_gwas_pred_lite$BP %in% (start:end))
  snp = IgG_gwas_pred_lite$SNP[cond]
  return(snp)
}

list_igg = list()

list_igg[[1]] = IgG_snp(chr = 1, bp = 25316924)
list_igg[[2]] = IgG_snp(chr = 2,	bp = 26148521, n = 350000)
list_igg[[3]] = IgG_snp(chr = 3,	bp = 186725887)
list_igg[[4]] = IgG_snp(chr = 6,	bp = 30000000, n = 2000000)
list_igg[[5]] = IgG_snp(chr = 6,	bp = 139629524)
list_igg[[6]] = IgG_snp(chr = 6,	bp = 143169723)
list_igg[[7]] = IgG_snp(chr = 7,	bp = 6502367)

list_igg[[8]] = IgG_snp(chr = 7,	bp = 50353144)
list_igg[[9]] = IgG_snp(chr = 7,	bp = 150942349)
list_igg[[10]] = IgG_snp(chr = 9,	bp = 33128021)
list_igg[[11]] = IgG_snp(chr = 17,	bp = 37912377)
list_igg[[12]] = IgG_snp(chr = 17,	bp = 43824360)
list_igg[[13]] = IgG_snp(chr = 17,	bp = 79217478)
list_igg[[14]] = IgG_snp(chr = 19,	bp = 19277296)

list_igg[[15]] = IgG_snp(chr = 20,	bp = 17832658)
list_igg[[16]] = IgG_snp(chr = 21,	bp = 36588180)
list_igg[[17]] = IgG_snp(chr = 22,	bp = 24170996)
list_igg[[18]] = IgG_snp(chr = 22,	bp = 39831986)

list_igg[[19]] = IgG_snp(chr = 3,	bp = 142902719)
list_igg[[20]] = IgG_snp(chr = 17,	bp = 16854480)

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

  
# Prepare X axis
axisdf <- don_plasma %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
ggplasma = ggplot(don_plasma, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey60", "gray73"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center,
                      expand = ggplot2::expansion(mult = 0.01),
                      position = "top" ) +
  ggplot2::labs(x = "", y = "Total blood plasma N-glycome") +
  #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  #geom_point(data=subset(don_plasma, is_highlight_other=="yes"), color="gray10", size=1.3) +
  geom_point(data=subset(don_plasma, is_highlight_igg=="yes"), color="limegreen", size=1.3) +
  geom_point(data=subset(don_plasma, is_highlight_tf=="yes"), color="magenta3", size=1.3) +
  
  # genime significanse line
  geom_hline(yintercept = -log10(5e-8/28), color = "red3", linetype = "dashed", size = 0.4) +
  
  # reverse MH plot
  scale_y_reverse(limits = c(70, 0),expand = ggplot2::expansion(mult = c(0, 0.02))) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(t = 0, b = 10, l = 10, r = 10)
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

# Prepare X axis
axisdf <- don_igg %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
ggigg = ggplot(don_igg, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey60", "gray73"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center,
                      expand = ggplot2::expansion(mult = 0.01),
                      guide = ggplot2::guide_axis(check.overlap =TRUE)) +
  scale_y_continuous(limits = c(0, 70),
                              expand = expansion(mult = c(0.02, 0))) +
  labs(x = "", y = "Synthetic IgG N-glycome") +    
  
  # Add highlighted points
  #geom_point(data=subset(don_plasma, is_highlight_other=="yes"), color="gray10", size=1.3) +
  geom_point(data=subset(don_igg, is_highlight_igg=="yes"), color="limegreen", size=1.3) +
  geom_point(data=subset(don_igg, is_highlight_tf=="yes"), color="magenta3", size=1.3) +
  
  # genime significanse line
  geom_hline(yintercept = -log10(5e-8/12), color = "red3", linetype = "dashed", size = 0.4) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(t = 10, b = 0, l = 10, r = 10)
  )

gridExtra::grid.arrange(ggigg, ggplasma, nrow = 2)


pdf("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/miami_plasma.pdf", width = 16, height = 9)
gridExtra::grid.arrange(ggigg, ggplasma, nrow = 2)
dev.off() 
