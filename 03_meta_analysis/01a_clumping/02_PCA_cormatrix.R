#This script makes PCA analysis on 16 synthetic IGG peaks (r>0.7)

load("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/MultiABEL/indep.snps.RData")
cormatrix <- sst$cor.pheno
cormatrix<-as.matrix(cormatrix)
df.cor <- princomp(x = cormatrix, covmat = cormatrix)
summary(df.cor)

###############################
Importance of components:
                          Comp.1    Comp.2    Comp.3     Comp.4     Comp.5    Comp.6      Comp.7
Standard deviation     2.3199029 2.1624780 1.7512351 1.18562579 0.85945457 0.6418337 0.358481212
Proportion of Variance 0.3363718 0.2922695 0.1916765 0.08785678 0.04616639 0.0257469 0.008031799
Cumulative Proportion  0.3363718 0.6286413 0.8203178 0.90817461 0.95434099 0.9800879 0.988119694
                            Comp.8      Comp.9     Comp.10     Comp.11     Comp.12      Comp.13
Standard deviation     0.259167469 0.229790266 0.167332899 0.164328488 0.078341914 0.0633888988
Proportion of Variance 0.004197986 0.003300223 0.001750019 0.001687741 0.000383591 0.0002511345
Cumulative Proportion  0.992317680 0.995617903 0.997367922 0.999055662 0.999439253 0.9996903879
                            Comp.14      Comp.15      Comp.16
Standard deviation     0.0507825949 0.0377992178 3.075941e-02
Proportion of Variance 0.0001611795 0.0000892988 5.913385e-05
Cumulative Proportion  0.9998515673 0.9999408662 1.000000e+00
