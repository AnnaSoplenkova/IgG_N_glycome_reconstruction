#power of replication
p_disc=5e-8/12
N_disc=12384

N_rep=5107
alpha=0.05/5

thr<-qchisq(alpha, df=1, lower.tail = FALSE)
NCP<-(qchisq(p_disc,df=1,lower.tail = FALSE)-1)*N_rep/N_disc
M<-pchisq(thr, df=1, ncp = NCP, lower.tail = FALSE)
M
