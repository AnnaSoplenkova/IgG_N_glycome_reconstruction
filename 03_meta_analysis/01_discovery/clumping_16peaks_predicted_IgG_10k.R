function_for_shlop_28_12_2017 <- function(locus_table,p_value="Pval",pos="bp",snp="rs_id",                                                                                                                         
                                          delta=2.5e5,chr="chr",thr=5e-8,trait=NULL){                                                                                                                                 
  
  locus_table[,p_value] <- as.numeric(locus_table[,p_value])                                                                                                                                                     
  
  if (!is.null(trait)){                                                                                                                                                                                          
    traits="traits"                                                                                                                                                                                            
    locus_table=cbind(locus_table,traits=locus_table[,trait])                                                                                                                                                  
    locus_table[,traits]=as.character(locus_table[,traits])                                                                                                                                                    
  }                                                                                                                                                                                                              
  
  out=locus_table[0,]                                                                                                                                                                                            
  
  locus_table=locus_table[locus_table[,p_value]<=thr,]                                                                                                                                                           
  
  i=1                                                                                                                                                                                                            
  if (nrow(locus_table)>0){                                                                                                                                                                                      
    locus_table[,pos]=as.numeric(locus_table[,pos])                                                                                                                                                            
    locus_table[,p_value]=as.numeric(locus_table[,p_value])                                                                                                                                                    
    Zx <-locus_table                                                                                                                                                                                           
    Zx=Zx[order(Zx[,p_value]),]                                                                                                                                                                                
    #n_traits=1                                                                                                                                                                                                
    #Zx=cbind(Zx,n_traits)                                                                                                                                                                                     
    i=1                                                                                                                                                                                                        
    while (nrow(Zx)>0){                                                                                                                                                                                        
      ind=which((abs(Zx[i,pos]-Zx[,pos])<=delta)&(Zx[i,chr]==Zx[,chr]))                                                                                                                                      
      
      if (!is.null(trait)){                                                                                                                                                                                  
        Zx[i,traits]=paste(unique(Zx[ind,trait]),collapse = ";")                                                                                                                                           
      }                                                                                                                                                                                                      
      
      out=rbind(out,Zx[i,])                                                                                                                                                                                  
      Zx=Zx[-ind,]                                                                                                                                                                                           
    }                                                                                                                                                                                                          
    rownames(out)=as.character(out[,snp])                                                                                                                                                                      
  }                                                                                                                                                                                                              
  
  if (!is.null(trait)){                                                                                                                                                                                          
    j=1                                                                                                                                                                                                        
    out=cbind(out,Ntraits=1)                                                                                                                                                                                   
    out[,"Ntraits"]=as.numeric(out[,"Ntraits"])                                                                                                                                                                
    for (j in 1:nrow(out)){                                                                                                                                                                                    
      trs=unique(unlist(strsplit(out[j,traits],split = ";")))                                                                                                                                                
      out[j,traits]=paste(trs,collapse = ";")                                                                                                                                                                
      out[j,"Ntraits"]=length(trs)                                                                                                                                                                           
    }                                                                                                                                                                                                          
  }                                                                                                                                                                                                              
  
  return(out) 
}






library(data.table)




#порог - пирсон 70%
peaks = c(1:15,20)

out <- NULL
thr <- 5e-8/8



for (i in peaks) {
  file_name = paste0("/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/unification_out/igp",i,"/igp",i,"_done.csv")
  sst <- fread(file_name, header=T, stringsAsFactors = F, data.table=F)
  sst <- sst[!is.na(sst$p),]
  sst$gwas_id = rep(paste0("igp_",i), nrow(sst))
  #sst_sm <- sst[(pmin(sst$eaf_ma,1-sst$eaf_ma)>=0.01),]
  #sst_sm <- sst_sm[sst_sm$n_ma>=(7500*0.9),]
  lt <- function_for_shlop_28_12_2017(sst,p_value="p",pos="bp",snp="rs_id", delta=2.5e5,chr="chr", thr=thr)
  if (nrow(lt)>0 ) {
    #lt=cbind(lt,trait)
    out=rbind(out,lt)
  }
}
head(out)
dim(out)
bt <- function_for_shlop_28_12_2017(out,trait="gwas_id",p_value="p",pos="bp",snp="rs_id",delta=6e5,chr="chr")
bt <- bt[order(bt$chr,bt$bp),]
colnames(bt)
str(bt)
bt$trait <- as.character(bt$trait)
str(bt)
library(openxlsx)
write.xlsx(bt, "/mnt/polyomica/projects/GPD/result_ma/IgG_OLS/Clumping_MA_loci_table_IgG_pred_10k_pearson07.xlsx")
