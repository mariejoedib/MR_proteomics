##load packages
library(data.table)
library(TwoSampleMR)
library(MendelianRandomization)
library(tidyverse)

##set working directory
setwd("/Users/marie-joe/Desktop/HF Proteomics/deCODE/")

#pre-allocate data frames
results<-data.frame(
  exposure=character(),outcome=character(),snps=numeric(),
  ivw.OR=numeric(),ivw.l=numeric(),ivw.u=numeric(),ivw.p=numeric(),
  wm.OR=numeric(),wm.l=numeric(),wm.u=numeric(),wm.p=numeric(),
  egger.OR=numeric(),egger.l=numeric(),egger.u=numeric(),egger.p=numeric(),
  Egger.intercept.p=numeric(),Q.p=numeric(),I2=numeric(),F.statistic=numeric(),
  stringsAsFactors = F
)

#generate file for storing results (in case first iteration has no significant SNPs - not really needed in this analysis)
write.table(results,file=paste0("/Users/marie-joe/Desktop/HF Proteomics/deCODE/results_deCODE_2SMR_HF_cis.tsv"),
            sep="\t",col.names=T,row.names=F,append=F,quote=F)

##load exposure
#This file has already been filtered to exclude SNPs with Pvalues > 5e-8 
pQTLs <- read.table("/Users/marie-joe/Desktop/HF Proteomics/deCODE/cis-pQTLs.tsv", header = T)

##select only the aptamers that were significant in the DHFA PHFS study
PHFS <- read.csv("/Users/marie-joe/Desktop/HF Proteomics/DHFA_significant_PHFS.csv", header = T)
PHFS <- PHFS[ ,-c(4:9)]
df<- pQTLs[pQTLs$aptamer_id.coords %in% PHFS$aptamer_id, ]
df$SNPID <- str_c(df$Chrom,df$Pos,sep="_")

##create a vector with the names of each exposure's Somascan ID
ids <- df$aptamer_id.coords
id <- unique(ids)

##load outcome GWAS: Heart Failure
outcome <- read.table("/Users/marie-joe/Desktop/HF Proteomics/HERMES_Jan2019_HeartFailure_summary_data.txt", header=TRUE)

##loop exposures: each aptamer_id

for (i in 1:length(id)){
  
  print(paste0("Starting iterations for ",id[i]))
  
  df<-pQTLs[pQTLs$aptamer_id.coords %in% id[i], ]
  
  df$SNPID <- str_c(df$Chrom,df$Pos,sep="_")
  
  #convert exposure to "TwoSampleMR" format
  exp<-format_data(
    df,
    type= "exposure",
    snps = NULL,
    snp_col = "rsids",
    beta_col = "Beta",
    se_col = "SE",
    effect_allele_col = "effectAllele",
    other_allele_col = "otherAllele",
    eaf_col = "effectAlleleFreq",
    pval_col = "Pval"
  )
  
  #name exposure correctly
  exp$exposure<-id[i]
  
  #keep unique SNPs
  outcome<- outcome[!duplicated(outcome$SNP),]
  #get SNPs from exposure data
  outc<-outcome[outcome$SNP %in% df$rsids, ]
  
  #If there are no SNPs available in outcome, move to next iteration
  if(nrow(outc)==0){
    
    print("no SNPs available")
    
    next
    
  }else{
    
    ##Convert outcome to "TwoSample" Format
    
    outc <- format_data(
      dat=outc,
      type="outcome",
      snps=outc$SNP,
      snp_col = "SNP",
      beta_col="b",
      se_col="se",
      eaf_col = "freq",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      pval_col = "p",
      chr_col = "CHR",
    )
    
    ##name outcome
    outc$outcome <- "HF"
    
    #Harmonise data
    dat<-harmonise_data(exp,outc,action=1)
    
    #Clump at r2=0.01
    dat <-clump_data(
      dat,
      clump_r2 = 0.01,
      pop = "EUR"
    )
    
    dat<- dat %>%
      filter(mr_keep == "TRUE")
    
    if(nrow(dat)==0){
      
      print("no SNPs available after harmonisation")
      
      next
      
    }else{
      
      #If there's only one independent SNP remaining
      if(nrow(dat)==1){
        
        #F statistic using chi-squared approximation
        dat$Fstat<-dat$beta.exposure^2/dat$se.exposure^2
        
        if(dat$mr_keep == FALSE){
          
          print("not amenable to MR")
          
          next
          
        }
        
        res<-mr(dat)
        
        #assign values in table
        results[1,1]<- id[i] #exposure name
        results[1,2]<-"HF" #outcome name
        results[1,3]<-res[1,6]#N SNPs
        
        results[1,4]<-exp(res[1,7])#Wald ratio OR
        results[1,5]<-exp(res[1,7]+qnorm(.025)*res[1,8])#Wald ratio lower 95%
        results[1,6]<-exp(res[1,7]+qnorm(.975)*res[1,8])#Wald ratio upper 95%
        results[1,7]<-res[1,9]#Wald ratio p-value
        
        
        results[1,8:18]<-NA
        
        results[1,19]<-mean(dat$Fstat)#mean F statistic
        
        #write results
        write.table(results,file=paste0("/Users/marie-joe/Desktop/HF Proteomics/deCODE/results_deCODE_2SMR_HF_cis.tsv"),
                    sep="\t",col.names=F,row.names=F,append=T,quote=F)
        
        
        print(paste("Analysis for HF Outcome with exposure",id[i],"finished"))
        
      }else{ #if there are two or more independent SNPs
        
        #Perform MR (IVW, WM, Egger)
        #F statistic using chi-squared approximation
        dat$Fstat<-dat$beta.exposure^2/dat$se.exposure^2
        
        #We want to account for LD so we need to use MendelianRandomization package
        #Convert to the MRInput object and obtain LD matrix for instruments
        dat2 <- dat_to_MRInput(dat, get_correlation = TRUE)
        res <- MendelianRandomization::mr_ivw(dat2[[1]], correl = TRUE)
        
        res2<-mr(dat,method_list=c("mr_egger_regression","mr_weighted_median"))
        
        #MR-Egger intercept
        pleio<-mr_pleiotropy_test(dat)
        
        #heterogeneity test
        het<-mr_heterogeneity(dat)
        
        het$I2<-ifelse(
          (het$Q-het$Q_df)/het$Q<0,
          0,
          100*(het$Q-het$Q_df)/het$Q
        )
        
        #assign values in table
        results[1,1]<-id[i] #exposure name
        results[1,2]<-"HF" #outcome name
        results[1,3]<-res@SNPs#N SNPs
        
        results[1,4]<-exp(res@Estimate)#IVW OR
        results[1,5]<-exp(res@CILower)#IVW lower 95%
        results[1,6]<-exp(res@CIUpper)#IVW upper 95%
        results[1,7]<-res@Pvalue#IVW p-value
        
        results[1,8]<-exp(res2[2,7])#WM OR
        results[1,9]<-exp(res2[2,7]+qnorm(.025)*res2[2,8])#WM lower 95%
        results[1,10]<-exp(res2[2,7]+qnorm(.975)*res2[2,8])#WM upper 95%
        results[1,11]<-res2[2,9]#WM p-value
        
        results[1,12]<-exp(res2[1,7])#Egger OR
        results[1,13]<-exp(res2[1,7]+qnorm(.025)*res2[1,8])#Egger lower 95%
        results[1,14]<-exp(res2[1,7]+qnorm(.975)*res2[1,8])#Egger upper 95%
        results[1,15]<-res2[1,9]#Egger p-value
        
        results[1,16]<-pleio[1,7]#Egger intercept p
        results[1,17]<-het[2,8]#IVW Q statistic p
        results[1,18]<-het[2,9]# I2 statistic
        results[1,19]<-mean(dat$Fstat)#mean F statistic
        
        #write results
        write.table(results,file=paste0("/Users/marie-joe/Desktop/HF Proteomics/deCODE/results_deCODE_2SMR_HF_cis.tsv"),
                    sep="\t",col.names=F,row.names=F,append=T,quote=F)
        
        print(paste("Analysis for HF Outcome",i,"finished"))
        
      }
    }
  }
}

##edit table
#load table
df<-read.table(file=paste0("/Users/marie-joe/Desktop/HF Proteomics/deCODE/results_deCODE_2SMR_HF_cis.tsv"),
               sep="\t",header=T,stringsAsFactors=F)


#add FDR correction
df$p.ajd <- p.adjust(df$ivw.p, method = "fdr")

#annotate
##load aptamer_ids of those significant for HF in PHFS
coords <- read.csv("/Users/marie-joe/Desktop/HF Proteomics/DHFA_significant_PHFS.csv", header = T)

##join dataframes
df1 <- right_join(
  coords,
  df,
  by = c("aptamer_id"="exposure"),
  copy = FALSE,
  suffix = c(".coords", ".gwas"),
  keep = FALSE,
  na_matches = c("na", "never")
)

#order by increasing p-value
df1<-df1[order(df1$ivw.p),]

#write results
write.table(df1,file=paste0("/Users/marie-joe/Desktop/HF Proteomics/deCODE/results_deCODE_2SMR_HF_cis.tsv"),
            sep="\t",row.names=F,quote=F)
