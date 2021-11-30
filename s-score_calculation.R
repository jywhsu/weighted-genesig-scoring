#=============================================
#STEP 0. LOAD LIBRARIES AND FORMAT DATA
#=============================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

#===========================================
#STEP 1. LOAD IN-HOUSE MADE FUNCTIONS
#===========================================
#load function for calculating weighted FDR score (based on paper) 
pm.score.calc = function(FDR){
  pm.score = c()
  p = c()
  for(i in 1:length(FDR)){
    if(-log(FDR[i],10) >= 4){
      p = 3
    }
    else if(-log(FDR[i],10) <= 1){
      p = 0
    }
    else(
      p = (-log(FDR[i],10)-1)
    )
    pm.score = rbind(pm.score,p)
  }
  return(pm.score)
}

SigScorePrep = function(casectl, testdata){
  m1 = data1 %>% 
    as.data.frame() %>%
    select(hgnc_symbol,FC,pm.score) %>% 
    distinct() %>% 
    filter(hgnc_symbol !="") %>% 
    na.omit()
  m2 = data2 %>% 
    as.data.frame() %>% 
    distinct() %>% 
    filter(hgnc_symbol!="") %>% 
    na.omit()
  m3 = left_join(m1,m2, by = "hgnc_symbol") %>% na.omit()
  data.output = as.data.frame(m3)
  return(m3)}

SigScoreCalc = function(data3){
  gene = data3$hgnc_symbol
  rho = data3$FC  
  rho = ifelse(rho<1, rho*-1,rho*1)
  p = data3$pm.score
  k = sum(abs(p))
  z = data3[,-(1:3)]
  s.score = c()
  for (i in 1:ncol(z)){
    s4 = 0
    for(j in 1:nrow(z)){
      s2 = rho[j]*p[j]*z[j,i]
      s4 = sum(s4,s2)
    }
    s1 = s4/k
    #print(paste(j,k,s4,s1))
    s.score = (rbind(s.score,s1))
  }
  data.output = cbind(colnames(z),s.score)
  colnames(data.output) = c("id","s.score")
  return(data.output)
}

#============================================
#STEP 2. IMPORT DATA
#============================================
#CASE-CTL DATA: 
#CASE = MYOTUBES, CTL = PAX7
#Note: this data frame must have  a column labeled "hgnc_symbol" with gene names
casectl.df = read.table("casectl_tube_pax7_edgeR_analysis.txt", header = T,sep = "\t") %>% 
  na.omit() 

#Z-SCORE CONVERTED CPM DATA: (example: SIX1 KD RNAseq CPM)
#Note: this data frame must have a column labeled "hgnc_symbol" with gene names
z.data = read.table("six1kd_rnaseq_z-score_cpm_data.txt", header = T, sep = "\t")

#====================================================================================
#STEP 3. CALCULATE FDR WEIGHT SCORE FOR CASE-CTL DATA;APPEND TO CASE-CTL DEG ANALYSIS
#====================================================================================
casectl.pm = casectl.df %>% 
  mutate(pm.score = pm.score.calc(FDR),FC = 2^logFC)

#===============================================
#STEP 4. CALCULATE S-SIGNATURE SCORE
#===============================================
x = SigScorePrep(casectl.pm,z.data)
s.scores = SigScoreCalc(x) 

