library(remotes)
library(sampling)
library(survey)
library(nonprobsvy)
library(ggplot2)
library(pps)
library(compositions)
library(robustbase)
library(MatchIt)
library(dummies)
library(MASS)
library(plyr)
library(KWML)
library(lme4)
library(depmixS4)
library(snowfall)
library(pROC)

############################################################Parameter Description
# B: number of replications, set as 3000 in this study             
# N: size of finite population     
# n_A: size of the nonprobability sample
# n_B: size of the reference sample
# R2: the correlation parameter ρ 
# n_LC: the number of latent classes, with possible values ranging from 1 to 5
# icc: the intraclass correlation coefficient 
# n_Z: the number of covariates whose effects differ in sign across subpopulations
# n_P: the number of subpopulations whose covariate-effect directions are opposite to those of the remaining subpopulations for the n_Z covariates
# Ytype: distribution of the outcome, with possible values "gaussian" and "binomial"
# dist: type of distribution parameters, with possible values 1 (Distributions of Z_i vary systematically across subpopulations) and 2 (Distributions of Z_i vary irregularly)
# maxclassnum: the upper limit for the number of latent classes
#######################################################################################

simucode <- function(B,N,n_A,n_B,R2,n_LC,icc,n_Z,n_P,Ytype,dist,maxclassnum=4){
  
  #Define the metrics that need to be calculated
  RB_naive<-RB_sp<-RB_ipw<-RB_dr<-RB_newsp<-RB_kwipw<-RB_kwdr<-NULL
  MSE_naive<-MSE_sp<-MSE_ipw<-MSE_dr<-MSE_newsp<-MSE_kwipw<-MSE_kwdr<-NULL
  CP_naive<-CP_sp<-CP_ipw<-CP_dr<-CP_newsp<-CP_kwipw<-CP_kwdr<-NULL
  WIDTH_naive<-WIDTH_sp<-WIDTH_ipw<-WIDTH_dr<-WIDTH_newsp<-WIDTH_kwipw<-WIDTH_kwdr<-NULL
  
for (i_B in 1:B) {

    #Set the type of distribution parameters
    if(dist==1){
      PV <- 
        rbind(    
          data.frame(p_LCi=0.1,l_LCi=0, u_LCi=3, lamda_LCi=0.5,df_LCi=1.0, m_LCi=1, sd_LCi=1),   #1
          data.frame(p_LCi=0.2,l_LCi=5, u_LCi=8, lamda_LCi=1.0,df_LCi=2.0, m_LCi=5, sd_LCi=1),   #2
          data.frame(p_LCi=0.3,l_LCi=10,u_LCi=13,lamda_LCi=1.5,df_LCi=3.0, m_LCi=9, sd_LCi=1),   #3
          data.frame(p_LCi=0.4,l_LCi=15,u_LCi=18,lamda_LCi=2.0,df_LCi=4.0, m_LCi=13,sd_LCi=1),   #4
          data.frame(p_LCi=0.5,l_LCi=20,u_LCi=23,lamda_LCi=2.5,df_LCi=5.0, m_LCi=17,sd_LCi=1))   #5
    }
    
    if(dist==2){
      PV <- 
        rbind(    
          data.frame(p_LCi=0.4,l_LCi=15,u_LCi=18,lamda_LCi=1.5,df_LCi=3, m_LCi=5, sd_LCi=1),   #1
          data.frame(p_LCi=0.1,l_LCi=0, u_LCi=3, lamda_LCi=0.5,df_LCi=4, m_LCi=17,sd_LCi=1),   #2
          data.frame(p_LCi=0.3,l_LCi=5, u_LCi=8, lamda_LCi=2.0,df_LCi=5, m_LCi=1, sd_LCi=1),   #3
          data.frame(p_LCi=0.2,l_LCi=10,u_LCi=13,lamda_LCi=1.0,df_LCi=1, m_LCi=13,sd_LCi=1),   #4
          data.frame(p_LCi=0.5,l_LCi=20,u_LCi=23,lamda_LCi=2.5,df_LCi=2, m_LCi=9, sd_LCi=1))   #5
    }

    
    #Calculate the number of subgroups for each heterogeneous scenario----
    P_LCgroup <- 1/n_LC
    if(n_LC==1){
      N1<-N;n_A1<-n_A}
    if(n_LC==2){
      N1<-round(N*P_LCgroup,0);N2<-N-N1
      n_A1<-round(n_A*P_LCgroup,0);n_A2<-n_A-n_A1}
    if(n_LC==3){
      N1<-N2<-round(N*P_LCgroup,0);N3<-N-N1-N2
      n_A1<-n_A2<-round(n_A*P_LCgroup,0);n_A3<-n_A-n_A1-n_A2}
    if(n_LC==4){
      N1<-N2<-N3<-round(N*P_LCgroup,0);N4<-N-N1-N2-N3
      n_A1<-n_A2<-n_A3<-round(n_A*P_LCgroup,0);n_A4<-n_A-n_A1-n_A2-n_A3}
    if(n_LC==5){
      N1<-N2<-N3<-N4<-round(N*P_LCgroup,0);N5<-N-N1-N2-N3-N4
      n_A1<-n_A2<-n_A3<-n_A4<-round(n_A*P_LCgroup,0);n_A5<-n_A-n_A1-n_A2-n_A3-n_A4}
    
    
    S1 <- S2 <- S3 <- S4 <- S5 <- 1:5
    S12 <- S22 <- S32 <- S42 <- S52 <- 1:5
    
    #The coefficients of Z1-Z5 show positive or negative values ​​in different subgroups
    if(n_LC==1){
      PS <- matrix(1,nrow=n_LC,ncol=5) }
    if(n_LC==2 & n_Z==1 & n_P==1){
      PS <-
        matrix(c(-1,1,1,1,1,  #1
                 1,1,1,1,1)   #2
               ,nrow=2,byrow = T)    }
    if(n_LC==2 & n_Z==2 & n_P==1){
      PS <-
        matrix(c(-1,1,1,1,1,  #1
                 1,-1,1,1,1)   #2
               ,nrow=2,byrow = T)    }
    if(n_LC==2 & n_Z==3 & n_P==1){
      PS <-
        matrix(c(-1,1,-1,1,1,  #1
                 1,-1,1,1,1)   #2
               ,nrow=2,byrow = T)    }
    if(n_LC==2 & n_Z==4 & n_P==1){
      PS <-
        matrix(c(-1,1,-1,1,1,  #1
                 1,-1,1,-1,1)   #2
               ,nrow=2,byrow = T)    }
    if(n_LC==2 & n_Z==5 & n_P==1){
      PS <-
        matrix(c(-1,1,-1,1,-1,  #1
                 1,-1,1,-1,1)   #2
               ,nrow=2,byrow = T)    }
    if(n_LC==3 & n_Z==2 & n_P==1){
      PS <-
        matrix(c(-1,1,1,1,1,   #1
                 1,-1,1,1,1,   #2
                 1,1,1,1,1)   #3
               ,nrow=3,byrow = T)    }
    if(n_LC==3 & n_Z==3 & n_P==1){
      PS <-
        matrix(c(-1,1,1,1,1,   #1
                 1,-1,1,1,1,   #2
                 1,1,-1,1,1)   #3
               ,nrow=3,byrow = T)    }
    if(n_LC==3 & n_Z==4 & n_P==1){
      PS <-
        matrix(c(-1,1,1,-1,1,   #1
                 1,-1,1,1,1,   #2
                 1,1,-1,1,1)   #3
               ,nrow=3,byrow = T)    }
    if(n_LC==3 & n_Z==5 & n_P==1){
      PS <-
        matrix(c(-1,1,1,-1,1,   #1
                 1,-1,1,1,-1,   #2
                 1,1,-1,1,1)   #3
               ,nrow=3,byrow = T)    }
    if(n_LC==4 & n_Z==3 & n_P==1){
      PS <-
        matrix(c(-1,1,1,1,1,   #1
                 1,-1,1,1,1,   #2
                 1,1,-1,1,1,   #3
                 1,1,1,1,1)    #4
               ,nrow=4,byrow = T)    }
    if(n_LC==4 & n_Z==4 & n_P==1){
      PS <-
        matrix(c(-1,1,1,1,1,   #1
                 1,-1,1,1,1,   #2
                 1,1,-1,1,1,   #3
                 1,1,1,-1,1)   #4
               ,nrow=4,byrow = T)    }
    if(n_LC==4 & n_Z==5 & n_P==1){
      PS <-
        matrix(c(-1,1,1,1,-1,   #1
                 1,-1,1,1,1,   #2
                 1,1,-1,1,1,   #3
                 1,1,1,-1,1)   #4
               ,nrow=4,byrow = T)    }
    if(n_LC==4 & n_Z==3 & n_P==2){
      PS <-
        matrix(c(-1,1,1,1,1,   #1
                 -1,-1,1,1,1,   #2
                 1,-1,-1,1,1,   #3
                 1,1,-1,1,1)    #4
               ,nrow=4,byrow = T)    }
    if(n_LC==4 & n_Z==4 & n_P==2){
      PS <-
        matrix(c(-1,1,1,-1,1,   #1
                 -1,-1,1,1,1,   #2
                 1,-1,-1,1,1,   #3
                 1,1,-1,-1,1)   #4
               ,nrow=4,byrow = T)    }
    if(n_LC==4 & n_Z==5 & n_P==2){
      PS <-
        matrix(c(-1,1,1,-1,-1,   #1
                 -1,-1,1,1,-1,   #2
                 1,-1,-1,1,1,   #3
                 1,1,-1,-1,1)   #4
               ,nrow=4,byrow = T)    }
    if(n_LC==5 & n_Z==4 & n_P==1){
      PS <-
        matrix(c(-1,1,1,1,1,   #1
                 1,-1,1,1,1,   #2
                 1,1,-1,1,1,   #3
                 1,1,1,-1,1,   #4
                 1,1,1,1,1)    #5
               ,nrow=5,byrow = T)    }
    if(n_LC==5 & n_Z==5 & n_P==1){
      PS <-
        matrix(c(-1,1,1,1,1,   #1
                 1,-1,1,1,1,   #2
                 1,1,-1,1,1,   #3
                 1,1,1,-1,1,   #4
                 1,1,1,1,-1)    #5
               ,nrow=5,byrow = T)    }
    if(n_LC==5 & n_Z==4 & n_P==2){
      PS <-
        matrix(c(-1,1,1,1,1,   #1
                 -1,-1,1,1,1,   #2
                 1,-1,-1,1,1,   #3
                 1,1,-1,-1,1,   #4
                 1,1,1,-1,1)    #5
               ,nrow=5,byrow = T)    }
    if(n_LC==5 & n_Z==5 & n_P==2){
      PS <-
        matrix(c(-1,1,1,1,-1,   #1
                 -1,-1,1,1,1,   #2
                 1,-1,-1,1,1,   #3
                 1,1,-1,-1,1,   #4
                 1,1,1,-1,-1)   #5
               ,nrow=5,byrow = T)    }
      
    #The number of heterogeneous groups is [1]----------
    if(n_LC %in% c(1,2,3,4,5)){
      #The independent variable of the first latent subclass
      CVdata_LC1 <- CVdatagenerate(N_LCi=N1,p_LCi=PV[S1[1],1],l_LCi=PV[S2[1],2],u_LCi=PV[S2[1],3],lamda_LCi=PV[S3[1],4],df_LCi=PV[S4[1],5],m_LCi=PV[S5[1],6],sd_LCi=PV[S5[1],7])
      Z11 <- CVdata_LC1$X1;var_Z11 <- CVdata_LC1$var_X1;m_Z11 <- CVdata_LC1$m_X1
      Z12 <- CVdata_LC1$X2;var_Z12 <- CVdata_LC1$var_X2;m_Z12 <- CVdata_LC1$m_X2
      Z13 <- CVdata_LC1$X3;var_Z13 <- CVdata_LC1$var_X3;m_Z13 <- CVdata_LC1$m_X3
      Z14 <- CVdata_LC1$X4;var_Z14 <- CVdata_LC1$var_X4;m_Z14 <- CVdata_LC1$m_X4
      Z15 <- CVdata_LC1$X5;var_Z15 <- CVdata_LC1$var_X5;m_Z15 <- CVdata_LC1$m_X5
      var_Z1 <- (PS[S12[1],1]^2)*var_Z11+(PS[S22[1],2]^2)*var_Z12+(PS[S32[1],3]^2)*var_Z13+(PS[S42[1],4]^2)*var_Z14+(PS[S52[1],5]^2)*var_Z15
      var_e1 <- var_Z1*(1-R2)/R2
      miu1 <- PS[S12[1],1]*m_Z11+PS[S22[1],2]*m_Z12+PS[S32[1],3]*m_Z13+PS[S42[1],4]*m_Z14+PS[S52[1],5]*m_Z15
    }
    if(n_LC == 1){
      if(Ytype=='binomial'){
        Y1s <- PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15+rnorm(n=N1,mean=0,sd=sqrt(var_e1))
        Y1c <- Y1s-min(Y1s)+2;pi_Y1 <- inclusionprobabilities(Y1c, N1/2)
        Y1 <- UPpivotal(pi_Y1)
      }
      
      if(Ytype=='gaussian'){
        Y1 <- PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15+rnorm(n=N1,mean=0,sd=sqrt(var_e1))
      }
      
      #Sampling scores and non-probability sampling indicator variables
      linear_A1 <- PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15+rnorm(n=N1,mean=0,sd=sqrt(var_e1))
      linear_A1_offset <- linear_A1-min(linear_A1)+2
      pi_A1 <- inclusionprobabilities(linear_A1_offset, n_A1)
      ind_A1 <- UPpoisson(pi_A1)
    }
    
    #The number of heterogeneous groups is [2]----------
    if(n_LC %in% c(2,3,4,5)){
      #The independent variable of the second latent subclass
      CVdata_LC2 <- CVdatagenerate(N_LCi=N2,p_LCi=PV[S1[2],1],l_LCi=PV[S2[2],2],u_LCi=PV[S2[2],3],lamda_LCi=PV[S3[2],4],df_LCi=PV[S4[2],5],m_LCi=PV[S5[2],6],sd_LCi=PV[S5[2],7])
      Z21 <- CVdata_LC2$X1;var_Z21 <- CVdata_LC2$var_X1;m_Z21 <- CVdata_LC2$m_X1
      Z22 <- CVdata_LC2$X2;var_Z22 <- CVdata_LC2$var_X2;m_Z22 <- CVdata_LC2$m_X2
      Z23 <- CVdata_LC2$X3;var_Z23 <- CVdata_LC2$var_X3;m_Z23 <- CVdata_LC2$m_X3
      Z24 <- CVdata_LC2$X4;var_Z24 <- CVdata_LC2$var_X4;m_Z24 <- CVdata_LC2$m_X4
      Z25 <- CVdata_LC2$X5;var_Z25 <- CVdata_LC2$var_X5;m_Z25 <- CVdata_LC2$m_X5
      var_Z2 <- (PS[S12[2],1]^2)*var_Z21+(PS[S22[2],2]^2)*var_Z22+(PS[S32[2],3]^2)*var_Z23+(PS[S42[2],4]^2)*var_Z24+(PS[S52[2],5]^2)*var_Z25
      var_e2 <- var_Z2*(1-R2)/R2
      miu2 <- PS[S12[2],1]*m_Z21+PS[S22[2],2]*m_Z22+PS[S32[2],3]*m_Z23+PS[S42[2],4]*m_Z24+PS[S52[2],5]*m_Z25
      
    }
    if(n_LC == 2){
      #Within-group variance
      var_I <- 0.5*(var_Z1+var_e1)+0.5*(var_Z2+var_e2)
      var_B <- var_I/(1-icc)*icc
      a <- sqrt(var_B)
      
      b10 <- -1*a-miu1
      b20 <- a-miu2
      
      #Outcome/Y
      if(Ytype=='binomial'){
        Y1s <- PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15+rnorm(n=N1,mean=0,sd=sqrt(var_e1))+b10+5
        Y1c <- Y1s-min(Y1s)+2;pi_Y1 <- inclusionprobabilities(Y1c, N1/2)
        Y1 <- UPpivotal(pi_Y1)
        
        Y2s <- PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25+rnorm(n=N2,mean=0,sd=sqrt(var_e2))+b20+5
        Y2c <- Y1s-min(Y2s)+2;pi_Y2 <- inclusionprobabilities(Y2c, N2/2)
        Y2 <- UPpivotal(pi_Y2)
      }
      
      if(Ytype=='gaussian'){
        Y1 <- PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15+rnorm(n=N1,mean=0,sd=sqrt(var_e1))+b10+5
        Y2 <- PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25+rnorm(n=N2,mean=0,sd=sqrt(var_e2))+b20+5
      }
      
      #Sampling scores and non-probability sampling indicator variables
      linear_A1 <- PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15+rnorm(n=N1,mean=0,sd=sqrt(var_e1))+b10
      linear_A1_offset <- linear_A1-min(linear_A1)+2
      pi_A1 <- inclusionprobabilities(linear_A1_offset, n_A1)
      ind_A1 <- UPpoisson(pi_A1)

      linear_A2 <- PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25+rnorm(n=N2,mean=0,sd=sqrt(var_e2))+b20
      linear_A2_offset <- linear_A2-min(linear_A2)+2
      pi_A2 <- inclusionprobabilities(linear_A2_offset, n_A2)
      ind_A2 <- UPpoisson(pi_A2)
    }
      
    #The number of heterogeneous groups is [3]----------
    if(n_LC %in% c(3,4,5)){
      #The independent variable of the third latent subclass
      CVdata_LC3 <- CVdatagenerate(N_LCi=N3,p_LCi=PV[S1[3],1],l_LCi=PV[S2[3],2],u_LCi=PV[S2[3],3],lamda_LCi=PV[S3[3],4],df_LCi=PV[S4[3],5],m_LCi=PV[S5[3],6],sd_LCi=PV[S5[3],7])
      Z31 <- CVdata_LC3$X1;var_Z31 <- CVdata_LC3$var_X1;m_Z31 <- CVdata_LC3$m_X1
      Z32 <- CVdata_LC3$X2;var_Z32 <- CVdata_LC3$var_X2;m_Z32 <- CVdata_LC3$m_X2
      Z33 <- CVdata_LC3$X3;var_Z33 <- CVdata_LC3$var_X3;m_Z33 <- CVdata_LC3$m_X3
      Z34 <- CVdata_LC3$X4;var_Z34 <- CVdata_LC3$var_X4;m_Z34 <- CVdata_LC3$m_X4
      Z35 <- CVdata_LC3$X5;var_Z35 <- CVdata_LC3$var_X5;m_Z35 <- CVdata_LC3$m_X5
      var_Z3 <- (PS[S12[3],1]^2)*var_Z31+(PS[S22[3],2]^2)*var_Z32+(PS[S32[3],3]^2)*var_Z33+(PS[S42[3],4]^2)*var_Z34+(PS[S52[3],5]^2)*var_Z35
      var_e3 <- var_Z3*(1-R2)/R2
      miu3 <- PS[S12[3],1]*m_Z31+PS[S22[3],2]*m_Z32+PS[S32[3],3]*m_Z33+PS[S42[3],4]*m_Z34+PS[S52[3],5]*m_Z35
    }
    if(n_LC == 3){
      #Within-group variance
      var_I <- 1/3*(var_Z1+var_e1)+1/3*(var_Z2+var_e2)+1/3*(var_Z3+var_e3)
      var_B <- var_I/(1-icc)*icc
      a <- sqrt(3/2*var_B)
      
      b10 <- -1*a-miu1
      b20 <- -miu2
      b30 <- a-miu3
      
      #Outcome/Y
      if(Ytype=='binomial'){
        Y1s <- PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15+rnorm(n=N1,mean=0,sd=sqrt(var_e1))+b10+5
        Y1c <- Y1s-min(Y1s)+2;pi_Y1 <- inclusionprobabilities(Y1c, N1/2)
        Y1 <- UPpivotal(pi_Y1)
        
        Y2s <- PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25+rnorm(n=N2,mean=0,sd=sqrt(var_e2))+b20+5
        Y2c <- Y1s-min(Y2s)+2;pi_Y2 <- inclusionprobabilities(Y2c, N2/2)
        Y2 <- UPpivotal(pi_Y2)
        
        Y3s <- PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35+rnorm(n=N3,mean=0,sd=sqrt(var_e3))+b30+5
        Y3c <- Y3s-min(Y3s)+2;pi_Y3 <- inclusionprobabilities(Y3c, N3/2)
        Y3 <- UPpivotal(pi_Y3)
      }
      
      if(Ytype=='gaussian'){
        Y1 <- PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15+rnorm(n=N1,mean=0,sd=sqrt(var_e1))+b10+5
        Y2 <- PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25+rnorm(n=N2,mean=0,sd=sqrt(var_e2))+b20+5
        Y3 <- PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35+rnorm(n=N3,mean=0,sd=sqrt(var_e3))+b30+5
      }
      
      #Sampling scores and non-probability sampling indicator variables
      linear_A1 <- PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15+rnorm(n=N1,mean=0,sd=sqrt(var_e1))+b10
      linear_A1_offset <- linear_A1-min(linear_A1)+2
      pi_A1 <- inclusionprobabilities(linear_A1_offset, n_A1)
      ind_A1 <- UPpoisson(pi_A1)
      
      linear_A2 <- PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25+rnorm(n=N2,mean=0,sd=sqrt(var_e2))+b20
      linear_A2_offset <- linear_A2-min(linear_A2)+2
      pi_A2 <- inclusionprobabilities(linear_A2_offset, n_A2)
      ind_A2 <- UPpoisson(pi_A2)
      
      linear_A3 <- PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35+rnorm(n=N3,mean=0,sd=sqrt(var_e3))+b30
      linear_A3_offset <- linear_A3-min(linear_A3)+2
      pi_A3 <- inclusionprobabilities(linear_A3_offset, n_A3)
      ind_A3 <- UPpoisson(pi_A3)

    }

    #The number of heterogeneous groups is [4]----------
    if(n_LC %in% c(4,5)){
      #The independent variable of the fourth latent subclass
      CVdata_LC4 <- CVdatagenerate(N_LCi=N4,p_LCi=PV[S1[4],1],l_LCi=PV[S2[4],2],u_LCi=PV[S2[4],3],lamda_LCi=PV[S3[4],4],df_LCi=PV[S4[4],5],m_LCi=PV[S5[4],6],sd_LCi=PV[S5[4],7])
      Z41 <- CVdata_LC4$X1;var_Z41 <- CVdata_LC4$var_X1;m_Z41 <- CVdata_LC4$m_X1
      Z42 <- CVdata_LC4$X2;var_Z42 <- CVdata_LC4$var_X2;m_Z42 <- CVdata_LC4$m_X2
      Z43 <- CVdata_LC4$X3;var_Z43 <- CVdata_LC4$var_X3;m_Z43 <- CVdata_LC4$m_X3
      Z44 <- CVdata_LC4$X4;var_Z44 <- CVdata_LC4$var_X4;m_Z44 <- CVdata_LC4$m_X4
      Z45 <- CVdata_LC4$X5;var_Z45 <- CVdata_LC4$var_X5;m_Z45 <- CVdata_LC4$m_X5
      var_Z4 <- (PS[S12[4],1]^2)*var_Z41+(PS[S22[4],2]^2)*var_Z42+(PS[S32[4],3]^2)*var_Z43+(PS[S42[4],4]^2)*var_Z44+(PS[S52[4],5]^2)*var_Z45
      var_e4 <- var_Z4*(1-R2)/R2
      miu4 <- PS[S12[4],1]*m_Z41+PS[S22[4],2]*m_Z42+PS[S32[4],3]*m_Z43+PS[S42[4],4]*m_Z44+PS[S52[4],5]*m_Z45
    }
    if(n_LC == 4){
      #Within-group variance
      var_I <- 1/4*(var_Z1+var_e1)+1/4*(var_Z2+var_e2)+1/4*(var_Z3+var_e3)+1/4*(var_Z4+var_e4)
      var_B <- var_I/(1-icc)*icc
      a <- sqrt(2/5*var_B)
      
      b10 <- -2*a-miu1
      b20 <- -1*a-miu2
      b30 <- a-miu3
      b40 <- 2*a-miu4

      
      #Outcome/Y
      if(Ytype=='binomial'){
        Y1s <- PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15+rnorm(n=N1,mean=0,sd=sqrt(var_e1))+b10+5
        Y1c <- Y1s-min(Y1s)+2;pi_Y1 <- inclusionprobabilities(Y1c, N1/2)
        Y1 <- UPpivotal(pi_Y1)
        
        Y2s <- PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25+rnorm(n=N2,mean=0,sd=sqrt(var_e2))+b20+5
        Y2c <- Y1s-min(Y2s)+2;pi_Y2 <- inclusionprobabilities(Y2c, N2/2)
        Y2 <- UPpivotal(pi_Y2)
        
        Y3s <- PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35+rnorm(n=N3,mean=0,sd=sqrt(var_e3))+b30+5
        Y3c <- Y3s-min(Y3s)+2;pi_Y3 <- inclusionprobabilities(Y3c, N3/2)
        Y3 <- UPpivotal(pi_Y3)
        
        Y4s <- PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45+rnorm(n=N4,mean=0,sd=sqrt(var_e4))+b40+5
        Y4c <- Y4s-min(Y4s)+2;pi_Y4 <- inclusionprobabilities(Y4c, N4/2)
        Y4 <- UPpivotal(pi_Y4)
      }
      
      if(Ytype=='gaussian'){
        Y1 <- PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15+rnorm(n=N1,mean=0,sd=sqrt(var_e1))+b10+5
        Y2 <- PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25+rnorm(n=N2,mean=0,sd=sqrt(var_e2))+b20+5
        Y3 <- PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35+rnorm(n=N3,mean=0,sd=sqrt(var_e3))+b30+5
        Y4 <- PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45+rnorm(n=N4,mean=0,sd=sqrt(var_e4))+b40+5
      }
      
      #Sampling scores and non-probability sampling indicator variables
      linear_A1 <- PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15+rnorm(n=N1,mean=0,sd=sqrt(var_e1))+b10
      linear_A1_offset <- linear_A1-min(linear_A1)+2
      pi_A1 <- inclusionprobabilities(linear_A1_offset, n_A1)
      ind_A1 <- UPpoisson(pi_A1)
      
      linear_A2 <- PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25+rnorm(n=N2,mean=0,sd=sqrt(var_e2))+b20
      linear_A2_offset <- linear_A2-min(linear_A2)+2
      pi_A2 <- inclusionprobabilities(linear_A2_offset, n_A2)
      ind_A2 <- UPpoisson(pi_A2)
      
      linear_A3 <- PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35+rnorm(n=N3,mean=0,sd=sqrt(var_e3))+b30
      linear_A3_offset <- linear_A3-min(linear_A3)+2
      pi_A3 <- inclusionprobabilities(linear_A3_offset, n_A3)
      ind_A3 <- UPpoisson(pi_A3)
      
      linear_A4 <- PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45+rnorm(n=N4,mean=0,sd=sqrt(var_e4))+b40
      linear_A4_offset <- linear_A4-min(linear_A4)+2
      pi_A4 <- inclusionprobabilities(linear_A4_offset, n_A4)
      ind_A4 <- UPpoisson(pi_A4)

    }

    
    #The number of heterogeneous groups is [5]----------
    if(n_LC == 5){
      #The independent variable of the fifth latent subclass
      CVdata_LC5 <- CVdatagenerate(N_LCi=N5,p_LCi=PV[S1[5],1],l_LCi=PV[S2[5],2],u_LCi=PV[S2[5],3],lamda_LCi=PV[S3[5],4],df_LCi=PV[S4[5],5],m_LCi=PV[S5[5],6],sd_LCi=PV[S5[5],7])
      Z51 <- CVdata_LC5$X1;var_Z51 <- CVdata_LC5$var_X1;m_Z51 <- CVdata_LC5$m_X1
      Z52 <- CVdata_LC5$X2;var_Z52 <- CVdata_LC5$var_X2;m_Z52 <- CVdata_LC5$m_X2
      Z53 <- CVdata_LC5$X3;var_Z53 <- CVdata_LC5$var_X3;m_Z53 <- CVdata_LC5$m_X3
      Z54 <- CVdata_LC5$X4;var_Z54 <- CVdata_LC5$var_X4;m_Z54 <- CVdata_LC5$m_X4
      Z55 <- CVdata_LC5$X5;var_Z55 <- CVdata_LC5$var_X5;m_Z55 <- CVdata_LC5$m_X5
      var_Z5 <- (PS[S12[5],1]^2)*var_Z51+(PS[S22[5],2]^2)*var_Z52+(PS[S32[5],3]^2)*var_Z53+(PS[S42[5],4]^2)*var_Z54+(PS[S52[5],5]^2)*var_Z55
      var_e5 <- var_Z5*(1-R2)/R2
      miu5 <- PS[S12[5],1]*m_Z51+PS[S22[5],2]*m_Z52+PS[S32[5],3]*m_Z53+PS[S42[5],4]*m_Z54+PS[S52[5],5]*m_Z55
      
      #Within-group variance
      var_I <- 1/5*(var_Z1+var_e1)+1/5*(var_Z2+var_e2)+1/5*(var_Z3+var_e3)+1/5*(var_Z4+var_e4)+1/5*(var_Z5+var_e5)
      var_B <- var_I/(1-icc)*icc
      a <- sqrt(1/2*var_B)
      
      b10 <- -2*a-miu1
      b20 <- -1*a-miu2
      b30 <- -miu3
      b40 <- a-miu4
      b50 <- 2*a-miu5
      
      #Outcome/Y
      if(Ytype=='binomial'){
      Y1s <- PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15+rnorm(n=N1,mean=0,sd=sqrt(var_e1))+b10+5
      Y1c <- Y1s-min(Y1s)+2;pi_Y1 <- inclusionprobabilities(Y1c, N1/2)
      Y1 <- UPpivotal(pi_Y1)

      Y2s <- PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25+rnorm(n=N2,mean=0,sd=sqrt(var_e2))+b20+5
      Y2c <- Y1s-min(Y2s)+2;pi_Y2 <- inclusionprobabilities(Y2c, N2/2)
      Y2 <- UPpivotal(pi_Y2)
      
      Y3s <- PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35+rnorm(n=N3,mean=0,sd=sqrt(var_e3))+b30+5
      Y3c <- Y3s-min(Y3s)+2;pi_Y3 <- inclusionprobabilities(Y3c, N3/2)
      Y3 <- UPpivotal(pi_Y3)
      
      Y4s <- PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45+rnorm(n=N4,mean=0,sd=sqrt(var_e4))+b40+5
      Y4c <- Y4s-min(Y4s)+2;pi_Y4 <- inclusionprobabilities(Y4c, N4/2)
      Y4 <- UPpivotal(pi_Y4)
      
      Y5s <- PS[S12[5],1]*Z51+PS[S22[5],2]*Z52+PS[S32[5],3]*Z53+PS[S42[5],4]*Z54+PS[S52[5],5]*Z55+rnorm(n=N5,mean=0,sd=sqrt(var_e5))+b50+5
      Y5c <- Y5s-min(Y5s)+2;pi_Y5 <- inclusionprobabilities(Y5c, N5/2)
      Y5 <- UPpivotal(pi_Y5)
      }
      
      if(Ytype=='gaussian'){
        Y1 <- PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15+rnorm(n=N1,mean=0,sd=sqrt(var_e1))+b10+5
        Y2 <- PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25+rnorm(n=N2,mean=0,sd=sqrt(var_e2))+b20+5
        Y3 <- PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35+rnorm(n=N3,mean=0,sd=sqrt(var_e3))+b30+5
        Y4 <- PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45+rnorm(n=N4,mean=0,sd=sqrt(var_e4))+b40+5
        Y5 <- PS[S12[5],1]*Z51+PS[S22[5],2]*Z52+PS[S32[5],3]*Z53+PS[S42[5],4]*Z54+PS[S52[5],5]*Z55+rnorm(n=N5,mean=0,sd=sqrt(var_e5))+b50+5
      }
      
      #入样分数和非概率入样指示变量
      linear_A1 <- PS[S12[1],1]*Z11+PS[S22[1],2]*Z12+PS[S32[1],3]*Z13+PS[S42[1],4]*Z14+PS[S52[1],5]*Z15+rnorm(n=N1,mean=0,sd=sqrt(var_e1))+b10
      linear_A1_offset <- linear_A1-min(linear_A1)+2
      pi_A1 <- inclusionprobabilities(linear_A1_offset, n_A1)
      ind_A1 <- UPpoisson(pi_A1)
      
      linear_A2 <- PS[S12[2],1]*Z21+PS[S22[2],2]*Z22+PS[S32[2],3]*Z23+PS[S42[2],4]*Z24+PS[S52[2],5]*Z25+rnorm(n=N2,mean=0,sd=sqrt(var_e2))+b20
      linear_A2_offset <- linear_A2-min(linear_A2)+2
      pi_A2 <- inclusionprobabilities(linear_A2_offset, n_A2)
      ind_A2 <- UPpoisson(pi_A2)
      
      linear_A3 <- PS[S12[3],1]*Z31+PS[S22[3],2]*Z32+PS[S32[3],3]*Z33+PS[S42[3],4]*Z34+PS[S52[3],5]*Z35+rnorm(n=N3,mean=0,sd=sqrt(var_e3))+b30
      linear_A3_offset <- linear_A3-min(linear_A3)+2
      pi_A3 <- inclusionprobabilities(linear_A3_offset, n_A3)
      ind_A3 <- UPpoisson(pi_A3)
      
      linear_A4 <- PS[S12[4],1]*Z41+PS[S22[4],2]*Z42+PS[S32[4],3]*Z43+PS[S42[4],4]*Z44+PS[S52[4],5]*Z45+rnorm(n=N4,mean=0,sd=sqrt(var_e4))+b40
      linear_A4_offset <- linear_A4-min(linear_A4)+2
      pi_A4 <- inclusionprobabilities(linear_A4_offset, n_A4)
      ind_A4 <- UPpoisson(pi_A4)
      
      linear_A5 <- PS[S12[5],1]*Z51+PS[S22[5],2]*Z52+PS[S32[5],3]*Z53+PS[S42[5],4]*Z54+PS[S52[5],5]*Z55+rnorm(n=N5,mean=0,sd=sqrt(var_e5))+b50
      linear_A5_offset <- linear_A5-min(linear_A5)+2
      pi_A5 <- inclusionprobabilities(linear_A5_offset, n_A5)
      ind_A5 <- UPpoisson(pi_A5)

    }
    

    if(n_LC==1){
      Z1 <- c(Z11)
      Z2 <- c(Z12)
      Z3 <- c(Z13)
      Z4 <- c(Z14)
      Z5 <- c(Z15)
      Y <- c(Y1)
      pi_A <- c(pi_A1)
      ind_A <- c(ind_A1)}
    
    if(n_LC==2){
      Z1 <- c(Z11,Z21)
      Z2 <- c(Z12,Z22)
      Z3 <- c(Z13,Z23)
      Z4 <- c(Z14,Z24)
      Z5 <- c(Z15,Z25)
      Y <- c(Y1,Y2)
      pi_A <- c(pi_A1,pi_A2)
      ind_A <- c(ind_A1,ind_A2)}
    
    if(n_LC==3){
      Z1 <- c(Z11,Z21,Z31)
      Z2 <- c(Z12,Z22,Z32)
      Z3 <- c(Z13,Z23,Z33)
      Z4 <- c(Z14,Z24,Z34)
      Z5 <- c(Z15,Z25,Z35)
      Y <- c(Y1,Y2,Y3)
      pi_A <- c(pi_A1,pi_A2,pi_A3)
      ind_A <- c(ind_A1,ind_A2,ind_A3)}
    
    if(n_LC==4){
      Z1 <- c(Z11,Z21,Z31,Z41)
      Z2 <- c(Z12,Z22,Z32,Z42)
      Z3 <- c(Z13,Z23,Z33,Z43)
      Z4 <- c(Z14,Z24,Z34,Z44)
      Z5 <- c(Z15,Z25,Z35,Z45)
      Y <- c(Y1,Y2,Y3,Y4)
      pi_A <- c(pi_A1,pi_A2,pi_A3,pi_A4)
      ind_A <- c(ind_A1,ind_A2,ind_A3,ind_A4)}
    
    if(n_LC==5){
      Z1 <- c(Z11,Z21,Z31,Z41,Z51)
      Z2 <- c(Z12,Z22,Z32,Z42,Z52)
      Z3 <- c(Z13,Z23,Z33,Z43,Z53)
      Z4 <- c(Z14,Z24,Z34,Z44,Z54)
      Z5 <- c(Z15,Z25,Z35,Z45,Z55)
      Y <- c(Y1,Y2,Y3,Y4,Y5)
      pi_A <- c(pi_A1,pi_A2,pi_A3,pi_A4,pi_A5)
      ind_A <- c(ind_A1,ind_A2,ind_A3,ind_A4,ind_A5)}

      ################################################################################### PPS sampling: Probability sampling
      pi_B <- rep(n_B/N,times=n_B)
      sample_B <- sample(1:N,n_B,replace = TRUE)
      ind_B <- rep(0,times=N)
      ind_B[sample_B] <- 1
      
      #Generate overall data
      pop_data <- data.frame(Z1,Z2,Z3,Z4,Z5,Y,pi_A,pi_B,ind_A,ind_B)
    
      #Separate the non-probability sample and the reference sample, and then combine them into one large sample
      data_A0 <- pop_data[ind_A==1,]
      data_B0 <- pop_data[ind_B==1,]
      data_AB0 <- rbind(cbind(data_A0,group=1),cbind(data_B0,group=0))
      
      #The adjusted sampling weight variable d*
      a <- 1-nrow(data_A0)/sum(1/data_B0$pi_B)
      sweight <- 1/(data_AB0$pi_B)*a
      sweight[data_AB0$group==1] <- 1
      
      data_AB <- cbind(data_AB0,sweight)
      data_A <- data_AB[data_AB0$group==1,]
      data_B <- data_AB[data_AB0$group==0,]
    
      ################################################################################ Exisitng methods: DR, IPW, and SP-----
      
      #First, define the sampling design for probability sampling
      sample_B_svy <- svydesign(ids = ~1, 
                                probs = ~ pi_B, 
                                data = data_B)
      
      ################################################################################ DR -----
      fit_dr <- 
        nonprob(selection = ~ Z1 + Z2 + Z3 + Z4 + Z5,
                outcome = Y ~ Z1 + Z2 + Z3 + Z4 + Z5,
                data = data_A,
                svydesign = sample_B_svy,
                method_selection = "logit",
                method_outcome = "glm",
                family_outcome = Ytype)
      
      #Calculate the mean and 95% confidence interval of Y
      point_dr <- fit_dr$output$mean
      CI_dr <- fit_dr$confidence_interval
      
      ################################################################################ IPW -----
      
      fit_ipw <- 
        nonprob(selection = ~ Z1 + Z2 + Z3 + Z4 + Z5 ,
                target = ~ Y,
                data = data_A,
                svydesign = sample_B_svy,
                method_selection = "logit",
                control_selection = controlSel(est_method_sel = "gee", h = 2))
      
      #Calculate the mean and 95% confidence interval of Y
      point_ipw <- fit_ipw$output$mean
      CI_ipw <- fit_ipw$confidence_interval
      
      ############################################################################## SP -----
      fit_sp <- 
        nonprob(outcome = Y ~ Z1 + Z2 + Z3 + Z4 + Z5, 
                data = data_A, 
                svydesign = sample_B_svy, 
                method_outcome = "glm", 
                family_outcome = Ytype)
      
      #Calculate the mean and 95% confidence interval of Y
      point_sp <- fit_sp$output$mean
      CI_sp   <- fit_sp$confidence_interval
      
      ######################################################################################################### the proposed NHDR method-----

      
      Ind_out <- try(out <-NPinfer_ALC(data_NP=data_A,data_P=data_B,
                                       CIV=c('Z1','Z2','Z3','Z4','Z5'),type_CIV=c('binomial','gaussian','gaussian','gaussian','gaussian'),
                                       DV_NP='Y',DV_type=Ytype,
                                       weigths_NP='sweight',weights_P='sweight',maxclassnum=maxclassnum),silent=TRUE)
      T1 <- 0
      #To prevent the above model from failing to run, a while loop is set up to ensure it runs successfully
      while (methods::is(Ind_out,"try-error")==TRUE & T1<=10) {
        Ind_out <- try(out <-NPinfer_ALC(data_NP=data_A,data_P=data_B,
                                         CIV=c('Z1','Z2','Z3','Z4','Z5'),type_CIV=c('binomial','gaussian','gaussian','gaussian','gaussian'),
                                         DV_NP='Y',DV_type=Ytype,
                                         weigths_NP='sweight',weights_P='sweight',maxclassnum=maxclassnum),silent=TRUE)
        T1 <- T1+1}
      
      if(methods::is(Ind_out,"try-error")==FALSE){
        point_kwdr <- out$point_kwdr
        CI_kwdr <- out$CI_kwdr
      }
      
      if(methods::is(Ind_out,"try-error")==TRUE){
        point_kwdr <- NA
        CI_kwdr <- data.frame(lower_bound=NA,upper_bound=NA)
      }
      ###############################################################################Calculate the original point estimates and confidence intervals
      point_naive <- mean(data_A[,'Y'])
      se_naive <- sqrt(var(data_A[,'Y']))/nrow(data_A)
      CI_naive <- data.frame('lower_bound'=point_naive-1.96*se_naive,'upper_bound'=point_naive+1.96*se_naive)
      
      ###############################################################################Compute the RB, MSE, CP, and WIDTH for each method
      
      y_true <- mean(pop_data$Y)  #真实值
      
      #######################################################RB
      RB_naive[i_B] <- (point_naive-y_true)/y_true*100
      
      RB_dr[i_B] <- (point_dr-y_true)/y_true*100
      RB_ipw[i_B] <- (point_ipw-y_true)/y_true*100
      RB_sp[i_B] <- (point_sp-y_true)/y_true*100
      
      RB_kwdr[i_B] <- (point_kwdr-y_true)/y_true*100
      #######################################################MAD
      MAD_naive[i_B] <- abs(point_naive-y_true)
      
      MAD_dr[i_B] <- abs(point_dr-y_true)
      MAD_ipw[i_B] <- abs(point_ipw-y_true)
      MAD_sp[i_B] <- abs(point_sp-y_true)
      
      MAD_kwdr[i_B] <- abs(point_kwdr-y_true)
      #######################################################MSE
      MSE_naive[i_B] <- (point_naive-y_true)^2
      
      MSE_dr[i_B] <- (point_dr-y_true)^2
      MSE_ipw[i_B] <- (point_ipw-y_true)^2
      MSE_sp[i_B] <- (point_sp-y_true)^2
      
      MSE_kwdr[i_B] <- (point_kwdr-y_true)^2
      #######################################################CP
      CP_naive[i_B] <- y_true>=CI_naive$lower_bound       & y_true<=CI_naive$upper_bound
      
      CP_dr[i_B] <- y_true>=CI_dr$lower_bound       & y_true<=CI_dr$upper_bound
      CP_ipw[i_B] <- y_true>=CI_ipw$lower_bound     & y_true<=CI_ipw$upper_bound
      CP_sp[i_B] <- y_true>=CI_sp$lower_bound       & y_true<=CI_sp$upper_bound
      
      CP_kwdr[i_B] <- y_true>=CI_kwdr$lower_bound     & y_true<=CI_kwdr$upper_bound
      #######################################################WIDTH
      WIDTH_naive[i_B] <- CI_naive$upper_bound-CI_naive$lower_bound
      
      WIDTH_dr[i_B] <- CI_dr$upper_bound-CI_dr$lower_bound
      WIDTH_ipw[i_B] <- CI_ipw$upper_bound-CI_ipw$lower_bound
      WIDTH_sp[i_B] <- CI_sp$upper_bound-CI_sp$lower_bound
      
      WIDTH_kwdr[i_B] <- CI_kwdr$upper_bound-CI_kwdr$lower_bound
      
}#i_B in 1:B
  ###############################################################################Output the RB, MSE, CP, and WIDTH for each method
  RB_mean <- c(mean(RB_naive,na.rm=TRUE),
               mean(RB_dr,na.rm=TRUE),mean(RB_ipw,na.rm=TRUE),mean(RB_sp,na.rm=TRUE),
               mean(RB_kwdr,na.rm=TRUE))
  
  
  MSE <- c(mean(MSE_naive,na.rm=TRUE),
           mean(MSE_dr,na.rm=TRUE),mean(MSE_ipw,na.rm=TRUE),mean(MSE_sp,na.rm=TRUE),
           mean(MSE_kwdr,na.rm=TRUE))
  
  CP <- c(mean(CP_naive,na.rm=TRUE),
          mean(CP_dr,na.rm=TRUE),mean(CP_ipw,na.rm=TRUE),mean(CP_sp,na.rm=TRUE),
          mean(CP_kwdr,na.rm=TRUE))
  
  WIDTH_mean <- c(mean(WIDTH_naive,na.rm=TRUE),
                  mean(WIDTH_dr,na.rm=TRUE),mean(WIDTH_ipw,na.rm=TRUE),mean(WIDTH_sp,na.rm=TRUE),
                  mean(WIDTH_kwdr,na.rm=TRUE))
  
  name <- c('naive','dr','ipw','sp','kwdr')
  
  output <- cbind('B'=B,'N'=N,'f_A'=f_A,'f_B'=f_B,'n_A'=n_A,'n_B'=n_B,'Ytype'=Ytype,'n_LC'=n_LC,'icc'=icc,'n_Z'=n_Z,'n_P'=n_P,'dist'=dist,
                  'R2'=R2,'maxclassnum'=maxclassnum,'name'=name,
                  'RB_mean'=RB_mean,'MSE'=MSE,'CP'=CP,'WIDTH_mean'=WIDTH_mean)
      
      return(output)
  
}#simucode
