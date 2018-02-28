###########################################################
# FUNCTIONS FOR COMPUTING GENERALIZED HOTELLING T-SQUARED TEST #
# ########################################################

# GHotelling_022618
# Read the data and prepare analysis files for GHotelling
# The data are n rows for individuals and the following 4 columns.
# (mean distance to HMP1 from group A sample, mean distance to HMP2 from group A sample, mean distance
# to HMP1 from group B sample, mean distance to HMP2 from group B sample)
# Use NA for unknowns.  For this program the only permitted unknown patterns are e.g. (1,2,NA,NA) and (NA,NA,1,2)
# The program allows for correlations induced by data like (1,2,3,4) where one individual contributes data to both
# group A and group B.  If all individuals contribute only to either group A or group B but not both, then all the 
# data look like (1,2,NA,NA)  or (NA,NA,1,2)


##################################################################################################################
GHotelling_function=function(data,nBoot){
##################             create separate data sets:
#Data12A, Data12B, Data12AB for bivariate to Group A only, bivariate to Group B only, bivariate for both A and B
#Data1A, Data1B, Data1AB univariate to stool for A only, B only or A and B
#Data2A, Data2B, Data2AB univariate to nasal for A only, B only or A and B
Data12AB=subset(data,!is.na(data[,1])&!is.na(data[,2])&!is.na(data[,3])&!is.na(data[,4]))
Data12A=subset(data,!is.na(data[,1])&!is.na(data[,2])&(is.na(data[,3])|is.na(data[,4])))
Data12A=Data12A[,1:2]
Data12B=subset(data,!is.na(data[,3])&!is.na(data[,4])&(is.na(data[,1])|is.na(data[,2])))
Data12B=Data12B[,3:4]

##################################################################################################################
###############  Analysis if some individuals contribute to group A and group B

if (nrow(Data12AB)>0) {          #There are DAB data
  out1=T2forAB(Data12A,Data12B,Data12AB,T2only=FALSE)
  #print(out1)
  T2_12=unlist(out1[17])
  #T2_12
  T2_1=unlist(out1[18])
  #T2_1
  T2_2=unlist(out1[19])
  #T2_2
  ############  begin bootstrap
  #center the data to put it at the null
  #centeroverall
  mDA=apply(Data12A,2,mean)
  nDA=nrow(Data12A)
  mDB=apply(Data12B,2,mean)
  nDB=nrow(Data12B)
  mDAB=apply(Data12AB,2,mean)
  nDAB=nrow(Data12AB)
  mA=(nDA*mDA+nDAB*mDAB[1:2])/(nDA+nDAB)
  mB=(nDB*mDB+nDAB*mDAB[3:4])/(nDB+nDAB)
  WA=Data12A
  WA[,1]=WA[,1]-mA[1]
  WA[,2]=WA[,2]-mA[2]
  WB=Data12B
  WB[,1]=WB[,1]-mB[1]
  WB[,2]=WB[,2]-mB[2]
  WAB=Data12AB
  WAB[,1]=WAB[,1]-mA[1]
  WAB[,2]=WAB[,2]-mA[2]
  WAB[,3]=WAB[,3]-mB[1]
  WAB[,4]=WAB[,4]-mB[2]
  out2=T2forAB(WA,WB,WAB,T2only=FALSE)
  
  nA=nrow(Data12A)
  nB=nrow(Data12B)
  nAB=nrow(Data12AB)
  Count_12=0
  Count_1=0
  Count_2=0
  ###########  begin bootstrap
  for (i in 1:nBoot){
    indexA=1:nA
    indexA=sample(indexA,nA,replace=TRUE)
    bWA=WA[indexA,]
    indexB=1:nB
    indexB=sample(indexB,nB,replace=TRUE)
    bWB=WB[indexB,]
    indexAB=1:nAB
    indexAB=sample(indexAB,nAB,replace=TRUE)
    bWAB=WAB[indexAB,]
    
    bT2=T2forAB(bWA,bWB,bWAB,T2only=TRUE)
    if (bT2[1]>=T2_12) Count_12=Count_12+1  #for T2_12
    if (bT2[2]>=T2_1) Count_1=Count_1+1  #for T2_1
    if (bT2[3]>=T2_2) Count_2=Count_2+1  #for T2_2
  }
}

################  Analysis if individuals contribute to only one group, either A or B

if (nrow(Data12AB)==0) {          #There are no DAB data
  out1=T2withnoAB(Data12A,Data12B,T2only=FALSE)
  #print(out1)
  T2_12=unlist(out1[15])
  #T2_12
  T2_1=unlist(out1[16])
  #T2_1
  T2_2=unlist(out1[17])
  #T2_2
  #############  begin bootstrap
  #center the data to put it at the null
  #centeroverall
  mDA=apply(Data12A,2,mean)
  nDA=nrow(Data12A)
  mDB=apply(Data12B,2,mean)
  nDB=nrow(Data12B)
  #mDAB=apply(Data12AB,2,mean)
  #nDAB=nrow(Data12AB)
  mA=mDA #(nDA*mDA+nDAB*mDAB[1:2])/(nDA+nDAB)
  mB=mDB #(nDB*mDB+nDAB*mDAB[3:4])/(nDB+nDAB)
  WA=Data12A
  WA[,1]=WA[,1]-mA[1]
  WA[,2]=WA[,2]-mA[2]
  WB=Data12B
  WB[,1]=WB[,1]-mB[1]
  WB[,2]=WB[,2]-mB[2]
  # WAB=Data12AB
  # WAB[,1]=WAB[,1]-mA[1]
  # WAB[,2]=WAB[,2]-mA[2]
  # WAB[,3]=WAB[,3]-mB[1]
  # WAB[,4]=WAB[,4]-mB[2]
  out2=T2withnoAB(WA,WB,T2only=FALSE)
  
  nA=nrow(Data12A)
  nB=nrow(Data12B)
  #nAB=nrow(Data12AB)
  Count_12=0
  Count_1=0
  Count_2=0
  for (i in 1:nBoot){
    indexA=1:nA
    indexA=sample(indexA,nA,replace=TRUE)
    bWA=WA[indexA,]
    indexB=1:nB
    indexB=sample(indexB,nB,replace=TRUE)
    bWB=WB[indexB,]
    # indexAB=1:nAB
    # indexAB=sample(indexAB,nAB,replace=TRUE)
    # bWAB=WAB[indexAB,]
    
    bT2=T2withnoAB(bWA,bWB,T2only=TRUE)
    if (bT2[1]>=T2_12) Count_12=Count_12+1  #for T2AB
    if (bT2[2]>=T2_1) Count_1=Count_1+1  #for T2AB
    if (bT2[3]>=T2_2) Count_2=Count_2+1  #for T2AB
  }
}

#################  output of esimated p values and original analysis####################################
phat_12=Count_12/nBoot
#phat_12
phat_1=Count_1/nBoot
#phat_1
phat_2=Count_2/nBoot
# phat_2
# T2_12
# T2_1
# T2_2
# nBoot
# out1
result=list(phast_12=phat_12,phat_1=phat_1,phat_2=phat_2,T2_12=T2_12,T2_1=T2_1,T2_2=T2_2,nBoot=nBoot,out1=out1)
return(result)
}

########################################################################################################
###########  Function to use if some individuals contribute to both groups A and B
T2forAB=function(DA,DB,DAB,T2only){
#DA is nrow(DA) x 2 matrix of pairs of distances from people with group A data only 
#DB is nrow(DB) x 2 matrix of pairs of distances from people with group B data only 
#DAB is nrwo(DAB) x 2 matrix of two pairs for distances from people with group A data and group B dat  
# If T2only==TRUE, then only T2 is returned.  Otherwise, the list including T2 is returned  
#requires all three types of data, DA, DB, DAB
  mDA=apply(DA,2,mean)
  nDA=nrow(DA)
  mDB=apply(DB,2,mean)
  nDB=nrow(DB)
  mDAB=apply(DAB,2,mean)
  nDAB=nrow(DAB)
  mA=(nDA*mDA+nDAB*mDAB[1:2])/(nDA+nDAB)
  mB=(nDB*mDB+nDAB*mDAB[3:4])/(nDB+nDAB)
  xmean=c(mA,mB)
  cA=cov(rbind(DA,DAB[,1:2])) #the diagonal elements are estimated by using the DA and DAB elements for A
  cB=cov(rbind(DB,DAB[,3:4])) #the diagonal elements are estimated by using the DB and DAB elements for B
  cAB=cov(DAB)
  Cov=matrix(0,nrow=4,ncol=4)
  Cov[1:2,3:4]=cAB[1:2,3:4] #the off-diagnonal elements are estimated from DAB alone
  Cov[3:4,1:2]=t(Cov[1:2,3:4])
  Cov[1:2,1:2]=cA
  Cov[3:4,3:4]=cB
  Covmeans=matrix(0,nrow=4,ncol=4)
  Covmeans[1:2,1:2]=Cov[1:2,1:2]/(nDA+nDAB)
  Covmeans[3:4,3:4]=Cov[3:4,3:4]/(nDB+nDAB)
  Covmeans[1:2,3:4]=Cov[1:2,3:4]*(nDAB/((nDB+nDAB)*(nDA+nDAB)))
  Covmeans[3:4,1:2]=t(Covmeans[1:2,3:4])
  Con12=matrix(c(1, 0,0,1, -1, 0, 0, -1),nrow=2,ncol=4)  #contrast matrix for bivariate
  U12=Con12%*%xmean
  T2_12=t(U12)%*%solve(Con12%*%Covmeans%*%t(Con12))%*%U12
  Con1=matrix(c(1,0,-1,0),1,4)  # for A contrast
  U1=Con1%*%xmean
  T2_1=t(U1)%*%solve(Con1%*%Covmeans%*%t(Con1))%*%U1
  Con2=matrix(c(0,1,0,-1),1,4)  # for B contrast
  U2=Con2%*%xmean
  T2_2=t(U2)%*%solve(Con2%*%Covmeans%*%t(Con2))%*%U2
  result=list(mDA=mDA,nDA=nDA,mDB=mDB,nDB=nDB,mDAB=mDAB,nDAB=nDAB,ma=mA,mB=mB,Cov=Cov,Covmeans=Covmeans,
              U12=U12,U1=U1,U2=U2,Contrast12=Con12,Contrast1=Con1,Contrast2=Con2,T2_12=T2_12,T2_1=T2_1,T2_2=T2_2)
  if (T2only) result=c(T2_12,T2_1,T2_2)  # in simulations we want only the T2 values
  return(result)
}

#########################################################################################################
#########   Function to use if individuals contribute only to group A or to group B but not both
T2withnoAB=function(DA,DB,T2only){
  #DA is nrow(DA) x 2 matrix of pairs of distances from people with group A data only 
  #DB is nrow(DB) x 2 matrix of pairs of distances from people with group B data only 
  
  #No DAB data are provided.  Two independent samples from groups A and B
  #DAB is nrwo(DAB) x 2 matrix of two pairs for distances from people with group A data and group B dat  
  # If T2only==TRUE, then only T2 is returned.  Otherwise, the list including T2 is returned  
  #requires all three types of data, DA, DB, DAB
  mDA=apply(DA,2,mean)
  nDA=nrow(DA)
  mDB=apply(DB,2,mean)
  nDB=nrow(DB)
  #mDAB=apply(DAB,2,mean)
  #nDAB=nrow(DAB)
  mA=mDA  # (nDA*mDA+nDAB*mDAB[1:2])/(nDA+nDAB)
  mB=mDB  #(nDB*mDB+nDAB*mDAB[3:4])/(nDB+nDAB)
  xmean=c(mA,mB)
  cA=cov(DA) #the diagonal elements are estimated by using the DA and DAB elements for A
  cB=cov(DB) #the diagonal elements are estimated by using the DB and DAB elements for B
  #cAB=cov(DAB)
  Cov=matrix(0,nrow=4,ncol=4)
  #Cov[1:2,3:4]=cAB[1:2,3:4] #the off-diagnonal elements are estimated from DAB alone
  #Cov[3:4,1:2]=cAB[3:4,1:2]
  Cov[1:2,1:2]=cA
  Cov[3:4,3:4]=cB
  Covmeans=matrix(0,nrow=4,ncol=4)
  Covmeans[1:2,1:2]=Cov[1:2,1:2]/nDA
  Covmeans[3:4,3:4]=Cov[3:4,3:4]/nDB
  Covmeans[1:2,3:4]=Cov[1:2,3:4]
  Covmeans[3:4,1:2]=t(Covmeans[1:2,3:4])
  Con12=matrix(c(1, 0,0,1, -1, 0, 0, -1),nrow=2,ncol=4)  #contrast matrix for bivariate
  U12=Con12%*%xmean
  T2_12=t(U12)%*%solve(Con12%*%Covmeans%*%t(Con12))%*%U12
  Con1=matrix(c(1,0,-1,0),1,4)  # for contrast to ref site 1
  U1=Con1%*%xmean
  T2_1=t(U1)%*%solve(Con1%*%Covmeans%*%t(Con1))%*%U1
  Con2=matrix(c(0,1,0,-1),1,4)  # for contrast to ref site 2
  U2=Con2%*%xmean
  T2_2=t(U2)%*%solve(Con2%*%Covmeans%*%t(Con2))%*%U2
  result=list(mDA=mDA,nDA=nDA,mDB=mDB,nDB=nDB,mA=mA,mB=mB,Cov=Cov,Covmeans=Covmeans,U12=U12,U1=U1,U2=U2,
              Contrast12=Con12,Contrast1=Con1,Contrast2=Con2,T2_12=T2_12,T2_1=T2_1,T2_2=T2_2)
  if (T2only) result=c(T2_12,T2_1,T2_2)  #return for simulations
  return(result)
}

