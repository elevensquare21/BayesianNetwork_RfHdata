##################################################################
##################################################################
#Apply Structural Equation Modeling from 'lavaan' package on RfH data
##################################################################
##################################################################

#load data
library(lavaan)
load("/Users/Selene/Desktop/BN/nv.rdata")
nv=nv[,-c(17,18)]

#run bnlearn first 
names(nv)=c("BMI","Age","CancerStage", "YearsDXRND","Neighborhood","Drink","Smoke","Insomnia","MB","Depression","Education","Sleep1","Sleep2","QOLp","QOLm","Arthritis","PA")
m=ncol(nv)-1
bl1=cbind.data.frame(names(nv)[-2],rep("Age",m))
names(bl1)=c('from','to')
bl2=cbind.data.frame(names(nv)[-11],rep("Education",m))
names(bl2)=c('from','to')
bl3=cbind.data.frame(names(nv)[-3],rep("CancerStage",m))
names(bl3)=c('from','to')
bl4=cbind.data.frame(names(nv)[-6],rep("Drink",m))
names(bl4)=c('from','to')
bl5=cbind.data.frame(names(nv)[-7],rep("Smoke",m))
names(bl5)=c('from','to')
bl6=cbind.data.frame(names(nv)[-4],rep("YearsDXRND",m))
names(bl6)=c('from','to')
bl7=cbind.data.frame(names(nv)[-5],rep("Neighborhood",m))
names(bl7)=c('from','to')
bl8=cbind.data.frame(rep("QOLp",m), names(nv)[-14])
names(bl8)=c('from','to')
bl9=cbind.data.frame(rep("QOLm",m), names(nv)[-15])
names(bl9)=c('from','to')
bl=rbind(bl1,bl2,bl3,bl4,bl5,bl6,bl7,bl8,bl9)

set.seed(3)
boot=boot.strength(data=nv, R=500, algorithm='hc',algorithm.args=list(score='aic-cg',blacklist=bl))
avg.boot=averaged.network(boot)
avg.boot$arcs


#build SEM model
#SEM has three components: measurement model, regressions, and residual correlations 
model <- '
  # measurement model
    Sleep.latent =~ Insomnia + Sleep1 + Sleep2
    Mental.latent =~ Depression + QOLm
    Physical.latent =~ Arthritis + QOLp + PA
  # regressions
    Mental.latent ~ Sleep.latent + BMI
    Physical.latent ~ Sleep.latent + BMI
  # residual correlations
    BMI ~~ QOLp + Insulin + CRP
    Age ~~ PA
    Smoke ~~ BMI
    Insomnia ~~ Depression + Sleep1
    Depression ~~ Sleep1 + Sleep2 + QOLm + Arthritis
    SLeep1 ~~ Sleep2
    Sleep2 ~~ QOLp + QOLm + PA
    Arthritis ~~ BMI + QOLp
    Insulin ~~ PA
    CancerStage ~~ CancerStage
    YearsDXRND ~~ YearsDXRND
    Neighborhood ~~ Neighborhood
    Drink ~~ Drink
    MB ~~ MB
    Education ~~ Education
'

#try only define the regressions 
model <- '
  # regressions
    BMI ~ QOLp + PA
    Age ~ PA
    Smoke ~ BMI
    Insomnia ~ Depression + Sleep1
    Depression ~ Sleep1 + Sleep2 + QOLm + Arthritis
    Sleep1 ~ Sleep2
    Sleep2 ~ QOLp + QOLm + PA
    Arthritis ~ BMI + QOLp
'
fit <- sem(model, data=nv)
summary(fit, standardized=TRUE, fit.measures=TRUE)


#move undirected arcs in equivalent class to corr 
directed.arcs(cpdag(avg.boot))
undirected.arcs(cpdag(avg.boot))
model <- '
  # regressions
    BMI ~ QOLp + PA
    Age ~ PA
    Smoke ~ BMI
    Sleep2 ~ QOLp + PA
    Arthritis ~ BMI + QOLp
  # residual correlations 
    Insomnia ~~ Depression + Sleep1
    Depression ~~ Sleep1 + Sleep2 + QOLm + Arthritis
    Sleep1 ~~ Sleep2
    Sleep2 ~~ QOLm
    
'
fit <- sem(model, data=nv)
summary(fit, standardized=TRUE, fit.measures=TRUE)

#try some latent variables
model <- '
  # measurement model
    Sleep.latent =~ Insomnia + Sleep1 + Sleep2
    Mental.latent =~ Depression + QOLm
    Physical.latent =~ Arthritis + QOLp + PA
  # regressions
    BMI ~ QOLp + PA
    Age ~ PA
    Smoke ~ BMI
    Insomnia ~ Depression + Sleep1
    Depression ~ Sleep1 + Sleep2 + QOLm + Arthritis
    Sleep1 ~ Sleep2
    Sleep2 ~ QOLp + QOLm + PA
    Arthritis ~ BMI + QOLp
  # residual correlations     
    Mental.latent ~~ Sleep.latent + BMI
    Physical.latent ~~ Sleep.latent + BMI
    
'
fit <- sem(model, data=nv)
summary(fit, standardized=TRUE, fit.measures=TRUE)






