##################################################################
#Data preparation
##################################################################
library(bnlearn)
load("/Users/Selene/Desktop/BN/agg.mvpa.rdata")
load("/Users/Selene/Desktop/BN/agg.sed.rdata")
nv=read.csv('/Users/Selene/Desktop/RFH biomarkers Dody.csv',h=T)
nv$marker=as.character(nv$marker)
nvi=nv[nv$time==0 & nv$marker=='Insulin',]
nvc=nv[nv$time==0 & nv$marker=='CRP',]
nv=nvi
nv$crp=nvc$value
nv=nv[,c(1,5,7)]
names(nv)=c('identifier','Insulin','CRP')
nv$identifier=as.character(nv$identifier)
nv$Insulin=log(nv$Insulin)
nv$CRP=log(nv$CRP)
nv.bio=nv

nv=read.csv("/Users/Selene/Desktop/variables.csv",h=T)
nv$insulin=nv.bio$Insulin
nv$CRP=nv.bio$CRP
nv=nv[,-1]
nv$mvpa=agg.id$mvpa
nv$sed=agg.id.sed$sed
#nv$ERSTATUS=as.factor(nv$ERSTATUS)
#nv$PRSTATUS=as.factor(nv$PRSTATUS)
nv$nei=as.numeric(nv$nei)
#nv$nei=as.factor(nv$nei)
nv$m_b=as.numeric(nv$m_b)
nv=nv[,-c(3,4)]
nv$depression=as.character(nv$depression)
nv$depression=ifelse(nv$depression=="Don't Know", NA,nv$depression)
nv$depression=as.factor(nv$depression)
nv$arthritis=as.character(nv$arthritis)
nv$arthritis=ifelse(nv$arthritis=="Don't Know", NA,nv$arthritis)
nv$arthritis=as.factor(nv$arthritis)
#nv$BL.BMI=ifelse(nv$BL.BMI<30,0,1)
#nv$BL.BMI=as.factor(nv$BL.BMI)
nv=na.omit(nv)


##################################################################
#run bnlearn
################################################################## 
names(nv)=c("BMI","Age","CancerStage", "YearsDXRND","Neighborhood","Drink","Smoke","Insomnia","MB","Depression","Education","Sleep1","Sleep2","QOLp","QOLm","Arthritis","Insulin","CRP","MVPA","Sedentary")

#create a blacklist disallowing arrows to go certain direction based on prior knowledge
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

#use boostrapped method to get average network 
set.seed(7)
boot=boot.strength(data=nv, R=500, algorithm='hc',algorithm.args=list(score='aic-cg',blacklist=bl))
avg.boot=averaged.network(boot)
sig=averaged.network(boot)$learning$args$threshold
boot[(boot$strength>sig)&(boot$direction>0.5),]
avg.boot$arcs
vstructs(avg.boot)
bn.fit(avg.boot,data=nv)
cpdag(avg.boot)
cpdag(avg.boot)$arcs

##################################################################
#plot using igraph
##################################################################
library(igraph)
library('RColorBrewer')
avg.boot$arcs
nodes=data.frame(names(nv))
names(nodes)=c('variable')
nodes$type=c(1,3,4,4,3,3,3,2,4,2,3,2,2,1,2,1,5,5,1,1)
links=as.data.frame(avg.boot$arcs)
links$strength=boot[(boot$strength>sig)&(boot$direction>0.5),]$strength
links$direction=boot[(boot$strength>sig)&(boot$direction>0.5),]$direction
links$direction=ifelse(links$direction<=0.57,0.6,links$direction)
links$colind=ceiling((links$direction-0.5)*100/(50/8))
#links$sign=c(-1,1,1,1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1)
links$sign=rep(1,nrow(links))

net <- graph.data.frame(links, nodes, directed=T) 
colrs=c("darkolivegreen3","tomato","gold","gray60","pink")
V(net)$color=colrs[V(net)$type]
V(net)$frame.color="white"
E(net)$width = E(net)$strength*5
E(net)$arrow.size=1

l <- layout.fruchterman.reingold(net)
l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
pa1=brewer.pal(8,"Blues")
pa2=brewer.pal(8,"Reds")
plot(net,rescale=F, layout=l*1.25, edge.color=ifelse(links$sign>0,pa1[links$colind],pa2[links$colind]),vertex.label.font=2,vertex.label.color="gray30",vertex.label.cex=0.7) 
library("SDMTools")
pnts = cbind(x =c(1.7,1.9,1.9,1.7)+0.2, y =c(1.4,1.4,1,1))
legend.gradient(pnts,col=pa1,c("Low","High"),title="Positive Direction",cex=0.7,text.col="gray30")
pnts2 = cbind(x =c(1.7,1.9,1.9,1.7)+0.2, y =c(0.9,0.9,0.5,0.5)-0.2)
legend.gradient(pnts2,col=pa2,c("Low","High"),title="Negative Direction",cex=0.7)
legend(-2.75,1.5, legend=c("Physical Health","Mental Health", "Demographics & Health Behavior","Clinical Factors","Biomarkers"),col=c("darkolivegreen3","tomato","gold","gray60","pink"),fill=colrs,ce=0.7,box.col="white")
legend(1.88,-0.8,legend=c("Low",rep("",5),"High"),lwd=1:7,title="Edge Strength",col="grey",box.col="white",cex=0.7)


#plot equivalence class network
nodes=data.frame(names(nv))
names(nodes)=c('variable')
nodes$type=c(1,3,4,4,3,3,3,2,4,2,3,2,2,1,2,1,5,5,1,1)
links=as.data.frame(cpdag(avg.boot)$arcs)
links$mode=c(0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,2,2,0,2)
links$col=ifelse(links$mode>0,2,1)
net <- graph.data.frame(links, nodes, directed=T) 
colrs=c("darkolivegreen3","tomato","gold","gray60","pink")
V(net)$color=colrs[V(net)$type]
V(net)$frame.color="white"
E(net)$arrow.mode=E(net)$mode
colrs.e=c("gray70","black")
plot(net,edge.color=colrs.e[links$col],rescale=F, layout=l*1.25, vertex.label.font=2,vertex.label.color="gray30",vertex.label.cex=0.7)
legend(-2.75,1.5, legend=c("Physical Health","Mental Health", "Demographics & Health Behavior","Clinical Factors","Biomarkers"),col=c("darkolivegreen3","tomato","gold","gray60","pink"),fill=colrs,ce=0.7,box.col="white")
legend(1.8,1.4,legend=c("Directed Edge","Ambiguous Edge"),col=c("black","grey70"),box.col="white",cex=0.7,lty=1)


##################################################################
#post analysis on the learned Bayesian network 
##################################################################

#linear models from the network 
summary(lm(BMI~Smoke,data=nv))
summary(lm(Sleep1~Insomnia+Depression,data=nv))
summary(lm(Sleep2~Depression+Sleep1,data=nv))
summary(lm(QOLp~BMI+Sleep2+Arthritis,data=nv))
summary(lm(QOLm~Depression+Sleep2,data=nv))
summary(lm(Insulin~BMI,data=nv))
summary(lm(CRP~BMI,data=nv))
summary(lm(Sedentary~Arthritis+MVPA,data=nv))

fit=glm(Depression~Insomnia, data=nv, family = binomial())
summary(fit)
fit=glm(Arthritis~Depression, data=nv, family = binomial())
summary(fit)

#make table for regression coef
setwd('/Users/Selene/Desktop')
sig=averaged.network(boot)$learning$args$threshold
links=as.data.frame(avg.boot$arcs)
links$strength=boot[(boot$strength>sig)&(boot$direction>0.5),]$strength
links$direction=boot[(boot$strength>sig)&(boot$direction>0.5),]$direction
links$direction=ifelse(links$direction<=0.57,0.6,links$direction)
links=links[order(links$to),]
write.csv(links,file='links2.csv')

#conditional indep for insulin and crp
mb(avg.boot,'Insulin')
mb(avg.boot,'CRP')

#analyze performance sensitivity (in terms of aic, bic, loglik) to changes in network
score(avg.boot,data=nv,type='aic-cg') #-15644.71
score(avg.boot,data=nv,type='bic-cg')  #-15644.71
score(avg.boot,data=nv,type='loglik-cg') #-15509.83

dag=drop.arc(avg.boot,from="MVPA",to="Sedentary")
score(dag,data=nv,type='aic-cg') #-15672.86
score(dag,data=nv,type='bic-cg')  #-15672.86+15644.71
score(dag,data=nv,type='loglik-cg') #-15543.71

dag=drop.arc(avg.boot,from="Arthritis",to="Sedentary")
score(dag,data=nv,type='aic-cg') #-15646.87
score(dag,data=nv,type='bic-cg')  #-15646.87+15644.71
score(dag,data=nv,type='loglik-cg') #-15517.72

#bn inference using cpquery (this tests how one part of network reacts to 
#disturbances from another part of the network)
fitted=bn.fit(avg.boot,nv)
cpquery(fitted,event=(BMI>30),evidence=(PA>335))

result=cpdist(fitted,nodes=c('BMI'),evidence=(MVPA<13 & MVPA >11))
mean(result$BMI)
result=cpdist(fitted,nodes=c('BMI'),evidence=(MVPA<17 & MVPA >15))
mean(result$BMI)

result=cpdist(fitted,nodes=c('BMI'),evidence=(MVPA >(150/7)))
mean(result$BMI)
result=cpdist(fitted,nodes=c('BMI'),evidence=(MVPA<(50/7)))
mean(result$BMI)

result=cpdist(fitted,nodes=c('Insulin'),evidence=(MVPA<13 & MVPA >11 & BMI<32 & BMI>30))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(MVPA<17 & MVPA >15 & BMI<30 & BMI>29))
mean(result$Insulin)

result=cpdist(fitted,nodes=c('Insulin'),evidence=(MVPA>(150/7) &  BMI<30))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(MVPA<(50/7) &  BMI>=30))
mean(result$Insulin)

result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sedentary<420 &   BMI<30))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sedentary>480 &   BMI>=30))
mean(result$Insulin)

result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sedentary<420 & MVPA>(150/7) &  BMI<30))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sedentary>480 & MVPA<(50/7) &   BMI>=30))
mean(result$Insulin)

result=cpdist(fitted,nodes=c('CRP'),evidence=(MVPA<13 & MVPA >11 & BMI<32 & BMI>30))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(MVPA<17 & MVPA >15 & BMI<30 & BMI>29))
mean(result$CRP)

result=cpdist(fitted,nodes=c('CRP'),evidence=(MVPA>(150/7) &  BMI<30))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(MVPA<(50/7) &  BMI>=30))
mean(result$CRP)

result=cpdist(fitted,nodes=c('CRP'),evidence=(Sedentary<420 &   BMI<30))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sedentary>480 &   BMI>=30))
mean(result$CRP)

result=cpdist(fitted,nodes=c('CRP'),evidence=(Sedentary<420 & MVPA>(150/7) &  BMI<30))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sedentary>480 & MVPA<(50/7) &   BMI>=30))
mean(result$CRP)

result=cpdist(fitted,nodes=c('Insulin'),evidence=(BMI>=30))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(BMI<30))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(BMI>30 & BMI<32))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(BMI<30 & BMI>28))
mean(result$Insulin)

result=cpdist(fitted,nodes=c('CRP'),evidence=(BMI>=30))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(BMI<30))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(BMI>30 & BMI<32))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(BMI<30 & BMI>28))
mean(result$CRP)

result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2>52.9))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2<41.4))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2>52.9))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2<41.4))
mean(result$CRP)

result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2>52.9 & BMI>=30))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2<41.4 & BMI<30))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2>52.9 & BMI>30 & BMI<32))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2<41.4 & BMI<30 & BMI>28))
mean(result$CRP)

result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2>52.9 & BMI>=30))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2<41.4 & BMI<30))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2>52.9 & BMI>30 & BMI<32))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2<41.4 & BMI<30 & BMI>28))
mean(result$Insulin)

result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2>52.9 & MVPA<13 & MVPA >11))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2<41.4 & MVPA<17 & MVPA >15))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2>52.9 & MVPA<(50/7)))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2<41.4 & MVPA>(150/7)))
mean(result$CRP)

result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2>52.9 & Sedentary>480))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2<41.4 & Sedentary<420))
mean(result$CRP)

result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2>52.9 & MVPA<(50/7) & Sedentary>480))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2<41.4 & MVPA>(150/7) & Sedentary<420))
mean(result$CRP)

result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2>52.9 & MVPA<13 & MVPA >11))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2<41.4 & MVPA<17 & MVPA >15))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2>52.9 & MVPA<(50/7)))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2<41.4 & MVPA>(150/7)))
mean(result$Insulin)

result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2>52.9 & Sedentary>480))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2<41.4 & Sedentary<420))
mean(result$Insulin)

result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2>52.9 & MVPA<(50/7) & Sedentary>480))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2<41.4 & MVPA>(150/7) & Sedentary<420))
mean(result$Insulin)

result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2>52.9 & MVPA<13 & MVPA >11 & BMI>30 & BMI<32))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2<41.4 & MVPA<17 & MVPA >15 & BMI<30 & BMI>28))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2>52.9 & MVPA<(50/7) & BMI>=30))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2<41.4 & MVPA>(150/7) & BMI<30))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2>52.9 & Sedentary>480 & BMI>=30))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2<41.4 & Sedentary<420 & BMI<30))
mean(result$CRP)

result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2>52.9 & MVPA<13 & MVPA >11 & BMI>30 & BMI<32))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2<41.4 & MVPA<17 & MVPA >15 & BMI<30 & BMI>28))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2>52.9 & MVPA<(50/7) & BMI>=30))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2<41.4 & MVPA>(150/7) & BMI<30))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2>52.9 & Sedentary>480 & BMI>=30))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2<41.4 & Sedentary<420 & BMI<30))
mean(result$Insulin)

result=cpdist(fitted,nodes=c('Insulin'),evidence=(MVPA<(50/7) ))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(MVPA>(150/7) ))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('CRP'),evidence=(MVPA<(50/7) ))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(MVPA>(150/7) ))
mean(result$CRP)


result=cpdist(fitted,nodes=c('Insulin'),evidence=( Sedentary>480 ))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=( Sedentary<420 ))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('CRP'),evidence=( Sedentary>480 ))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=( Sedentary<420 ))
mean(result$CRP)

result=cpdist(fitted,nodes=c('Insulin'),evidence=( MVPA<(50/7) & Sedentary>480 ))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=( MVPA>(150/7) & Sedentary<420 ))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('CRP'),evidence=( MVPA<(50/7) & Sedentary>480 ))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=( MVPA>(150/7) & Sedentary<420 ))
mean(result$CRP)

result=cpdist(fitted,nodes=c('QOLm'),evidence=(Sleep2>52.9))
mean(result$QOLm)
result=cpdist(fitted,nodes=c('QOLm'),evidence=(Sleep2<41.4))
mean(result$QOLm)

result=cpdist(fitted,nodes=c('QOLp'),evidence=(Sleep2>52.9))
mean(result$QOLp)
result=cpdist(fitted,nodes=c('QOLp'),evidence=(Sleep2<41.4))
mean(result$QOLp)

#for sed to exceed 3rd quartile: 492, look at bmi, mvpa,age, sleep2
result=cpdist(fitted,nodes=c('Sedentary'),evidence=( MVPA<2.4 & BMI>30 & Age>63))
mean(result$Sedentary)

#for grant
result=cpdist(fitted,nodes=c('BMI'),evidence=(MVPA <=12))
mean(result$BMI)
result=cpdist(fitted,nodes=c('BMI'),evidence=(MVPA>=21))
mean(result$BMI)


result=cpdist(fitted,nodes=c('BMI'),evidence=(MVPA <=12 & Sleep2>52.9))
mean(result$BMI)
result=cpdist(fitted,nodes=c('BMI'),evidence=(MVPA>=21 &Sleep2<41.4))
mean(result$BMI)










