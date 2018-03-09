##################################################################
##################################################################
#Apply Bayesian Network model from 'bnlearn' package on RfH data
##################################################################
##################################################################

#get appropriate physical activity data (mvpa, sedentary time, and average activity)
library(bnlearn)
load("/Users/Selene/Desktop/Updated Missing Data/extended data with activity and date.rdata")
#name is edata
edata=data.e
load("/Users/Selene/Desktop/MVPA/aggregate data.rdata")
#name is agg
agg.d=agg
agg.d$wt=1440-agg.d$miss
agg.d$mvpa.adj=agg.d$mvpa/agg.d$wt*720
agg.id=aggregate(agg.d$mvpa.adj,list(agg.d$identifier),mean,na.rm=TRUE)
names(agg.id)=c('identifier','mvpa')
save(agg.id,file='agg.mvpa.rdata')

sed.0=list()
for(i in 1:length(unique(agg$identifier))){
	x=edata[edata$identifier==unique(agg$identifier)[i],]
	s=rep(0,length(unique(x$dt)))
	for(j in 1:length(unique(x$dt))){
		s[j]=length(x[x$dt==unique(x$dt)[j] & x$activity<100,3])
	}
	sed.0[[i]]=s
}
sed.0=unlist(sed.0)
agg$sed.0=sed.0
agg$sed=agg$sed.0-agg$miss
agg$wt=1440-agg$miss
agg$sed.adj=agg$sed/agg$wt*720
agg.id.sed=aggregate(agg$sed.adj,list(agg$identifier),mean,na.rm=TRUE)
names(agg.id.sed)=c('identifier','sed')
save(agg.id.sed,file='agg.sed.rdata')

agg=aggregate(edata$activity,list(edata$identifier),mean,na.rm=TRUE)
names(agg)=c('identifier','PA')

#get biomarkers at baseline
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

#get other variables
nv=read.csv("/Users/Selene/Desktop/variables.csv",h=T)
nv$insulin=nv.bio$Insulin
nv$CRP=nv.bio$CRP

#get treatment type and estrogen etc.
nv=read.csv("/Users/Selene/Desktop/Reach for Health.csv",h=T)
nv=nv[,c(11,18,19)]
nv[[1]]=as.character(nv[[1]])
li=list()
for (i in 1:nrow(nv)){
	a=nv[i,1]
	li[[i]]=as.numeric(unlist(strsplit(gsub('\\|','',a),'')))
}
lis=unlist(li)
sum(lis==1)/333
sum(lis==2)/333
sum(lis==3)/333
sum(lis==8)/333
sum(lis==9)/333
sum(nv[[2]]==1)/333
sum(nv[[3]]==1)/333

nv=nv[,-1]
nv$PA=agg$PA
nv$ERSTATUS=as.factor(nv$ERSTATUS)
nv$PRSTATUS=as.factor(nv$PRSTATUS)
nv$nei=as.numeric(nv$nei)
nv$m_b=as.numeric(nv$m_b)
nv=nv[,-c(3,4)]
nv$depression=as.character(nv$depression)
nv$depression=ifelse(nv$depression=="Don't Know", NA,nv$depression)
nv$depression=as.factor(nv$depression)
nv$arthritis=as.character(nv$arthritis)
nv$arthritis=ifelse(nv$arthritis=="Don't Know", NA,nv$arthritis)
nv$arthritis=as.factor(nv$arthritis)
nv=na.omit(nv)
setwd('/Users/Selene/Desktop/BN')
save(nv,file='nv.rdata')


##################################################################
#run bnlearn
################################################################## 
load("/Users/Selene/Desktop/BN/nv.rdata")
names(nv)=c("BMI","Age","CancerStage", "YearsDXRND","Neighborhood","Drink","Smoke","Insomnia","MB","Depression","Education","Sleep1","Sleep2","QOLp","QOLm","Arthritis","Insulin","CRP","PA")

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
set.seed(3)
boot=boot.strength(data=nv, R=500, algorithm='hc',algorithm.args=list(score='aic-cg',blacklist=bl))
avg.boot=averaged.network(boot)
sig=averaged.network(boot)$learning$args$threshold
boot[(boot$strength>sig)&(boot$direction>0.5),]
avg.boot$arcs
vstructs(avg.boot)
bn.fit(avg.boot,data=nv)

#can adjust threshold for averaged network if needed 
boot[(boot$strength>0.85)&(boot$direction>0.5),]
avg.boot=averaged.network(boot,threshold=0.85)

#find equivalence class network
cpdag(avg.boot)
cpdag(avg.boot)$arcs
directed.arcs(cpdag(avg.boot))


##################################################################
#analyze performance sensitivity (in terms of aic, bic, loglik) to changes in network
##################################################################
score(avg.boot,data=nv,type='aic-cg') #-14483.52
score(avg.boot,data=nv,type='bic-cg')  #-14483.52
score(avg.boot,data=nv,type='loglik-cg') #-14345.76

#isolate sleep1, sleep2 or both, compare BIC scores
mat=as.data.frame(avg.boot$arcs)
mat1=mat[mat$from!="Sleep1",]
mat1=mat1[mat1$to!="Sleep1",]
dag = empty.graph(names(nv))
arcs(dag) = as.matrix(mat1)
score(dag,data=nv,type='aic-cg') #-14637.4
score(dag,data=nv,type='bic-cg')  #-14637.4+14483.52
score(dag,data=nv,type='loglik-cg') #-14513.99

mat=as.data.frame(avg.boot$arcs)
mat1=mat[mat$from!="Sleep2",]
mat1=mat1[mat1$to!="Sleep2",]
dag = empty.graph(names(nv))
arcs(dag) = as.matrix(mat1)
score(dag,data=nv,type='aic-cg') #-14673.86
score(dag,data=nv,type='bic-cg')  #-14673.86+14483.52
score(dag,data=nv,type='loglik-cg') #-14559.07

mat=as.data.frame(avg.boot$arcs)
mat1=mat[mat$from!="Sleep2",]
mat1=mat1[mat1$to!="Sleep2",]
mat1=mat1[mat1$to!="Sleep1",]
mat1=mat1[mat1$from!="Sleep1",]
dag = empty.graph(names(nv))
arcs(dag) = as.matrix(mat1)
score(dag,data=nv,type='aic-cg') #-14739.17
score(dag,data=nv,type='bic-cg')  #-14739.17+14483.52
score(dag,data=nv,type='loglik-cg') #-14632.98

#drop edge from bmi to biomarkers & CRP
dag=drop.arc(avg.boot,from="BMI",to="Insulin")
score(dag,data=nv,type='aic-cg') #-14503.2
score(dag,data=nv,type='bic-cg')  #-14503.2+14483.52
score(dag,data=nv,type='loglik-cg') #-14368.31
dag=drop.arc(avg.boot,from="BMI",to="CRP")
score(dag,data=nv,type='aic-cg') #-14509.8
score(dag,data=nv,type='bic-cg')  #-14509.8+14483.52
score(dag,data=nv,type='loglik-cg') #-14374.92

#more
dag=drop.arc(avg.boot,from="Depression",to="Arthritis")
score(dag,data=nv,type='aic-cg') #-14485.08
score(dag,data=nv,type='bic-cg')  #-14485.08+14483.52=-1.56
score(dag,data=nv,type='loglik-cg') #-14350.19

dag=drop.arc(avg.boot,from="Arthritis",to="BMI")
score(dag,data=nv,type='aic-cg') #-14483.61
score(dag,data=nv,type='bic-cg')  #-14483.61+14483.52=-0.09
score(dag,data=nv,type='loglik-cg') #-14351.6

dag=drop.arc(avg.boot,from="Smoke",to="BMI")
score(dag,data=nv,type='aic-cg') #-14488.52
score(dag,data=nv,type='bic-cg')  #-14488.52+14483.52=-5
score(dag,data=nv,type='loglik-cg') #-14356.5

#try random start
nodes=names(nv)
start=random.graph(nodes=nodes, method='ic-dag',num=500)
netlist=lapply(start,function(net) {hc(nv,score='aic-cg',start=net,blacklist=bl)})
rnd=custom.strength(netlist,nodes=nodes)


##################################################################
#plot using igraph
##################################################################
library(igraph)
library('RColorBrewer')
avg.boot$arcs
nodes=data.frame(names(nv))
names(nodes)=c('variable')
nodes$type=c(1,3,4,4,3,3,3,2,4,2,3,2,2,1,2,1,5,5,1)
links=as.data.frame(avg.boot$arcs)
links$strength=boot[(boot$strength>sig)&(boot$direction>0.5),]$strength
links$direction=boot[(boot$strength>sig)&(boot$direction>0.5),]$direction
links$direction=ifelse(links$direction<=0.57,0.6,links$direction)
links$colind=ceiling((links$direction-0.5)*100/(50/8))
links$sign=c(-1,1,1,-1,1,    1,1,1,1,-1,     1,1,-1,-1,-1,     1,-1,-1)

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
nodes$type=c(1,3,4,4,3,3,3,2,4,2,3,2,2,1,2,1,5,5,1)
links=as.data.frame(cpdag(avg.boot)$arcs)
links$mode=c(2,2,2,2,2,   0,0,0,0,0,    0,0,0,0,0,   0,0,2,0,2,   0,0,2,0,2,2)
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
summary(lm(BMI~Smoke+Arthritis,data=nv))
summary(lm(Sleep1~Insomnia+Depression,data=nv))
summary(lm(Sleep2~Depression+Sleep1,data=nv))
summary(lm(QOLp~BMI+Sleep2+Arthritis,data=nv))
summary(lm(QOLm~Depression+Sleep2,data=nv))
summary(lm(PA~Age+Sleep2+Insulin,data=nv))
summary(lm(Insulin~BMI,data=nv))
summary(lm(CRP~BMI,data=nv))

fit=glm(Depression~Insomnia, data=nv, family = binomial())
summary(fit)
exp(cbind(OR = coef(fit), confint(fit)))
fit=glm(Arthritis~Depression, data=nv, family = binomial())
summary(fit)
exp(cbind(OR = coef(fit), confint(fit)))

#make table for regression coef
setwd('/Users/Selene/Desktop')
sig=averaged.network(boot)$learning$args$threshold
links=as.data.frame(avg.boot$arcs)
links$strength=boot[(boot$strength>sig)&(boot$direction>0.5),]$strength
links$direction=boot[(boot$strength>sig)&(boot$direction>0.5),]$direction
links$direction=ifelse(links$direction<=0.57,0.6,links$direction)
links=links[order(links$to),]
write.csv(links,file='links.csv')


#conditional indep analysis through markov blanket
mb(avg.boot,'PA')
mb(avg.boot,'BMI')
mb(avg.boot,'QOLp')
mb(avg.boot,'QOLm')
mb(avg.boot,'Insulin')
mb(avg.boot,'CRP')


#analysis on sub-network after discretization
padf=data.frame(nv$PA, nv$Age, nv$BMI, nv$Sleep2)
names(padf)=c("PA","Age","BMI","Sleep")
padf$Age=ifelse(padf$Age<65,0,1)
padf$BMI=ifelse(padf$BMI<30,0,1)
padf$Sleep=ifelse(padf$Sleep<median(padf$Sleep),0,1)
d=aggregate(padf$PA,list(padf$Age,padf$BMI,padf$Sleep),mean)


#bn inference using cpquery (this tests how one part of network reacts to 
#disturbances from another part of the network)
fitted=bn.fit(avg.boot,nv)

#study how mb(Insulin) (excpet age)  -->  insulin(raw)>870 (log(870)=6.768493)
mb(avg.boot,'Insulin') #BMI, Sleep2, PA
result=cpdist(fitted,nodes=c('Insulin'),evidence=( PA<155  & BMI>=40 ))
mean(result$Insulin) #need BMI to exceed 95 %-tile, PA to be less than 10 %-tile

cpquery(fitted,event=(BMI>30),evidence=(PA>335))

result=cpdist(fitted,nodes=c('PA'),evidence=(Arthritis=='Yes'))
mean(result$PA)
result=cpdist(fitted,nodes=c('PA'),evidence=(Arthritis=='No'))
mean(result$PA)

result=cpdist(fitted,nodes=c('BMI'),evidence=(PA<270))
mean(result$BMI)
result=cpdist(fitted,nodes=c('BMI'),evidence=(PA >=380))
mean(result$BMI)

result=cpdist(fitted,nodes=c('Insulin'),evidence=(PA<260 & PA >240 & BMI<32 & BMI>30))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(PA<310 & PA >290 & BMI<30 & BMI>29))
mean(result$Insulin)

result=cpdist(fitted,nodes=c('CRP'),evidence=(PA<260 & PA >240 & BMI<32 & BMI>30))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(PA<310 & PA >290 & BMI<30 & BMI>29))
mean(result$CRP)

result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2>52.9 & PA<260 & PA >240))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2<41.4 & PA<310 & PA >290))
mean(result$CRP)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2>52.9 & PA<260 & PA >240))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2<41.4 & PA<310 & PA >290))
mean(result$Insulin)

result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2>52.9 & PA<260 & PA >240 & BMI<32 & BMI>30))
mean(result$CRP)
result=cpdist(fitted,nodes=c('CRP'),evidence=(Sleep2<41.4 & PA<310 & PA >290 & BMI<30 & BMI>29))
mean(result$CRP)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2>52.9 & PA<260 & PA >240 & BMI<32 & BMI>30))
mean(result$Insulin)
result=cpdist(fitted,nodes=c('Insulin'),evidence=(Sleep2<41.4 & PA<310 & PA >290 & BMI<30 & BMI>29))
mean(result$Insulin)

result=cpdist(fitted,nodes=c('QOLp'),evidence=(BMI>=30))
mean(result$QOLp)
result=cpdist(fitted,nodes=c('QOLp'),evidence=(BMI<30))
mean(result$QOLp)

result=cpdist(fitted,nodes=c('QOLm'),evidence=(BMI>=30))
mean(result$QOLm)
result=cpdist(fitted,nodes=c('QOLm'),evidence=(BMI<30))
mean(result$QOLm)

result=cpdist(fitted,nodes=c('QOLp'),evidence=(Sleep2>52.9))
mean(result$QOLp)
result=cpdist(fitted,nodes=c('QOLp'),evidence=(Sleep2<41.4))
mean(result$QOLp)

result=cpdist(fitted,nodes=c('QOLm'),evidence=(Sleep2>52.9))
mean(result$QOLm)
result=cpdist(fitted,nodes=c('QOLm'),evidence=(Sleep2<41.4))
mean(result$QOLm)

result=cpdist(fitted,nodes=c('QOLp'),evidence=(Sleep2>52.9 & BMI>=30))
mean(result$QOLp)
result=cpdist(fitted,nodes=c('QOLp'),evidence=(Sleep2<41.4 & BMI<30))
mean(result$QOLp)

result=cpdist(fitted,nodes=c('QOLm'),evidence=(Sleep2>52.9 & BMI>=30))
mean(result$QOLm)
result=cpdist(fitted,nodes=c('QOLm'),evidence=(Sleep2<41.4 & BMI<30))
mean(result$QOLm)

result=cpdist(fitted,nodes=c('BMI'),evidence=(PA<270 &Sleep2>52.9))
mean(result$BMI)
result=cpdist(fitted,nodes=c('BMI'),evidence=(PA >=380 & Sleep2<41.4))
mean(result$BMI)









