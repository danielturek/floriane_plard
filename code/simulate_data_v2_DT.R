inisize=100 ##number of individuals in year 1
nbr.an=15 ###number of years
nn <- 40 # nn of classes within the model



modsr=  3
pc=   0.85


###fonctions
library(boot)

#parameters
Sur=matrix(c(-1,1,0.5,0,0,
             0.5,1,0.5,-0.18,0,
             1,1,0.5,-0.55,0.25),nrow=3,byrow=T)
Rec=matrix(c(0.5,0,0.17,0,0,
             0.6,0,0.17,-0.08,0,
             0.7,0,0.15, -0.5,0.17),nrow=3,byrow=T)
Inh=c(0,0.3,1,0,0)



minsize <- -5 # minimum possible size 
maxsize <- 6 # maximum possible size 
me<-0
sde<-1
n.age <- 3 # number of ages -- senescents pooled into final class
L<-1.1*minsize; U<-1.1*maxsize
# boundary points b and mesh points y
b<-L+c(0:nn)*(U-L)/nn
meshpoints<-0.5*(b[1:nn]+b[2:(nn+1)])

COV=c(1:nbr.an)
Covs=scale(COV)



##########################

#     WHICH MODEL?
##########################

PARAM=matrix(c(Sur[modsr,], Rec[modsr,],Inh),nrow=3, byrow=T)
prs=pinh=pc    


##Simulation
 sim=1
year=id=PdsMere=X=count=Mark=S=R=inh=agef=rep(0,nbr.an*inisize*1000)
Tpop=data.frame(id,year,PdsMere,X,count,Mark,S,R,inh,agef)
siz=rep(0,nbr.an+1)
  
  
  ##First year
  siz[1]=inisize
  Tpop$year[1:inisize]=rep(1,inisize)
  Tpop$id[1:inisize]=c(1:inisize)
  Tpop$PdsMere[1:inisize]=rep(NA,inisize)
  Tpop$X[1:inisize]=rnorm(inisize,me,sde)
  Tpop$agef[1:inisize]=rep(1,inisize)
  
  Tpop$S[1:inisize]=rbinom(length(Tpop$X[Tpop$year==1]),1,inv.logit(PARAM[1,1]+PARAM[1,3]*Tpop$X[Tpop$year==1]+PARAM[1,2]*Tpop$agef[Tpop$year==1]+PARAM[1,4]*COV[1]+PARAM[1,5]*COV[1]*Tpop$X[Tpop$year==1]))
  Tpop$R[1:inisize]=rpois(length(Tpop$X[Tpop$year==1]),exp(PARAM[2,1]+PARAM[2,3]*Tpop$X[Tpop$year==1]+PARAM[2,4]*COV[1]+PARAM[2,5]*COV[1]*Tpop$X[Tpop$year==1]))
  Tpop$R[Tpop$S==0]=0
  siz[2]=sum(Tpop$S[Tpop$year==1]+Tpop$R[Tpop$year==1])+siz[1]
  
  for (t in c(1:(nbr.an-1))){
    Tpop$year[(siz[t]+1):(siz[t+1])]=t+1
    Tpop$agef[(siz[t]+1):(siz[t]+sum(Tpop$S[Tpop$year==t]))]=1
    Tpop$agef[(siz[t]+1+sum(Tpop$S[Tpop$year==t])):siz[t+1]]=0
    
    ##Survivors
    Tpop[Tpop$year==(t+1)&Tpop$agef==1,c(1,3,4)]=Tpop[Tpop$year==t&Tpop$S==1,c(1,3,4)]
    
    ###Recruits
    Tpop$id[Tpop$year==(t+1)&Tpop$agef==0]=c((siz[t]+1+sum(Tpop$S[Tpop$year==t])):siz[t+1])
    Tpop$PdsMere[Tpop$year==(t+1)&Tpop$agef==0]=rep(Tpop$X[Tpop$year==t&Tpop$R>0],times=Tpop$R[Tpop$year==t&Tpop$R>0])
    
    ##Inheritance
    Tpop$X[Tpop$year==(t+1)&Tpop$agef==0]=rnorm(length(Tpop$X[Tpop$year==(t+1)&Tpop$agef==0]),PARAM[3,1]+PARAM[3,2]*Tpop$PdsMere[Tpop$year==(t+1)&Tpop$agef==0],sd=PARAM[3,3])
    
    ### next year survival and recruitment
    Tpop$S[Tpop$year==(t+1)]=rbinom(length(Tpop$X[Tpop$year==(t+1)]),1,inv.logit(PARAM[1,1]+PARAM[1,3]*Tpop$X[Tpop$year==(t+1)]+PARAM[1,2]*Tpop$agef[Tpop$year==(t+1)]+PARAM[1,4]*COV[(t+1)]+PARAM[1,5]*COV[(t+1)]*Tpop$X[Tpop$year==(t+1)]))
    Tpop$R[Tpop$year==(t+1)]=rpois(length(Tpop$X[Tpop$year==(t+1)]),exp(PARAM[2,1]+PARAM[2,3]*Tpop$X[Tpop$year==(t+1)]+PARAM[2,4]*COV[(t+1)]+PARAM[2,5]*COV[(t+1)]*Tpop$X[Tpop$year==(t+1)]))
    Tpop$R[Tpop$year==(t+1)&Tpop$S==0]=0
    siz[t+2]=sum(Tpop$S[Tpop$year==(t+1)]+Tpop$R[Tpop$year==(t+1)])+siz[t+1]
  }
  
  
  ###################################################
  
  # DIFFERENT DATA COLLECTION
  
  #################################################
  
   
  Tpop=subset(Tpop,Tpop$year>0)
  
  
  Tpop$Mark[1:inisize]=rbinom(inisize,1,prs) 
  for (t in c(1:(nbr.an-1))){
    Tpop$Mark[Tpop$year==(t+1)&Tpop$agef==1]=Tpop$Mark[Tpop$year==t&Tpop$S==1]
    
    Markj=Tpop$Mark[Tpop$year==t&Tpop$R>0]
    Markj[Tpop$Mark[Tpop$year==t&Tpop$R>0]==1]=rbinom(length(Markj[Tpop$Mark[Tpop$year==t&Tpop$R>0]==1]),1,pinh)##chance d'etre marqué si mere marquée
    Markj[Tpop$Mark[Tpop$year==t&Tpop$R>0]==0]=rbinom(length(Markj[Tpop$Mark[Tpop$year==t&Tpop$R>0]==0]),1,prs*(1-pinh))##sinon
    
    Tpop$Mark[Tpop$year==(t+1)&Tpop$agef==0]=rep(Markj,times=Tpop$R[Tpop$year==t&Tpop$R>0])
    Tpop$inh[Tpop$year==(t+1)&Tpop$agef==0]=rep(Tpop$Mark[Tpop$year==t&Tpop$R>0],times=Tpop$R[Tpop$year==t&Tpop$R>0])
  }
  Tpop$count[Tpop$agef==0]=0
  Tpop$count[Tpop$agef==1]=rbinom(length(Tpop$count[Tpop$agef==1]),1,pc)
  
  Tpop$Xs=scale(Tpop$X)
  Tpop$PdsMeres=(Tpop$PdsMere-mean(Tpop$X))/sd(Tpop$X)
  Tpop$Covs=Covs[Tpop$year]
  
  Spop=subset(Tpop, Tpop$Mark==1)
  Count=subset(Tpop, Tpop$count==1,select=c('year','Xs'))

###############initial pop
  iniPop2=matrix(rep(0,nn*(n.age-1)),nrow=n.age-1,byrow=T)
  iniPop2[2,]=table(cut(Count$Xs[Count$year==1],breaks=b))
  ypop=as.numeric(table(Count$year))
  

if (modsr==3){as=glm(S~agef+Xs*Covs, family='binomial', data=Tpop)
              ar=glm(R~Xs*Covs, family='poisson', data=Tpop[Tpop$S==1,])}
ai=lm(Xs~PdsMeres, data=Tpop)
air=lm(ai$residuals^2~1, data=Tpop)
