library(nimble)

ipop1Code <- nimbleCode({
  nad ~ dpois(y1) 
  Ny[1] <- 0
  Nj[1] <- 0
  Nad[1] <- nad
  Ntot[1] <- Nad[1]+Ny[1]
  
  sj[1] <- 0
  logit(sa[1]) <- Sp[2]+Sepsa[1]
  log(f[1]) <- Rp+Reps[1]
  Seps[1]<- 0
  Reps[1]~dnorm(0,tauR)
  Sepsa[1]~dnorm(0,tauSa)
  # Survival and recapture probabilities, as well as productivity
  for (t in 2:(nyears-1)){
    logit(sj[t]) <- Sp[1]+Seps[t]
    logit(sa[t]) <- Sp[1]+Sp[2]+Sepsa[t]
    log(f[t]) <- Rp+Reps[t]
    
    Seps[t]~dnorm(0,tauS)
    Sepsa[t]~dnorm(0,tauSa)
    Reps[t]~dnorm(0,tauR)
  }
  for (u in 1:2){
    Sp[u]~ dnorm(0, 0.01)
  }       
  Rp~ dnorm(0, 0.01)
  tauS <- pow(sigmaS, -2)
  sigmaS ~ dunif(0, 10)
  tauSa <- pow(sigmaSa, -2)
  sigmaSa ~ dunif(0, 10)    
  tauR <- pow(sigmaR, -2)
  sigmaR ~ dunif(0, 10)
  
  #-------------------------------------------------
  # 2. Derived parameters
  #-------------------------------------------------
  # Population growth rate
  for (t in 1:(nyears-1)){
   l[t]<-Ntot[t+1] / (Ntot[t]+0.00001)
  }
  
  #-------------------------------------------------
  # 3. The likelihoods of the single data sets
  #-------------------------------------------------
  # 3.1. Likelihood for population population count data (state-space model)
  # 3.1.1 System process
  for (t in 2:nyears){
    Nad[t] ~ dbin(sa[t-1], Ntot[t-1])
    Ny[t] ~ dbin(sj[t-1], Nj[t-1])
    Ntot[t] <- Nad[t] + Ny[t]
    Nj[t] ~ dpois(f[t-1]*Ntot[t])
  }

  # 3.1.2 Observation process
  for (t in 1:nyears){
    y[t] ~ dpois(Ntot[t])
  }
  
  # 3.2 Likelihood for survival (2 age classes)
  for (t in 1:(nyears-1)){
    NSurv[1,t] ~ dbin(sj[t],NInd[1,t])
    NSurv[2,t] ~ dbin(sa[t],NInd[2,t])
  }
  
  # 3.3. Likelihood for productivity data: Poisson regression
  for (t in 1:(nyears-1)){
    Rec[t] ~ dpois(rho[t])
    rho[t] <- Nest[t]*f[t]
  }
}
)

parameters <- c('Rp','Sp','sigmaS','sigmaR','sigmaSa',"sj", "sa", "f", "Ny", "Nad", "l")

Consts <- list(nyears =nbr.an)
data <- list(y = as.numeric(table(Count$year)),y1=as.numeric(table(Count$year))[1], Nest = as.numeric(table(Spop$year[Spop$S==1]))[1:(nbr.an-1)] , Rec = aggregate(x =Spop$R[Spop$S==1], by = list(Spop$year[Spop$S==1]), FUN = "sum")$x[1:(nbr.an-1)], NSurv=table(Spop$agef[Spop$S==1],Spop$year[Spop$S==1])[,1:(nbr.an-1)], NInd=table(Spop$agef,Spop$year)[,1:(nbr.an-1)])
Inits <- list(Rp = 0, Sp = c(0.5,0.5), sigmaS=0.1, sigmaSa=0.1, sigmaR=0.1)
projp <- nimbleModel(code = ipop1Code , name = 'projp' , constants = Consts,
                    data = data, inits = Inits)
Cproj <- compileNimble(projp)
spec <- configureMCMC(projp, monitors=c('Rp','Sp','sigmaS','sigmaR','sigmaSa',"sj", "sa", "f", "Ny", "Nad", "l"),thin=5)
projpMCMC <- buildMCMC(spec)
CprojpMCMC <- compileNimble(projpMCMC,project=projp)
CprojpMCMC$run(2000)
MCMCsamples <- as.matrix(CprojpMCMC$mvSamples)
apply(MCMCsamples[101:400,],2, mean)





ipm2_1Code <- nimbleCode({
    
    # Initial population density
    for (a in 1:2){
    for (trait in 1:nn){
    vf[a,trait,1] <- inipop[a,trait]
    }
    }
    #-------------------------------------------------
    # 1. Define the priors for the parameters
    #-------------------------------------------------
    
    # Survival, recruitment and inheritance parameters
    for (u in 1:3){
    Sp[u] ~ dnorm(0,0.01)
    }  
    for (u in 1:2){
    Rp[u] ~ dnorm(0,0.01)
    }
    for (u in 1:2){
    Ip[u] ~ dnorm(0,0.01)
    }
    Ipv~ dunif(0, 10)
    
    
    
    
    #-------------------------------------------------
    # 2. Derived parameters
    #-------------------------------------------------
    for(trait in 1:nn){
      for(traitY in 1:nn){
        I[traitY,trait] <-exp(-((X[traitY]-(Ip[1]+Ip[2]*X[trait]))^2)/(2*Ipv*Ipv))/(sqrt(2*pi)*Ipv)
      }
      SI[trait]<-sum(I[1:nn,trait])+0.000001
    }
    
    for(trait in 1:nn){
      for(traitY in 1:nn){
        Inorm[traitY,trait] <- I[traitY,trait]/SI[trait]
      }
    }
    
    for (it in 1:(nbr.an-1)){
    for(trait in 1:nn){
    logit(Sj[trait,it]) <- Sp[1]+Sp[3]*X[trait]
    logit(Sa[trait,it]) <- Sp[2]+Sp[3]*X[trait]
    log(R [trait,it]) <-   Rp[1]+Rp[2]*X[trait]
    
    vf[2,trait,it+1] <- vf[1,trait,it]*Sj[trait,it]+vf[2,trait,it]*Sa[trait,it]
    for(traitY in 1:nn){
    v1junk[traitY,trait,(it+1)] <- Inorm[traitY,trait]*(R[trait,it]*vf[2,trait,(it+1)])
    }
    }
    for(trait in 1:nn){
      vf[1,trait,(it+1)] <- sum(v1junk[trait,1:nn,(it+1)])
    }
    
    # Population growth rate
    l[it] <- sum(vf[1:2,1:nn,(it+1)])/sum(vf[1:2,1:nn,(it)])
    }
    
    #-------------------------------------------------
    # 3. The likelihoods of the single data sets
    #-------------------------------------------------
    # 3.1. Likelihood for population count data
    for (t in 1:nbr.an){
    for (trait in 1:nn){
    y[trait,t] ~ dpois(vf[2,trait,t])
    }
    }
    
    #-------------------------------------------------
    # 3.2. Survival
    for (i in 1:Ndatasur){
    Sdata[i] ~ dbern(psur[i])
    logit(psur[i]) <- Sp[Sage[i]+1]+Sp[3]*SX[i]
    }
    
    # 3.3. Recruitment
    for (i in 1:Ndatarec){
    Rdata[i] ~ dpois(Rlambda[i])
    log(Rlambda[i]) <- Rp[1]+Rp[2]*RX[i]
    }
    
    # 3.4. Inheritance
    for (i in 1:Ndatainh){
    Idata[i] ~ dnorm(mu[i],tau_I)
    mu[i] <- Ip[1]+Ip[2]*IX[i]
    }
    
    }
)

data2 <- list( y = table(cut(Count$X,breaks=b),Count$year), inipop=iniPop2 ,  Sdata=Spop$S, SX=Spop$X,  Rdata=Spop$R[Spop$S==1], RX=Spop$X[Spop$S==1], Idata=Spop$X[Spop$inh==1], IX=Spop$PdsMere[Spop$inh==1])
Consts2 <- list(X=meshpoints,pi=pi, nbr.an=nbr.an, nn=nn, Ndatasur=length(Spop$id), Ndatarec=length(Spop$id[Spop$S==1]), Ndatainh=length(Spop$X[Spop$inh==1]),Sage=Spop$agef)
Inits2 <- list(Rp = c(0,0),Ip = c(0,0),Ipv=0, Sp = c(0,0.5,0.5))
mod2 <- nimbleModel(code = ipm2_1Code , name = 'mod2' , constants = Consts2,
                    data = data2, inits = Inits2)
Cmod2 <- compileNimble(mod2)
# spec <- configureMCMC(mod2, monitors=c('Rp','Sp','Ip','Ipv','l','vf'),thin=5)
spec <- configureMCMC(mod2, monitors=c('Rp','Sp','Ip','Ipv'),thin=5)
mod2MCMC <- buildMCMC(spec)
Cmod2MCMC <- compileNimble(mod2MCMC,project=mod2)
Cmod2MCMC$run(5000)
MCMCsamplesmod2 <- as.matrix(Cmod2MCMC$mvSamples)


apply(MCMCsamples[401:1000,],2, mean)




