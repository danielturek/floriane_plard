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
    Ipv~ dunif(0, 100)
    tau_I <- 1 / (Ipv^2)
    
    
    
    
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
    logit(Sa[trait,it]) <- Sp[1]+Sp[2]+Sp[3]*X[trait]
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
    logit(psur[i]) <- Sp[1]+Sp[2]*Sage[i]+Sp[3]*SX[i]
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
Inits2 <- list(Rp = c(0,0),Ip = c(0,0),Ipv=1, Sp = c(0,0.5,0.5))
mod2 <- nimbleModel(code = ipm2_1Code , name = 'mod2' , constants = Consts2,
                    data = data2, inits = Inits2)
Cmod2 <- compileNimble(mod2)
# spec <- configureMCMC(mod2, monitors=c('Rp','Sp','Ip','Ipv','l','vf'),thin=5)
spec <- configureMCMC(mod2, monitors=c('Rp','Sp','Ip','Ipv'),thin=5, useConjugacy=FALSE)
spec$removeSamplers('Rp', print = FALSE)
spec$addSampler(target = 'Rp[1:2]', type = 'RW_block')
spec$removeSamplers('Sp', print = FALSE)
spec$addSampler(target = 'Sp[1:3]', type = 'RW_block')
spec$removeSamplers('Ip', print = FALSE)
spec$addSampler(target = 'Ip[1:2]', type = 'RW_block')


mod2MCMC <- buildMCMC(spec)
Cmod2MCMC <- compileNimble(mod2MCMC,project=mod2,resetFunctions = TRUE)
Cmod2MCMC$run(5000)
MCMCsamplesmod2 <- as.matrix(Cmod2MCMC$mvSamples)[401:1000,]

par(mfrow=c(4,2))
plot(MCMCsamplesmod2[ , 1], type = 'l', xlab = 'iteration', ylab = 'Ip1')
plot(MCMCsamplesmod2[ , 2], type = 'l', xlab = 'iteration', ylab = 'Ip2')
plot(MCMCsamplesmod2[ , 3], type = 'l', xlab = 'iteration', ylab = 'Ipv')
plot(MCMCsamplesmod2[ , 4], type = 'l', xlab = 'iteration', ylab = 'Rp1')
plot(MCMCsamplesmod2[ , 5], type = 'l', xlab = 'iteration', ylab = 'Rp2')
plot(MCMCsamplesmod2[ , 6], type = 'l', xlab = 'iteration', ylab = 'Sp1')
plot(MCMCsamplesmod2[ , 7], type = 'l', xlab = 'iteration', ylab = 'Sp2')
plot(MCMCsamplesmod2[ , 8], type = 'l', xlab = 'iteration', ylab = 'Sp3')

apply(MCMCsamplesmod2,2, mean)


################## Third MODEL

ipm2_3Code <- nimbleCode({
  
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
  for (u in 1:5){
    Sp[u] ~ dnorm(0,0.01)
  }  
  for (u in 1:4){
    Rp[u] ~ dnorm(0,0.01)
  }
  for (u in 1:2){
    Ip[u] ~ dnorm(0,0.01)
  }
  Ipv~ dunif(0, 100)
  tau_I <- 1 / (Ipv^2)
  sdY~ dunif(0, 100)
  tauY <- 1 / (sdY^2)
  
  
  
  
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
      logit(Sj[trait,it]) <- Sp[1]+Sp[3]*X[trait]+Sp[4]*COV[it]+Sp[5]*X[trait]*COV[it]
      logit(Sa[trait,it]) <- Sp[1]+Sp[2]+Sp[3]*X[trait]+Sp[4]*COV[it]+Sp[5]*X[trait]*COV[it]
      log(R [trait,it]) <-   Rp[1]+Rp[2]*X[trait]+Rp[3]*COV[it]+Rp[4]*X[trait]*COV[it]
      
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
      y[trait,t] ~ dnorm(vf[2,trait,t],tauY)
    }
  }
  
  #-------------------------------------------------
  # 3.2. Survival
  for (i in 1:Ndatasur){
    Sdata[i] ~ dbern(psur[i])
    logit(psur[i]) <- Sp[1]+Sp[2]*Sage[i]+Sp[3]*SX[i]+Sp[4]*SCOV[i]+Sp[5]*SX[i]*SCOV[i]
  }
  
  # 3.3. Recruitment
  for (i in 1:Ndatarec){
    Rdata[i] ~ dpois(Rlambda[i])
    log(Rlambda[i]) <- Rp[1]+Rp[2]*RX[i]+Rp[3]*RCOV[i]+Rp[4]*RX[i]*RCOV[i]
  }
  
  # 3.4. Inheritance
  for (i in 1:Ndatainh){
    Idata[i] ~ dnorm(mu[i],tau_I)
    mu[i] <- Ip[1]+Ip[2]*IX[i]
  }
  
}
)


data3 <- list( y = table(cut(Count$X,breaks=b),Count$year), COV=Covs[1:(nbr.an-1)], inipop=iniPop2 ,  Sdata=Spop$S, SX=Spop$Xs[,1], SCOV=Spop$Covs,Sage=Spop$agef,  Rdata=Spop$R[Spop$S==1], RX=Spop$Xs[Spop$S==1], RCOV=Spop$Covs[Spop$S==1], Idata=Spop$Xs[Spop$inh==1], IX=Spop$PdsMeres[Spop$inh==1])
Consts3 <- list(X=meshpoints,pi=pi, nbr.an=nbr.an, nn=nn, Ndatasur=length(Spop$id), Ndatarec=length(Spop$id[Spop$S==1]), Ndatainh=length(Spop$X[Spop$inh==1]))
# Inits3 <- list(Rp = c(0,0,0,0),Ip = c(0,0),Ipv=1,sdY=1, Sp = c(0.5,0.5,0,0,0))
inits3 <-list(Rp= as.numeric(ar$coefficients[1:4]), Sp= as.numeric(as$coefficients[1:5]), Ip=as.numeric(ai$coefficients[1:2]), Ipv= as.numeric(air$coefficients[1]), sdY=1)
mod3 <- nimbleModel(code = ipm2_3Code , name = 'mod3' , constants = Consts3,
                    data = data3, inits = inits3)
Cmod3 <- compileNimble(mod3)
#  spec <- configureMCMC(mod3, monitors=c('Rp','Sp','Ip','Ipv','sdY','l','vf'),thin=5, useConjugacy=FALSE)
 spec <- configureMCMC(mod3, monitors=c('Rp','Sp','Ip','Ipv','sdY'),thin=5, useConjugacy=FALSE)
spec$removeSamplers('Rp', print = FALSE)
spec$addSampler(target = 'Rp[1:4]', type = 'RW_block')
spec$removeSamplers('Sp', print = FALSE)
spec$addSampler(target = 'Sp[1:5]', type = 'RW_block')
spec$removeSamplers('Ip', print = FALSE)
spec$addSampler(target = 'Ip[1:2]', type = 'RW_block')


mod3MCMC <- buildMCMC(spec)
Cmod3MCMC <- compileNimble(mod3MCMC,project=mod3,resetFunctions = TRUE)
Cmod3MCMC$run(50000)
MCMCsamplesmod3 <- as.matrix(Cmod3MCMC$mvSamples)[1:10000,]

par(mfrow=c(3,2))
plot(MCMCsamplesmod3[ , 1], type = 'l', xlab = 'iteration', ylab = 'Ip1')
plot(MCMCsamplesmod3[ , 2], type = 'l', xlab = 'iteration', ylab = 'Ip2')
plot(MCMCsamplesmod3[ , 3], type = 'l', xlab = 'iteration', ylab = 'Ipv')
plot(MCMCsamplesmod3[ , 4], type = 'l', xlab = 'iteration', ylab = 'Rp1')
plot(MCMCsamplesmod3[ , 5], type = 'l', xlab = 'iteration', ylab = 'Rp2')
plot(MCMCsamplesmod3[ , 6], type = 'l', xlab = 'iteration', ylab = 'Rp3')
plot(MCMCsamplesmod3[ , 7], type = 'l', xlab = 'iteration', ylab = 'Rp4')
plot(MCMCsamplesmod3[ , 8], type = 'l', xlab = 'iteration', ylab = 'Sp1')
plot(MCMCsamplesmod3[ , 9], type = 'l', xlab = 'iteration', ylab = 'Sp2')
plot(MCMCsamplesmod3[ , 10], type = 'l', xlab = 'iteration', ylab = 'Sp3')
plot(MCMCsamplesmod3[ , 11], type = 'l', xlab = 'iteration', ylab = 'Sp4')
plot(MCMCsamplesmod3[ , 12], type = 'l', xlab = 'iteration', ylab = 'Sp5')

apply(MCMCsamplesmod3,2, mean)


