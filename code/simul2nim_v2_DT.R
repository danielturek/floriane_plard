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

################## Floraine, start here:


library(nimble)
load(file = '../data/data_v2.RData')  ## I've saved the simulated data


## TO CHANGE DATA: you need to define all these functions,
## only once.  This defines 3 custom distributions, and registers
## them for use with NIMBLE

dSdataBern <- nimbleFunction(
    run = function(x = double(1), prob = double(1), length = integer(), log.p = double()) {
        ll <- 0
        for(i in 1:length) {
            ll <- ll + dbinom(x[i], prob = prob[i], size = 1, log=TRUE)
        }
        returnType(double())
        if(log.p) return(ll) else return(exp(ll))
    }
)
rSdataBern <- nimbleFunction(
    run = function(n = integer(), prob = double(1), length = integer()) {
        print('this should never run')
        ##x <- numeric(length)
        declare(x, double(1, length))
        returnType(double(1))
        return(x)
    }
)
dRdataPois <- nimbleFunction(
    run = function(x = double(1), lambda = double(1), length = integer(), log.p = double()) {
        ll <- 0
        for(i in 1:length) {
            ll <- ll + dpois(x[i], lambda[i], log=TRUE)
        }
        returnType(double())
        if(log.p) return(ll) else return(exp(ll))
    }
)
rRdataPois <- nimbleFunction(
    run = function(n = integer(), lambda = double(1), length = integer()) {
        print('this should never run')
        ##x <- numeric(length)
        declare(x, double(1, length))
        returnType(double(1))
        return(x)
    }
)
dIdataNorm <- nimbleFunction(
    run = function(x = double(1), mean = double(1), tau = double(), length = integer(), log.p = double()) {
        sd <- 1/sqrt(tau)
        ll <- 0
        for(i in 1:length) {
            ll <- ll + dnorm(x[i], mean[i], sd=sd, log=TRUE)
        }
        returnType(double())
        if(log.p) return(ll) else return(exp(ll))
    }
)
rIdataNorm <- nimbleFunction(
    run = function(n = integer(), mean = double(1), tau = double(), length = integer()) {
        print('this should never run')
        ##x <- numeric(length)
        declare(x, double(1, length))
        returnType(double(1))
        return(x)
    }
)
registerDistributions(list(
    dSdataBern = list(
        BUGSdist = 'dSdataBern(prob, length)',
        types    = c('value = double(1)', 'prob = double(1)', 'length = integer()')
    ),
    dRdataPois = list(
        BUGSdist = 'dRdataPois(lambda, length)',
        types    = c('value = double(1)', 'lambda = double(1)', 'length = integer()')
    ),
    dIdataNorm = list(
        BUGSdist = 'dIdataNorm(mean, tau, length)',
        types    = c('value = double(1)', 'mean = double(1)', 'tau = double()', 'length = integer()')
    )
))


## there are a number of changes I made to the model, so speed it up.
## these are denoted by CHANGE1, CHANGE2, and CHANGE3.
## some timings are shown later, so show the speedup.

ipm2_3Code <- nimbleCode({
    ## Initial population density
    ## CHANGE1: vectorize this initial vf declaration:
    ##for (a in 1:2){
    ##    for (trait in 1:nn){
    ##        vf[a,trait,1] <- inipop[a,trait]
    ##    }
    ##}
    vf[1:2, 1:nn, 1] <- inipop[1:2, 1:nn]
    ##-------------------------------------------------
    ## 1. Define the priors for the parameters
    ##-------------------------------------------------
    ## Survival, recruitment and inheritance parameters
    for (u in 1:5)
        Sp[u] ~ dnorm(0,0.01)
    for (u in 1:4)
        Rp[u] ~ dnorm(0,0.01)
    for (u in 1:2)
        Ip[u] ~ dnorm(0,0.01)
    Ipv ~ dunif(0, 100)
    tau_I <- 1 / (Ipv^2)
    sdY~ dunif(0, 100)
    tauY <- 1 / (sdY^2)
    ##-------------------------------------------------
    ## 2. Derived parameters
    ##-------------------------------------------------
    for(trait in 1:nn){
        ## CHANGE2: vectorize I[,] and Inorm[,] over traitY
        ##for(traitY in 1:nn){
        ##    I[traitY,trait] <- exp(-((X[traitY]-(Ip[1]+Ip[2]*X[trait]))^2)/(2*Ipv*Ipv))/(sqrt(2*pi)*Ipv)
        ##    Inorm[traitY,trait] <- I[traitY,trait]/SI[trait]
        ##}
        I[1:nn, trait] <- exp(-((X[1:nn] - (Ip[1]+Ip[2]*X[trait]))^2)/(2*Ipv*Ipv))/(sqrt(2*pi)*Ipv)
        SI[trait] <- sum(I[1:nn, trait]) + 0.000001
        Inorm[1:nn, trait] <- I[1:nn, trait] / SI[trait]
    }
    for (it in 1:(nbr.an-1)){
        for(trait in 1:nn){
            logit(Sj[trait,it]) <- Sp[1]+Sp[3]*X[trait]+Sp[4]*COV[it]+Sp[5]*X[trait]*COV[it]
            logit(Sa[trait,it]) <- Sp[1]+Sp[2]+Sp[3]*X[trait]+Sp[4]*COV[it]+Sp[5]*X[trait]*COV[it]
            log(R [trait,it]) <-   Rp[1]+Rp[2]*X[trait]+Rp[3]*COV[it]+Rp[4]*X[trait]*COV[it]
            vf[2,trait,it+1] <- vf[1,trait,it]*Sj[trait,it]+vf[2,trait,it]*Sa[trait,it]
            ## CHANGE2: vectorize I[,] and Inorm[,] over traitY
            ##for(traitY in 1:nn){
            ##    v1junk[traitY,trait,(it+1)] <- Inorm[traitY,trait]*(R[trait,it]*vf[2,trait,(it+1)])
            ##}
            v1junk[1:nn, trait, (it+1)] <- Inorm[1:nn, trait] * (R[trait,it] * vf[2,trait,(it+1)])
        }
        for(trait in 1:nn){
            vf[1,trait,(it+1)] <- sum(v1junk[trait,1:nn,(it+1)])
        }
        ## Population growth rate
        l[it] <- sum(vf[1:2,1:nn,(it+1)])/sum(vf[1:2,1:nn,(it)])
    }
    ##-------------------------------------------------
    ## 3. The likelihoods of the single data sets
    ##-------------------------------------------------
    ## 3.1. Likelihood for population count data
    for (t in 1:nbr.an){
        for (trait in 1:nn){
            y[trait,t] ~ dnorm(vf[2,trait,t],tauY)
        }
    }
    ##-------------------------------------------------
    ##CHANGE3: using custom distributions here, which allows for variable data.
    #### 3.2. Survival
    ##for (i in 1:Ndatasur){
    ##    Sdata[i] ~ dbern(psur[i])
    ##    logit(psur[i]) <- Sp[1]+Sp[2]*Sage[i]+Sp[3]*SX[i]+Sp[4]*SCOV[i]+Sp[5]*SX[i]*SCOV[i]
    ##}
    #### 3.3. Recruitment
    ##for (i in 1:Ndatarec){
    ##    Rdata[i] ~ dpois(Rlambda[i])
    ##    log(Rlambda[i]) <- Rp[1]+Rp[2]*RX[i]+Rp[3]*RCOV[i]+Rp[4]*RX[i]*RCOV[i]
    ##}
    #### 3.4. Inheritance
    ##for (i in 1:Ndatainh){
    ##    Idata[i] ~ dnorm(mu[i],tau_I)
    ##    mu[i] <- Ip[1]+Ip[2]*IX[i]
    ##}
    ## 3.2. Survival
    logit(psur[1:Ndatasur]) <- Sp[1]+Sp[2]*Sage[1:Ndatasur]+Sp[3]*SX[1:Ndatasur]+Sp[4]*SCOV[1:Ndatasur]+Sp[5]*SX[1:Ndatasur]*SCOV[1:Ndatasur]
    Sdata[1:Ndatasur] ~ dSdataBern(psur[1:Ndatasur], SdataLength)
    ## 3.3. Recruitment
    log(Rlambda[1:Ndatarec]) <- Rp[1]+Rp[2]*RX[1:Ndatarec]+Rp[3]*RCOV[1:Ndatarec]+Rp[4]*RX[1:Ndatarec]*RCOV[1:Ndatarec]
    Rdata[1:Ndatarec] ~ dRdataPois(Rlambda[1:Ndatarec], RdataLength)
    ## 3.4. Inheritance
    mu[1:Ndatainh] <- Ip[1]+Ip[2]*IX[1:Ndatainh]
    Idata[1:Ndatainh] ~ dIdataNorm(mu[1:Ndatainh], tau_I, IdataLength)
})

##CHANGE3: using custom distributions here, which allows for variable data.
## these variables define the length of our actual (first) data set
Ndatasur <- length(Spop$id)
Ndatarec <- length(Spop$id[Spop$S==1])
Ndatainh <- length(Spop$X[Spop$inh==1])

## padding data with NA's, to be the correct length.
data3 <- list(y = table(cut(Count$X,breaks=b),Count$year),
              COV=Covs[1:(nbr.an-1)],
              inipop=iniPop2,
              ##Sdata=Spop$S,
              ##SX=Spop$Xs[,1],
              ##SCOV=Spop$Covs,
              ##Sage=Spop$agef,
              ##Rdata=Spop$R[Spop$S==1],
              ##RX=Spop$Xs[Spop$S==1],
              ##RCOV=Spop$Covs[Spop$S==1],
              ##Idata=Spop$Xs[Spop$inh==1],
              ##IX=Spop$PdsMeres[Spop$inh==1]
              Sdata = c(Spop$S,               rep(as.numeric(NA), 2*Ndatasur-length(Spop$S))),
              Rdata = c(Spop$R[Spop$S==1],    rep(as.numeric(NA), 2*Ndatarec-length(Spop$R[Spop$S==1]))),
              Idata = c(Spop$Xs[Spop$inh==1], rep(as.numeric(NA), 2*Ndatainh-length(Spop$Xs[Spop$inh==1])))
              )

## here we define the size of the data structures to be twice as large.
Consts3 <- list(X=meshpoints,
                pi=pi,
                nbr.an=nbr.an,
                nn=nn,
                ##CHANGE3: using custom distributions here, which allows for variable data.
                ##doubling the constants, to allow for extra room in data structures, 
                ## you can increase these more if necessary.
                Ndatasur=Ndatasur * 2,
                Ndatarec=Ndatarec * 2,
                Ndatainh=Ndatainh * 2
                )

inits3 <-list(Rp= as.numeric(ar$coefficients[1:4]),
              Sp= as.numeric(as$coefficients[1:5]),
              Ip=as.numeric(ai$coefficients[1:2]),
              Ipv= as.numeric(air$coefficients[1]),
              sdY=1,
              ## CHANGE (temp) new inits:
              ## this are for the length of the actual (non-NA) data:
              SdataLength = Ndatasur,
              RdataLength = Ndatarec,
              IdataLength = Ndatainh,
              ## the covariates are also padded with NA's:
              SX = c(Spop$Xs[,1], rep(as.numeric(NA), 2*Ndatasur-length(Spop$S))),
              SCOV = c(Spop$Covs, rep(as.numeric(NA), 2*Ndatasur-length(Spop$S))),
              Sage = c(Spop$agef, rep(as.numeric(NA), 2*Ndatasur-length(Spop$S))),
              RX = c(Spop$Xs[Spop$S==1], rep(as.numeric(NA), 2*Ndatarec-length(Spop$R[Spop$S==1]))),
              RCOV = c(Spop$Covs[Spop$S==1], rep(as.numeric(NA), 2*Ndatarec-length(Spop$R[Spop$S==1]))),
              IX = c(Spop$PdsMeres[Spop$inh==1], rep(as.numeric(NA), 2*Ndatainh-length(Spop$Xs[Spop$inh==1])))
              )


## timings (seconds)
## original: 32
## change1: 
## change2: 14
## change3: 7
system.time(mod3 <- nimbleModel(code = ipm2_3Code, name = 'mod3', constants = Consts3, data = data3, inits = inits3, check = TRUE))

## timings (seconds)
## original: 86
## change1: 83
## change2: 39
## change3: 22
system.time(Cmod3 <- compileNimble(mod3))

spec <- configureMCMC(mod3, monitors=c('Rp','Sp','Ip','Ipv','sdY'),thin=5, useConjugacy=FALSE)
spec$printSamplers()

spec$removeSamplers(c('Sdata', 'Rdata', 'Idata'))
spec$removeSamplers('Rp', print = FALSE)
spec$addSampler(target = 'Rp[1:4]', type = 'RW_block')
spec$removeSamplers('Sp', print = FALSE)
spec$addSampler(target = 'Sp[1:5]', type = 'RW_block')
spec$removeSamplers('Ip', print = FALSE)
spec$addSampler(target = 'Ip[1:2]', type = 'RW_block')
spec$printSamplers()

## timings (seconds)
## origina: 12
## change1:
## change2:
## change3: 9
system.time(mod3MCMC <- buildMCMC(spec))


## timings (seconds)
## original: 46
## change1: 
## change2: 19
## change3: 11
system.time(Cmod3MCMC <- compileNimble(mod3MCMC, project=mod3, resetFunctions = TRUE))


## timings (seconds)
## original: 100 seconds (10,000 samples)
## change1: 
## change2: 36
## change3: 25
set.seed(0)
system.time(Cmod3MCMC$run(10000))   ## my shorter
##system.time(Cmod3MCMC$run(50000))   ## his original

## if you want to remove the warning from running the MCMC,
## then remove all NA's from the covariates (right now there
## are lot's of NA values in them).



MCMCsamplesmod3 <- as.matrix(Cmod3MCMC$mvSamples)[1:2000,]   ## my shorter
##MCMCsamplesmod3 <- as.matrix(Cmod3MCMC$mvSamples)[1:10000,]   ## his original
##save(MCMCsamplesmod3, file = '../data/MCMCsamplesmod3_v2.RData')


## change1: same results
## change2: 
apply(MCMCsamplesmod3,2, mean)
##     Ip[1]      Ip[2]        Ipv      Rp[1]      Rp[2]      Rp[3]      Rp[4] 
##-0.8733089  0.2634541  0.9672236 -1.4506731  1.6371766 -1.3791696  0.8503266 
##     Sp[1]      Sp[2]      Sp[3]      Sp[4]      Sp[5]        sdY 
##-0.3667458  1.1877480  2.8614412 -1.1191852  1.2038239  1.8905793 


## to change the dataset (without rebuilding model or re-compiling),
## just do the things below, with a new data set, then re-run the
## compiled MCMC:

## TO CHANGE DATA: set the length of each data set, each time you change the data
Cmod3$SdataLength <- Ndatasur
Cmod3$RdataLength <- Ndatarec
Cmod3$IdataLength <- Ndatainh

## TO CHANGE DATA: create the new data sets
thisSdata <- c(Spop$S,               rep(as.numeric(NA), 2*Ndatasur-length(Spop$S)))
thisRdata <- c(Spop$R[Spop$S==1],    rep(as.numeric(NA), 2*Ndatarec-length(Spop$R[Spop$S==1])))
thisIdata <- c(Spop$Xs[Spop$inh==1], rep(as.numeric(NA), 2*Ndatainh-length(Spop$Xs[Spop$inh==1])))

## TO CHANGE DATA: set the new data, padded with NA's for setData()
Cmod3$setData(list(Sdata = thisSdata, Rdata = thisRdata, Idata = thisIdata))
unique(removeIndexing(Cmod3$getNodeNames(dataOnly = TRUE)))

## TO CHANGE DATA: set the new covaraites, don't need to pad with NA's for this
Cmod3$SX[1:Ndatasur] <- Spop$Xs[,1]
Cmod3$SCOV[1:Ndatasur] <- Spop$Covs
Cmod3$Sage[1:Ndatasur] <- Spop$agef
Cmod3$RX[1:Ndatarec] <- Spop$Xs[Spop$S==1]
Cmod3$RCOV[1:Ndatarec] <- Spop$Covs[Spop$S==1]
Cmod3$IX[1:Ndatainh] <- Spop$PdsMeres[Spop$inh==1]



Cmod3MCMC$run(10000)





dev.new(); par(mfrow=c(3,2))
plot(MCMCsamplesmod3[ , 1], type = 'l', xlab = 'iteration', ylab = 'Ip1')
plot(MCMCsamplesmod3[ , 2], type = 'l', xlab = 'iteration', ylab = 'Ip2')
plot(MCMCsamplesmod3[ , 3], type = 'l', xlab = 'iteration', ylab = 'Ipv')
plot(MCMCsamplesmod3[ , 4], type = 'l', xlab = 'iteration', ylab = 'Rp1')
plot(MCMCsamplesmod3[ , 5], type = 'l', xlab = 'iteration', ylab = 'Rp2')
plot(MCMCsamplesmod3[ , 6], type = 'l', xlab = 'iteration', ylab = 'Rp3')
dev.new(); par(mfrow=c(3,2))
plot(MCMCsamplesmod3[ , 7], type = 'l', xlab = 'iteration', ylab = 'Rp4')
plot(MCMCsamplesmod3[ , 8], type = 'l', xlab = 'iteration', ylab = 'Sp1')
plot(MCMCsamplesmod3[ , 9], type = 'l', xlab = 'iteration', ylab = 'Sp2')
plot(MCMCsamplesmod3[ , 10], type = 'l', xlab = 'iteration', ylab = 'Sp3')
plot(MCMCsamplesmod3[ , 11], type = 'l', xlab = 'iteration', ylab = 'Sp4')
plot(MCMCsamplesmod3[ , 12], type = 'l', xlab = 'iteration', ylab = 'Sp5')

