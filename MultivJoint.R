# INLA install:
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA) # Bayesian inference with INLA
library(rstanarm) # MCMC with STAN
library(PermAlgo) # Permutation algorithm to generate survival times dependent on time-varying covariates
library(mvtnorm)
inla.setOption(inla.mode="experimental")
set.seed(1)

# inla.setOption(pardiso.license = "~/pardiso.lic")
# pardiso provides a high performance computing environment with 
# parallel computing support using OpenMP (not avail. on windows)
# see https://pardiso-project.org/r-inla/

# data simulation
nsujet=500 # number of indivduals

# Y1 (continuous)
b1_0=0.2 # intercept
b1_1=-0.1 # slope
b1_2=0.1 # continuous covariate
b1_3=-0.2 # binary covariate
b1_e=0.4 # residual error
# Y2 (counts)
b2_0=3 # intercept
b2_1=-0.1 # slope
b2_2=0.1 # continuous covariate
b2_3=-0.2 # binary covariate
# Y3 (binary)
b3_0=1 # intercept
b3_1=-1 # slope
b3_2=1 # continuous covariate
b3_3=-1 # binary covariate
# Survival
b1_s=0.5 # effect of Y1's linear predictor on the risk of event
b2_s=-0.2 # effect of Y2's linear predictor on the risk of event
b3_s=0.3 # effect of Y3's linear predictor on the risk of event

gapLongi=0.1 # gap between longi measurements
gap=0.01 # used to generate a lot of time points because the permutation
# algorithm chooses among those time points to define survival times
followup=5 # follow-up time
mestime=seq(0,followup,gap) # measurement times
timesLongi=mestime[which(round(mestime-round(mestime/gapLongi,0)*gapLongi,6)==0)] # visit times
time=rep(mestime, nsujet) # time column
nmesindiv=followup/gap+1 # max. number of individual measurements
nmesy= nmesindiv*nsujet # max. total number of longi measurements
# the number is reduced because of censoring by the terminal event
idY<-rep(1:nsujet, each=nmesindiv) # individual id

# random effects variance and covariance matrix
Sigma <- matrix(c(0.16, 0.03, 0.02, 0.04, 0.00,
                  0.03, 0.09, 0.03, 0.00, -0.06, 
                  0.02, 0.03, 0.25, 0.08, 0.05,
                  0.04, 0.00, 0.08, 0.16, 0.04, 
                  0.00, -0.06, 0.05, 0.04, 0.25),ncol=5,nrow=5)
MVnorm <- rmvnorm(nsujet, rep(0, 5), Sigma)
b1_int <- rep(MVnorm[,1], each=nmesindiv) # random intercept Y1
b1_slo <- rep(MVnorm[,2], each=nmesindiv) # random slope Y1
b2_int <- rep(MVnorm[,3], each=nmesindiv) # random intercept Y2
b2_slo <- rep(MVnorm[,4], each=nmesindiv) # random slope Y2
b3_int <- rep(MVnorm[,5], each=nmesindiv) # random intercept Y3

ctsX=rep(rnorm(nsujet,1, 0.5) , each=nmesindiv) # continuous covariate
binX=rep(rbinom(nsujet,1, 0.5), each=nmesindiv) # binary covariate

# linear predictors
linPredY1 <- (b1_0+b1_int) + (b1_1+b1_slo)*time + b1_2*ctsX + b1_3*binX
linPredY2 <- (b2_0+b2_int) + (b2_1+b2_slo)*time + b2_2*ctsX + b2_3*binX
linPredY3 <- (b3_0+b3_int) + b3_1*time + b3_2*ctsX + b3_3*binX
# continuous outcome Y1
Y1 <- rnorm(nmesy, linPredY1, b1_e)
# count outcome Y2
Y2 <- rpois(nmesy, exp(linPredY2))
# binary outcome Y3
Y3 <- rbinom(nmesy,1, exp(linPredY3)/(1+exp(linPredY3)))

# Permutation algorithm to generate survival times that depends on the linear predictors
DatTmp <- permalgorithm(nsujet,nmesindiv,Xmat=matrix(c(linPredY1, linPredY2, linPredY3), nrow=nsujet*nmesindiv),
                        eventRandom = round(rexp(nsujet, 0.003)+1,0), # ~40% death
                       censorRandom=runif(nsujet,1,nmesindiv), # uniform random censoring
                       XmatNames=c("linPredY1", "linPredY2", "linPredY3"), # association
                       betas=c(b1_s,b2_s,b3_s)) # association parameters

# extract last line for each Id (= death/censoring time)
DatTmp2=DatTmp[c(which(diff(DatTmp[,"Id"])==1), dim(DatTmp)[1]), c("Id","Event","Stop")]
DatTmp2$deathTimes <- mestime[DatTmp2$Stop+1] # deathtimes
survDat <- DatTmp2[, c("Id", "deathTimes", "Event")]
DatTmp$time <- mestime[DatTmp$Start+1] # measurements times of the biomarker
DatTmp$Uid <- paste(DatTmp$Id, DatTmp$time) # unique identifier to match covariates and observed biomarker values
longDat3 <- merge(DatTmp[,c("Uid", "Id", "time")], cbind("Uid"=paste(idY,time), ctsX, binX, Y1, Y2, Y3), by=c("Uid"))
longDat <- sapply(longDat3[longDat3$time%in%timesLongi,-1], as.numeric)
longDat <- as.data.frame(longDat[order(longDat[, "Id"], longDat[, "time"]),])
summary(survDat) # survival dataset
summary(longDat) # longitudinal dataset

# Model fit
NL <- nrow(longDat)
NS <- nrow(survDat)
# Prepare data for INLA
# Cox model structure for survival (with Bayesian smooth splines for the baseline hazard, i.e., "rw2")
cox_ext = inla.coxph(YS ~  -1+Intercept, control.hazard=list(model="rw2", scale.model=TRUE, diagonal=1e-2, 
                  constr=TRUE, hyper=list(prec=list(prior="pc.prec", param=c(0.5,0.01)))),
                  data = c(list(YS = inla.surv(time = c(survDat$deathTimes), event = c(survDat$Event)),
                  Intercept = rep(1, NS), eta1 = survDat$Id, eta2 = NS+survDat$Id, eta3 = NS*2+survDat$Id)))

# time weight for time-dependent covariates in survival (i.e., fixed and random slope)
t.weight = cox_ext$data$baseline.hazard.time +  0.5 * cox_ext$data$baseline.hazard.length
ns_cox = nrow(cox_ext$data) # number of intervals for survival
NY = NL + ns_cox
IDcox = cox_ext$data$expand..coxph # random effects id for survival
cox_ext$data$eta1 = 1:ns_cox # unique id for shared linear predictor Y1
cox_ext$data$eta2 = ns_cox+1:ns_cox # unique id for shared linear predictor Y2
cox_ext$data$eta3 = ns_cox*2+1:ns_cox # unique id for shared linear predictor Y3

# match binary and continuous covariates with survival time intervals
IDBC = merge(cox_ext$data, unique(cbind("expand..coxph"=longDat$Id, "ctsX"=longDat$ctsX, "binX"=longDat$binX)), by="expand..coxph")

covariates <- list(
InteY1 = c(rep(1,NY), rep(NA, NY*2)), # intercept Y1
TIMEY1 = c(longDat$time,t.weight, rep(NA, NY*2)), # time Y1
ctsXY1 = c(longDat$ctsX, IDBC$ctsX, rep(NA, NY*2)), # ctsX Y1
binXY1 = c(longDat$binX, IDBC$binX, rep(NA, NY*2)), # binX Y1
InteY2 = c(rep(NA,NY), rep(1, NY), rep(NA, NY)), # intercept Y2
TIMEY2 = c(rep(NA,NY), longDat$time, t.weight, rep(NA, NY)), # time Y2
ctsXY2 = c(rep(NA,NY), longDat$ctsX, IDBC$ctsX, rep(NA, NY)), # ctsX Y2
binXY2 = c(rep(NA,NY), longDat$binX, IDBC$binX, rep(NA, NY)), # binX Y2
InteY3 = c(rep(NA,NY*2), rep(1, NY)), # intercept Y3
TIMEY3 = c(rep(NA,NY*2), longDat$time, t.weight), # time Y3
ctsXY3 = c(rep(NA,NY*2), longDat$ctsX, IDBC$ctsX), # ctsX Y3
binXY3 = c(rep(NA,NY*2), longDat$binX, IDBC$binX), # binX Y3
IDY1 = c(longDat$Id, IDcox, rep(NA, NY*2)), # random intercept Y1
IDY1_s = c(NS+longDat$Id, NS+IDcox, rep(NA, NY*2)), # random slope Y1
WY1_s = c(longDat$time, t.weight, rep(NA, NY*2)), # weight random slope Y1
IDY2 = c(rep(NA, NY), NS*2+longDat$Id, NS*2+IDcox, rep(NA, NY)), # random intercept Y2
IDY2_s = c(rep(NA, NY), NS*3+longDat$Id, NS*3+IDcox, rep(NA, NY)), # random slope Y2
WY2_s = c(rep(NA, NY), longDat$time, t.weight, rep(NA, NY)), # weight random slope Y2
IDY3 = c(rep(NA, NY*2), NS*4+longDat$Id,  NS*4+IDcox), # random intercept Y3
u1 = c(rep(NA, NL), cox_ext$data$eta1, rep(NA, NY*2)), # Y1 association survival
w1 = c(rep(NA, NL), rep(-1, ns_cox), rep(NA, NY*2)), # Y1 weight association survival
u2 = c(rep(NA, NL+NY), cox_ext$data$eta2, rep(NA, NY)), # Y2 association survival
w2 = c(rep(NA, NL+NY), rep(-1, ns_cox), rep(NA, NY)), # Y2 weight association survival
u3 = c(rep(NA, NL+NY*2), cox_ext$data$eta3), # Y3 association survival
w3 = c(rep(NA, NL+NY*2), rep(-1, ns_cox))) # Y3 weight association survival

Y.joint <- list(Y1 = c(longDat$Y1, rep(NA, NY*2+ns_cox)), # Y1
                Y.eta1 = c(rep(NA, NL), rep(0, ns_cox),rep(NA, NY*2)), # Y1 association with survival
                Y2 = c(rep(NA, NY), longDat$Y2,rep(NA, NY+ns_cox)), # Y2
                Y.eta2 = c(rep(NA, NL+NY), rep(0, ns_cox),rep(NA, NY)), # Y2 association with survival
                Y3 = c(rep(NA, NY*2), longDat$Y3,rep(NA, ns_cox)), # Y3
                Y.eta3 = c(rep(NA, NL+NY*2), rep(0, ns_cox))) # Y3 association with survival

jointdf = data.frame(covariates, Y.joint)
joint.data_cox <- c(as.list(inla.rbind.data.frames(jointdf, cox_ext$data)), cox_ext$data.list)
joint.data_cox$E..coxph[(NY+1):(NL*2+ns_cox)] <- 1 # for Poisson distribution (Y2)
Yjoint = joint.data_cox[c("Y1", "Y.eta1", "Y2", "Y.eta2", "Y3", "Y.eta3", "y..coxph")] # outcomes (longi and survival)
joint.data_cox$Y <- Yjoint

# update formula from the cox model structure to add longitudinal part
formulaJ= update(cox_ext$formula, Yjoint ~ . -1 + InteY1 + TIMEY1 + ctsXY1 + binXY1 +
                   InteY2 + TIMEY2 + ctsXY2 + binXY2 +
                   InteY3 + TIMEY3 + ctsXY3 + binXY3 +
                   f(IDY1, model="iidkd", order=5, n=NS*5,constr=FALSE, 
                     hyper = list(theta1 = list(param = c(10,rep(1,5),rep(0,10)))))+
                   f(IDY1_s, WY1_s, copy="IDY1") + f(IDY2, copy="IDY1") +
                   f(IDY2_s, WY2_s, copy="IDY1") + f(IDY3, copy="IDY1") +
                   f(u1, w1, model="iid", hyper = list(prec = list(initial = -6, fixed=TRUE)), constr=F)+
                   f(eta1, copy="u1", hyper = list(beta = list(fixed = FALSE, param = c(0, 0.16), initial = 1)))+
                   f(u2, w2, model="iid", hyper = list(prec = list(initial = -6, fixed=TRUE)), constr=F)+
                   f(eta2, copy="u2", hyper = list(beta = list(fixed = FALSE, param = c(0, 0.16), initial = 1)))+
                   f(u3, w3, model="iid", hyper = list(prec = list(initial = -6, fixed=TRUE)), constr=F)+
                   f(eta3, copy="u3", hyper = list(beta = list(fixed = FALSE, param = c(0, 0.16), initial = 1))))

# inla function call
JMinla <- inla(formulaJ,family = c("gaussian", "gaussian", "poisson", "gaussian", "binomial", "gaussian", cox_ext$family),
               data=joint.data_cox,
               control.fixed = list(mean=0, prec=0.16, mean.intercept=0, prec.intercept=0.16),
               control.family = list(
                 list(), list(hyper = list(prec = list(initial = 12, fixed=TRUE))),
                 list(), list(hyper = list(prec = list(initial = 12, fixed=TRUE))),
                 list(), list(hyper = list(prec = list(initial = 12, fixed=TRUE))),
                 list()),
               E = joint.data_cox$E..coxph,
               # int.strategy="eb" for empirical Bayes, default is full Bayesian (remove next line to switch to default).
               control.inla = list(int.strategy="eb"))
print(summary(JMinla))

# Get random effects covariance from Cholesky
MC_samples <- inla.iidkd.sample(10^4, JMinla, "IDY1", return.cov=TRUE) 
VarCov <- matrix(unlist(MC_samples), nrow = 5^2)
VarCovMeans <- matrix(rowMeans(VarCov),5,5);round(VarCovMeans,2)


# Model fit with rstanarm
if(F){ # Long computation time!
JMrs <- stan_jm(
  formulaLong = list(
    Y1 ~ time+ctsX+binX + (time | Id),
    Y2 ~ time+ctsX+binX + (time | Id),
    Y3 ~ time+ctsX+binX + (1 | Id)),
  formulaEvent = survival::Surv(deathTimes, Event) ~ 1,
  dataLong = longDat, dataEvent = survDat,
  family = list(gaussian, poisson(link = "log"), binomial(link="logit")),
  time_var = "time",
  control = list(max_treedepth = 20),
  priorLong_intercept = normal(0, 2.5),
  priorLong = normal(0, 2.5),
  priorEvent_assoc = normal(0, 2.5),
  # default is chains = 4 and iter = 2000
  chains = 1, iter = 1000, seed = 12345) 
summary(JMrs)
}


