library(INLAjoint)
library(JM) # This package contains the dataset
set.seed(1)
data(pbc2) # dataset
# extract some variable of interest without missing values
Longi <- na.omit(pbc2[, c("id", "years", "status","drug","age",
                          "sex","year","serBilir","SGOT", "albumin", "edema",
                          "platelets", "alkaline","spiders", "ascites")])
Surv <- Longi[c(which(diff(as.numeric(Longi[,which(colnames(Longi)=="id")]))==1),
                length(Longi[,which(colnames(Longi)=="id")])),-c(7:10, 12:16)]
Surv$death <- ifelse(Surv$status=="dead",1,0) # competing event 1
Surv$trans <- ifelse(Surv$status=="transplanted",1,0) # competing event 2

DTH <- inla.surv(time = Surv$years, event = Surv$death) # survival outcomes
TSP <- inla.surv(time = Surv$years, event = Surv$trans)

# set up natural cubic splines for longitudinal markers's trajectories
Nsplines <- ns(Longi$year, knots=c(1,4))
f1 <- function(x) predict(Nsplines, x)[,1]
f2 <- function(x) predict(Nsplines, x)[,2]
f3 <- function(x) predict(Nsplines, x)[,3]

JM_INLA <- joint(formSurv = list(DTH ~ drug, TSP ~ drug),
                 formLong = list(serBilir ~ (1 + f1(year) + f2(year) + f3(year)) * drug +
                                            (1 + f1(year) + f2(year) + f3(year) | id),
                                 SGOT ~ (1 + f1(year) + f2(year) + f3(year)) * drug +
                                        (1 + f1(year) + f2(year) + f3(year) | id),
                                 albumin ~ (1 + year) * drug + (1 + year | id),
                                 platelets ~ (1 + f1(year) + f2(year) + f3(year)) * drug +
                                             (1 + f1(year) + f2(year) + f3(year) | id),
                                 spiders ~ (1 + year) * drug + (1 + year | id)),
                 dataLong = Longi, id = "id", timeVar = "year",
                 family = c("lognormal", "lognormal", "gaussian", "poisson", "binomial"),
                 basRisk = c("rw2", "rw1"), NbasRisk = 15, assoc = list(c("CV_CS", "CV"),
                                 c("CV", ""), c("CV", "CV"), c("CV", "CV"), c("CV", "")),
                 control=list(int.strategy="eb", cfg=T))
summary(JM_INLA)
# standard deviations and correlations for residual error and random effects (instead of variance-covariance)
summary(JM_INLA, sdcor=T)

