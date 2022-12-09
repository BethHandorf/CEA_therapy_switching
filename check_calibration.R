# Microsimulation model for metastatic prostate cancer

#check calibration of survival curves from simulated model against true trial data
#Required inputs: Fitted microsimulation model from "mPC_microsimulation.R"
#Overall survival datasets from target trials


library(survival)

#Load the fitted microsimulation model



setwd("./")
load("Fitted_model.Rdata")


#Read in inferred overall survival data

OS.DCT<-read.csv("DCT_OS_inferred.csv")
OS.AB<-read.csv("AB_OS_inferred.csv")
#OS.L2.AB<-read.csv("AB_L2_OS_inferred.csv")
#OS.L2.DCT<-read.csv("DCT_L2_OS_inferred.csv")
#OS.PD<-read.csv("PD_OS_inferred.csv")


###################################

#Model based results
#3 year survival rate
1-tr.AB_DCT[52,8]
#5 year survival rate
1-tr.AB_DCT[87,8]

#Trial based
s.AA<-survfit(Surv(Time*(30.4/21),Event)~1, data=OS.AB)
s.est.AA<-summary(s.AA, times=c(52,87))
s.est.AA$surv


#3 year death rate
1-tr.DCT_AB[52,8]
#5 year death rate
1-tr.DCT_AB[87,8]

#Trial-based
s.DCT<-survfit(Surv(Time*(30.4/21),Event)~1, data=OS.DCT)
s.est.DCT<-summary(s.DCT, times=c(52,87))
s.est.DCT$surv


#Plot both fitted curves agianst trial data
par(mfrow=c(1,2))
plot(survfit(Surv(Time*(30.4/21),Event)~1, data=OS.DCT), conf.int=TRUE,xlim=c(0,87), main="DCT", xlab="Cycle (3 weeks)",ylab="Overall survival")
lines(tr.DCT_AB$cycle,1-tr.DCT_AB$Death, col="firebrick4", lwd=2)
legend("bottomleft",c("Trial data","Model-based"),lty=1, col=c("black","firebrick4"))

plot(survfit(Surv(Time*(30.4/21),Event)~1, data=OS.AB), conf.int=TRUE,xlim=c(0,87), main="AA", xlab="Cycle (3 weeks)",ylab="Overall survival")
lines(tr.AB_DCT$cycle,1-tr.AB_DCT$Death, col="firebrick4", lwd=2)
legend("bottomleft",c("Trial data","Model-based"),lty=1, col=c("black","firebrick4"))

#Komogorov-Smirnov test
#DCT
cdf.trial<-1-summary(s.DCT, times = c(0:87))[['surv']]
ks.test(cdf.trial, tr.DCT_AB$Death)

cdf.trial<-1-summary(s.AA, times = c(0:87))[['surv']]
ks.test(cdf.trial, tr.AB_DCT$Death[1:59])


