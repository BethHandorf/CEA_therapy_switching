

############################################################################################
################# mPCA microsimulation model ###################
############################################################################################

# Microsimulation model for metastatic prostate cancer

# Based on code from:


# 'Microsimulation modeling for health decision sciences using R: a tutorial' 
# Authors: Eline Krijkamp, Fernando Alarid-Escudero, 
#          Eva Enns, Hawre Jalal, Myriam Hunink and  Petros Pechlivanoglou
#
# See GitHub for more information or code updates
# https://github.com/DARTH-git/Microsimulation-tutorial
#


rm(list = ls())  # remove any variables in R's memory 

library(flexsurv)
library(ggplot2)
library(reshape2)

inputdir<-""
resdir<-""



##################################### Model input #########################################
# Model input
n.i   <- 100000                # number of simulated individuals
n.t   <- 87                    # time horizon, 87 3 week cycles (5 years)

#8 states:
#LN1: castrate sensitivie M1 PCA on line 1 therapy
#AEL1: No progression, Adverse event (line 1) 
#PostL1: After line 1
#LN2: castrate resistant M1 PCA on line 2 therapy
#AEL2: No progression, Adverse event (line 2) 
#PostL2: After line 2
#PD: Progression
#D: Death

v.n   <- c("LN1","AEL1","PostL1","LN2","AEL2","PostL2", "PD", "Death")  
n.s   <- length(v.n)           # the number of health states
v.M_1 <- rep("LN1", n.i)         # everyone begins in the healthy state 
d.c   <- d.e <- 0.03/17.3           # equal discounting of costs and QALYs 
#d.c   <- d.e <- 0  


#2 strategies:
# 1. DCT -> AB: Docetaxel followed by Abiraterone
# 2. AB -> DCT: Aberaterone followed by Docetaxel

v.Trt <- c("DCT_AB", "AB_DCT") # store the strategy names

###########
# Get Transition probabilities (by cycle number)
###########

setwd(inputdir)


#Read in survival transition probabilities
#Note: these are pre-computed conditional probabilities (i.e. the probability of transitioning from one )
surv.trans.pr<-read.csv("transition_probs_prog_death.csv")

#remove the first row - this is all zeros (starting point)
#surv.trans.pr<-surv.trans.pr[2:88,]

attach(surv.trans.pr)
#Apply a hazard ratio to line 2 death/progression probabilities based on the time spent in line 1
#HR.max: Maximum hazard ratio applied to line 2 and PD

#t1.lim: Time at which HR decreases to 1.


#Optimal DCT correction 
 HR.max.DCT = HR.max.AB =2.21
 t1.lim.DCT = t1.lim.AB =87

#Optimal AB correction 
# HR.max.DCT = HR.max.AB = 5.07
# t1.lim.DCT = t1.lim.AB = 36
 
#Use optimal correction for each treatment
#
# HR.max.DCT = 2.21
# HR.max.AB = 5.07
#
# t1.lim.DCT = 87
# t1.lim.AB = 36
 
# No correction
# HR.max.DCT = HR.max.AB=1
# t1.lim.DCT = t1.lim.AB=1

#######
# Other transition probabilities (per cycle)
#######

#stopping probabilities
p.stop.DCT.L1=c(0.031, 0.019,0.019, 0.035, 0.043,1,rep(1,81)) #conditional probability of stopping
#16% overall discontinuation rate
p.stop.DCT.L2=c(rep(0.016,10),rep(1,77)) #conditional probability of stopping

p.stop.AB.L1= c(rep(0.12/9,9),rep(0,78)) #Assume 12% probability of stopping over the first 9 cycles
p.stop.AB.L2= c(rep(0.1/9,9),rep(0,78)) #Assume 10% probability of stopping over the first 9 cycles

#Adverse events
p.fatigue.AB = c(rep(0.02/9,9),rep(0,78)) # Assume 2% probability of fatigue, first 6 months only


p.NP_FNP.DCT.L1 = c(rep((0.120+0.061)/6,6),rep(0,81))
p.NP_FNP.DCT.L2 = c(rep((0.163+0.044)/6,6),rep(0,81)) #Although they get DCT for 10 cycles, assume neutropenia would happen in the first 6

#probability of going directly to PD after progression

 p.PD.DCTL1<-0.1
 p.PD.ABL1<-0.1


#p.PD.DCTL1<-0
#p.PD.ABL1<-0


########
#Costs
########

 
cost.pred<-2.8*21/30
 
#cost.AB = (3420*(21/30)+cost.pred+38.62) #brand name 2021

cost.AB = (420*(21/30)+cost.pred+38.62) #generic
#cost.AB = (9368*(21/30)+cost.pred+38.62) #on patent
 
 
cost.DCT=2386+287.27  #Medication plus administration
cost.post.DCT =38.62
cost.ADT = 1145.06+7.63 #Medication plus administration - per cycle cost


#Costs of AEs
cost.NP = 8142.94 #Neutropenia
cost.FNP = 19675.38 #Febrile neutropenia
cost.FAT = 0 #Fatigue
cost.PD = 5477 #Being in PD

########
#Utilities
########
#Line 1
u.L1.AB=0.83
u.L1.ADT=0.83
u.L1.DCT=0.78
u.AE.L1.FAT = 0.78
u.AE.L1.NP.cycle1 = (0.06/(0.06+0.12))*0.47+ (1-(0.06/(0.06+0.12)))*u.L1.DCT #Weighted prob of being in neutropenia/febrile neutropenia
u.AE.L1.NP = 0.78


#Line 2 - use uniform utility decrements of 0.105 for now (0.83-0.725)
u.L2.AB=((0.62+0.83)/2)             #Assume 50% asymptomatic ((0.62+0.83)/2)
u.L2.ADT=((0.62+0.83)/2)  
u.L2.DCT=u.L1.DCT -0.105
u.AE.L2.FAT = u.L2.AB-.05
u.AE.L2.NP.cycle1 = (0.044/(0.044+0.163))*(0.47-.105)+ (1-(0.044/(0.044+0.163)))*u.L2.DCT
u.AE.L2.NP = u.L1.DCT -0.105

u.PD = 0.62 #Distant hormone resistant symptomatic

## To get estimate in LY, use these utilities

# u.L1.AB=1
# u.L1.DCT=1
# u.AE.L1.FAT = 1
# u.AE.L1.NP.cycle1 = 1 #Weighted prob of being in neutropenia/febrile neutropenia
# u.AE.L1.NP = 1
# 

# #Line 2 - 
# u.L2.AB=1        
# u.L2.DCT=1
# u.AE.L2.FAT = 1
# u.AE.L2.NP.cycle1 = 1
# u.AE.L2.NP = 1
# 
# u.PD = 1 #Distant hormone resistant symptomatic


##################################### Functions ###########################################

MicroSim <- function(v.M_1, n.i, n.t, v.n, d.c, d.e, TR.out = TRUE, TS.out = TRUE, Trt = 5, seed = 1) {
# Arguments:  
  # v.M_1:   vector of initial states for individuals 
  # n.i:     number of individuals
  # n.t:     total number of cycles to run the model
  # v.n:     vector of health state names
  # d.c:     discount rate for costs
  # d.e:     discount rate for health outcome (QALYs)
  # TR.out:  should the output include a Microsimulation trace? (default is TRUE)
  # TS.out:  should the output include a transition array between states? (default is TRUE)
  # Trt:     are the n.i individuals receiving treatment? (scalar with a Boolean value, default is FALSE)
  # seed:    starting seed number for random number generator (default is 1)
# Makes use of:
  # Probs:   function for the estimation of transition probabilities
  # Costs:   function for the estimation of cost state values
  # Effs:    function for the estimation of state specific health outcomes (QALYs)
 
  v.dwc <- 1 / (1 + d.c) ^ (0:n.t)   # calculate the cost discount weight based on the discount rate d.c 
  v.dwe <- 1 / (1 + d.e) ^ (0:n.t)   # calculate the QALY discount weight based on the discount rate d.e
  
 # create the matrix capturing the state name/costs/health outcomes for all individuals at each time point
 # Added matricies to keep track of time spent in the current state, and time spent in line2
  m.M <- m.C <- m.E <- m.stime<- m.l2time<-m.pdtime<- matrix(nrow = n.i, ncol = n.t + 1, 
                               dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                               paste("cycle", 0:n.t, sep = " ")))  
  
  m.M[, 1] <- v.M_1                     # indicate the initial health state   
  m.stime[, 1]<-1    #Initial time in cycle is 1
  m.l2time[, 1]<-0    #Initial time in line 2 is 0
  m.pdtime[, 1]<-0    #Initial time in PD is 0
  
  v.L1.tot.time<-rep(0,n.i)
  
  for (i in 1:n.i) {
    set.seed(seed + i)                  # set the seed for every individual for the random number generator
    m.C[i, 1] <- Costs(m.M[i, 1], Trt)  # estimate costs per individual of the initial health state conditional on treatment
    m.E[i, 1] <- Effs (m.M[i, 1], Trt)  # estimate QALYs per individual of the initial health state conditional on treatment
    
    for (t in 1:n.t) {
      v.p <- Probs(m.M[i, t],stime=m.stime[i,t],mtime=t, l2time=m.l2time[i,t], pdtime=m.pdtime[i,t], Trt,L1.tot.time=v.L1.tot.time[i])           # calculate the transition probabilities at cycle t 
      
      m.M[i, t + 1] <- sample(v.n, prob = v.p, size = 1)  # sample the next health state and store that state in matrix m.M 
      m.C[i, t + 1] <- Costs(m.M[i, t + 1], Trt)   # estimate costs per individual during cycle t + 1 conditional on treatment
      m.E[i, t + 1] <- Effs( m.M[i, t + 1], Trt)   # estimate QALYs per individual during cycle t + 1 conditional on treatment
      
      #Update state time, line 2 time
      m.stime[i, t+1] <- as.numeric(m.M[i,t]==m.M[i, t + 1])* (m.stime[i, t]+1) + #Add 1 cycle if cycle doesn't change
                         as.numeric(m.M[i,t]!=m.M[i, t + 1])*1                    #reset if new cycle
      
      #Time spent in L1
      v.L1.tot.time[i]<-v.L1.tot.time[i]+as.numeric(m.M[i, t + 1] %in% c("LN1","AEL1","PostL1"))
      
      #increment L2 counter for every cycle in LN2, AEL2, PostL2
      m.l2time[i, t+1] <- m.l2time[i, t] + as.numeric(m.M[i, t + 1] %in% c("LN2","AEL2","PostL2")) 
      
      #increment PD counter for every cycle in PD
      m.pdtime[i, t+1] <- m.pdtime[i, t] + as.numeric(m.M[i, t + 1] %in% c("PD")) 
      
      
    } # close the loop for the time points 
    if(i/100 == round(i/100,0)) {          # display the progress of the simulation
      cat('\r', paste(i/n.i * 100, "% done", sep = " "))
    }
  } # close the loop for the individuals 
  
  tc <- m.C %*% v.dwc       # total discounted cost per individual
  te <- m.E %*% v.dwe       # total discounted QALYs per individual 
  
  tc_hat <- mean(tc)        # average discounted cost 
  te_hat <- mean(te)        # average discounted QALYs

  if (TS.out == TRUE) {  # create a  matrix of transitions across states
    TS <- paste(m.M, cbind(m.M[, -1], NA), sep = "->") # transitions from one state to the other
    TS <- matrix(TS, nrow = n.i)
    rownames(TS) <- paste("Ind",   1:n.i, sep = " ")   # name the rows 
    colnames(TS) <- paste("Cycle", 0:n.t, sep = " ")   # name the columns 
  } else {
    TS <- NULL
  }
  
  if (TR.out == TRUE) { # create a trace from the individual trajectories
    TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE))))
    TR <- TR / n.i                                       # create a distribution trace
    rownames(TR) <- paste("Cycle", 0:n.t, sep = " ")     # name the rows 
    colnames(TR) <- v.n                                  # name the columns 
  } else {
    TR <- NULL
  }
  
  results <- list(m.M = m.M, m.C = m.C, m.E = m.E, tc = tc, te = te, tc_hat = tc_hat, te_hat = te_hat, TS = TS, TR = TR) # store the results from the simulation in a list  
  return(results)  # return the results
}  # end of the MicroSim function  


#Function to calculate the probability of an event at t2 based on the CI curve 
#applying a hazard ratio based on the observed pfs time in line 1
#HR = HR.max at t0
#decreases linearly to 1 at t1.lim
#
#ci = referece cumulative incidence curve
#t2 = current time in line2
#t1.obs = observed time 
calc.p.event.HR<-function(ci,t2,t1.obs,HR.max,t1.lim) {
  prob<-0
  HR.subject<-1
  HR.subject[t1.obs<=t1.lim]<-(HR.max-t1.obs*(HR.max-1)/t1.lim)
  ci.HR<-1-(1-ci)^HR.subject
  prob[t2==1]<-ci.HR[1]
  prob[t2>1]<-(ci.HR[t2]-ci.HR[t2-1])/(1-ci.HR[t2-1])
  return(prob)
}


#### Probability function
# The Probs function that updates the transition probabilities of every cycle is shown below.


Probs <- function(M_it, stime=1, mtime=1, l2time=0, pdtime=0, Trt=5, L1.tot.time=0) { 
  # M_it:    health state occupied by individual i at cycle t (character variable)
  # stime: state time
  # mtime: model time
  # l2time: line2 time
  # pdtime: time in PD
  # Trt: treatment code

  v.p.it <- rep(0, n.s)     # create vector of state transition probabilities - initialize to zero
  names(v.p.it) <- v.n       # name the vector
  
  #v.p.it is a vector of length 8, each representing the probability of
  #transitioning to (or remaining in) a given state 
  #Elements of v.p.it are
  #v.p.it[8] LN1 (line 1)
  #v.p.it[7] AEL1 (adverse event from line 1 treatment)
  #v.p.it[6] PostL1 (not on line 1 therapy, but no progression)
  #v.p.it[5] LN2 (line 2)
  #v.p.it[4] AEL2 (adverse event from line 2 treatment)
  #v.p.it[3] PostL2 (not on line 2 therapy, but no progression)
  #v.p.it[2] PD (Progressive disease, or "extensive disease" state)
  #v.p.it[1] Death
  

  #Probabilities vary by cycle time, model time, line 2 time
  #Note: probabilities of progressing/dying are already conditional
  
  #Note - assigning values to zero unnecessary b/c vector was initialized to zero
  #Commented out zeros included only for clarity
  
  
  #Allow HR.max and t1.lim to vary by treatment
  
  HR.max<-HR.max.DCT*(Trt==1)+HR.max.AB*(Trt==2)
  t1.lim<-t1.lim.DCT*(Trt==1)+t1.lim.AB*(Trt==2)

  if(M_it=="LN1") {
    #Probability of death
    v.p.it[8]<-(p.DCTL1.death[mtime] * as.numeric(Trt==1) + 
                p.ABL1.death[mtime] * as.numeric(Trt==2) )
    #PD state
    v.p.it[7]<-
      (p.PD.DCTL1)*p.DCTL1.prog[mtime] * as.numeric(Trt==1) + 
         (p.PD.ABL1)*p.ABL1.prog[mtime] * as.numeric(Trt==2)  
    #v.p.it[6]<-0 #Can't get to line 2 AE/post trt from here
    #v.p.it[5]<-0
    #L2
    v.p.it[4]<-
      ((1-p.PD.DCTL1)*p.DCTL1.prog[mtime] * as.numeric(Trt==1) + 
         (1-p.PD.ABL1)*p.ABL1.prog[mtime] * as.numeric(Trt==2))
    #Post L1
    v.p.it[3]<-(1-sum(v.p.it))*
      (p.stop.DCT.L1[mtime] * as.numeric(Trt==1) +
      p.stop.AB.L1[mtime] * as.numeric(Trt==2))
    #AE
    v.p.it[2]<-(1-sum(v.p.it))*
      (p.NP_FNP.DCT.L1[mtime] * as.numeric(Trt==1) + 
         p.fatigue.AB[mtime] * as.numeric(Trt==2) )
    #LN1
    v.p.it[1]<-1-sum(v.p.it)
  }


  if(M_it=="AEL1") {
    v.p.it[8]<-(p.DCTL1.death[mtime] * as.numeric(Trt==1) + 
                  p.ABL1.death[mtime] * as.numeric(Trt==2) )
    v.p.it[7]<-
      ((p.PD.DCTL1)*p.DCTL1.prog[mtime] * as.numeric(Trt==1) + 
         (p.PD.ABL1)*p.ABL1.prog[mtime] * as.numeric(Trt==2)      )  
    #v.p.it[6]<-0 #Can't get to line 2 AE/post trt from here
    #v.p.it[5]<-0
    v.p.it[4]<-
      ((1-p.PD.DCTL1)*p.DCTL1.prog[mtime] * as.numeric(Trt==1) + 
         (1-p.PD.ABL1)*p.ABL1.prog[mtime] * as.numeric(Trt==2))
    v.p.it[3]<-(1-sum(v.p.it))*
      (p.stop.DCT.L1[mtime] * as.numeric(Trt==1) +
         p.stop.AB.L1[mtime] * as.numeric(Trt==2))
    v.p.it[2]<-1-sum(v.p.it)
    #v.p.it[1]<-0
  }
  
  if(M_it=="PostL1") {
    v.p.it[8]<-(p.DCTL1.death[mtime] * as.numeric(Trt==1) + 
                  p.ABL1.death[mtime] * as.numeric(Trt==2))
    v.p.it[7]<-
      ((p.PD.DCTL1)*p.DCTL1.prog[mtime] * as.numeric(Trt==1) + 
         (p.PD.ABL1)*p.ABL1.prog[mtime] * as.numeric(Trt==2)     )   
    #v.p.it[6]<-0 #Can't get to line 2 AE/post trt from here
    #v.p.it[5]<-0
    v.p.it[4]<-
      ((1-p.PD.DCTL1)*p.DCTL1.prog[mtime] * as.numeric(Trt==1) + 
         (1-p.PD.ABL1)*p.ABL1.prog[mtime] * as.numeric(Trt==2) )
    v.p.it[3]<-1-sum(v.p.it)
    #v.p.it[2]<-0
    #v.p.it[1]<-0
  }
  
  

  if(M_it=="LN2") {
    #Probability of death - DCT curves are too high to represent current care - always use AB
    v.p.it[8]<-calc.p.event.HR(cif.ABL2.death,t2=stime,t1.obs=L1.tot.time,HR.max=HR.max,t1.lim=t1.lim) 
    
    #PD
    v.p.it[7]<-
      (calc.p.event.HR(cif.ABL2.prog,t2=stime,t1.obs=L1.tot.time,HR.max=HR.max,t1.lim=t1.lim) * 
         as.numeric(Trt %in% c(1,3)) + 
         calc.p.event.HR(cif.DCTL2.prog,t2=stime,t1.obs=L1.tot.time,HR.max=HR.max,t1.lim=t1.lim) * 
         as.numeric(Trt %in% c(2,4)))
    #PostL2
    v.p.it[6]<-(1-sum(v.p.it))*
      (p.stop.AB.L2[stime] * as.numeric(Trt %in% c(1,3)) +
        p.stop.DCT.L2[stime] * as.numeric(Trt %in% c(2,4)))
    #AEL2
    v.p.it[5]<-(1-sum(v.p.it))*
      (p.fatigue.AB[stime] * as.numeric(Trt %in% c(1,3)) + 
         p.NP_FNP.DCT.L2[stime] * as.numeric(Trt %in% c(2,4)) ) 
        
    #LN2
    v.p.it[4]<-1-sum(v.p.it)
    #v.p.it[3]<-0
    #v.p.it[2]<-0
    #v.p.it[1]<-0
  }
  
  
  if(M_it=="AEL2") {
    #Probability of death
    v.p.it[8]<-calc.p.event.HR(cif.ABL2.death,t2=l2time,t1.obs=L1.tot.time,HR.max=HR.max,t1.lim=t1.lim) 

    #PD
    v.p.it[7]<-
      (calc.p.event.HR(cif.ABL2.prog,t2=l2time,t1.obs=L1.tot.time,HR.max=HR.max,t1.lim=t1.lim) * 
          as.numeric(Trt %in% c(1,3)) + 
          calc.p.event.HR(cif.DCTL2.prog,t2=l2time,t1.obs=L1.tot.time,HR.max=HR.max,t1.lim=t1.lim) * 
          as.numeric(Trt %in% c(2,4)))
    #PostL2
    v.p.it[6]<-(1-sum(v.p.it))*
      (p.stop.AB.L2[l2time] * as.numeric(Trt %in% c(1,3)) +
         p.stop.DCT.L2[l2time] * as.numeric(Trt %in% c(2,4)))
    #AEL2
    v.p.it[5]<-1-sum(v.p.it)
    
    #LN2
    #v.p.it[4]<-0
    #v.p.it[3]<-0
    #v.p.it[2]<-0
    #v.p.it[1]<-0
  }
  
  
  if(M_it=="PostL2") {
    #Probability of death
    v.p.it[8]<-calc.p.event.HR(cif.ABL2.death,t2=l2time,t1.obs=L1.tot.time,HR.max=HR.max,t1.lim=t1.lim)
    
    #PD
    v.p.it[7]<-
      (calc.p.event.HR(cif.ABL2.prog,t2=l2time,t1.obs=L1.tot.time,HR.max=HR.max,t1.lim=t1.lim) * 
         as.numeric(Trt %in% c(1,3)) + 
         calc.p.event.HR(cif.DCTL2.prog,t2=l2time,t1.obs=L1.tot.time,HR.max=HR.max,t1.lim=t1.lim) * 
         as.numeric(Trt %in% c(2,4)))
    #PostL2
    v.p.it[6]<-1-sum(v.p.it)
    #v.p.it[5]<-0
    #v.p.it[4]<-0
    #v.p.it[3]<-0
    #v.p.it[2]<-0
    #v.p.it[1]<-0
  }
  
  if(M_it=="PD") {
    #Probability of death
    #v.p.it[8]<-p.PD.death[stime]
    v.p.it[8]<-calc.p.event.HR(cif.PD.death,t2=pdtime,t1.obs=L1.tot.time,HR.max=HR.max,t1.lim=t1.lim)
    #PD
    v.p.it[7]<-1-sum(v.p.it)
    #v.p.it[6]<-0
    #v.p.it[5]<-0
    #v.p.it[4]<-0
    #v.p.it[3]<-0
    #v.p.it[2]<-0
    #v.p.it[1]<-0
  }
  
  if(M_it=="Death") {
    #Probability of death
    v.p.it[8]<-1
    #v.p.it[7]<-0
    #v.p.it[6]<-0
    #v.p.it[5]<-0
    #v.p.it[4]<-0
    #v.p.it[3]<-0
    #v.p.it[2]<-0
    #v.p.it[1]<-0
  }
  
# Use all.equal instead of == to prevent numeric errors
 ifelse( isTRUE(all.equal(sum(v.p.it),1)) & isTRUE(all.equal(sum(v.p.it<0),0)) , return(v.p.it), print("Probabilities do not sum to 1 or negative probs present")) # return the transition probabilities or produce an error
}       


### Costs function
# The Costs function estimates the costs at every cycle.

#c("LN1","AEL1","PostL1","LN2","AEL2","PostL2", "PD", "Death")



Costs <- function (M_it, Trt = 5, stime=0) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  
  c.it <- 0                                  # by default the cost for everyone is zero 

  #Everyone gets ADT
  c.it[M_it == "LN1"]  <- cost.ADT + 
                      cost.DCT*as.numeric(Trt==1)  +
                      cost.AB* as.numeric(Trt==2)                  
  #No extra cost of fatigue, 
  c.it[M_it == "AEL1"] <- cost.ADT +
                      cost.NP * as.numeric(Trt==1) +  #Cost of treating uncomplicated neutropenial
                      (cost.FNP-cost.NP) * (0.06/(0.06+0.12)) * as.numeric(Trt==1 & stime==1)  #Based on trial data - proportion of with neutropenia starting as febrile
  c.it[M_it == "PostL1"] <- cost.ADT  + cost.post.DCT 
  c.it[M_it == "LN2"]<- cost.ADT +
                        cost.AB *as.numeric (Trt %in% c(1,3)) +
                        cost.DCT *as.numeric (Trt %in% c(2,4))
  c.it[M_it == "AEL2"] <- cost.ADT  +
                      cost.NP * as.numeric(Trt==1) +  #Cost of treating uncomplicated neutropenial
                      (cost.FNP-cost.NP) * (0.044/(0.044+0.163)) * as.numeric(Trt==1 & stime==1)  #Based on trial data - proportion of with neutropenia starting as febrile
  c.it[M_it == "PostL2"] <- cost.ADT  
  c.it[M_it == "PD"] <- cost.PD
  return(c.it)        		                   # return the costs
}


### Health outcome function 
# The Effs function to update the utilities at every cycle.

Effs <- function (M_it, Trt = 5, cl = 0.057534, stime=0) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Trt:  what is the treatment 
  # cl:   cycle length (3 week cycles in years)
  
  u.it <- 0                                  # by default the utility for everyone is zero 
  u.it[M_it == "LN1"]  <- u.L1.AB * (Trt==2) + u.L1.DCT * (Trt==1)                 
  u.it[M_it == "AEL1"] <- u.AE.L1.FAT * (Trt %in% c(2)) + 
                          u.AE.L1.NP *(Trt==1 & stime!=1)+
                          u.AE.L1.NP.cycle1 *(Trt==1 & stime==1) #Weighted utility of FNP and NP
  #Utility accounts for 9 cycle recovery phase after DCT 
  u.it[M_it == "PostL1"] <- u.L1.ADT +
                           (-1)*(u.L1.ADT-u.L1.DCT)/9*(9-stime)*as.numeric(Trt==1 & stime<= 9)   
  #Line 2 utilties
  u.it[M_it == "LN2"]<- u.L2.AB*as.numeric(Trt %in% c(1)) + 
                        u.L2.DCT * as.numeric(Trt %in% c(2))
  u.it[M_it == "AEL2"] <- u.AE.L2.FAT * (Trt %in% c(1)) + 
                        u.AE.L2.NP *(Trt %in% c(2) & stime!=1)+
                        u.AE.L2.NP.cycle1 *(Trt %in% c(2) & stime==1) #Weighted utility of FNP and NP
  #Utility accounts for 9 cycle recovery phase after DCT 
  u.it[M_it == "PostL2"] <- u.L2.ADT +
                        (-1)*(u.L2.ADT-u.L2.DCT)/9*(9-stime)*as.numeric(Trt %in% c(2) & stime<= 9)  
  u.it[M_it == "PD"] <- u.PD
  
  
  QALYs <-  u.it * cl            # calculate the QALYs during cycle t
  return(QALYs)                  # return the QALYs
}


##################################### Run the simulation ##################################

sim_DCT_AB  <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, Trt = 1)  # 1. DCT -> AB

sim_AB_DCT  <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, Trt = 2)  # 2. AB -> DCT


################################# Plot state proportions by cycle (for each strategy)
#Tables summarizing counts by state
sim_DCT_AB$TR
sim_AB_DCT$TR

tr.DCT_AB<-data.frame(sim_DCT_AB$TR)
tr.AB_DCT<-data.frame(sim_AB_DCT$TR)

L1DCT<-tr.DCT_AB$LN1+tr.DCT_AB$AEL1+tr.DCT_AB$PostL1
L1AB<-tr.AB_DCT$LN1+tr.AB_DCT$AEL1+tr.AB_DCT$PostL1

plot(c(0:87), L1DCT, type="l", col="red")
lines(c(0:87), L1AB, col="blue")

# DCT -> AB
tr.DCT_AB<-data.frame(sim_DCT_AB$TR)
tr.DCT_AB$cycle<-c(0:87)

mdata <- melt(tr.DCT_AB, id=c("cycle")) 
names(mdata)<-c("cycle","state","proportion")

ggplot(data=mdata, aes(x=cycle, y=proportion, color=state,group=state)) +
  geom_line()+
  geom_point()+
  scale_color_brewer(palette="Dark2")+
  ggtitle("1. DCT -> AB")


# AB->DCT
tr.AB_DCT<-data.frame(sim_AB_DCT$TR)
tr.AB_DCT$cycle<-c(0:87)

mdata <- melt(tr.AB_DCT, id=c("cycle")) 
names(mdata)<-c("cycle","state","proportion")

ggplot(data=mdata, aes(x=cycle, y=proportion, color=state,group=state)) +
  geom_line()+
  geom_point()+
  scale_color_brewer(palette="Dark2")+
  ggtitle("2. AB -> DCT")







################################# Cost-effectiveness analysis #############################



# store the mean costs (and the MCSE) of each strategy in a new variable v.C (vector costs)
v.C  <- c(sim_DCT_AB$tc_hat, sim_AB_DCT$tc_hat)
sd.C <- c(sd(sim_DCT_AB$tc), sd(sim_AB_DCT$tc))/  sqrt(n.i)

# store the mean QALYs (and the MCSE) of each strategy in a new variable v.E (vector health outcomes)

v.E  <- c(sim_DCT_AB$te_hat, sim_AB_DCT$te_hat)
sd.E <- c(sd(sim_DCT_AB$te), sd(sim_AB_DCT$te))/ sqrt(n.i)


# Summarize results for each 
res<-cbind.data.frame(v.C,sd.C,v.E,sd.E)
rownames(res)<-c("DCT_AB", "AB_DCT")
res<-res[order(res$v.C),]



sum(1-tr.DCT_AB$Death)*3/52

res$delta.C<-c(NA, res$v.C[2:dim(res)[1]]-res$v.C[1:(dim(res)[1]-1)])
res$delta.E<-c(NA, res$v.E[2:dim(res)[1]]-res$v.E[1:(dim(res)[1]-1)])
res$ICER<-res$delta.C/res$delta.E


print(res)

########Number of months given AB 
#First line - proportion of time in L1 * total number of months in study 87*(21/30.4)=60
mean(tr.AB_DCT$LN1+tr.AB_DCT$AEL1)*87*(21/30.4)


#Save results as workspace

setwd(resdir)
save.image("Fitted_model.Rdata")

