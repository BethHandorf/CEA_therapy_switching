# CEA_therapy_switching

This repository contains code to support the manuscript:
Cost-effectiveness analysis for therapy sequence in advanced cancer: A microsimulation approach with application to metastatic prostate cancer
https://arxiv.org/abs/2210.05086

Acknowledgement: this project uses code adapted from
 'Microsimulation modeling for health decision sciences using R: a tutorial' 
 Authors: Eline Krijkamp, Fernando Alarid-Escudero, 
          Eva Enns, Hawre Jalal, Myriam Hunink and  Petros Pechlivanoglou

 See GitHub for more information or code updates
 https://github.com/DARTH-git/Microsimulation-tutorial

# Contents

## 1. Transition probabilities from survival curve
These are probabilities of transitioning between health states based on progression events or death events.  Note that they are pre-computed conditional probabilities based off the cumulative incidence functions (i.e. the probaiblity of transitioning at time i given that the simulated subject survived to time i-1)

It has the follwing columns:  
cycle	- cycle number (discrete)  
p.DCTL1.prog	- Conditional probability of progressing on docetaxel in line 1 therapy  
p.DCTL1.death	- Conditional probability of dying on docetaxel in line 1 therapy  
p.ABL1.prog	- Conditional probability of progressing on abiraterone in line 1 therapy  
p.ABL1.death - Conditional probability of dying on abiraterone in line 1 therapy	  
p.DCTL2.prog	- Conditional probability of progressing on docetaxel in line 2 therapy  
p.DCTL2.death	- Conditional probability of dying on docetaxel in line 2 therapy  
p.ABL2.prog	- Conditional probability of progressing on abiraterone in line 2 therapy  
p.ABL2.death	- Conditional probability of dying on abiraterone in line 2 therapy	  
p.PD.death	- Conditional probability of dying in the extensive disease state therapy	  
cif.ABL2.prog	- cumulative incidence of progression events from abiraterone in line 2 therapy  
cif.ABL2.death	- cumulative incidence of death events from abiraterone in line 2 therapy  
cif.DCTL2.prog	- cumulative incidence of progression events from docetaxel in line 2 therapy  
cif.DCTL2.death	- cumulative incidence of death events from docetaxel in line 2 therapy  
cif.PD.death - cumulative incidence of death events from the extensive disease state  


## 2. Code to run the cost-effectiveness model
mPC_microsimulation.R

Note: you will need to set a working directory so you can read in the transition probability file.

This code runs the microsimulation model with pre-specified inputs (as discussed in the manuscript).  It produces estimates of costs, effects, their differences, and the corresponding ICER.  It also produces several diagnostic plots.

The calibration parameters can be changed (line 77-101)
For Life Years Saved (LYS) instead of QALYS, utilities can be set to zero (lines 156-194)



