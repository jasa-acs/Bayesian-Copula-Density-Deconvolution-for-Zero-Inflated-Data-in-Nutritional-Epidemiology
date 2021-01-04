
#####################################
### Bayesian Copula Deconvolution ###
#####################################

# Codes accompanying "Bayesian Copula Density Deconvolution for Zero-Inflated Data with Applications in Nutritional Epidemiology" by Sarkar, Pati, Mallick and Carroll.
# Codes written by Abhra Sarkar (abhra.sarkar@utexas.edu), last modified on Dec 15, 2019, in Austin, TX


rm(list=ls())

################
### Set Path ###
################

setwd("//Users//as85655//Dropbox//Projects\ with\ Dr.\ C//X\ ME\ Density\ Deconvolution\ for\ Zero\ Inflated\ Data//JASA\ Submission//Codes")



###################
### Set Seed No ###
###################

seedno <- 1
set.seed(seedno)  



######################################
### Add Libraries and Source Files ###
######################################

# Load libraries
library(foreach)
library(doParallel)
library(ks)
library(mvtnorm) 
library(MCMCpack) 
library(msm)
library(corpcor)
library(RColorBrewer)

# Detect cores to run univariate methods in parallel 
numCores <- detectCores()-1

# Source the following R files
source("Bayes_Copula_Decon_Functs.R")
source("Bayes_Copula_Decon_Univariate_Regular.R")
source("Bayes_Copula_Decon_Univariate_Episodic.R")
source("Bayes_Copula_Decon_MVT.R")



#####################
### Run Main File ###
#####################

# Some arguments are shared by functions nested with the main function Bayes_Copula_Decon_MVT() called here. 
# To facilitate passing of the shared arguments to the nested functions, set their values outside in the global environment.
# To run the function on a different data set (with default choices for other parameters), 
# simply change the csv filename below to one containg that data set. 
# To conduct a simulation study, run this file repeatedly on multiple synthetic data sets. 
# For the real data analysis, change the input csv file to one containing the real data set. 
# Descriptions of the plots produced by Bayes_Copula_Decon_MVT() function can be found in its body. 


filename = "Simulated_Data.csv"
Run_Univariate_Models = TRUE
Plot_Results = TRUE
Save_Results = TRUE
simsize = 5000
burnin = 3000
simsize_univ = 500
burnin_univ = 300
K.t = 10
z.x.max = 10
z.u.max = 10
z.x.max.univ  =  5
z.u.max.univ = 5

Bayes_Copula_Decon_MVT(filename,Run_Univariate_Models,Plot_Results,Save_Results,simsize,burnin,simsize_univ,burnin_univ,K.t,z.x.max,z.u.max,z.x.max.univ,z.u.max.univ)

