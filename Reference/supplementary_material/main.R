#This program is for Data generation and parameter estimation for the Sequential Response Model
#under the 1000 simple size, short sequence and informative prior condition.
#The details are in the paper:  A Sequential Response Model for 
#Analyzing Process Data on Technology-Based Problem-solving Tasks. 


#One needs to change the following line according to where the file is stored on your compute.
source('functions.R')
#####################################################################
#                            Data simulation                        #
#####################################################################
#Problem structure setting
stateSet = c('A','B','C','D','E','F','G')
transSet = c('A#','AB',      
             'BA','B#','BC', 
             'CA','CB','CD', 
             'DB','DC','D#','DE', 
             'EC','E#','EF', 
             'FC','FD','FG') 
u <- length(transSet) #number of state transitions
z <- length(stateSet) - 1 #number of states minus 1
n <- 1000 #number of respondents

#Generate true value
set.seed(1)
theta <- rnorm(n)
#One needs to change the following line according to where the file is stored on your compute.
lambda <- as.matrix(read.csv('lambda19.csv')) #data simulation
dataList <- SimData0(lambda, theta, transSet)
scv <- dataList$scv
stl <- dataList$stl
#save data
write.csv(cbind(scv, theta), 'True Value of Ability and Response Sequence.csv')

#####################################################################
#                        Bayesian estimation                        #
#####################################################################
sd.lam <-  rep(0.1, z)   #sd of proposal distribution
mod1 <- updatePara(stl, 1, 2000, 1, sd.lam) #1 chain, 2000 sampling
sd.lam <-  c(0.3,0.06,0.03,0.05,0.08,0.1)  #adjust sd of proposal distribution
mod1 <- updatePara(stl, 1, 3000, 1, sd.lam) #1 chain, 3000 sampling

M = 5     # number of chains
L = 10000 # length of each chain
B = 3001  # burn-in + 1
mod <- updatePara(stl, M, L, B, sd.lam) #sampling parameters
theMat <- mod$theMat
lamMat <- mod$lamMat

#####################################################################
#                               Results                             #
#####################################################################
library(coda)
library(MCMCvis)
library(ggmcmc)
library(stringr)

#Convert sampling matrix into MCMC list  
mcmclist <- mat.to.mcmc(B,theMat,lamMat)
the.mcmc.list <- mcmclist$the.mcmc.list
lam.mcmc.list <- mcmclist$lam.mcmc.list

#Summary of the Marginal Posterior Distributions
(lam.summary <- cbind( MCMCsummary(lam.mcmc.list), MCMCsummary(lam.mcmc.list,HPD = T)[,c(3,4)]))
the.summary <- cbind( MCMCsummary(the.mcmc.list), MCMCsummary(the.mcmc.list,HPD = T)[,c(3,4)])
write.csv(the.summary,'Summary for Abilities.csv')
write.csv(lam.summary,'Summary for State Transition Parameters.csv')

#Estimation accuracy 
CorMse(theta,lambda,the.summary[,1],lam.summary[,1],scv)

#Plot
plotMCMC(lamMat,lam.summary[,1],transSet) # trace plot
densplot(lam.mcmc.list)
#ggs_density(ggs(lam.mcmc.list)) # enough plot space required

#Model fit
modFit(theMat,lamMat,the.summary[,1],lam.summary[,1])
ppmcPlot(transSet,theMat,lamMat)


