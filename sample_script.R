############# Bayesian Multivariate Matched Proportions #############
# Sample script for both half-Cauchy and vMCMC BMMaPs	
#	from Meyer and Knutson 2018 									
#																	
# Created:  10/26/2018												
# Modified: 10/26/2018												
#																	
# By: Mark J Meyer													
#																	
#####################################################################

### load bmmaps.R ###
source('~/Documents/Code/BMMaPs/bmmaps.R')

### set data directory ###
setwd('~/Documents/Code/BMMaPs/Data')

### load data examples ##

## 50% Discordant Pairs ##
#  Each .RData file contains and array
# 	with three dimensions. The first is
#	the pairs, second is outcome, and third
#	is varying types of discordancy
#	
#  Data below is from the 50% concordant
#	setting (as % of original data)
X1ash	<- readRDS('X1ash.RData')
X2ash	<- readRDS('X2ash.RData')

## analysis where both outcomes have lower left discordancy of 1 ##
d	<- 1			# change to 2 or 3 to vary level of discordancy
X1	<- X1ash[,,d]	# matrix of outcomes for first treatment
X2	<- X2ash[,,d]	# matrix of outcomes for second treatment

## additional simulated data ##
#	all artificial datasets are included in the data folder
#	X1asf.RData and X2asf.RData	for 100% concordant (as % of original data)
#	X1as50.RData and X2as50.RData for 150% concordant (as % of original data)
#	X1asD.RData and X2asD.RData	for 200% concordant (as % of original data)

## set priors for half-Cauchy BMMaPs ##
pvSpec	<- list(lambda = 300, omega = 100,
				tTune = 0.04, nuTau = 5^2,
				sTune = c(0.0015, 0.002), nuSig = 3^2)

## run model ##
B		<- 100	# 1000
burnin	<- 100	# 1000 
modelHC	<- bmmaps(X1, X2, link = 'logit', B = B, burn = burnin, up = 100, priorVar = "HC",
				pvSpec = pvSpec)

## set priors for vMCMC BMMaPs ##
pvSpec	<- list(a_lambda = 0.001, a_gamma = 0.001, b_gamma = 0.001,
				a_omega = 0.001, a_eta = 0.001, b_eta = 0.001,
				tol = 1e-10, maxIter = 10)
	
## run model ##
B		<- 100	# 1000
burnin	<- 100
modelVB	<- bmmaps(X1, X2, link = 'logit', B = B, burn = burnin, up = 100, priorVar = "VB",
				pvSpec = pvSpec)

## chain length ##
#	for illustration, we set B and burnin to 100, but for an acutal analysis
#	we recommend B = 1000 for both models and burnin = 1000 for the half-Cauchy

## changing link ##
# 	the link argument can be set to logit, probit, or cloglog

## compare model results ##
summary(modelHC)
summary(modelVB)
