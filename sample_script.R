############# Bayesian Multivariate Matched Proportions #############
# Sample script for both half-Cauchy and vMCMC BMMaPs	
#	from Meyer and Knutson 2018 									
#																	
# Created:  10/26/2018												
# Modified: 02/14/2019												
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

## set priors for multivariate item response ##
pvSpec	<- list(a_lambda = 0.001, a_gamma = 0.001, b_gamma = 0.001,
				a_omega = 0.001, a_eta = 0.001, b_eta = 0.001,
				tol = 1e-10, maxIter = 10, mut = rep(0, ncol(X1)),
				mub = rep(0, ncol(X1)), mua1 = rep(0, ncol(X1)),
				mua2 = rep(0, ncol(X1)))
	
## run multivariate item response model ##
B		<- 100	# 1000
burnin	<- 100
modelVB	<- bmir(X1, X2, link = 'logit', B = B, burn = burnin, up = 100, priorVar = "VB",
				pvSpec = pvSpec)

## chain length ##
#	for illustration, we set B and burnin to 100, but for an acutal analysis
#	we recommend B = 1000 for both models and burnin = 1000 for the half-Cauchy

## changing link ##
# 	the link argument can be set to logit, probit, or cloglog

## compare model results ##
summary(modelHC)
summary(modelVB)

## set priors for independent multinomials model ##
prior		<- list(alpha = matrix(1/2, nrow = ncol(X1), ncol = 4))

## run independent multinomials model ##
B			<- 1000
burnin		<- 1000
modelIMM	<- bim(X1, X2, B = B, burn = burnin, prior = prior)

summary(modelIMM)


prior		<- list(alpha = matrix(c(rep(1/2, 4*3), rep(1, 4*2)), nrow = ncol(X1), ncol = 4, byrow = TRUE))

## run independent multinomials model ##
B			<- 3000
burnin		<- 1000
prior		<- list(alpha = matrix(c(rep(1/2, 4*3)), nrow = ncol(X1), ncol = 4, byrow = TRUE))
modelIMM	<- bim(X1, X2, B = B, burn = burnin, prior = prior)

summary(modelIMM)

dim(modelIMM$theta)

apply(modelIMM$theta, c(2,3), mean)

exp(cmhExact(X1, X2)$est)

gbmmps(X1, X2)

exp(gbmmps(X1, X2)$est)


setwd('/Users/mjm556/Dropbox/Research/Drafts/Matched Proportions/Data')
library(readstata13)
socw		<- read.dta13('knutson_data_wide.dta')

pc		<- socw[,c('PC_school', 'PC_CPS', 'PC_MH', 'PC_JJ', 'PC_DD')]
colnames(pc)	<- c('SS', 'CPS', 'MH', 'JJ', 'DD')

sc		<- socw[,c('Psych_school', 'Psych_CPS', 'Psych_MH', 'Psych_JJ', 'Psych_DD')]
colnames(sc)	<- c('SS', 'CPS', 'MH', 'JJ', 'DD')

pc		<- socw[,c('PC_DD', 'PC_MH', 'PC_JJ', 'PC_CPS', 'PC_school')]
colnames(pc)	<- c('DD', 'MH', 'JJ', 'CPS', 'SS')

sc		<- socw[,c('Psych_DD', 'Psych_MH', 'Psych_JJ', 'Psych_CPS', 'Psych_school')]
colnames(sc)	<- c('DD', 'MH', 'JJ', 'CPS', 'SS')


getCounts(pc, sc, 5)

cmhExact(pc, sc)

pc2		<- socw[,c('PC_school', 'PC_MH')]
colnames(pc2)	<- c('SS', 'MH')

sc2		<- socw[,c('Psych_school', 'Psych_MH')]
colnames(sc2)	<- c('SS', 'MH')

X1			<- pc
X2			<- sc
B			<- 5000
burnin		<- 5000

modelIMM	<- bim(X1, X2, B = B, burn = burnin)

summary(modelIMM)
apply(modelIMM$theta, c(2,3), mean)

sparseTables1	<- apply(modelIMM$MCMCsp$prior$sparsityCheck, 2, sum)

paste(paste(which(sparseTables1 > 0), collapse = ", "), '.', sep = '')

## set up for new BMIR sampler ##
X1			<- pc
X2			<- sc
B			<- 5000
burnin		<- 5000

n		<- nrow(X1)
K		<- ncol(X1)
Y		<- cbind(X1, X2)
gDat	<- data.frame(y	= c(t(Y)), s = rep(1:(2*K), n), id = rep(1:n, each = 2*K))

X		<- model.matrix(y ~ factor(s) - 1, data = gDat)
U		<- model.matrix(y ~ factor(id) - 1, data = gDat)

XsCols	<- c(which(apply(matrix(unlist(lapply(getCounts(X1, X2, K), function(x) x == 0)), nrow = 4, byrow = TRUE), 2, sum) > 0), which(apply(matrix(unlist(lapply(getCounts(X1, X2, K), function(x) x == 0)), nrow = 4, byrow = TRUE), 2, sum) > 0) + K)

Xn	<- X[,-XsCols]
Xs	<- X[,XsCols]


