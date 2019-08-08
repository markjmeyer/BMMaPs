### load bmmaps.R ###
source('~/Documents/Code/BMMaPs/bmmaps.R')

setwd('/Users/mjm556/Dropbox/Research/Drafts/Matched Proportions/Data')
library(readstata13)
socw		<- read.dta13('knutson_data_wide.dta')

pc		<- socw[,c('PC_DD', 'PC_MH', 'PC_JJ', 'PC_CPS', 'PC_school')]
colnames(pc)	<- c('DD', 'MH', 'JJ', 'CPS', 'SS')

sc		<- socw[,c('Psych_DD', 'Psych_MH', 'Psych_JJ', 'Psych_CPS', 'Psych_school')]
colnames(sc)	<- c('DD', 'MH', 'JJ', 'CPS', 'SS')


getCounts(pc, sc, 5)


## set up for new BMIR sampler ##
library(msm)
library(mvtnorm)
# library(MCMCpack)

X1			<- pc
X2			<- sc
B			<- 5000
burnin		<- 5000
v			<- 8
A			<- 5
B			<- 5

n		<- nrow(X1)
K		<- ncol(X1)
Y		<- cbind(X1, X2)
gDat	<- data.frame(y	= c(t(Y)), s = rep(1:(2*K), n), id = rep(1:n, each = 2*K))

X		<- model.matrix(y ~ factor(s) - 1, data = gDat)
U		<- model.matrix(y ~ factor(id) - 1, data = gDat)

XsCols	<- c(which(apply(matrix(unlist(lapply(getCounts(X1, X2, K), function(x) x == 0)), nrow = 4, byrow = TRUE), 2, sum) > 0), which(apply(matrix(unlist(lapply(getCounts(X1, X2, K), function(x) x == 0)), nrow = 4, byrow = TRUE), 2, sum) > 0) + K)

Xn	<- X[,-XsCols]
Xs	<- X[,XsCols]

Yv	<- gDat$y
N	<- length(Yv)
N1	<- sum(Yv)
N0	<- N - N1

## parameter storage ##

# mean components #
betan	<- matrix(0, nrow = B + burnin, ncol = ncol(Xn))
betas	<- matrix(0, nrow = B + burnin, ncol = ncol(Xs))
theta	<- matrix(0, nrow = B + burnin, ncol = n)

# variance components #
lambda	<- matrix(1, nrow = B + burnin, ncol = N)
sig2i	<- matrix(1, nrow = B + burnin, ncol = n)
omega	<- vector('numeric', length = B + burnin)
gamma2	<- matrix(1, nrow = B + burnin, ncol = ncol(Xs))
alpha	<- vector('numeric', length = B + burnin)

# latent variable #
Z		<- rep(0, N)

# set starting values #
omega[1]	<- 1
alpha[1]	<- 1

# prior parameters #
a_lam	<- v/2 + 1/2
a_sig	<- 1
a_ome	<- n/2 + 1/2
a_gam	<- 1
a_alp	<- ncol(Xs)/2 + 1/2

for(b in 2:(B + burnin)){
	## update non-sparse betas ##
	Zn		<- Z - Xs%*%betas[b-1,] - U%*%theta[b-1,]
	Lambda	<- diag(lambda[b-1,])
	var_bn	<- solve(t(Xn)%*%Lambda%*%Xn)
	mu_bnz	<- t(Xn)%*%Lambda%*%Zn
	mu_bn	<- var_bn%*%mu_bnz
	
	betan[b,]	<- rmvnorm(1, mu_bn, var_bn)
	
	## update sparse betas ##
	Zs		<- Z - Xn%*%betan[b-1,] - U%*%theta[b-1,]
	sig_gam	<- diag(gamma2[b-1,])
	prec_bs	<- t(Xs)%*%Lambda%*%Xs + solve(sig_gam)
	var_bs	<- solve(prec_bs)
	mu_bsz	<- t(Xs)%*%Lambda%*%Zs
	mu_bs	<- var_bs%*%mu_bsz
	
	betas[b,]	<- rmvnorm(1, mu_bs, var_bs)
	
	## update thetas ##
	Zt		<- Z - Xn%*%betan[b-1,] - Xs%*%betas[b-1,]
	sig_th	<- diag(sig2i[b-1,])
	prec_th	<- t(U)%*%Lambda%*%U + solve(sig_th)
	var_th	<- solve(prec_th)
	mu_thz	<- t(U)%*%Lambda%*%Zt
	mu_th	<- var_th%*%mu_thz
	
	theta[b,]	<- rmvnorm(1, mu_th, var_th)
	
	## update lambdas ##
	b_lam		<- v/2 + (1/2)*(Z - (Xn%*%betan[b-1,] + Xs%*%betas[b-1,] + U%*%theta[b-1,]))^2
	lambda[b,]	<- rgamma(N, a_lam, b_lam)
	
	## update sigmas ##
	b_sig		<- 1/omega[b-1] + (1/2)*(theta[b-1,]^2)
	sig2i[b,]	<- rinvgamma(n, a_sig, b_sig)
	
	## update omega ##
	b_ome		<- 1/(A^2) + sum(1/sig2i[b-1,])
	omega[b]	<- rinvgamma(1, a_ome, b_ome)
	
	## update gammas ##
	b_gam		<- 1/alpha[b-1] + (1/2)*(betas[b-1,]^2)
	gamma2[b,]	<- rinvgamma(1, a_gam, b_gam)
	
	## update alpha ##
	b_alp		<- 1/(B^2)	+ sum(1/gamma2[b,])
	alpha[b]	<- rinvgamma(1, a_alp, b_alp)
	
	## update Z's ##
	XB	<- Xn%*%betan[b-1,] + Xs%*%betas[b-1,] + U%*%theta[b-1,]
	
	# yi == 0 #
	Z[which(Yv == 0)]	<- rtnorm(N0, mean = XB[which(Yv == 0)], sd = sqrt(lambda[b-1, which(Yv == 0)]), upper = 0)
	
	# yi == 1 #
	Z[which(Yv == 1)]	<- rtnorm(N1, mean = XB[which(Yv == 1)], sd = sqrt(lambda[b-1, which(Yv == 1)]), lower = 1)
	
}




