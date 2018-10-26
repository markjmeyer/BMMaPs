############# Bayesian Multivariate Matched Proportions #############
# R functions to run Bayesian Multivariate Matched Proportions		
#	from Meyer and Knutson 2018 									
#																	
# Created:  07/25/2018												
# Modified: 10/26/2018												
#																	
# By: Mark J Meyer													
#																	
#####################################################################

## require libraries ##
suppressPackageStartupMessages(require(MASS))
suppressPackageStartupMessages(require(MCMCpack))
suppressPackageStartupMessages(require(numbers))
suppressPackageStartupMessages(require(msm))

## functions for varying links ##

exbit	<- function(x){exp(x)/(1 + exp(x))}

dexbit	<- function(x){exp(x)/((1 + exp(x))^2)}

domexb	<- function(x){-exp(x)/((1 + exp(x))^2)}

evCDF	<- function(x){1-exp(-exp(x))}

evPDF	<- function(x){exp(x)*exp(-exp(x))}

f	<- function(x, link){
	switch(link,
		logit	= exbit(x),
		probit	= pnorm(x),
		cloglog	= evCDF(x)
	)
}

dF	<- function(x, link){
	switch(link,
		logit	= dexbit(x),
		probit	= dnorm(x),
		cloglog	= evPDF(x)
	)
}

dFb	<- function(x, link){
	switch(link,
		logit	= domexb(x),
		probit	= -dnorm(x),
		cloglog	= -evPDF(x)
	)
}

## convergence diagnostics function ##

getGelman	<- function(mcmcOut, chains = 4){
	if(chains != 4 & chains != 2){
		stop('Split chains in half or quarters only')
	}
	K		<- ncol(mcmcOut)
	B		<- nrow(mcmcOut)
	GRR		<- matrix(0, ncol = 2, nrow = K)
	cuts		<- split(1:B, rep(1:chains, each = B/chains))
	
	for(k in 1:K){
		if(chains == 4){
			mcmc1	<- mcmc(mcmcOut[cuts[[1]],k])
			mcmc2	<- mcmc(mcmcOut[cuts[[2]],k])
			mcmc3	<- mcmc(mcmcOut[cuts[[3]],k])
			mcmc4	<- mcmc(mcmcOut[cuts[[4]],k])
			gdb		<- gelman.diag(mcmc.list(mcmc1, mcmc2, mcmc3, mcmc4))
		
			GRR[k,]	<- gdb$psrf
		} else {
			mcmc1	<- mcmc(mcmcOut[cuts[[1]],k])
			mcmc2	<- mcmc(mcmcOut[cuts[[2]],k])
			gdb		<- gelman.diag(mcmc.list(mcmc1, mcmc2))
		
			GRR[k,]	<- gdb$psrf		
		}
	}
	colnames(GRR)	<- c('Point est.', 'Upper C.I.')
	rownames(GRR)	<- 1:K

	return(GRR)
	
}

## functions for ARS ##

lptb	<- function(c1, x1, x2, c2, a1, a2, s, link, mutb = 0){
	in1		<- x1*log(f(c1 + c2 + a1, link)) + (1-x1)*log(1 - f(c1 + c2 + a1, link))
	in2		<- x2*log(f(c1 + c2 + a2, link)) + (1-x2)*log(1 - f(c1 + c2 + a2, link))
	out		<- sum(in1 + in2) - (1/(2*s))*((c1 - mutb)^2)
	return(out)
}

dlptb	<- function(c1, x1, x2, c2, a1, a2, s, link, mutb = 0){
	in1		<- (x1/f(c1 + c2 + a1, link))*dF(c1 + c2 + a1, link) + ((1-x1)/(1-f(c1 + c2 + a1, link)))*dFb(c1 + c2 + a1, link)
	in2		<- (x2/f(c1 + c2 + a2, link))*dF(c1 + c2 + a2, link) + ((1-x2)/(1-f(c1 + c2 + a2, link)))*dFb(c1 + c2 + a2, link)
	out		<- sum(in1 + in2) - (1/s)*(c1 - mutb)
	return(out)
}

lpalp	<- function(aj, th, be, x, link, s = 1e500, mua = 0){
	in1		<- x*log(f(th + be + aj, link)) + (1 - x)*log(1 - f(th + be + aj, link))
	out		<- sum(in1) - (1/(2*s)*(aj - mua)^2)
	return(out)
}

dlpalp	<- function(aj, th, be, x, link, s = 1e500, mua = 0){
	in1		<- (x/f(th + be + aj, link))*dF(th + be + aj, link) + ((1 - x)/(1 - f(th + be + aj, link)))*dFb(th + be + aj, link)
	out		<- sum(in1) - (1/(s)*(aj - mua))
	return(out)
}

lpsigk	<- function(s, n, K, th, nu){
	out	<- -((n*K)/2)*log(s) - (1/(2*s))*sum(th^2) - log(s + nu)
	return(out)
}

lptau	<- function(tau, n, K, be, nu){
	out	<- -((n*K)/2)*log(tau) - (n/(2*tau))*sum(be^2) - log(tau + nu)
	return(out)
}

zftb	<- function(yfixed, x1, x2, c2, a1, a2, s, link){
	yf0f		<- head(yfixed, n=-1)
	yf1f		<- tail(yfixed, n=-1)
	zfixed	<- vector('numeric', length = length(yfixed)-1)
	for(zl in 1:length(zfixed)){
		yf0			<- yf0f[zl]
		yf1			<- yf1f[zl]
		zfixed[zl] <- yf0 + (lptb(yf0, x1, x2, c2, a1, a2, s, link) - lptb(yf1, x1, x2, c2, a1, a2, s, link) + (yf1 - yf0)*dlptb(yf1, x1, x2, c2, a1, a2, s, link)) / (dlptb(yf1, x1, x2, c2, a1, a2, s, link) - dlptb(yf0, x1, x2, c2, a1, a2, s, link))
	}
	return(zfixed)	
}

zfalp	<- function(yfixed, th, be, x, link){
	yf0f		<- head(yfixed, n=-1)
	yf1f		<- tail(yfixed, n=-1)
	zfixed	<- vector('numeric', length = length(yfixed)-1)
	for(zl in 1:length(zfixed)){
		yf0			<- yf0f[zl]
		yf1			<- yf1f[zl]
		zfixed[zl]	<- yf0 + (lpalp(yf0, th, be, x, link) - lpalp(yf1, th, be, x, link) + (yf1 - yf0)*dlpalp(yf1, th, be, x, link)) / (dlpalp(yf1, th, be, x, link) - dlpalp(yf0, th, be, x, link))
	}
	return(zfixed)	
}

utb <- function(y, yfixed, x1, x2, c2, a1, a2, s, link, ymin = -Inf, ymax = Inf){
	res				<- rep(0, length(y))
	zfixed			<- zftb(yfixed, x1, x2, c2, a1, a2, s, link)

	piecewise.idx	<- findInterval(y, c(ymin, zfixed, ymax))
	npieces 		<- length(zfixed) + 2
	for(pidx in 1:npieces){
		yp 			<- y[piecewise.idx == pidx]
		xx 			<- lptb(yfixed[pidx], x1, x2, c2, a1, a2, s, link) + (yp - yfixed[pidx])*dlptb(yfixed[pidx], x1, x2, c2, a1, a2, s, link)
		res[piecewise.idx == pidx]	<- xx
	}
	return(res)
}

ualp <- function(y, yfixed, th, be, x, link, ymin = -Inf, ymax = Inf){
	res				<- rep(0, length(y))
	zfixed			<- zfalp(yfixed, th, be, x, link)

	piecewise.idx	<- findInterval(y, c(ymin, zfixed, ymax))
	npieces 		<- length(zfixed) + 2
	for(pidx in 1:npieces){
		yp 			<- y[piecewise.idx == pidx]
		xx 			<- lpalp(yfixed[pidx], th, be, x, link) + (yp - yfixed[pidx])*dlpalp(yfixed[pidx], th, be, x, link)
		res[piecewise.idx == pidx]	<- xx
	}
	return(res)
}

stb		<- function(vals, yfixed, x1, x2, c2, a1, a2, s, link){
	zfixed		<- zftb(yfixed, x1, x2, c2, a1, a2, s, link)

	zlen 		<- length(zfixed)
	pct 			<- numeric(length(vals))
	norm.const	<- 0
	for(zi in 0:zlen){
		if(zi == 0){
			zm	<- -Inf
		} else {
			zm	<- zfixed[zi]
		}

		if(zi == zlen){
			zp	<- Inf
		} else {
			zp	<- zfixed[zi+1]
		}

		yp			<- yfixed[zi+1]
		ds			<- exp(lptb(yp, x1, x2, c2, a1, a2, s, link))/dlptb(yp, x1, x2, c2, a1, a2, s, link) * ( exp((zp - yp)*dlptb(yp, x1, x2, c2, a1, a2, s, link)) - exp((zm - yp)*dlptb(yp, x1, x2, c2, a1, a2, s, link)) )

		cidx		<- zm < vals & vals <= zp
		hidx		<- vals > zp

		pct[cidx]	<- pct[cidx] + exp(lptb(yp, x1, x2, c2, a1, a2, s, link))/dlptb(yp, x1, x2, c2, a1, a2, s, link) * ( exp((vals[cidx] - yp)*dlptb(yp, x1, x2, c2, a1, a2, s, link)) - exp((zm - yp)*dlptb(yp, x1, x2, c2, a1, a2, s, link)) )
		pct[hidx]	<- pct[hidx] + ds

		norm.const <- norm.const + ds
	}

	l <- list(pct = pct / norm.const, norm.const = norm.const)
	return(l)
}

salp		<- function(vals, yfixed, th, be, x, link){
	zfixed		<- zfalp(yfixed, th, be, x, link)

	zlen 		<- length(zfixed)
	pct 			<- numeric(length(vals))
	norm.const	<- 0
	for(zi in 0:zlen){
		if(zi == 0){
			zm	<- -Inf
		} else {
			zm	<- zfixed[zi]
		}

		if(zi == zlen){
			zp	<- Inf
		} else {
			zp	<- zfixed[zi+1]
		}

		yp			<- yfixed[zi+1]
		ds			<- exp(lpalp(yp, th, be, x, link))/dlpalp(yp, th, be, x, link) * ( exp((zp - yp)*dlpalp(yp, th, be, x, link)) - exp((zm - yp)*dlpalp(yp, th, be, x, link)) )

		cidx		<- zm < vals & vals <= zp
		hidx		<- vals > zp

		pct[cidx]	<- pct[cidx] + exp(lpalp(yp, th, be, x, link))/dlpalp(yp, th, be, x, link) * ( exp((vals[cidx] - yp)*dlpalp(yp, th, be, x, link)) - exp((zm - yp)*dlpalp(yp, th, be, x, link)) )
		pct[hidx]	<- pct[hidx] + ds

		norm.const <- norm.const + ds
	}

	l <- list(pct = pct / norm.const, norm.const = norm.const)
	return(l)
}

stb.sample <- function(samp.size, yfixed, x1, x2, c2, a1, a2, s, link, ymin){
	zfixed			<- zftb(yfixed, x1, x2, c2, a1, a2, s, link)
	gp				<- stb(zfixed, yfixed, x1, x2, c2, a1, a2, s, link)
	zpct			<- gp$pct		
	norm.const		<- gp$norm.const
	ub				<- c(0, zpct, 1)

	unif.samp		<- runif(samp.size)

	fidx				<- findInterval(unif.samp, ub)
	num.intervals	<- length(ub) - 1
	zlow				<- c(ymin, zfixed)
	res				<- rep(NaN, length(unif.samp))
	for(ii in 1:num.intervals){
		ui	<- unif.samp[ fidx == ii ]

		if(length(ui) == 0){
			next
		}

		yp	<- yfixed[ii]
		zm	<- zlow[ii]
		tmp	<- (ui - ub[ii]) * dlptb(yp, x1, x2, c2, a1, a2, s, link) * norm.const / exp(lptb(yp, x1, x2, c2, a1, a2, s, link)) + exp( (zm - yp)*dlptb(yp, x1, x2, c2, a1, a2, s, link) )
		tmp	<- yp + log(tmp) / dlptb(yp, x1, x2, c2, a1, a2, s, link)
		res[ fidx == ii ] <- tmp
	}
	return(res)
}

salp.sample <- function(samp.size, yfixed, th, be, x, link, ymin){
	zfixed			<- zfalp(yfixed, th, be, x, link)
	gp				<- salp(zfixed, yfixed, th, be, x, link)
	zpct				<- gp$pct		
	norm.const		<- gp$norm.const
	ub				<- c(0, zpct, 1)

	unif.samp		<- runif(samp.size)

	fidx				<- findInterval(unif.samp, ub)
	num.intervals	<- length(ub) - 1
	zlow				<- c(ymin, zfixed)
	res				<- rep(NaN, length(unif.samp))
	for(ii in 1:num.intervals){
		ui	<- unif.samp[ fidx == ii ]

		if(length(ui) == 0){
			next
		}

		yp	<- yfixed[ii]
		zm	<- zlow[ii]
		tmp	<- (ui - ub[ii]) * dlpalp(yp, th, be, x, link) * norm.const / exp(lpalp(yp, th, be, x, link)) + exp( (zm - yp)*dlpalp(yp, th, be, x, link) )
		tmp	<- yp + log(tmp) / dlpalp(yp, th, be, x, link)
		res[ fidx == ii ] <- tmp
	}
	return(res)
}


tbUpdate	<- function(draws, yfixed, x1, x2, c2, a1, a2, s, link, ymin, ymax, upFix = TRUE, dlpThresh){
	thb	<- vector('numeric', length = draws)

	count	<- 0
	tt		<- 0
	while(tt < draws){
		count	<- count + 1

		# draw Ys from sthe #
		# Ys	<- sthe.sample(samp.size, yfixed, x1, x2, a1, a2, s, link, ymin)
		Ys	<- stb.sample(1, yfixed, x1, x2, c2, a1, a2, s, link, ymin)

		# draw U(0,1) #
		U	<- runif(1)

		test <- U < exp(lptb(Ys, x1, x2, c2, a1, a2, s, link) - utb(Ys, yfixed, x1, x2, c2, a1, a2, s, link, ymin, ymax))

		if(test){
			tt			<- tt + 1
			thb[tt]		<- Ys
		} else{
			dlpEval		<- dlptb(Ys, x1, x2, c2, a1, a2, s, link)
			addFix		<- abs(dlpEval) < abs(dlpThresh)
			if(addFix & upFix){
				yb 			<- c(yfixed, Ys)
				yfixed		<- sort(yb)
			}
		}
	}
	l	<- list(thetab = thb, yfixedb = yfixed, count = count)
	return(l)		
}

alphaUpdate	<- function(samp.size, yfixed, th, be, x, link, ymin, ymax, upFix = TRUE, dlpThresh){
	alb	<- vector('numeric', length = samp.size)

	count	<- 0
	tt		<- 0
	while(tt < samp.size){
		count	<- count + 1

		# draw Ys from sthe #
		Ys		<- salp.sample(samp.size, yfixed, th, be, x, link, ymin)

		# draw U(0,1) #
		U		<- runif(1)

		test	<- U < exp(lpalp(Ys, th, be, x, link) - ualp(Ys, yfixed, th, be, x, link, ymin, ymax))

		if(test){
			tt			<- tt + 1
			alb[tt]		<- Ys
		} else{
			dlpEval		<- dlpalp(Ys, th, be, x, link)
			addFix		<- abs(dlpEval) < abs(dlpThresh)
			if(addFix & upFix){
				yb 			<- c(yfixed, Ys)
				yfixed		<- sort(yb)
			}
		}
	}
	l	<- list(alphab = alb, yfixedb = yfixed, count = count)
	return(l)		
}

yfTBSelect	<- function(tInt = c(-2,2), inc = 0.005, cent = 0.15, x1, x2, c2, a1, a2, s, link){
	tseq		<- matrix(seq(tInt[1], tInt[2], by = inc), ncol = 1)
	lp		<- apply(tseq, 1, function(z) lptb(z, x1, x2, c2, a1, a2, s, link))

	cselp	<- cumsum(exp(lp))
	if(min(cselp/sum(exp(lp))) > cent){
		stop('percentile too small (thetas or betas), consider increasing cent or searchInt')
	}
	low		<- tseq[max(which(cselp/sum(exp(lp)) <= cent))]
	mid		<- tseq[max(which(cselp/sum(exp(lp)) <= 0.5))]
	high		<- tseq[max(which(cselp/sum(exp(lp)) <= 1-cent))]
	dlpl		<- dlptb(low, x1, x2, c2, a1, a2, s, link)
	dlpm		<- dlptb(mid, x1, x2, c2, a1, a2, s, link)
	dlph		<- dlptb(high, x1, x2, c2, a1, a2, s, link)
	
	l	<- list(yfix = c(low, mid, high), deriv = c(dlpl, dlpm, dlph))
	
	return(l)

}

yfAlphaSelect	<- function(aInt = c(-10, 10), inc = 0.025, cent = 0.15, th, be, x, link){
	aseq		<- matrix(seq(aInt[1], aInt[2], by = inc), ncol = 1)
	lp		<- apply(aseq, 1, function(z) lpalp(z, th, be, x, link))
	
	cselp	<- cumsum(exp(lp))
	if(min(cselp/sum(exp(lp))) > cent){
		stop('percentile too small (alphas), consider increasing cent or searchInt')
	}
	low		<- aseq[max(which(cselp/sum(exp(lp)) <= cent))]
	mid		<- aseq[max(which(cselp/sum(exp(lp)) <= 0.5))]
	high		<- aseq[max(which(cselp/sum(exp(lp)) <= 1-cent))]
	dlpl		<- dlpalp(low, th, be, x, link)
	dlpm		<- dlpalp(mid, th, be, x, link)
	dlph		<- dlpalp(high, th, be, x, link)

	l	<- list(yfix = c(low, mid, high), deriv = c(dlpl, dlpm, dlph))

	return(l)	
	
}

## functions for updating variance components ##

sigma2UpdateIG	<- function(samp.size, n, K, th, as, bs){
	ab	<- (n*K)/2 + as
	bb	<- bs + (K/2)*sum(th^2)
	out	<- rinvgamma(samp.size, ab, bb)
	return(out)
}

tau2UpdateIG	<- function(samp.size, n, K, be, at, bt){
	ab	<- (n*K)/2 + at
	bb	<- bt + (n/2)*sum(be^2)
	out	<- rinvgamma(samp.size, ab, bb)
	return(out)
}

sigma2UpdateGA	<- function(samp.size, n, K, th, as, bs){
	ab	<- (n*K)/2 + as
	bb	<- bs + (K/2)*sum(th^2)
	out	<- rgamma(samp.size, ab, bb)
	return(out)
}

tau2UpdateGA	<- function(samp.size, n, K, be, at, bt){
	ab	<- (n*K)/2 + at
	bb	<- bt + (n/2)*sum(be^2)
	out	<- rgamma(samp.size, ab, bb)
	return(out)
}

sigma2UpdateHC	<- function(samp.size, sig, sMean, sTune, n, K, th, nuSig){
	sstar	<- rtnorm(samp.size, mean = sMean, sd = sTune, lower = 0)
	gstar	<- dtnorm(sstar, mean = sig, sd = sTune, lower = 0, log = TRUE)
	gsig	<- dtnorm(sig, mean = sstar, sd = sTune, lower = 0, log = TRUE)
	lpsstar	<- lpsigk(sstar, n, K, th, nuSig)
	lpssig	<- lpsigk(sig, n, K, th, nuSig)
				
	sprob	<- exp(lpsstar - lpssig - gstar + gsig)
	arsig	<- 0
	Us		<- runif(1)
	if(Us < min(1, sprob)){
		sig		<- sstar
		arsig	<- 1
	}
	l	<- list(est = sig, ar = arsig)
	return(l)
}

tau2UpdateHC	<- function(samp.size, tau, tMean, tTune, n, K, be, nuTau){
	tstar	<- rtnorm(samp.size, mean = tMean, sd = tTune, lower = 0)
	gstar	<- dtnorm(tstar, mean = tau, sd = tTune, lower = 0, log = TRUE)
	gtau	<- dtnorm(tau, mean = tstar, sd = tTune, lower = 0, log = TRUE)
	lptstar	<- lptau(tstar, n, K, be, nuTau)
	lpttau	<- lptau(tau, n, K, be, nuTau)
	
	tprob	<- exp(lptstar - lpttau - gstar + gtau)
	artau	<- 0
	Ut		<- runif(1)
	if(Ut < min(1, tprob)){
		tau		<- tstar
		artau	<- 1
	}
	l	<- list(est = tau, ar = artau)
	return(l)
}

sigma2VB	<- function(n, K, th, a_lambda, a_gamma, b_gamma, B_gamma = 0.01, tol = 1e-5, maxIter = 10){
	log_delta	<- 0.5
	iter		<- 0
	log_prev	<- 0
	post		<- matrix(0, maxIter, K+2)
	B_lambda	<- vector('numeric', length = K)
	colnames(post)	<- c(paste('E(lambda',1:K,')', sep = ''), 'E(gamma)', 'KL')
	
	while(log_delta > tol & iter <= maxIter){
		iter		<- iter + 1
		for(k in 1:K){
			B_lambda[k]	<- a_gamma/B_gamma + (1/2)*sum(th[,k]^2)
		}
		E_lambda	<- (n/2 + a_lambda)/B_lambda
		B_gamma		<- b_gamma + sum(E_lambda)
		E_gamma		<- a_gamma/B_gamma
		
		log_p		<- -E_gamma*sum(E_lambda)
		log_delta	<- abs(log_p - log_prev)
		log_prev	<- log_p
		post[iter,]	<- c(E_lambda, E_gamma, log_p)	
	}
	
	if(iter >= maxIter){
		warning('max iteration reached in VB step, consider increasing maxIter')
	}
	
	out	<- list(est = 1/E_lambda, E_lambda = E_lambda, B_lambda = B_lambda, ar = iter)
	return(out)
	
}

tau2VB		<- function(n, K, be, a_omega, a_eta, b_eta, B_eta = 0.01, tol = 1e-5, maxIter = 10){
	log_delta	<- 0.5
	iter		<- 0
	log_prev	<- 0
	post		<- matrix(0, maxIter, 3)
	colnames(post)	<- c('E(omega)', 'E(eta)', 'KL')
	
	while(log_delta > tol & iter <= maxIter){
		iter		<- iter + 1
		B_omega		<- a_eta/B_eta + (n/2)*sum(be^2)
		E_omega		<- ((n*K)/2 + a_omega)/B_omega
		B_eta		<- b_eta + E_omega
		E_eta		<- a_eta/B_eta
		
		log_p		<- -E_eta*E_omega
		log_delta	<- abs(log_p - log_prev)
		log_prev	<- log_p
		post[iter,]	<- c(E_omega, E_eta, log_p)
	}

	if(iter >= maxIter){
		warning('max iteration reached in VB step, consider increasing maxIter')
	}
	
	out	<- list(est = 1/E_omega, E_omega = E_omega, B_omega = B_omega, ar = iter)
	return(out)
	
}

## extract cell counts from pop averaged tables ##

getCounts	<- function(X1, X2, K){
	n21		<- vector('numeric', length = K)
	n12		<- vector('numeric', length = K)
	n11		<- vector('numeric', length = K)
	n22		<- vector('numeric', length = K)

	for(j in 1:K){
		rawtab	<- table(X2[,j], X1[,j])
		temp	<- matrix(0, nrow = 2, ncol = 2)
		temp[1:dim(rawtab)[1], 1:dim(rawtab)[2]]	<- rawtab
	
		patab	<- matrix(c(temp[2,2], temp[1,2], temp[2,1], temp[1,1]), nrow = 2)

		n21[j]	<- patab[2,1]
		n12[j]	<- patab[1,2]
	
		n11[j]	<- patab[1,1]
		n22[j]	<- patab[2,2]
	}
	
	l	<- list(n21 = n21, n12 = n12, n11 = n11, n22 = n22)
	return(l)
}

## primary function for BMMaPs ##

bmmaps	<- function(X1, X2, link = 'logit', B, burn = 100, up = 1000, searchIntT = NULL, incT = 0.005, centT = 0.1, searchIntB = NULL, incB = 0.005, centB = 0.1, searchIntA = NULL, incA = 0.01, centA = 0.1, priorVar = 'VB', pvSpec = NULL, ymin = -Inf, ymax = Inf, upFix = TRUE, dlpThresh = 100){

	if(ncol(X1) != ncol(X2)){
		stop('column dimensions of X1 and X2 must agree')
	}

	if(nrow(X1) != nrow(X2)){
		stop('row dimensions of X1 and X2 must agree')
	}

	if(link != 'logit' & link != 'probit' & link != 'cloglog'){
		stop('improper link specified use either logit, probit, or cloglog')
	}

	if(priorVar != 'VB' & priorVar != 'FX' & priorVar != 'GA' & priorVar != 'HC'){
		stop('improper prior variance specified use either VB, HC, GA, or FX')
	}
		
	if(B < 1){
		stop('B must be at least 1')
	}
	
	if(is.null(searchIntA)){
		if(link == 'probit'){
			searchIntA	<- c(-7,7)
		} else if(link == 'cloglog'){
			searchIntA	<- c(-13, 3)
		} else {
			searchIntA	<- c(-10, 10)
		}
	}
	
	if(is.null(searchIntB)){
		if(link == 'cloglog'){
			searchIntB	<- c(-5, 3)
		} else {
			searchIntB	<- c(-5, 5)
		}
	}

	if(is.null(searchIntT)){
		if(link == 'cloglog'){
			searchIntT	<- c(-5, 3)
		} else {
			searchIntT	<- c(-5, 5)
		}
	}
	
	if(is.null(pvSpec)){
		if(priorVar == 'VB'){
			a_lambda	<- 0.001
			a_gamma		<- 0.001
			b_gamma		<- 0.001
			a_omega		<- 0.001
			a_eta		<- 0.001
			b_eta		<- 0.001
			tol			<- 1e-5
			maxIter		<- 10
			
			pvSpec	<- list(a_lambda = a_lambda, a_gamma = a_gamma, b_gamma = b_gamma,
							a_omega = a_omega, a_eta = a_eta, b_eta = b_eta,
							tol = tol, maxIter = maxIter)
		} else if(priorVar == 'GA'){
			lambda		<- 100
			omega		<- 100
			at			<- 1
			bt			<- 1
			as			<- 1
			bs			<- 1
			
			pvSpec		<- list(lambda = lambda, omega = omega,
								at = at, bt = bt, as = as, bs = bs)
		} else if(priorVar == 'HC'){
			lambda		<- 100
			omega		<- 100
			tTune		<- 0.01
			nuTau		<- 2
			sTune		<- rep(0.01, ncol(X1))
			nuSig		<- 2
			
			pvSpec		<- list(lambda = lambda, omega = omega,
								tTune = tTune, nuTau = nuTau,
								sTune = sTune, nuSig = nuSig)
		} else {
			lambda		<- 100
			omega		<- 100
			
			pvSpec		<- list(lambda = lambda, omega = omega)
		}
	} else {
		if(priorVar == 'VB'){
			a_lambda	<- pvSpec$a_lambda
			a_gamma		<- pvSpec$a_gamma
			b_gamma		<- pvSpec$b_gamma
			a_omega		<- pvSpec$a_omega
			a_eta		<- pvSpec$a_eta
			b_eta		<- pvSpec$b_eta
			tol			<- pvSpec$tol
			maxIter		<- pvSpec$maxIter
		} else if(priorVar == 'GA'){
			lambda		<- pvSpec$lambda
			omega		<- pvSpec$omega
			at			<- pvSpec$at
			bt			<- pvSpec$bt
			as			<- pvSpec$as
			bs			<- pvSpec$bs
		} else if(priorVar == 'HC'){
			lambda		<- pvSpec$lambda
			omega		<- pvSpec$omega
			tTune		<- pvSpec$tTune
			nuTau		<- pvSpec$nuTau
			sTune		<- pvSpec$sTune
			nuSig		<- pvSpec$nuSig
		} else {
			lambda		<- pvSpec$lambda
			omega		<- pvSpec$omega
		}		
	}
		
	K		<- ncol(X1)	# number of outcomes
	n		<- nrow(X1)	# number of pairs

	id00		<- rep(list(NULL), K)
	id10		<- rep(list(NULL), K)
	id01		<- rep(list(NULL), K)
	id11		<- rep(list(NULL), K)
	
	for(k in 1:K){
		id00[[k]]		<- which(X1[,k] == 0 & X2[,k] == 0)
		id10[[k]]		<- which(X1[,k] == 1 & X2[,k] == 0)
		id01[[k]]		<- which(X1[,k] == 0 & X2[,k] == 1)
		id11[[k]]		<- which(X1[,k] == 1 & X2[,k] == 1)
		
		if(length(id00[[k]]) == 0 | length(id10[[k]]) == 0 | length(id01[[k]]) == 0 | length(id11[[k]]) == 0){
			stop('Population average tables must have non-zero counts')
		}
	}

	a1k		<- matrix(0, nrow = burn+B, ncol = K)
	a2k		<- matrix(0, nrow = burn+B, ncol = K)
	
	theta	<- array(0, dim = c(n, K, burn+B))
	beta	<- matrix(0, nrow = burn+B, ncol = K)

	a1k[1,]		<- rep(0, K)
	a2k[1,]		<- rep(0, K)
	theta[,,1]	<- matrix(rep(runif(n, min = -0.1, max = 0.1)), nrow = n, ncol = K)
	beta[1,]	<- rep(0.1, K)
	
	ar1			<- matrix(0, nrow = burn+B, ncol = K)
	ar2			<- matrix(0, nrow = burn+B, ncol = K)
	art			<- array(0, dim = c(4, K, burn+B))
	arb			<- matrix(0, nrow = burn+B, ncol = K)

	if(priorVar == 'VB'){
		# Sample \sigma_k^2 #
		sigbk	<- sigma2VB(n = n, K = K, th = theta[,,1], a_lambda, a_gamma, b_gamma, B_gamma = 0.01, tol, maxIter)
		sigk		<- matrix(rep(sigbk$est, each = burn + B), nrow = burn + B, ncol = K)

		# Sample \tau^2 #
		tau2b		<- tau2VB(n = n, K = K, be = beta[1,], a_omega, a_eta, b_eta, B_eta = 0.01, tol, maxIter)
		tau2		<- rep(tau2b$est, burn + B)
	} else if(priorVar == 'HC') {
		sigbk		<- sigma2VB(n = n, K = K, th = theta[,,1], a_lambda = 0.001, a_gamma = 0.001, b_gamma = 0.001, B_gamma = 0.01, tol = 1e-10, maxIter = 10)
		sVB			<- sigbk$est
		
		tau2b		<- tau2VB(n = n, K = K, be = beta[1,], a_omega = 0.001, a_eta = 0.001, b_eta = 0.001, B_eta = 0.01, tol = 1e-5, maxIter = 10)
		tVB			<- tau2b$est

		sigk		<- matrix(1/lambda, nrow = burn+B, ncol = K)
		tau2		<- rep(1/omega, burn+B)
	} else {
		sigk		<- matrix(1/lambda, nrow = burn+B, ncol = K)
		tau2		<- rep(1/omega, burn+B)
	}

	arTau		<- vector('numeric', length = burn+B)
	arSig		<- matrix(0, nrow = burn+B, ncol = K)
	
	# sampler #
	for(b in 2:(burn+B)){
		
		if(priorVar == 'HC'){
			# Sample \tau^2 #
			tau2b	<- tau2UpdateHC(1, tau = tau2[b-1], tMean = tVB, tTune = tTune, n = n, K = K, be = beta[b-1,], nuTau = nuTau)
			tau2[b]		<- tau2b$est
			arTau[b]		<- tau2b$ar			
			
			# Sample \sigma_k^2 #
			# sigbk	<- sigma2UpdateHC(1, sig = sigk[b-1,1], sMean = sVB, sTune = sTune, n = n, K = K, th = theta[,,b-1], nuSig = nuSig)
			# sigk[b,]	<- rep(sigbk$est, K)
			# arSig[b,]	<- rep(sigbk$ar, K)	
		} else if(priorVar == 'GA'){
			# Sample \tau^2 #
			omegab		<- tau2UpdateGA(1, n, K, beta[b-1,], at, bt)
			tau2[b]		<- 1/omegab

			# Sample \sigma_k^2 #
			lambdab	<- apply(theta[,,b-1], 2, sigma2UpdateGA, samp.size = 1, n = n, K = 1, as = as, bs = bs)
			sigk[b,]	<- 1/lambdab
		}
				
		# loop over outcomes for  #
		for(k in 1:K){
			if(priorVar == 'HC'){
				# Sample \sigma_k^2 #
				sigbk	<- sigma2UpdateHC(1, sig = sigk[b-1,k], sMean = sVB[k], sTune = sTune[k], n = n, K = 1, th = theta[,k,b-1], nuSig = nuSig)
				sigk[b,k]	<- sigbk$est
				arSig[b,k]	<- sigbk$ar
			}
			
			# sample \theta_ik #
			id00k	<- id00[[k]]
			id10k	<- id10[[k]]
			id01k	<- id01[[k]]
			id11k	<- id11[[k]]

			n00k		<- length(id00k)
			n10k		<- length(id10k)
			n01k		<- length(id01k)
			n11k		<- length(id11k)
			
			yfix00	<- yfTBSelect(tInt = searchIntT, inc = incT, cent = centT, x1 = 0, x2 = 0, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link)			
			yfix10	<- yfTBSelect(tInt = searchIntT, inc = incT, cent = centT, x1 = 1, x2 = 0, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link)
			yfix01	<- yfTBSelect(tInt = searchIntT, inc = incT, cent = centT, x1 = 0, x2 = 1, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link)
			yfix11	<- yfTBSelect(tInt = searchIntT, inc = incT, cent = centT, x1 = 1, x2 = 1, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link)

			thUp00k	<- tbUpdate(n00k, yfix00$yfix, x1 = 0, x2 = 0, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, ymin, ymax, upFix, dlpThresh)
			thUp10k	<- tbUpdate(n10k, yfix10$yfix, x1 = 1, x2 = 0, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, ymin, ymax, upFix, dlpThresh)
			thUp01k	<- tbUpdate(n01k, yfix01$yfix, x1 = 0, x2 = 1, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, ymin, ymax, upFix, dlpThresh)
			thUp11k	<- tbUpdate(n11k, yfix11$yfix, x1 = 1, x2 = 1, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, ymin, ymax, upFix, dlpThresh)

			theta[id00k,k,b]	<- thUp00k$thetab
			theta[id10k,k,b]	<- thUp10k$thetab
			theta[id01k,k,b]	<- thUp01k$thetab
			theta[id11k,k,b]	<- thUp11k$thetab

			art[1,k,b]			<- n00k/thUp00k$count
			art[2,k,b]			<- n10k/thUp10k$count
			art[3,k,b]			<- n01k/thUp01k$count
			art[4,k,b]			<- n11k/thUp11k$count			

			# \alpha_1k
			yfixa1k		<- yfAlphaSelect(aInt = searchIntA, inc = incA, cent = centA, theta[,k,b-1], beta[b-1,k], X1[,k], link)
			yfika1		<- yfixa1k$yfix
			a1Upk		<- alphaUpdate(1, yfika1, theta[,k,b-1], beta[b-1,k], X1[,k], link, ymin, ymax, upFix, dlpThresh)
			a1k[b,k]		<- a1Upk$alphab
			ar1[b,k]		<- a1Upk$count
			
			# \alpha_2k
			yfixa2k		<- yfAlphaSelect(aInt = searchIntA, inc = incA, cent = centA, theta[,k,b-1], beta[b-1,k], X2[,k], link)
			yfika2		<- yfixa2k$yfix
			a2Upk		<- alphaUpdate(1, yfika2, theta[,k,b-1], beta[b-1,k], X2[,k], link, ymin, ymax, upFix, dlpThresh)
			a2k[b,k]		<- a2Upk$alphab
			ar2[b,k]		<- a2Upk$count
			
			# Sample \beta_k #
			yfixbk		<- yfTBSelect(tInt = searchIntB, inc = incB, cent = centB, X1[,k], X2[,k], theta[,k,b-1], a1k[b-1,k], a2k[b-1,k], tau2[b-1], link)
			yfbk		<- yfixbk$yfix
			bkUp		<- tbUpdate(1, yfbk, X1[,k], X2[,k], theta[,k,b-1], a1k[b-1,k], a2k[b-1,k], tau2[b-1], link, ymin, ymax, upFix, dlpThresh)
			beta[b,k]	<- bkUp$thetab
			arb[b,k]	<- bkUp$count
		}
		
		if(priorVar == 'VB'){
			indv		<- vector(length = K)
			for(k in 1:K){
				indv[k]	<- sum(theta[,k,b]^2) > 0.08 & sum(theta[,k,b]^2) < 0.3
			}
			if(sum(indv) == K){
				# Sample \sigma_k^2 #
				sigbk					<- sigma2VB(n = n, K = K, th = theta[,,b], a_lambda, a_gamma, b_gamma, B_gamma = 0.01, tol, maxIter)
				sigk[b:(B+burn),]	<- sigbk$est
			}
			
			if(sum(beta[b,]^2) > 0.02 & sum(beta[b,]^2) < 0.2){
				# Sample \tau^2 #
				tau2b					<- tau2VB(n = n, K = K, be = beta[b,], a_omega, a_eta, b_eta, B_eta = 0.01, tol = 1e-5, maxIter = 10)
				tau2[b:(B+burn)]	<- rep(tau2b$est)
			}
		}


		if(b == burn){
			cat("burn-in complete, begin sampling\n")
		}
		
		if(b > burn){
			if(mod(b, up) == 0){
				cat(paste(b - burn,"samples completed\n"))
			}
		}
	}

	# calculate difference in treatment effects
	rho		<- a1k - a2k

	# name output matrices
	colnames(rho)	<- paste("rho", 1:K, sep = "")
	colnames(beta)	<- paste("beta", 1:K, sep = "")
	colnames(sigk)	<- paste("sigma", 1:K, sep = "")
	colnames(a1k)	<- colnames(X1)
	colnames(a2k)	<- colnames(X2)
	colnames(ar1)	<- colnames(X1)
	colnames(ar2)	<- colnames(X2)
	
	alphas	<- list(alpha1 = a1k, alpha2 = a2k)
	aRate	<- list(arTheta = art, arBeta = arb, arAlpha1 = ar1, arAlpha2 = ar2, arTau = arTau, arSigma = arSig)
	if(priorVar == 'VB'){
		varComp	<- list(sigk2 = sigbk$est, E_lambda = sigbk$E_lambda, B_lambda = sigbk$B_lambda, tau2 = tau2b$est, E_omega = tau2b$E_omega, B_omega = tau2b$B_omega, fsig = sigk, ftau = tau2)
	} else if(priorVar == 'FX'){
		varComp	<- list(sigk2 = rep(1/lambda, K), tau2 = 1/omega)
	} else {
		varComp	<- list(sigk2 = sigk, tau2 = tau2)
	}
	MCMCsp	<- list(B = B, burn = burn, priorVar = priorVar, pvSpec = pvSpec)		
	data	<- list(X1 = X1, X2 = X2)
	
	l		<- list(rho = rho, alphas = alphas, vars = varComp, theta = theta, beta = beta, accept = aRate, MCMCsp = MCMCsp, link = link, data	= data)
	
	class(l)	<- 'bmmaps'

	return(l)
	
}

## class related function for bmmaps objects ##

print.bmmaps		<- function(mod){
	K		<- ncol(mod$rho)
	n		<- dim(mod$theta)[1]
	B		<- mod$MCMCsp$B
	burn	<- mod$MCMCsp$burn
	cat("Bayesian Multivariate Matched Proportions\n")
	cat(paste("K =", K, "outcomes for n =", n, "pairs,", mod$link), "link\n")
	
	# rho results
	meds			<- matrix(apply(mod$rho[-c(1:burn),], 2, median), ncol = K, nrow = 1)
	colnames(meds)	<- 1:K
	rownames(meds)	<- ""
	cat("\nMedian Difference in Treatment Effects:\n")
	print(round(meds, digits = 3))
	
}

summary.bmmaps	<- function(mod, alpha = 0.05, rho0 = 0, chains = 4){
	K			<- ncol(mod$rho)
	n			<- dim(mod$theta)[1]
	B			<- mod$MCMCsp$B
	burn		<- mod$MCMCsp$burn
	priorVar	<- mod$MCMCsp$priorVar
	cat("Bayesian Multivariate Matched Proportions\n")
	cat(paste("K =", K, "outcomes for n =", n, "pairs,", mod$link), "link\n")

	# rho summaries
	low				<- alpha/2
	high			<- 1-(alpha/2)
	tab				<- t(apply(mod$rho[-c(1:burn),], 2, quantile, probs = c(0.5, low, high)))
	pprob			<- apply(mod$rho[-c(1:burn),] > rho0, 2, mean)
	GR				<- getGelman(mod$rho[-c(1:burn),], chains)
	tab1				<- cbind(tab, pprob, GR)	
	colnames(tab1)	<- c('Median', paste(100*low,"%", sep = ""), paste(100*high,"%", sep = ""), paste("P(>",rho0,")", sep = ""), 'GR Est.', 'Upper GR')
	rownames(tab1)	<- 1:K
	cat("\nDifference in Treatment Effects\n\n")
	print(round(tab1, 3))
	cat("Note: GR denotes Gelman-Rubin, < 1.1 suggests convergence\n")
	
	# MCMC specs
	acr1	<- B/apply(mod$accept$arAlpha1[-c(1:burn),], 2, sum)
	acr2	<- B/apply(mod$accept$arAlpha2[-c(1:burn),], 2, sum)
	minAr	<- min(c(acr1, acr2))
	maxAr	<- max(c(acr1, acr2))
	arThe	<- apply(mod$accept$arTheta[,,-c(1:burn)], c(1,2), mean)
	minArt	<- min(arThe)
	maxArt	<- max(arThe)
	arBet	<- B/apply(mod$accept$arBeta[-c(1:burn),], 2, sum)
	minBrt	<- min(arBet)
	maxBrt	<- max(arBet)

	cat(paste("\nPosterior samples based on", B, "samples after a burn-in of\n"))
	cat(paste(burn, " (discarded). ARS acceptance rates range from ", round(100*minAr,2), "%\n", sep = ""))
	cat(paste("to ", round(100*maxAr,2), "% (alphas), from ", round(100*minArt,2), "% to ", round(100*maxArt,2), "% (thetas), and\n", sep = ""))
	cat(paste("from ", round(100*minBrt,2), "% to ", round(100*maxBrt, 2), "% (betas).\n", sep = ""))
	
	# transformed rho
	if(mod$link == 'logit'){
		tab2		<- exp(tab)
		colnames(tab2)	<- c('Odds Ratio', paste(100*low,"%", sep = ""), paste(100*high,"%", sep = ""))
		rownames(tab2)	<- 1:K
		cat("\nExponentiated Difference in Treatment Effects\n\n")
		print(round(tab2, 3))
	}
	
	# variance components, if sampled
	if(priorVar == 'HC' | priorVar == 'GA'){
		sigTab		<- apply(mod$vars$sigk2[-c(1:burn),], 2, quantile, probs = c(0.5, low, high))
		tauTab		<- quantile(mod$vars$tau2[-c(1:burn)], probs = c(0.5, low, high))
		tab			<- rbind(t(sigTab), tauTab)
		GRs			<- getGelman(mod$vars$sigk2[-c(1:burn),], chains)
		GRt			<- getGelman(matrix(mod$vars$tau2[-c(1:burn)], ncol = 1), chains)
		GR			<- rbind(GRs, GRt)
		tab3		<- cbind(tab, GR)	
		colnames(tab3)	<- c('Median', paste(100*low,"%", sep = ""), paste(100*high,"%", sep = ""),  'GR Est.', 'Upper GR')
		# if(priorVar == 'HC'){
			# tab3	<- tab3[-2,]
			# rownames(tab3)	<- c('Subject', 'Outcome')
		# } else {
			rownames(tab3)	<- c(paste('Subject', 1:K), 'Outcome')
		# }
		switch(priorVar,
			HC = cat("\nVariance Components with Half-Cauchy Prior\n\n"),
			GA = cat("\nVariance Components with Gamma Prior\n\n")
		)
		print(round(tab3,4))
		if(priorVar == 'HC'){	
			arSig		<- apply(mod$accept$arSigma[-c(1:burn),], 2, mean)
			arTau		<- mean(mod$accept$arTau[-c(1:burn)])
			cat(paste("\nM-H acceptance rates are ", round(100*arTau,2), "% (tau) and ", paste(round(100*apply(mod$accept$arSigma[-c(1:burn),], 2, mean), 2), sep = "", collapse = "%, "), "% (sigmas).\n", sep = ""))
		}
	} else if(priorVar == 'VB'){
		sigk	<- mod$vars$sigk2
		E_lam	<- mod$vars$E_lambda
		B_lam	<- mod$vars$B_lambda
		a_lam	<- mod$MCMCsp$pvSpec$a_lambda
		vb1		<- cbind(sigk, E_lam, qgamma(low, n/2 + a_lam, B_lam), qgamma(high, n/2 + a_lam, B_lam))

		tau2	<- mod$vars$tau2
		E_ome	<- mod$vars$E_omega
		B_ome	<- mod$vars$B_omega
		a_ome	<- mod$MCMCsp$pvSpec$a_omega
		vb2		<- c(tau2, E_ome, qgamma(low, (n*K)/2 + a_ome, B_ome), qgamma(high, (n*K)/2 + a_ome, B_ome))
		
		vb				<- rbind(vb1, vb2)
		colnames(vb)	<- c('Estimate', 'Tuning Parameter', paste(100*low,"%", sep = ""), paste(100*high,"%", sep = ""))
		rownames(vb)	<- c(paste('Subject', 1:K), 'Outcome')
		
		cat("\nVariational Bayes Estimates of Variance Components\n\n")
		print(round(vb,4))
		cat("Note: Tuning parmeter intervals based on product density\n")
		cat("variational approximation")
		
	} else {
		lambda	<- mod$MCMCsp$pvSpec$lambda
		omega	<- mod$MCMCsp$pvSpec$omega
		cat(paste("\nTuning parameters fixed at ", 1/lambda, "(sigmas) and ", 1/omega, "(tau)"))
	}
}

coef.bmmaps		<- function(mod){
	burn		<- mod$MCMCsp$burn
	rhoMed		<- apply(mod$rho[-c(1:burn),], 2, median)
	
	return(rhoMed)
}

sd.bmmaps		<- function(mod){
	burn		<- mod$MCMCsp$burn
	rhoSd		<- apply(mod$rho[-c(1:burn),], 2, sd)
	
	return(rhoSd)
}

ci		<- function(mod, alpha = 0.05){
	burn		<- mod$MCMCsp$burn
	low			<- alpha/2
	high		<- 1-(alpha/2)
	tab			<- t(apply(mod$rho[-c(1:burn),], 2, quantile, probs = c(low, high)))
	
	return(tab)
	
}

pprob	<- function(mod, rho0 = 0){
	burn	<- mod$MCMCsp$burn
	probs	<- apply(mod$rho[-c(1:burn),] > rho0, 2, mean)
	
	return(probs)
}

getVar	<- function(mod){
	burn			<- mod$MCMCsp$burn
	priorVar		<- mod$MCMCsp$priorVar
	if(priorVar == 'HC' | priorVar == 'GA'){
		sigk		<- apply(mod$vars$sigk2[-c(1:burn),], 2, median)
		tau2		<- median(mod$vars$tau2[-c(1:burn)])
	} else if(priorVar == 'VB'){
		sigk		<- mod$vars$sigk2
		tau2		<- mod$vars$tau2
	}
	
	l	<- list(sigk = sigk, tau2 = tau2)
	
	return(l)
}








