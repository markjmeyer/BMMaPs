############# Bayesian Multivariate Matched Proportions #############
# R functions to run Bayesian Multivariate Matched Proportions		
#	from Meyer and Knutson 2019 									
#																	
# Created:  07/25/2018												
# Modified: 07/15/2019												
#																	
# By: Mark J Meyer													
#																	
#####################################################################

## require libraries ##
suppressPackageStartupMessages(require(MASS))
suppressPackageStartupMessages(require(MCMCpack))
suppressPackageStartupMessages(require(mvtnorm))
suppressPackageStartupMessages(require(numbers))
suppressPackageStartupMessages(require(msm))
suppressPackageStartupMessages(require(geepack))

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

lptb	<- function(c1, x1, x2, c2, a1, a2, s, link, mutb){
	in1		<- x1*log(f(c1 + c2 + a1, link)) + (1-x1)*log(1 - f(c1 + c2 + a1, link))
	in2		<- x2*log(f(c1 + c2 + a2, link)) + (1-x2)*log(1 - f(c1 + c2 + a2, link))
	out		<- sum(in1 + in2) - (1/(2*s))*((c1 - mutb)^2)
	return(out)
}

dlptb	<- function(c1, x1, x2, c2, a1, a2, s, link, mutb){
	in1		<- (x1/f(c1 + c2 + a1, link))*dF(c1 + c2 + a1, link) + ((1-x1)/(1-f(c1 + c2 + a1, link)))*dFb(c1 + c2 + a1, link)
	in2		<- (x2/f(c1 + c2 + a2, link))*dF(c1 + c2 + a2, link) + ((1-x2)/(1-f(c1 + c2 + a2, link)))*dFb(c1 + c2 + a2, link)
	out		<- sum(in1 + in2) - (1/s)*(c1 - mutb)
	return(out)
}

lpalp	<- function(aj, th, be, x, link, mua, siga = 1e500){
	in1		<- x*log(f(th + be + aj, link)) + (1 - x)*log(1 - f(th + be + aj, link))
	out		<- sum(in1) - (1/(2*siga)*(aj - mua)^2)
	return(out)
}

dlpalp	<- function(aj, th, be, x, link, mua, siga = 1e500){
	in1		<- (x/f(th + be + aj, link))*dF(th + be + aj, link) + ((1 - x)/(1 - f(th + be + aj, link)))*dFb(th + be + aj, link)
	out		<- sum(in1) - (1/(siga)*(aj - mua))
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

zftb	<- function(yfixed, x1, x2, c2, a1, a2, s, link, mutb){
	yf0f		<- head(yfixed, n=-1)
	yf1f		<- tail(yfixed, n=-1)
	zfixed	<- vector('numeric', length = length(yfixed)-1)
	for(zl in 1:length(zfixed)){
		yf0			<- yf0f[zl]
		yf1			<- yf1f[zl]
		zfixed[zl] <- yf0 + (lptb(yf0, x1, x2, c2, a1, a2, s, link, mutb) - lptb(yf1, x1, x2, c2, a1, a2, s, link, mutb) + (yf1 - yf0)*dlptb(yf1, x1, x2, c2, a1, a2, s, link, mutb)) / (dlptb(yf1, x1, x2, c2, a1, a2, s, link, mutb) - dlptb(yf0, x1, x2, c2, a1, a2, s, link, mutb))
	}
	return(zfixed)	
}

zfalp	<- function(yfixed, th, be, x, link, mua, siga){
	yf0f		<- head(yfixed, n=-1)
	yf1f		<- tail(yfixed, n=-1)
	zfixed	<- vector('numeric', length = length(yfixed)-1)
	for(zl in 1:length(zfixed)){
		yf0			<- yf0f[zl]
		yf1			<- yf1f[zl]
		zfixed[zl]	<- yf0 + (lpalp(yf0, th, be, x, link, mua, siga) - lpalp(yf1, th, be, x, link, mua, siga) + (yf1 - yf0)*dlpalp(yf1, th, be, x, link, mua, siga)) / (dlpalp(yf1, th, be, x, link, mua, siga) - dlpalp(yf0, th, be, x, link, mua, siga))
	}
	return(zfixed)	
}

utb <- function(y, yfixed, x1, x2, c2, a1, a2, s, link, mutb, ymin = -Inf, ymax = Inf){
	res				<- rep(0, length(y))
	zfixed			<- zftb(yfixed, x1, x2, c2, a1, a2, s, link, mutb)

	piecewise.idx	<- findInterval(y, c(ymin, zfixed, ymax))
	npieces 		<- length(zfixed) + 2
	for(pidx in 1:npieces){
		yp 			<- y[piecewise.idx == pidx]
		xx 			<- lptb(yfixed[pidx], x1, x2, c2, a1, a2, s, link, mutb) + (yp - yfixed[pidx])*dlptb(yfixed[pidx], x1, x2, c2, a1, a2, s, link, mutb)
		res[piecewise.idx == pidx]	<- xx
	}
	return(res)
}

ualp <- function(y, yfixed, th, be, x, link, mua, siga, ymin = -Inf, ymax = Inf){
	res				<- rep(0, length(y))
	zfixed			<- zfalp(yfixed, th, be, x, link, mua, siga)

	piecewise.idx	<- findInterval(y, c(ymin, zfixed, ymax))
	npieces 		<- length(zfixed) + 2
	for(pidx in 1:npieces){
		yp 			<- y[piecewise.idx == pidx]
		xx 			<- lpalp(yfixed[pidx], th, be, x, link, mua, siga) + (yp - yfixed[pidx])*dlpalp(yfixed[pidx], th, be, x, link, mua, siga)
		res[piecewise.idx == pidx]	<- xx
	}
	return(res)
}

stb		<- function(vals, yfixed, x1, x2, c2, a1, a2, s, link, mutb){
	zfixed		<- zftb(yfixed, x1, x2, c2, a1, a2, s, link, mutb)

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
		ds			<- exp(lptb(yp, x1, x2, c2, a1, a2, s, link, mutb))/dlptb(yp, x1, x2, c2, a1, a2, s, link, mutb) * ( exp((zp - yp)*dlptb(yp, x1, x2, c2, a1, a2, s, link, mutb)) - exp((zm - yp)*dlptb(yp, x1, x2, c2, a1, a2, s, link, mutb)) )

		cidx		<- zm < vals & vals <= zp
		hidx		<- vals > zp

		pct[cidx]	<- pct[cidx] + exp(lptb(yp, x1, x2, c2, a1, a2, s, link, mutb))/dlptb(yp, x1, x2, c2, a1, a2, s, link, mutb) * ( exp((vals[cidx] - yp)*dlptb(yp, x1, x2, c2, a1, a2, s, link, mutb)) - exp((zm - yp)*dlptb(yp, x1, x2, c2, a1, a2, s, link, mutb)) )
		pct[hidx]	<- pct[hidx] + ds

		norm.const <- norm.const + ds
	}

	l <- list(pct = pct / norm.const, norm.const = norm.const)
	return(l)
}

salp		<- function(vals, yfixed, th, be, x, link, mua, siga){
	zfixed		<- zfalp(yfixed, th, be, x, link, mua, siga)

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
		ds			<- exp(lpalp(yp, th, be, x, link, mua, siga))/dlpalp(yp, th, be, x, link, mua, siga) * ( exp((zp - yp)*dlpalp(yp, th, be, x, link, mua, siga)) - exp((zm - yp)*dlpalp(yp, th, be, x, link, mua, siga)) )

		cidx		<- zm < vals & vals <= zp
		hidx		<- vals > zp

		pct[cidx]	<- pct[cidx] + exp(lpalp(yp, th, be, x, link, mua, siga))/dlpalp(yp, th, be, x, link, mua, siga) * ( exp((vals[cidx] - yp)*dlpalp(yp, th, be, x, link, mua, siga)) - exp((zm - yp)*dlpalp(yp, th, be, x, link, mua, siga)) )
		pct[hidx]	<- pct[hidx] + ds

		norm.const <- norm.const + ds
	}

	l <- list(pct = pct / norm.const, norm.const = norm.const)
	return(l)
}

stb.sample <- function(samp.size, yfixed, x1, x2, c2, a1, a2, s, link, mutb, ymin){
	zfixed			<- zftb(yfixed, x1, x2, c2, a1, a2, s, link, mutb)
	gp				<- stb(zfixed, yfixed, x1, x2, c2, a1, a2, s, link, mutb)
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
		tmp	<- (ui - ub[ii]) * dlptb(yp, x1, x2, c2, a1, a2, s, link, mutb) * norm.const / exp(lptb(yp, x1, x2, c2, a1, a2, s, link, mutb)) + exp( (zm - yp)*dlptb(yp, x1, x2, c2, a1, a2, s, link, mutb) )
		tmp	<- yp + log(tmp) / dlptb(yp, x1, x2, c2, a1, a2, s, link, mutb)
		res[ fidx == ii ] <- tmp
	}
	return(res)
}

salp.sample <- function(samp.size, yfixed, th, be, x, link, mua, siga, ymin){
	zfixed			<- zfalp(yfixed, th, be, x, link, mua, siga)
	gp				<- salp(zfixed, yfixed, th, be, x, link, mua, siga)
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
		tmp	<- (ui - ub[ii]) * dlpalp(yp, th, be, x, link, mua, siga) * norm.const / exp(lpalp(yp, th, be, x, link, mua, siga)) + exp( (zm - yp)*dlpalp(yp, th, be, x, link, mua, siga) )
		tmp	<- yp + log(tmp) / dlpalp(yp, th, be, x, link, mua, siga)
		res[ fidx == ii ] <- tmp
	}
	return(res)
}


tbUpdate	<- function(draws, yfixed, x1, x2, c2, a1, a2, s, link, mutb, ymin, ymax, upFix = TRUE, dlpThresh){
	thb	<- vector('numeric', length = draws)

	count	<- 0
	tt		<- 0
	while(tt < draws){
		count	<- count + 1

		# draw Ys from sthe #
		# Ys	<- sthe.sample(samp.size, yfixed, x1, x2, a1, a2, s, link, ymin)
		Ys	<- stb.sample(1, yfixed, x1, x2, c2, a1, a2, s, link, mutb, ymin)

		# draw U(0,1) #
		U	<- runif(1)

		test <- U < exp(lptb(Ys, x1, x2, c2, a1, a2, s, link, mutb) - utb(Ys, yfixed, x1, x2, c2, a1, a2, s, link, mutb, ymin, ymax))

		if(test){
			tt			<- tt + 1
			thb[tt]		<- Ys
		} else{
			dlpEval		<- dlptb(Ys, x1, x2, c2, a1, a2, s, link, mutb)
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

alphaUpdate	<- function(samp.size, yfixed, th, be, x, link, mua, siga, ymin, ymax, upFix = TRUE, dlpThresh){
	alb	<- vector('numeric', length = samp.size)

	count	<- 0
	tt		<- 0
	while(tt < samp.size){
		count	<- count + 1

		# draw Ys from sthe #
		Ys		<- salp.sample(samp.size, yfixed, th, be, x, link, mua, siga, ymin)

		# draw U(0,1) #
		U		<- runif(1)

		test	<- U < exp(lpalp(Ys, th, be, x, link, mua, siga) - ualp(Ys, yfixed, th, be, x, link, mua, siga, ymin, ymax))

		if(test){
			tt			<- tt + 1
			alb[tt]		<- Ys
		} else{
			dlpEval		<- dlpalp(Ys, th, be, x, link, mua, siga)
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

yfTBSelect	<- function(tInt = c(-2,2), inc = 0.005, cent = 0.15, x1, x2, c2, a1, a2, s, link, mutb){
	tseq		<- matrix(seq(tInt[1], tInt[2], by = inc), ncol = 1)
	lp			<- apply(tseq, 1, function(z) lptb(z, x1, x2, c2, a1, a2, s, link, mutb))

	cselp	<- cumsum(exp(lp))
	if(min(cselp/sum(exp(lp))) > cent){
		stop('percentile too small (thetas or betas), consider increasing cent or searchInt')
	}
	low		<- tseq[max(which(cselp/sum(exp(lp)) <= cent))]
	mid		<- tseq[max(which(cselp/sum(exp(lp)) <= 0.5))]
	high		<- tseq[max(which(cselp/sum(exp(lp)) <= 1-cent))]
	dlpl		<- dlptb(low, x1, x2, c2, a1, a2, s, link, mutb)
	dlpm		<- dlptb(mid, x1, x2, c2, a1, a2, s, link, mutb)
	dlph		<- dlptb(high, x1, x2, c2, a1, a2, s, link, mutb)
	
	l	<- list(yfix = c(low, mid, high), deriv = c(dlpl, dlpm, dlph))
	
	return(l)

}

yfAlphaSelect	<- function(aInt = c(-10, 10), inc = 0.025, cent = 0.15, th, be, x, link, mua, siga){
	aseq		<- matrix(seq(aInt[1], aInt[2], by = inc), ncol = 1)
	lp		<- apply(aseq, 1, function(z) lpalp(z, th, be, x, link, mua, siga))
	
	cselp	<- cumsum(exp(lp))
	if(min(cselp/sum(exp(lp))) > cent){
		stop('percentile too small (alphas), consider increasing cent or searchInt')
	}
	low		<- aseq[max(which(cselp/sum(exp(lp)) <= cent))]
	mid		<- aseq[max(which(cselp/sum(exp(lp)) <= 0.5))]
	high		<- aseq[max(which(cselp/sum(exp(lp)) <= 1-cent))]
	dlpl		<- dlpalp(low, th, be, x, link, mua, siga)
	dlpm		<- dlpalp(mid, th, be, x, link, mua, siga)
	dlph		<- dlpalp(high, th, be, x, link, mua, siga)

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

## primary function for bmir ##
bmir	<- function(X1, X2, link = 'logit', B, burn = 100, up = 1000, searchIntT = NULL, incT = 0.005, centT = 0.1, searchIntB = NULL, incB = 0.005, centB = 0.1, searchIntA = NULL, incA = 0.01, centA = 0.1, priorVar = 'VB', pvSpec = NULL, ymin = -Inf, ymax = Inf, upFix = TRUE, dlpThresh = 100){

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

	K		<- ncol(X1)	# number of outcomes
	n		<- nrow(X1)	# number of pairs
	
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
			
			mut			<- rep(0, K)
			mub			<- rep(0, K)
			mua1		<- rep(0, K)
			mua2		<- rep(0, K)
			
			siga1		<- rep(10, K)	# shrink more? shrink by sparsity?
			siga2		<- rep(10, K)	# shrink more?
			
			pvSpec	<- list(a_lambda = a_lambda, a_gamma = a_gamma, b_gamma = b_gamma,
							a_omega = a_omega, a_eta = a_eta, b_eta = b_eta,
							tol = tol, maxIter = maxIter, mut = mut, mub = mub,
							mua1 = mua1, mua2 = mua2, siga1 = siga1, siga2 = siga2)
		} else if(priorVar == 'GA'){
			lambda		<- 100
			omega		<- 100
			at			<- 1
			bt			<- 1
			as			<- 1
			bs			<- 1
			
			mut			<- rep(0, K)
			mub			<- rep(0, K)
			mua1		<- rep(0, K)
			mua2		<- rep(0, K)

			siga1		<- rep(10, K)
			siga2		<- rep(10, K)

			pvSpec		<- list(lambda = lambda, omega = omega,
								at = at, bt = bt, as = as, bs = bs,
								mut = mut, mub = mub, mua1 = mua1, mua2 = mua2,
								siga1 = siga1, siga2 = siga2)
		} else if(priorVar == 'HC'){
			lambda		<- 100
			omega		<- 100
			tTune		<- 0.01
			nuTau		<- 2
			sTune		<- rep(0.01, ncol(X1))
			nuSig		<- 2
			
			mut			<- rep(0, K)
			mub			<- rep(0, K)
			mua1		<- rep(0, K)
			mua2		<- rep(0, K)

			siga1		<- rep(10, K)
			siga2		<- rep(10, K)

			pvSpec		<- list(lambda = lambda, omega = omega,
								tTune = tTune, nuTau = nuTau,
								sTune = sTune, nuSig = nuSig,
								mut = mut, mub = mub, mua1 = mua1, mua2 = mua2,
								siga1 = siga1, siga2 = siga2)
		} else {
			lambda		<- 100
			omega		<- 100

			mut			<- rep(0, K)
			mub			<- rep(0, K)
			mua1		<- rep(0, K)
			mua2		<- rep(0, K)

			siga1		<- rep(10, K)
			siga2		<- rep(10, K)
			
			pvSpec		<- list(lambda = lambda, omega = omega, mut = mut, mub = mub,
								mua1 = mua1, mua2 = mua2, siga1 = siga1, siga2 = siga2)
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
			
			if(is.null(pvSpec$mut)){
					mut			<- rep(0, K)
			} else {
				mut			<- pvSpec$mut			
			}
			if(is.null(pvSpec$mub)){
				mub			<- rep(0, K)
			} else {
				mub			<- pvSpec$mub			
			}
			if(is.null(pvSpec$mua1)){
					mua1		<- rep(0, K)
			} else {
				mua1		<- pvSpec$mua1
			}
			if(is.null(pvSpec$mua2)){
					mua2		<- rep(0, K)		
			} else {
				mua2		<- pvSpec$mua2
			}
			if(is.null(pvSpec$siga1)){
					siga1		<- rep(10, K)
			} else {
				siga1		<- pvSpec$siga1
			}
			if(is.null(pvSpec$siga2)){
					siga2		<- rep(10, K)		
			} else {
				siga2		<- pvSpec$siga2
			}
		} else if(priorVar == 'GA'){
			lambda		<- pvSpec$lambda
			omega		<- pvSpec$omega
			at			<- pvSpec$at
			bt			<- pvSpec$bt
			as			<- pvSpec$as
			bs			<- pvSpec$bs

			if(is.null(pvSpec$mut)){
					mut			<- rep(0, K)
			} else {
				mut			<- pvSpec$mut			
			}
			if(is.null(pvSpec$mub)){
				mub			<- rep(0, K)
			} else {
				mub			<- pvSpec$mub			
			}
			if(is.null(pvSpec$mua1)){
					mua1		<- rep(0, K)
			} else {
				mua1		<- pvSpec$mua1
			}
			if(is.null(pvSpec$mua2)){
					mua2		<- rep(0, K)		
			} else {
				mua2		<- pvSpec$mua2
			}
			if(is.null(pvSpec$siga1)){
					siga1		<- rep(10, K)
			} else {
				siga1		<- pvSpec$siga1
			}
			if(is.null(pvSpec$siga2)){
					siga2		<- rep(10, K)		
			} else {
				siga2		<- pvSpec$siga2
			}
		} else if(priorVar == 'HC'){
			lambda		<- pvSpec$lambda
			omega		<- pvSpec$omega
			tTune		<- pvSpec$tTune
			nuTau		<- pvSpec$nuTau
			sTune		<- pvSpec$sTune
			nuSig		<- pvSpec$nuSig

			if(is.null(pvSpec$mut)){
					mut			<- rep(0, K)
			} else {
				mut			<- pvSpec$mut			
			}
			if(is.null(pvSpec$mub)){
				mub			<- rep(0, K)
			} else {
				mub			<- pvSpec$mub			
			}
			if(is.null(pvSpec$mua1)){
					mua1		<- rep(0, K)
			} else {
				mua1		<- pvSpec$mua1
			}
			if(is.null(pvSpec$mua2)){
					mua2		<- rep(0, K)		
			} else {
				mua2		<- pvSpec$mua2
			}
			if(is.null(pvSpec$siga1)){
					siga1		<- rep(10, K)
			} else {
				siga1		<- pvSpec$siga1
			}
			if(is.null(pvSpec$siga2)){
					siga2		<- rep(10, K)		
			} else {
				siga2		<- pvSpec$siga2
			}
		} else {
			lambda		<- pvSpec$lambda
			omega		<- pvSpec$omega

			if(is.null(pvSpec$mut)){
					mut			<- rep(0, K)
			} else {
				mut			<- pvSpec$mut			
			}
			if(is.null(pvSpec$mub)){
				mub			<- rep(0, K)
			} else {
				mub			<- pvSpec$mub			
			}
			if(is.null(pvSpec$mua1)){
					mua1		<- rep(0, K)
			} else {
				mua1		<- pvSpec$mua1
			}
			if(is.null(pvSpec$mua2)){
					mua2		<- rep(0, K)		
			} else {
				mua2		<- pvSpec$mua2
			}
			if(is.null(pvSpec$siga1)){
					siga1		<- rep(10, K)
			} else {
				siga1		<- pvSpec$siga1
			}
			if(is.null(pvSpec$siga2)){
					siga2		<- rep(10, K)		
			} else {
				siga2		<- pvSpec$siga2
			}

		}		
	}
	
	id00		<- rep(list(NULL), K)
	id10		<- rep(list(NULL), K)
	id01		<- rep(list(NULL), K)
	id11		<- rep(list(NULL), K)

	for(k in 1:K){
		id00[[k]]		<- which(X1[,k] == 0 & X2[,k] == 0)
		id10[[k]]		<- which(X1[,k] == 1 & X2[,k] == 0)
		id01[[k]]		<- which(X1[,k] == 0 & X2[,k] == 1)
		id11[[k]]		<- which(X1[,k] == 1 & X2[,k] == 1)
		
		# if(length(id00[[k]]) == 0 | length(id10[[k]]) == 0 | length(id01[[k]]) == 0 | length(id11[[k]]) == 0){
			# stop('Population average tables must have non-zero counts')
		# }
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
	sigma1		<- vector('numeric', length = burn+B)
	sigma2		<- vector('numeric', length = burn+B)
	
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
			
			if(n00k > 0){
				yfix00	<- yfTBSelect(tInt = searchIntT, inc = incT, cent = centT, x1 = 0, x2 = 0, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, mut[k])
				thUp00k	<- tbUpdate(n00k, yfix00$yfix, x1 = 0, x2 = 0, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, mut[k], ymin, ymax, upFix, dlpThresh)
				theta[id00k,k,b]	<- thUp00k$thetab
				art[1,k,b]			<- n00k/thUp00k$count
			}
			
			if(n10k > 0){
				yfix10	<- yfTBSelect(tInt = searchIntT, inc = incT, cent = centT, x1 = 1, x2 = 0, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, mut[k])
				thUp10k	<- tbUpdate(n10k, yfix10$yfix, x1 = 1, x2 = 0, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, mut[k], ymin, ymax, upFix, dlpThresh)
				theta[id10k,k,b]	<- thUp10k$thetab
				art[2,k,b]			<- n10k/thUp10k$count
			}
			
			if(n01k > 0){
				yfix01	<- yfTBSelect(tInt = searchIntT, inc = incT, cent = centT, x1 = 0, x2 = 1, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, mut[k])
				thUp01k	<- tbUpdate(n01k, yfix01$yfix, x1 = 0, x2 = 1, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, mut[k], ymin, ymax, upFix, dlpThresh)
				theta[id01k,k,b]	<- thUp01k$thetab
				art[3,k,b]			<- n01k/thUp01k$count
			}
			
			if(n11k > 0){
				yfix11	<- yfTBSelect(tInt = searchIntT, inc = incT, cent = centT, x1 = 1, x2 = 1, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, mut[k])
				thUp11k	<- tbUpdate(n11k, yfix11$yfix, x1 = 1, x2 = 1, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, mut[k], ymin, ymax, upFix, dlpThresh)
				theta[id11k,k,b]	<- thUp11k$thetab
				art[4,k,b]			<- n11k/thUp11k$count			
			}
			
			# \alpha_1k
			yfixa1k		<- yfAlphaSelect(aInt = searchIntA, inc = incA, cent = centA, theta[,k,b-1], beta[b-1,k], X1[,k], link, mua1[k], siga1[k])
			yfika1		<- yfixa1k$yfix
			a1Upk		<- alphaUpdate(1, yfika1, theta[,k,b-1], beta[b-1,k], X1[,k], link, mua1[k], siga1[k], ymin, ymax, upFix, dlpThresh)
			a1k[b,k]		<- a1Upk$alphab
			ar1[b,k]		<- a1Upk$count
			
			# \alpha_2k
			yfixa2k		<- yfAlphaSelect(aInt = searchIntA, inc = incA, cent = centA, theta[,k,b-1], beta[b-1,k], X2[,k], link, mua2[k], siga2[k])
			yfika2		<- yfixa2k$yfix
			a2Upk		<- alphaUpdate(1, yfika2, theta[,k,b-1], beta[b-1,k], X2[,k], link, mua2[k], siga2[k], ymin, ymax, upFix, dlpThresh)
			a2k[b,k]		<- a2Upk$alphab
			ar2[b,k]		<- a2Upk$count
			
			# Sample \beta_k #
			yfixbk		<- yfTBSelect(tInt = searchIntB, inc = incB, cent = centB, X1[,k], X2[,k], theta[,k,b-1], a1k[b-1,k], a2k[b-1,k], tau2[b-1], link, mub[k])
			yfbk		<- yfixbk$yfix
			bkUp		<- tbUpdate(1, yfbk, X1[,k], X2[,k], theta[,k,b-1], a1k[b-1,k], a2k[b-1,k], tau2[b-1], link, mub[k], ymin, ymax, upFix, dlpThresh)
			beta[b,k]	<- bkUp$thetab
			arb[b,k]	<- bkUp$count
		}

		# sample \sigma_a's #
		# if(link == 'logit'){
			sigma1[b]	<- rinvgamma(1, (K/2) + 1, (1/2)*t(a1k[b-1,] - mua1)%*%(a1k[b-1,] - mua1) + 1)
			sigma2[b]	<- rinvgamma(1, (K/2) + 1, (1/2)*t(a2k[b-1,] - mua2)%*%(a2k[b-1,] - mua2) + 1)
			siga1	<- rep(sigma1[b], K)
			siga2	<- rep(sigma2[b], K)
		# }
		
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
				tau2b					<- tau2VB(n = 1, K = K, be = beta[b,], a_omega, a_eta, b_eta, B_eta = 0.01, tol = 1e-5, maxIter = 10)
				# tau2b					<- tau2VB(n = n, K = K, be = beta[b,], a_omega, a_eta, b_eta, B_eta = 0.01, tol = 1e-5, maxIter = 10)
				tau2[b:(B+burn)]	<- rep(tau2b$est)
			}
		}

		if(mod(b, 10) == 0){
			cat('.')
		}
		
		if(mod(b, up) == 0){
			cat(paste("\n",b,"samples completed\n"))
		}

		# if(b == burn){
			# cat("\nburn-in complete, begin sampling\n")
		# }

		# if(b > burn){
			# if(mod(b, up) == 0){
				# cat(paste("\n",b - burn,"samples completed\n"))
			# }
		# }
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
	varComp$tune	<- list(sigma1 = sigma1, sigma2 = sigma2)
	MCMCsp	<- list(B = B, burn = burn, priorVar = priorVar, pvSpec = pvSpec)		
	data	<- list(X1 = X1, X2 = X2)
	
	l		<- list(rho = rho, alphas = alphas, vars = varComp, theta = theta, beta = beta, accept = aRate, MCMCsp = MCMCsp, link = link, data	= data)
	
	class(l)	<- 'bmir'

	return(l)
	
}

# bmir_old	<- function(X1, X2, link = 'logit', B, burn = 100, up = 1000, searchIntT = NULL, incT = 0.005, centT = 0.1, searchIntB = NULL, incB = 0.005, centB = 0.1, searchIntA = NULL, incA = 0.01, centA = 0.1, priorVar = 'VB', pvSpec = NULL, ymin = -Inf, ymax = Inf, upFix = TRUE, dlpThresh = 100){

	# if(ncol(X1) != ncol(X2)){
		# stop('column dimensions of X1 and X2 must agree')
	# }

	# if(nrow(X1) != nrow(X2)){
		# stop('row dimensions of X1 and X2 must agree')
	# }

	# if(link != 'logit' & link != 'probit' & link != 'cloglog'){
		# stop('improper link specified use either logit, probit, or cloglog')
	# }

	# if(priorVar != 'VB' & priorVar != 'FX' & priorVar != 'GA' & priorVar != 'HC'){
		# stop('improper prior variance specified use either VB, HC, GA, or FX')
	# }
		
	# if(B < 1){
		# stop('B must be at least 1')
	# }
	
	# if(is.null(searchIntA)){
		# if(link == 'probit'){
			# searchIntA	<- c(-7,7)
		# } else if(link == 'cloglog'){
			# searchIntA	<- c(-13, 3)
		# } else {
			# searchIntA	<- c(-10, 10)
		# }
	# }
	
	# if(is.null(searchIntB)){
		# if(link == 'cloglog'){
			# searchIntB	<- c(-5, 3)
		# } else {
			# searchIntB	<- c(-5, 5)
		# }
	# }

	# if(is.null(searchIntT)){
		# if(link == 'cloglog'){
			# searchIntT	<- c(-5, 3)
		# } else {
			# searchIntT	<- c(-5, 5)
		# }
	# }

	# K		<- ncol(X1)	# number of outcomes
	# n		<- nrow(X1)	# number of pairs
	
	# if(is.null(pvSpec)){
		# if(priorVar == 'VB'){
			# a_lambda	<- 0.001
			# a_gamma		<- 0.001
			# b_gamma		<- 0.001
			# a_omega		<- 0.001
			# a_eta		<- 0.001
			# b_eta		<- 0.001
			# tol			<- 1e-5
			# maxIter		<- 10
			
			# mut			<- rep(0, K)
			# mub			<- rep(0, K)
			# mua1		<- rep(0, K)
			# mua2		<- rep(0, K)
			
			# siga1		<- rep(10, K)
			# siga2		<- rep(10, K)
			
			# pvSpec	<- list(a_lambda = a_lambda, a_gamma = a_gamma, b_gamma = b_gamma,
							# a_omega = a_omega, a_eta = a_eta, b_eta = b_eta,
							# tol = tol, maxIter = maxIter, mut = mut, mub = mub,
							# mua1 = mua1, mua2 = mua2, siga1 = siga1, siga2 = siga2)
		# } else if(priorVar == 'GA'){
			# lambda		<- 100
			# omega		<- 100
			# at			<- 1
			# bt			<- 1
			# as			<- 1
			# bs			<- 1
			
			# mut			<- rep(0, K)
			# mub			<- rep(0, K)
			# mua1		<- rep(0, K)
			# mua2		<- rep(0, K)

			# siga1		<- rep(10, K)
			# siga2		<- rep(10, K)

			# pvSpec		<- list(lambda = lambda, omega = omega,
								# at = at, bt = bt, as = as, bs = bs,
								# mut = mut, mub = mub, mua1 = mua1, mua2 = mua2,
								# siga1 = siga1, siga2 = siga2)
		# } else if(priorVar == 'HC'){
			# lambda		<- 100
			# omega		<- 100
			# tTune		<- 0.01
			# nuTau		<- 2
			# sTune		<- rep(0.01, ncol(X1))
			# nuSig		<- 2
			
			# mut			<- rep(0, K)
			# mub			<- rep(0, K)
			# mua1		<- rep(0, K)
			# mua2		<- rep(0, K)

			# siga1		<- rep(10, K)
			# siga2		<- rep(10, K)

			# pvSpec		<- list(lambda = lambda, omega = omega,
								# tTune = tTune, nuTau = nuTau,
								# sTune = sTune, nuSig = nuSig,
								# mut = mut, mub = mub, mua1 = mua1, mua2 = mua2,
								# siga1 = siga1, siga2 = siga2)
		# } else {
			# lambda		<- 100
			# omega		<- 100

			# mut			<- rep(0, K)
			# mub			<- rep(0, K)
			# mua1		<- rep(0, K)
			# mua2		<- rep(0, K)

			# siga1		<- rep(10, K)
			# siga2		<- rep(10, K)
			
			# pvSpec		<- list(lambda = lambda, omega = omega, mut = mut, mub = mub,
								# mua1 = mua1, mua2 = mua2, siga1 = siga1, siga2 = siga2)
		# }
	# } else {
		# if(priorVar == 'VB'){
			# a_lambda	<- pvSpec$a_lambda
			# a_gamma		<- pvSpec$a_gamma
			# b_gamma		<- pvSpec$b_gamma
			# a_omega		<- pvSpec$a_omega
			# a_eta		<- pvSpec$a_eta
			# b_eta		<- pvSpec$b_eta
			# tol			<- pvSpec$tol
			# maxIter		<- pvSpec$maxIter
			
			# if(is.null(pvSpec$mut)){
					# mut			<- rep(0, K)
			# } else {
				# mut			<- pvSpec$mut			
			# }
			# if(is.null(pvSpec$mub)){
				# mub			<- rep(0, K)
			# } else {
				# mub			<- pvSpec$mub			
			# }
			# if(is.null(pvSpec$mua1)){
					# mua1		<- rep(0, K)
			# } else {
				# mua1		<- pvSpec$mua1
			# }
			# if(is.null(pvSpec$mua2)){
					# mua2		<- rep(0, K)		
			# } else {
				# mua2		<- pvSpec$mua2
			# }
			# if(is.null(pvSpec$siga1)){
					# siga1		<- rep(10, K)
			# } else {
				# siga1		<- pvSpec$siga1
			# }
			# if(is.null(pvSpec$siga2)){
					# siga2		<- rep(10, K)		
			# } else {
				# siga2		<- pvSpec$siga2
			# }
		# } else if(priorVar == 'GA'){
			# lambda		<- pvSpec$lambda
			# omega		<- pvSpec$omega
			# at			<- pvSpec$at
			# bt			<- pvSpec$bt
			# as			<- pvSpec$as
			# bs			<- pvSpec$bs

			# if(is.null(pvSpec$mut)){
					# mut			<- rep(0, K)
			# } else {
				# mut			<- pvSpec$mut			
			# }
			# if(is.null(pvSpec$mub)){
				# mub			<- rep(0, K)
			# } else {
				# mub			<- pvSpec$mub			
			# }
			# if(is.null(pvSpec$mua1)){
					# mua1		<- rep(0, K)
			# } else {
				# mua1		<- pvSpec$mua1
			# }
			# if(is.null(pvSpec$mua2)){
					# mua2		<- rep(0, K)		
			# } else {
				# mua2		<- pvSpec$mua2
			# }
			# if(is.null(pvSpec$siga1)){
					# siga1		<- rep(10, K)
			# } else {
				# siga1		<- pvSpec$siga1
			# }
			# if(is.null(pvSpec$siga2)){
					# siga2		<- rep(10, K)		
			# } else {
				# siga2		<- pvSpec$siga2
			# }
		# } else if(priorVar == 'HC'){
			# lambda		<- pvSpec$lambda
			# omega		<- pvSpec$omega
			# tTune		<- pvSpec$tTune
			# nuTau		<- pvSpec$nuTau
			# sTune		<- pvSpec$sTune
			# nuSig		<- pvSpec$nuSig

			# if(is.null(pvSpec$mut)){
					# mut			<- rep(0, K)
			# } else {
				# mut			<- pvSpec$mut			
			# }
			# if(is.null(pvSpec$mub)){
				# mub			<- rep(0, K)
			# } else {
				# mub			<- pvSpec$mub			
			# }
			# if(is.null(pvSpec$mua1)){
					# mua1		<- rep(0, K)
			# } else {
				# mua1		<- pvSpec$mua1
			# }
			# if(is.null(pvSpec$mua2)){
					# mua2		<- rep(0, K)		
			# } else {
				# mua2		<- pvSpec$mua2
			# }
			# if(is.null(pvSpec$siga1)){
					# siga1		<- rep(10, K)
			# } else {
				# siga1		<- pvSpec$siga1
			# }
			# if(is.null(pvSpec$siga2)){
					# siga2		<- rep(10, K)		
			# } else {
				# siga2		<- pvSpec$siga2
			# }
		# } else {
			# lambda		<- pvSpec$lambda
			# omega		<- pvSpec$omega

			# if(is.null(pvSpec$mut)){
					# mut			<- rep(0, K)
			# } else {
				# mut			<- pvSpec$mut			
			# }
			# if(is.null(pvSpec$mub)){
				# mub			<- rep(0, K)
			# } else {
				# mub			<- pvSpec$mub			
			# }
			# if(is.null(pvSpec$mua1)){
					# mua1		<- rep(0, K)
			# } else {
				# mua1		<- pvSpec$mua1
			# }
			# if(is.null(pvSpec$mua2)){
					# mua2		<- rep(0, K)		
			# } else {
				# mua2		<- pvSpec$mua2
			# }
			# if(is.null(pvSpec$siga1)){
					# siga1		<- rep(10, K)
			# } else {
				# siga1		<- pvSpec$siga1
			# }
			# if(is.null(pvSpec$siga2)){
					# siga2		<- rep(10, K)		
			# } else {
				# siga2		<- pvSpec$siga2
			# }

		# }		
	# }
	
	# id00		<- rep(list(NULL), K)
	# id10		<- rep(list(NULL), K)
	# id01		<- rep(list(NULL), K)
	# id11		<- rep(list(NULL), K)

	# for(k in 1:K){
		# id00[[k]]		<- which(X1[,k] == 0 & X2[,k] == 0)
		# id10[[k]]		<- which(X1[,k] == 1 & X2[,k] == 0)
		# id01[[k]]		<- which(X1[,k] == 0 & X2[,k] == 1)
		# id11[[k]]		<- which(X1[,k] == 1 & X2[,k] == 1)
		
		# if(length(id00[[k]]) == 0 | length(id10[[k]]) == 0 | length(id01[[k]]) == 0 | length(id11[[k]]) == 0){
			# stop('Population average tables must have non-zero counts')
		# }
	# }

	# a1k		<- matrix(0, nrow = burn+B, ncol = K)
	# a2k		<- matrix(0, nrow = burn+B, ncol = K)
	
	# theta	<- array(0, dim = c(n, K, burn+B))
	# beta	<- matrix(0, nrow = burn+B, ncol = K)

	# a1k[1,]		<- rep(0, K)
	# a2k[1,]		<- rep(0, K)
	# theta[,,1]	<- matrix(rep(runif(n, min = -0.1, max = 0.1)), nrow = n, ncol = K)
	# beta[1,]	<- rep(0.1, K)
	
	# ar1			<- matrix(0, nrow = burn+B, ncol = K)
	# ar2			<- matrix(0, nrow = burn+B, ncol = K)
	# art			<- array(0, dim = c(4, K, burn+B))
	# arb			<- matrix(0, nrow = burn+B, ncol = K)

	# if(priorVar == 'VB'){
		# # Sample \sigma_k^2 #
		# sigbk	<- sigma2VB(n = n, K = K, th = theta[,,1], a_lambda, a_gamma, b_gamma, B_gamma = 0.01, tol, maxIter)
		# sigk		<- matrix(rep(sigbk$est, each = burn + B), nrow = burn + B, ncol = K)

		# # Sample \tau^2 #
		# tau2b		<- tau2VB(n = n, K = K, be = beta[1,], a_omega, a_eta, b_eta, B_eta = 0.01, tol, maxIter)
		# tau2		<- rep(tau2b$est, burn + B)
	# } else if(priorVar == 'HC') {
		# sigbk		<- sigma2VB(n = n, K = K, th = theta[,,1], a_lambda = 0.001, a_gamma = 0.001, b_gamma = 0.001, B_gamma = 0.01, tol = 1e-10, maxIter = 10)
		# sVB			<- sigbk$est
		
		# tau2b		<- tau2VB(n = n, K = K, be = beta[1,], a_omega = 0.001, a_eta = 0.001, b_eta = 0.001, B_eta = 0.01, tol = 1e-5, maxIter = 10)
		# tVB			<- tau2b$est

		# sigk		<- matrix(1/lambda, nrow = burn+B, ncol = K)
		# tau2		<- rep(1/omega, burn+B)
	# } else {
		# sigk		<- matrix(1/lambda, nrow = burn+B, ncol = K)
		# tau2		<- rep(1/omega, burn+B)
	# }

	# arTau		<- vector('numeric', length = burn+B)
	# arSig		<- matrix(0, nrow = burn+B, ncol = K)
	
	# # sampler #
	# for(b in 2:(burn+B)){
		
		# if(priorVar == 'HC'){
			# # Sample \tau^2 #
			# tau2b	<- tau2UpdateHC(1, tau = tau2[b-1], tMean = tVB, tTune = tTune, n = n, K = K, be = beta[b-1,], nuTau = nuTau)
			# tau2[b]		<- tau2b$est
			# arTau[b]		<- tau2b$ar			
			
			# # Sample \sigma_k^2 #
			# # sigbk	<- sigma2UpdateHC(1, sig = sigk[b-1,1], sMean = sVB, sTune = sTune, n = n, K = K, th = theta[,,b-1], nuSig = nuSig)
			# # sigk[b,]	<- rep(sigbk$est, K)
			# # arSig[b,]	<- rep(sigbk$ar, K)	
		# } else if(priorVar == 'GA'){
			# # Sample \tau^2 #
			# omegab		<- tau2UpdateGA(1, n, K, beta[b-1,], at, bt)
			# tau2[b]		<- 1/omegab

			# # Sample \sigma_k^2 #
			# lambdab	<- apply(theta[,,b-1], 2, sigma2UpdateGA, samp.size = 1, n = n, K = 1, as = as, bs = bs)
			# sigk[b,]	<- 1/lambdab
		# }
				
		# # loop over outcomes for  #
		# for(k in 1:K){
			# if(priorVar == 'HC'){
				# # Sample \sigma_k^2 #
				# sigbk	<- sigma2UpdateHC(1, sig = sigk[b-1,k], sMean = sVB[k], sTune = sTune[k], n = n, K = 1, th = theta[,k,b-1], nuSig = nuSig)
				# sigk[b,k]	<- sigbk$est
				# arSig[b,k]	<- sigbk$ar
			# }
			
			# # sample \theta_ik #
			# id00k	<- id00[[k]]
			# id10k	<- id10[[k]]
			# id01k	<- id01[[k]]
			# id11k	<- id11[[k]]

			# n00k		<- length(id00k)
			# n10k		<- length(id10k)
			# n01k		<- length(id01k)
			# n11k		<- length(id11k)
			
			# yfix00	<- yfTBSelect(tInt = searchIntT, inc = incT, cent = centT, x1 = 0, x2 = 0, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, mut[k])
			# yfix10	<- yfTBSelect(tInt = searchIntT, inc = incT, cent = centT, x1 = 1, x2 = 0, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, mut[k])
			# yfix01	<- yfTBSelect(tInt = searchIntT, inc = incT, cent = centT, x1 = 0, x2 = 1, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, mut[k])
			# yfix11	<- yfTBSelect(tInt = searchIntT, inc = incT, cent = centT, x1 = 1, x2 = 1, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, mut[k])

			# thUp00k	<- tbUpdate(n00k, yfix00$yfix, x1 = 0, x2 = 0, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, mut[k], ymin, ymax, upFix, dlpThresh)
			# thUp10k	<- tbUpdate(n10k, yfix10$yfix, x1 = 1, x2 = 0, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, mut[k], ymin, ymax, upFix, dlpThresh)
			# thUp01k	<- tbUpdate(n01k, yfix01$yfix, x1 = 0, x2 = 1, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, mut[k], ymin, ymax, upFix, dlpThresh)
			# thUp11k	<- tbUpdate(n11k, yfix11$yfix, x1 = 1, x2 = 1, beta[b-1,k], a1k[b-1,k], a2k[b-1,k], sigk[b-1,k], link, mut[k], ymin, ymax, upFix, dlpThresh)

			# theta[id00k,k,b]	<- thUp00k$thetab
			# theta[id10k,k,b]	<- thUp10k$thetab
			# theta[id01k,k,b]	<- thUp01k$thetab
			# theta[id11k,k,b]	<- thUp11k$thetab

			# art[1,k,b]			<- n00k/thUp00k$count
			# art[2,k,b]			<- n10k/thUp10k$count
			# art[3,k,b]			<- n01k/thUp01k$count
			# art[4,k,b]			<- n11k/thUp11k$count			

			# # \alpha_1k
			# yfixa1k		<- yfAlphaSelect(aInt = searchIntA, inc = incA, cent = centA, theta[,k,b-1], beta[b-1,k], X1[,k], link, mua1[k], siga1[k])
			# yfika1		<- yfixa1k$yfix
			# a1Upk		<- alphaUpdate(1, yfika1, theta[,k,b-1], beta[b-1,k], X1[,k], link, mua1[k], siga1[k], ymin, ymax, upFix, dlpThresh)
			# a1k[b,k]		<- a1Upk$alphab
			# ar1[b,k]		<- a1Upk$count
			
			# # \alpha_2k
			# yfixa2k		<- yfAlphaSelect(aInt = searchIntA, inc = incA, cent = centA, theta[,k,b-1], beta[b-1,k], X2[,k], link, mua2[k], siga2[k])
			# yfika2		<- yfixa2k$yfix
			# a2Upk		<- alphaUpdate(1, yfika2, theta[,k,b-1], beta[b-1,k], X2[,k], link, mua2[k], siga2[k], ymin, ymax, upFix, dlpThresh)
			# a2k[b,k]		<- a2Upk$alphab
			# ar2[b,k]		<- a2Upk$count
			
			# # Sample \beta_k #
			# yfixbk		<- yfTBSelect(tInt = searchIntB, inc = incB, cent = centB, X1[,k], X2[,k], theta[,k,b-1], a1k[b-1,k], a2k[b-1,k], tau2[b-1], link, mub[k])
			# yfbk		<- yfixbk$yfix
			# bkUp		<- tbUpdate(1, yfbk, X1[,k], X2[,k], theta[,k,b-1], a1k[b-1,k], a2k[b-1,k], tau2[b-1], link, mub[k], ymin, ymax, upFix, dlpThresh)
			# beta[b,k]	<- bkUp$thetab
			# arb[b,k]	<- bkUp$count
		# }
		
		# if(priorVar == 'VB'){
			# indv		<- vector(length = K)
			# for(k in 1:K){
				# indv[k]	<- sum(theta[,k,b]^2) > 0.08 & sum(theta[,k,b]^2) < 0.3
			# }
			# if(sum(indv) == K){
				# # Sample \sigma_k^2 #
				# sigbk					<- sigma2VB(n = n, K = K, th = theta[,,b], a_lambda, a_gamma, b_gamma, B_gamma = 0.01, tol, maxIter)
				# sigk[b:(B+burn),]	<- sigbk$est
			# }
			
			# if(sum(beta[b,]^2) > 0.02 & sum(beta[b,]^2) < 0.2){
				# # Sample \tau^2 #
				# tau2b					<- tau2VB(n = 1, K = K, be = beta[b,], a_omega, a_eta, b_eta, B_eta = 0.01, tol = 1e-5, maxIter = 10)
				# # tau2b					<- tau2VB(n = n, K = K, be = beta[b,], a_omega, a_eta, b_eta, B_eta = 0.01, tol = 1e-5, maxIter = 10)
				# tau2[b:(B+burn)]	<- rep(tau2b$est)
			# }
		# }

		# if(mod(b, 10) == 0){
			# cat('.')
		# }
		
		# if(mod(b, up) == 0){
			# cat(paste("\n",b,"samples completed\n"))
		# }

		# # if(b == burn){
			# # cat("\nburn-in complete, begin sampling\n")
		# # }

		# # if(b > burn){
			# # if(mod(b, up) == 0){
				# # cat(paste("\n",b - burn,"samples completed\n"))
			# # }
		# # }
	# }

	# # calculate difference in treatment effects
	# rho		<- a1k - a2k

	# # name output matrices
	# colnames(rho)	<- paste("rho", 1:K, sep = "")
	# colnames(beta)	<- paste("beta", 1:K, sep = "")
	# colnames(sigk)	<- paste("sigma", 1:K, sep = "")
	# colnames(a1k)	<- colnames(X1)
	# colnames(a2k)	<- colnames(X2)
	# colnames(ar1)	<- colnames(X1)
	# colnames(ar2)	<- colnames(X2)
	
	# alphas	<- list(alpha1 = a1k, alpha2 = a2k)
	# aRate	<- list(arTheta = art, arBeta = arb, arAlpha1 = ar1, arAlpha2 = ar2, arTau = arTau, arSigma = arSig)
	# if(priorVar == 'VB'){
		# varComp	<- list(sigk2 = sigbk$est, E_lambda = sigbk$E_lambda, B_lambda = sigbk$B_lambda, tau2 = tau2b$est, E_omega = tau2b$E_omega, B_omega = tau2b$B_omega, fsig = sigk, ftau = tau2)
	# } else if(priorVar == 'FX'){
		# varComp	<- list(sigk2 = rep(1/lambda, K), tau2 = 1/omega)
	# } else {
		# varComp	<- list(sigk2 = sigk, tau2 = tau2)
	# }
	# MCMCsp	<- list(B = B, burn = burn, priorVar = priorVar, pvSpec = pvSpec)		
	# data	<- list(X1 = X1, X2 = X2)
	
	# l		<- list(rho = rho, alphas = alphas, vars = varComp, theta = theta, beta = beta, accept = aRate, MCMCsp = MCMCsp, link = link, data	= data)
	
	# class(l)	<- 'bmir'

	# return(l)
	
# }

## class related function for bmir objects ##

print.bmir		<- function(mod){
	K		<- ncol(mod$rho)
	n		<- dim(mod$theta)[1]
	B		<- mod$MCMCsp$B
	burn	<- mod$MCMCsp$burn
	cat("Bayesian Multivariate Matched Proportions\n")
	cat(paste("K =", K, "outcomes for n =", n, "pairs,", mod$link), "link\n")
	cat("Multivariate Item Response\n")
	
	# rho results
	meds			<- matrix(apply(mod$rho[-c(1:burn),], 2, median), ncol = K, nrow = 1)
	colnames(meds)	<- 1:K
	rownames(meds)	<- ""
	cat("\nMedian Difference in Treatment Effects:\n")
	print(round(meds, digits = 3))
	
}

summary.bmir	<- function(mod, alpha = 0.05, rho0 = 0, chains = 4){
	K			<- ncol(mod$rho)
	n			<- dim(mod$theta)[1]
	B			<- mod$MCMCsp$B
	burn		<- mod$MCMCsp$burn
	priorVar	<- mod$MCMCsp$priorVar
	cat("Bayesian Multivariate Matched Proportions\n")
	cat(paste("K =", K, "outcomes for n =", n, "pairs,", mod$link), "link\n")
	cat("Multivariate Item Response\n")

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

coef.bmir		<- function(mod, CI = FALSE, alpha = 0.05){
	burn	<- mod$MCMCsp$burn

	orTable	<- apply(mod$rho[-c(1:burn),], 2, quantile, probs = c(0.5, alpha/2, (1 - alpha/2)))
	
	if(CI){
		or	<- orTable
	} else {
		or	<- orTable[1, ]
	}
	
	return(or)
}

bsd		<- function(mod, ...){
	UseMethod("bsd")
}

bsd.bmir		<- function(mod){
	burn		<- mod$MCMCsp$burn
	rhoSd		<- apply(mod$rho[-c(1:burn),], 2, sd)
	
	return(rhoSd)
}

pprob	<- function(mod, ...){
	UseMethod("pprob")
}

pprob.bmir	<- function(mod, rho0 = 0){
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


### Independent Multinomials ###
# Bayesian Independent Dirichlet-Multinomial Model: bim #

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

setPrior	<- function(X1, X2){
	if(ncol(X1) != ncol(X2)){
		stop('column dimensions of X1 and X2 must agree')
	}

	if(nrow(X1) != nrow(X2)){
		stop('row dimensions of X1 and X2 must agree')
	}

	K	<- ncol(X1)
	
	sparsityCheck	<- matrix(unlist(lapply(getCounts(X1, X2, K), function(x) x == 0)), nrow = 4, byrow = TRUE)
	
	if(sum(sparsityCheck) > 0){
		jeff		<- rep(1/2, ncol = 4)
		cube		<- rep(1, ncol = 4)
		
		alpha		<- matrix(0, nrow = K, ncol = 4)
		for(k in 1:K){
			if(sum(sparsityCheck[,k]) == 0){
				alpha[k,]	<- jeff
			} else {
				alpha[k,]	<- cube				
			}
		}
	} else {
		alpha		<- matrix(1/2, nrow = K, ncol = 4)
	}
	
	l	<- list(alpha = alpha, sparsityCheck = sparsityCheck)
	
	return(l)
	
}

bim	<- function(X1, X2, B, burnin = NULL, prior = NULL,...){
	if(ncol(X1) != ncol(X2)){
		stop('column dimensions of X1 and X2 must agree')
	}

	if(nrow(X1) != nrow(X2)){
		stop('row dimensions of X1 and X2 must agree')
	}

	if(B < 1){
		stop('B must be at least 1')
	}

	if(is.null(burnin)){
		burnin	<- B
	}
	
	K	<- ncol(X1)
	
	if(is.null(prior)){
		prior		<- setPrior(X1, X2)
		alpha		<- prior$alpha
	} else {
		sparsityCheck	<- matrix(unlist(lapply(getCounts(X1, X2, K), function(x) x == 0)), nrow = 4, byrow = TRUE)
		alpha		<- prior$alpha
		prior$sparsityCheck	<- sparsityCheck
	}
	
	counts	<- getCounts(X1, X2, K)
	
	theta	<- array(0, dim = c(burnin + B, 4, K))
	for(b in 1:(burnin + B)){
		for(k in 1:K){
			xk				    <- c(counts$n21[k], counts$n12[k], counts$n11[k], counts$n22[k])
			theta[b, , k]	<- rdirichlet(1, alpha[k,] + xk)
		}
	}
	
	MCMCsp	<- list(B = B, burnin = burnin, prior = prior)
	data	<- list(X1 = X1, X2 = X2)
	
	l		<- list(theta = theta, MCMCsp = MCMCsp, counts = counts, data = data)
	
	class(l)	<- 'bim'
	
	return(l)
}

print.bim		<- function(mod){
	K		<- dim(mod$theta)[3]
	n		<- dim(mod$data$X1)[1]
	B		<- mod$MCMCsp$B
	burnin	<- mod$MCMCsp$burnin
	cat("Bayesian independent Dirichlet-multinomial model\n")
	cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
	# cat("Independent Multinomials\n")
	
	## odds ratio ##
	or				<- matrix(apply(mod$theta[-(1:burnin),1,]/mod$theta[-(1:burnin),2,], 2, median), ncol = K, nrow = 1)
	colnames(or)	<- 1:K
	rownames(or)	<- ""
	cat("\nMedian Odds Ratio\n")
	print(round(or, digits = 3))

}

summary.bim	<- function(mod, ci_alpha = 0.05, chains = 4){
	K		<- dim(mod$theta)[3]
	n		<- dim(mod$data$X1)[1]
	B		<- mod$MCMCsp$B
	burnin	<- mod$MCMCsp$burnin
	cat("Bayesian independent Dirichlet-multinomial model\n")
	cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
	# cat("Independent Multinomials\n")

	## odds ratio ##
	or		<- mod$theta[-c(1:burnin),1,]/mod$theta[-c(1:burnin),2,]
	low		<- ci_alpha/2
	high	<- 1-(ci_alpha/2)
	tab		<- t(apply(or, 2, quantile, probs = c(0.5, low, high)))
	prob	<- apply(mod$theta[-c(1:burnin),1,] > mod$theta[-c(1:burnin),2,], 2, mean)
	GR		<- getGelman(or, chains)
	tab1	<- cbind(tab, prob, GR)
	colnames(tab1)	<- c('Median', paste(100*low,"%", sep = ""), paste(100*high,"%", sep = ""), "P(OR > 1)", 'GR Est.', 'Upper GR')
	rownames(tab1)	<- 1:K
	cat("\nOdds Ratio Estimates\n\n")
	print(round(tab1, 3))
	cat('Notes: P(OR > 1) corresponds to P(theta21 > theta12).\nGR denotes Gelman-Rubin, < 1.1 suggests convergence.\n')

	cat(paste("\nPosterior samples based on", B, "samples after a burn-in of\n"))
	cat(paste(burnin, "(discarded). "))
	
	sparseTables	<- apply(mod$MCMCsp$prior$sparsityCheck, 2, sum)
	if(sum(sparseTables) > 0){
		cat(paste('Sparse prior used for tables ', paste(which(sparseTables > 0), collapse = ", "), '.', sep = '')) 
	}
}

coef.bim		<- function(mod, CI = FALSE, alpha = 0.05){
	burnin	<- mod$MCMCsp$burnin

	orTable		<- apply(mod$theta[-(1:burnin),1,]/mod$theta[-(1:burnin),2,], 2, quantile, probs = c(0.5, alpha/2, (1 - alpha/2)))
	
	if(CI){
		or	<- orTable
	} else {
		or	<- orTable[1, ]
	}
	
	return(or)
}

bsd.bim		<- function(mod, logOR = FALSE){
	burnin	<- mod$MCMCsp$burnin

	if(logOR){
		orSD	<- apply(log(mod$theta[-(1:burnin),1,]/mod$theta[-(1:burnin),2,]), 2, sd)		
	} else {
		orSD	<- apply(mod$theta[-(1:burnin),1,]/mod$theta[-(1:burnin),2,], 2, sd)
	}
	return(orSD)
}

pprob.bim	<- function(mod){
	burnin	<- mod$MCMCsp$burnin
	probs	<- apply(mod$theta[-c(1:burnin),1,] > mod$theta[-c(1:burnin),2,], 2, mean)
	
	return(probs)
}


### Penalized Marginal Probability Models ###
# Bayesian SParsity Adjusted Matched-proportions model: BSPAM #

bspam	<- function(X1, X2, B, burnin = NULL, penalty = c('Global', 'Local', 't'), prior = NULL, mu = NULL,
                  ptype = c('Normal', 'pe', 'Group'), covstr = 'Independence', up = 100, dots = up/10){
	if(ncol(X1) != ncol(X2)){
		stop('column dimensions of X1 and X2 must agree')
	}

	if(nrow(X1) != nrow(X2)){
		stop('row dimensions of X1 and X2 must agree')
	}

	if(B < 1){
		stop('B must be at least 1')
	}
	
	if(length(penalty) > 1){
		penalty <- penalty[1]
	}
	
	if(penalty != 'Global' & penalty != 'Local' & penalty != 't'){
		stop('penalty must be either Global, Local, or t')
	}

	if(is.null(burnin)){
		burnin    <- B
	}
  
  if(is.null(mu)){
    mu   <- 1/(n + 2)
  }
	
	if(is.null(prior)){
		if(penalty == 'Global'){
			prior	<- list()
			at		<- prior$at	<- 1
			bt		<- prior$bt <- 1
			as		<- prior$as <- 1
			bs		<- prior$bs <- 1
		} else if(penalty == 'Local'){
		  ## CHECK PRIORS ##
			prior	<- list()
			aw		<- prior$aw <- 1
			bw		<- prior$bw <- 1
			aa		<- prior$aa <- 1
			ba		<- prior$ba <- 1
		} else {
			prior	<- list()
			nu		<- prior$nu <- 5
		}
	} else {
		if(penalty == 'Global'){
			at		<- prior$at
			bt		<- prior$bt
			as		<- prior$as
			bs		<- prior$bs
		} else if(penalty == 'Local'){
			aw		<- prior$aw
			bw		<- prior$bw
			aa		<- prior$aa
			ba		<- prior$ba
		} else {
			al		<- prior$al
			bl		<- prior$bl
		}		
	}

	if(length(ptype) > 1){
		ptype <- ptype[1]
	}
	
	if(penalty == 't'){
		prior$type	<- 't'
	} else {
		prior$type	<- ptype
	}
	
	n		<- nrow(X1)
	K		<- ncol(X1)
	Ym		<- as.matrix(cbind(X1, X2))
	ids		<- rep(1:n, each = 2*K)
	outs	<- rep(rep(1:K, 2), n)
	gDat	<- data.frame(y	= c(t(Ym)), s = rep(1:(2*K), n), id = ids, outcome = outs, idk = paste(ids, outs, sep = ''))
	L		<- cbind(diag(K), -diag(K))

	Y		<- gDat$y
	N		<- length(Y)
	X		<- model.matrix(y ~ factor(s) - 1, data = gDat)

	sumX	<- c(apply(X1, 2, sum), apply(X2, 2, sum))
	XsCols	<- c(which(apply(matrix(unlist(lapply(getCounts(X1, X2, K), function(x) x == 0)), nrow = 4, byrow = TRUE), 2, sum) > 0), which(apply(matrix(unlist(lapply(getCounts(X1, X2, K), function(x) x == 0)), nrow = 4, byrow = TRUE), 2, sum) > 0) + K)
	
	if(covstr == 'Independence'){
		if(penalty == 'Global'){
			xtx			                <- t(X)%*%X
			beb                     <- rep(mu, 2*K)
			beb[which(sumX == 0)]   <- mu
			bebf			              <- (sumX + 1)/(n + 2)
			prior$beb               <- beb
			
			beta		<- matrix(0, nrow = B + burnin, ncol = ncol(X))
			lambda		<- vector('numeric', length = B + burnin)
			tau			<- vector('numeric', length = B + burnin)
			sigk		<- matrix(1, nrow = B + burnin, ncol = ncol(X))
			
			lambda[1]	<- N/(sum((Y-X%*%bebf)^2))
			tau[1]		<- (2*K)/(t(bebf)%*%(bebf))
			
			for(b in 2:(B + burnin)){
				## sample betas ##
				precb		<- lambda[b-1]*xtx + tau[b-1]*diag(sigk[b-1,])
				sigb		<- 1/diag(precb)
				mub			<- sigb*(lambda[b-1]*t(X)%*%Y + tau[b-1]*diag(sigk[b-1,])%*%beb)
				beta[b,]	<- rnorm(ncol(X), mean = mub, sd = sqrt(sigb))
				
				## sample lambda ##
				b_lam		<- (1/2)*(t(Y - X%*%beta[b-1,])%*%(Y - X%*%beta[b-1,]))
				lambda[b]	<- rgamma(1, N/2, b_lam)
				
				## sample tau ##
				b_tau		<- bt + (1/2)*sum(sigk[b-1,]*((beta[b-1,] - beb)^2))
				tau[b]		<- rgamma(1, K + at, b_tau)
				
				## sample sigmas ##
				if(ptype == 'pe'){
					b_sigk		<- bs + (tau[b]/2)*(beta[b,] - beb)^2
					sigk[b,]	<- rgamma(length(b_sigk), (1/2) + as, b_sigk)
				}

				if(mod(b, dots) == 0){
					cat('.')
				}
		
				if(mod(b, up) == 0){
					cat(paste("\n",b,"samples completed\n"))
				}
			}
	
			rho			    <- t(L%*%t(beta))
			var.comp	  <- list(lambda = lambda, tau = tau, sigk = sigk)
			mcmc.specs	<- list(B = B, burnin = burnin, penalty = penalty, prior = prior, covstr = covstr)
			data.list	  <- list(Y = Y, X = X, SparseCols = XsCols)
			
			out			<- list(rho = rho, beta = beta, var.comp = var.comp, mcmc.specs = mcmc.specs, data = data.list)
			
		} else if(penalty == 'Local'){
			xtx			      <- t(X)%*%X
			beb           <- rep(0, 2*K)
			# beb[which(sumX == 0)]		<- 1/(n + 2)
			beb[XsCols]   <- mu
			bebf			    <- sumX/n + beb
			prior$beb     <- beb
			
			ps			    <- length(XsCols)
			pn			    <- ncol(X) - length(XsCols)
			
			beta		    <- matrix(0, nrow = B + burnin, ncol = ncol(X))
			lambda      <- vector('numeric', length = B + burnin)
			omegas      <- matrix(1, nrow = B + burnin, ncol = 2)
			gammas      <- matrix(1, nrow = B + burnin, ncol = ncol(X))
			
			lambda[1]   <- N/(sum((Y-X%*%bebf)^2))
			omegas[1,]	<- c(length(XsCols)/(t(bebf[XsCols])%*%(bebf[XsCols])), (ncol(X) - length(XsCols))/(t(bebf[-XsCols])%*%(bebf[-XsCols])))
			beta[1,]    <- bebf

			for(b in 2:(B + burnin)){
				## sample betas ##
				was				    <- rep(0, ncol(X))
				was[XsCols]		<- omegas[b-1,1]*gammas[b-1, XsCols]
				if(ptype == 'Group'){
  				was[-XsCols]	<- omegas[b-1,2]*gammas[b-1, -XsCols]
				}
				precb			    <- lambda[b-1]*xtx + diag(was)
				sigb			    <- 1/diag(precb)
				mub				    <- sigb*(lambda[b-1]*t(X)%*%Y + diag(was)%*%beb)
				beta[b,]		  <- rnorm(ncol(X), mean = mub, sd = sqrt(sigb))
				
				## sample lambda ##
				b_lam			<- (1/2)*(t(Y - X%*%beta[b-1,])%*%(Y - X%*%beta[b-1,]))
				lambda[b]		<- rgamma(1, N/2, b_lam)
			
				## sample omegas ##
				b_omegas		<- bw + (1/2)*sum(gammas[b-1, XsCols]*((beta[b-1, XsCols] - beb[XsCols])^2))
				if(ptype == 'Group'){
				  b_omegan		<- bw + (1/2)*sum(gammas[b-1, -XsCols]*((beta[b-1, -XsCols] - beb[-XsCols])^2))
				}
				omegas[b,1]		<- rgamma(1, ps/2 + aw, b_omegas)
				if(ptype == 'Group'){
				  omegas[b,2]		<- rgamma(1, pn/2 + aw, b_omegan)
				}
				
				## sample gammas ## -- remove for local?
				if(ptype == 'pe'){
					b_gammas			<- ba + (omegas[b,1]/2)*(beta[b, XsCols] - beb[XsCols])^2
					if(ptype == 'Group'){
					  b_gamman			<- ba + (omegas[b,2]/2)*(beta[b, -XsCols] - beb[-XsCols])^2
					}
					gammas[b, XsCols]	<- rgamma(length(b_gammas), (1/2) + aa, b_gammas)
					if(ptype == 'Group'){
					  gammas[b, -XsCols]	<- rgamma(length(b_gamman), (1/2) + aa, b_gamman)
					}
				}
					
				if(mod(b, dots) == 0){
					cat('.')
				}
		
				if(mod(b, up) == 0){
					cat(paste("\n",b,"samples completed\n"))
				}
			}
			
			rho			    <- t(L%*%t(beta))
			var.comp	  <- list(lambda = lambda, omegas = omegas, gammas = gammas)
			mcmc.specs	<- list(B = B, burnin = burnin, penalty = penalty, prior = prior, covstr = covstr)
			data.list	  <- list(Y = Y, X = X, SparseCols = XsCols)
			
			out			    <- list(rho = rho, beta = beta, var.comp = var.comp, mcmc.specs = mcmc.specs, data = data.list)
			
		} else {
		  xtx			                <- t(X)%*%X
		  beb                     <- rep(mu, 2*K)
		  # beb[which(sumX == 0)]   <- mu
		  bebf			              <- (sumX + 1)/(n + 2)
		  prior$beb               <- beb
		  
		  lambda	<- vector('numeric', length = B + burnin)
			beta		<- matrix(0, nrow = B + burnin, ncol = ncol(X))
			tau			<- vector('numeric', length = B + burnin)
			sigk		<- matrix(1, nrow = B + burnin, ncol = ncol(X))

			lambda[1]	<- 1
			tau[1]		<- 1
			beta[1,]	<- bebf
			
			for (b in 2:(B + burnin)) { 
			  ## sample betas ##
			  precb		<- lambda[b-1]*xtx + diag(1/sigk[b-1,])
			  sigb		<- 1/diag(precb)
			  mub			<- sigb*(lambda[b-1]*t(X)%*%Y + diag(1/sigk[b-1,])%*%beb)
			  beta[b,]	<- rnorm(ncol(X), mean = mub, sd = sqrt(sigb))
			  
			  ## sample lambda ##
			  b_lam		<- (1/2)*(t(Y - X%*%beta[b-1,])%*%(Y - X%*%beta[b-1,]))
			  lambda[b]	<- rgamma(1, N/2, b_lam)

			  ## sample sigmas ##
			  b_sigk		<- (nu*tau[b-1]/2) + (1/2)*(beta[b-1,] - beb)^2
			  sigk[b,]	<- rgamma(length(b_sigk), (nu+1)/2, b_sigk)

			  ## sample tau ##
			  b_tau		  <- (nu/2)*sum(1/sigk[b,])
			  tau[b]		<- rgamma(1, (2*K*nu)/2, b_tau)
			  
			  if(mod(b, dots) == 0){
			    cat('.')
			  }
			  
			  if(mod(b, up) == 0){
			    cat(paste("\n",b,"samples completed\n"))
			  }
			}
			
			rho			    <- t(L%*%t(beta))
			var.comp	  <- list(lambda = lambda, sigk = sigk, tau = tau)
			mcmc.specs	<- list(B = B, burnin = burnin, penalty = penalty, prior = prior, covstr = covstr)
			data.list	  <- list(Y = Y, X = X, SparseCols = XsCols)
			
			out			    <- list(rho = rho, beta = beta, var.comp = var.comp, mcmc.specs = mcmc.specs, data = data.list)
	
		}
	} else if(covstr == 'Unstructured'){
		if(penalty == 'Global'){
			xtx		<- t(X)%*%X
			beb                     <- rep(0, 2*K)
			beb[which(sumX == 0)]		<- mu
			bebf			              <- (sumX + 1)/(n + 2)
			ones	                  <- matrix(1, nrow = n, ncol = 1)

			beta	<- matrix(0, nrow = B + burnin, ncol = ncol(X))
			lambda	<- array(0, dim = c(ncol(X), ncol(X), B + burnin))
			tau		<- vector('numeric', length = B + burnin)
			sigk	<- matrix(1, nrow = B + burnin, ncol = ncol(X))
			
			lambda[,,1]	<- diag(ncol(X))
			tau[1]		<- (2*K)/(t(bebf)%*%(bebf))
			beta[1,]	<- beb
			
			for(b in 2:(B + burnin)){
				## sample betas ##
				Lam		<- kronecker(lambda[,,b-1], diag(n))
				xtLx	<- t(X)%*%Lam%*%X
				xtLY	<- t(X)%*%Lam%*%Y
				precb	<- xtLx + tau[b-1]*diag(sigk[b-1,])
				sigb	<- solve(precb)
				mub			<- sigb%*%(xtLY + tau[b-1]*diag(sigk[b-1,])%*%beb)
				beta[b,]	<- rmvnorm(1, mean = mub, sigma = sigb)
				
				## sample lambda ##
				LS			<- solve(t(Ym - ones%*%beta[b-1,])%*%(Ym - ones%*%beta[b-1,]) + diag(ncol(X)))
				lambda[,,b]	<- rwish(n + 2*K + 1, LS)
				
				## sample tau ##
				b_tau		<- bt + (1/2)*sum(sigk[b-1,]*((beta[b-1,] - beb)^2))
				tau[b]		<- rgamma(1, K + at, b_tau)
				
				if(ptype == 'pe'){
					## sample sigmas ##
					b_sigk		<- bs + (tau[b]/2)*(beta[b,] - beb)^2
					sigk[b,]	<- rgamma(length(b_sigk), (1/2) + as, b_sigk)
				}

				if(mod(b, dots) == 0){
					cat('.')
				}
		
				if(mod(b, up) == 0){
					cat(paste("\n",b,"samples completed\n"))
				}
			}

			rho			<- t(L%*%t(beta))
			var.comp	<- list(lambda = lambda, tau = tau, sigk = sigk)
			mcmc.specs	<- list(B = B, burnin = burnin, penalty = penalty, prior = prior, covstr = covstr)
			data.list	<- list(Y = Y, X = X, SparseCols = XsCols)
			
			out			<- list(rho = rho, beta = beta, var.comp = var.comp, mcmc.specs = mcmc.specs, data = data.list)
			
		} else if(penalty == 'Local'){
		  xtx			      <- t(X)%*%X
			beb           <- rep(0, 2*K)
			# beb[which(sumX == 0)]		<- 1/(n + 2)
			beb[XsCols]   <- mu
			bebf			    <- sumX/n + beb
			ones		      <- matrix(1, nrow = n, ncol = 1)
			prior$beb	    <- beb
		
			ps			    <- length(XsCols)
			pn			    <- ncol(X) - length(XsCols)
			
			beta		    <- matrix(0, nrow = B + burnin, ncol = ncol(X))
			lambda		  <- array(0, dim = c(ncol(X), ncol(X), B + burnin))
			omegas		  <- matrix(1, nrow = B + burnin, ncol = 2)
			gammas		  <- matrix(1, nrow = B + burnin, ncol = ncol(X))
			
			lambda[,,1]	<- diag(ncol(X))
			omegas[1,]	<- c(length(XsCols)/(t(bebf[XsCols])%*%(bebf[XsCols])), (ncol(X) - length(XsCols))/(t(bebf[-XsCols])%*%(bebf[-XsCols])))
			beta[1,]	  <- bebf
			
			for(b in 2:(B + burnin)){
				## sample betas ##
				was				    <- rep(0, ncol(X))
				was[XsCols]		<- omegas[b-1,1]*gammas[b-1, XsCols]
				if(ptype == 'Group'){
				  was[-XsCols]	<- omegas[b-1,2]*gammas[b-1, -XsCols]
				}
				Lam				<- kronecker(lambda[,,b-1], diag(n))
				xtLx			<- t(X)%*%Lam%*%X
				xtLY			<- t(X)%*%Lam%*%Y				
				precb			<- xtLx + diag(was)
				sigb			<- solve(precb)
				mub				<- sigb%*%(xtLY + diag(was)%*%beb)
				beta[b,]  <- rmvnorm(1, mean = mub, sigma = sigb)
				
				## sample lambda ##
				LS			<- solve(t(Ym - ones%*%beta[b-1,])%*%(Ym - ones%*%beta[b-1,]) + diag(ncol(X)))
				lambda[,,b]	<- rwish(n + 2*K + 1, LS)
			
				## sample omegas ##
				b_omegas		<- bw + (1/2)*sum(gammas[b-1, XsCols]*((beta[b-1, XsCols] - beb[XsCols])^2))
				if(ptype == 'Group'){
				  b_omegan		<- bw + (1/2)*sum(gammas[b-1, -XsCols]*((beta[b-1, -XsCols] - beb[-XsCols])^2))
				}
				omegas[b,1]		<- rgamma(1, ps/2 + aw, b_omegas)
				if(ptype == 'Group'){
				  omegas[b,2]		<- rgamma(1, pn/2 + aw, b_omegan)
				}
				
				if(ptype == 'pe'){
					## sample gammas ##
					b_gammas			<- ba + (omegas[b,1]/2)*(beta[b, XsCols] - beb[XsCols])^2
					if(ptype == 'Group'){
					  b_gamman			<- ba + (omegas[b,2]/2)*(beta[b, -XsCols] - beb[-XsCols])^2
					}
					gammas[b, XsCols]	<- rgamma(length(b_gammas), (1/2) + aa, b_gammas)
					if(ptype == 'Group'){
					  gammas[b, -XsCols]	<- rgamma(length(b_gamman), (1/2) + aa, b_gamman)
					}
				}
					
				if(mod(b, dots) == 0){
					cat('.')
				}
		
				if(mod(b, up) == 0){
					cat(paste("\n",b,"samples completed\n"))
				}
			}
			
			rho			    <- t(L%*%t(beta))
			var.comp	  <- list(lambda = lambda, omegas = omegas, gammas = gammas)
			mcmc.specs	<- list(B = B, burnin = burnin, penalty = penalty, prior = prior, covstr = covstr)
			data.list	  <- list(Y = Y, X = X, SparseCols = XsCols)
			
			out			    <- list(rho = rho, beta = beta, var.comp = var.comp, mcmc.specs = mcmc.specs, data = data.list)
			
		} else {
			lambda		<- array(0, dim = c(ncol(X), ncol(X), B + burnin))
			beta		<- matrix(0, nrow = B + burnin, ncol = ncol(X))
			xtxi		<- solve(t(X)%*%X) 
			bols		<- xtxi%*%t(X)%*%Y
			ones		<- matrix(1, nrow = n, ncol = 1)
			
			lambda[,,1]	<- diag(ncol(X))
			beta[1,]	<- bols
			
			for (b in 2:(B + burnin)) { 
				# Sample Beta #
				Lam			<- kronecker(lambda[,,b-1], diag(n))
				xtLx		<- t(X)%*%Lam%*%X
				sigb		<- solve(xtLx)				
				beta[b,]	<- c(rmvnorm(1, bols, sigb))
			
				# Sample lambda #
				LS			<- solve(t(Ym - ones%*%beta[b-1,])%*%(Ym - ones%*%beta[b-1,]) + diag(ncol(X)))
				lambda[,,b]	<- rwish(n + 2*K + 1, LS)

				if(mod(b, dots) == 0){
					cat('.')
				}
		
				if(mod(b, up) == 0){
					cat(paste("\n",b,"samples completed\n"))
				}
			}
			
			rho			<- t(L%*%t(beta))
			var.comp	<- list(lambda = lambda)
			mcmc.specs	<- list(B = B, burnin = burnin, penalty = penalty, prior = prior, covstr = covstr)
			data.list	<- list(Y = Y, X = X, SparseCols = XsCols)
			
			out			<- list(rho = rho, beta = beta, var.comp = var.comp, mcmc.specs = mcmc.specs, data = data.list)
			
		}
	} else {
		stop('covstr must be either Independence or Unstructured')
	}
	class(out)	<- 'bspam'
	
	return(out)

}

print.bspam		<- function(mod){
	K		<- ncol(mod$data$X)/2
	n		<- nrow(mod$data$X)/(2*K)
	B		<- mod$mcmc.specs$B
	burnin	<- mod$mcmc.specs$burnin
	cat("Bayesian sparsity adjusted matched-proportions model\n")
	cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
	cat(paste("Covariance structure:", mod$mcmc.specs$covstr,"\n"))
	cat(paste("Penalty:", mod$mcmc.specs$penalty,"\n"))
	
	## difference in marginal probabilities ##
	dmp				<- matrix(apply(mod$rho[-(1:burnin),], 2, median), ncol = K, nrow = 1)
	colnames(dmp)	<- 1:K
	rownames(dmp)	<- ""
	cat("\nMedian difference in marginal probability\n")
	print(round(dmp, digits = 3))

}

summary.bspam	<- function(mod, ci.alpha = 0.05, chains = 4){
	K		<- ncol(mod$data$X)/2
	n		<- nrow(mod$data$X)/(2*K)
	B		<- mod$mcmc.specs$B
	burnin	<- mod$mcmc.specs$burnin
	cat("Bayesian sparsity adjusted matched-proportions model\n")
	cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
	cat(paste("Covariance structure:", mod$mcmc.specs$covstr,"\n"))
	cat(paste("Penalty:", mod$mcmc.specs$penalty,"\n"))

	## difference in marginal probability ##
	rho		<- mod$rho
	low		<- ci.alpha/2
	high	<- 1-(ci.alpha/2)
	tab		<- t(apply(rho[-c(1:burnin),], 2, quantile, probs = c(0.5, low, high)))
	prob	<- apply(rho[-c(1:burnin),] > 0, 2, mean)
	GR		<- getGelman(rho[-c(1:burnin),], chains)
	tab1	<- cbind(tab, prob, GR)
	colnames(tab1)	<- c('Median', paste(100*low,"%", sep = ""), paste(100*high,"%", sep = ""), "P(d > 0)", 'GR Est.', 'Upper GR')
	rownames(tab1)	<- 1:K
	cat("\nDifference in marginal probability\n\n")
	print(round(tab1, 3))
	cat('Notes: P(r > 0) is posterior probability difference > 0.\nGR denotes Gelman-Rubin, < 1.1 suggests convergence.\n')

	cat(paste("\nPosterior samples based on", B, "samples after a burn-in of\n"))
	cat(paste(burnin, "(discarded). "))
	
	# sparseTables	<- apply(mod$MCMCsp$prior$sparsityCheck, 2, sum)
	if(length(mod$data$SparseCols) > 0){
		cat(paste('Sparse tables detected: ', paste(mod$data$SparseCols[1:(length(mod$data$SparseCols)/2)], collapse = ", "), '.', sep = '')) 
	}
}

coef.bspam		<- function(mod, CI = FALSE, ci.alpha = 0.05){
	burnin		<- mod$mcmc.specs$burnin

	dmpTable	<- apply(mod$rho[-(1:burnin),], 2, quantile, probs = c(0.5, ci.alpha/2, (1 - ci.alpha/2)))
	
	if(CI){
		dmpt	<- dmpTable
	} else {
		dmpt	<- dmpTable[1, ]
	}
	
	return(dmpt)
}

bsd.bspam	<- function(mod){
	burnin	<- mod$mcmc.specs$burnin

	dmpSD	<- apply(mod$rho[-(1:burnin),], 2, sd)
	
	return(dmpSD)
}

pprob.bspam	<- function(mod, rho0 = 0){
	burnin	<- mod$mcmc.specs$burnin
	probs	<- apply(mod$rho[-c(1:burnin),] > rho0, 2, mean)
	
	return(probs)
}


### CMH, Liu and Chang 2016 ###
cmhlc	<- function(X1, X2, exact = TRUE, alpha = 0.05, direction = "two.sided"){
	K		<- ncol(X1)
	counts	<- getCounts(X1, X2, K)
	
	est		<- log(counts$n21/counts$n12)
	sdv		<- sqrt(1/counts$n21 + 1/counts$n12)
	
	pval	<- vector('numeric', length = K)
	int		<- matrix(0, nrow = K, ncol = 2)
	for(k in 1:K){
		if(exact){
			n21k	<- counts$n21[k]
			n12k	<- counts$n12[k]
			temp	<- binom.test(n21k, n21k + n12k, p = 0.5, alternative = direction, conf.level = 1 - alpha)
			pval[k]	<- temp$p.value
			int[k,]	<- log(temp$conf.int/(1-temp$conf.int))
		} else {
			int[k,]	<- c(est[k] - qnorm(1-alpha/2)*sdv[k], est[k] + qnorm(1-alpha/2)*sdv[k])
			if(direction == "two.sided"){
				pval[k]	<- 1 - pchisq((est[k]/sdv[k])^2, df = 1)
			} else if(direction == "less"){
				pval[k]	<- pnorm(est[k]/sdv[k])
			} else {
				pval[k]	<- 1 - pnorm(est[k]/sdv[k])				
			}
		}	
	}
		
	l		<- list(est = est, sdv = sdv, int = int, pval = pval, data = list(X1 = X1, X2 = X2), args = list(exact = exact, alpha = alpha, direction = direction))
	
	class(l)	<- 'cmhlc'
	
	return(l)
}

print.cmhlc	<- function(mod){
	K		<- ncol(mod$data$X1)
	n		<- nrow(mod$data$X1)
	cat("CMH-based multivariate matched proportions\n")
	cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))

	## difference in log odds ratios probabilities ##
	dmp				<- matrix(c(mod$est), ncol = K, nrow = 1)
	colnames(dmp)	<- 1:K
	rownames(dmp)	<- ""
	cat("\nLog odds ratio\n")
	print(round(dmp, digits = 3))
	
}

summary.cmhlc	<- function(mod){
	K		<- ncol(mod$data$X1)
	n		<- nrow(mod$data$X1)
	exact	<- mod$args$exact
	cat("CMH-based multivariate matched proportions\n")
	cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))

	## odds ratio ##
	if(exact){
		rho		<- mod$est
		sdv		<- mod$sdv
		pval	<- mod$pval
		tab1	<- cbind(rho, sdv, pval)
		colnames(tab1)	<- c('Estimate', 'Std.err', 'Exact P')
		rownames(tab1)	<- 1:K
		cat("\nLog odds ratio\n\n")
		print(round(tab1, 3))
		cat('Notes: Estimates based on Lui and Chang (2016).\nStd.err estimated using aysmptotic estimator.\n')	
	} else {
		rho		<- mod$est
		sdv		<- mod$sdv
		z		<- rho/sdv
		pval	<- mod$pval
		tab1	<- cbind(rho, sdv, z, pval)
		colnames(tab1)	<- c('Estimate', 'Std.err', 'Z', 'Pr(>|Z|)')
		rownames(tab1)	<- 1:K
		cat("\nLog odds ratio\n\n")
		print(round(tab1, 3))
		cat('Notes: Estimates based on Lui and Chang (2016).\nStd.err estimated using aysmptotic estimator.\n')	
	}
	

}

coef.cmhlc		<- function(mod){
	rho	<- c(mod$est)
		
	return(rho)
}


### Klingenberg and Agresti (2006) GEE ###
gbmp	<- function(X1, X2, family = gaussian(link = 'identity'), corstr = 'independence'){
	n		<- nrow(X1)
	K		<- ncol(X1)
	Y		<- as.matrix(cbind(X1, X2))
	gDat	<- data.frame(y	= c(t(Y)), s = rep(1:(2*K), n), id = rep(1:n, each = 2*K))
	
	model	<- geeglm(y ~ factor(s) - 1, family = family, data = gDat, id = id, corstr = corstr)
	
	beta	<- coef(model)
	covb	<- summary(model)$cov.scaled
	# covb	<- (t(Y - beta)%*%(Y - beta))/(n^2)


	# tests of differences #
	L	      <- cbind(diag(K), -diag(K))
	rho		  <- L%*%beta
	seRho   <- sqrt(diag(L%*%covb%*%t(L)))
	lower	  <- rho - qnorm(0.975)*seRho
	upper	  <- rho + qnorm(0.975)*seRho
	wald	  <- (rho/seRho)^2
	pval	  <- 1 - pchisq(wald, 1)

	## Simultaneous Marginal Homogeneity ##
	# KA Generalized Score Test #
	zse     <- which(diag(covb) == 0)
# 	if(length(zse) == 0){
#   	Tgs		<- (n^2)*t(rho)%*%solve(L%*%(t(Y)%*%Y)%*%t(L))%*%rho
#   	pTgs	<- 1 - pchisq(Tgs, nrow(L))
# 	} else {
	  Tgs		<- NaN
	  pTgs	<- NaN
  # }

	## if estimating 
	if(corstr == 'unstructured'){
		cory	<- matrix(0, 2*K, 2*K)
		cory[lower.tri(cory, diag = FALSE)]	<- summary(model)$corr[,1]
		cory[upper.tri(cory, diag = FALSE)]	<- t(cory)[upper.tri(t(cory), diag = FALSE)]		
		diag(cory)	<- rep(1, 2*K)
	} else {
		cory	<- NULL
	}
	
	l		<- list(rho = rho, beta = beta, sdv = seRho, cov = list(covb = covb, cory = cory, zse = zse), lower = lower, upper = upper, wald = wald, pval = pval, model = model, smh = list(Tgs = Tgs, pTgs = pTgs), data = list(X1 = X1, X2 = X2))

	class(l)	<- 'gbmp'
	
# 	if(zse > 0){
#   	warning('Sparse columns detected, variance estimates may be singular')
# 	}
	
	return(l)
}

print.gbmp	<- function(mod){
	K		<- ncol(mod$data$X1)
	n		<- nrow(mod$data$X1)
	cat("GEE-based multivariate matched proportions model\n")
	cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))

	## difference in marginal probabilities ##
	dmp				<- matrix(c(mod$rho), ncol = K, nrow = 1)
	colnames(dmp)	<- 1:K
	rownames(dmp)	<- ""
	cat("\nDifference in marginal probability\n")
	print(round(dmp, digits = 3))
	
}

summary.gbmp	<- function(mod){
	K		<- ncol(mod$data$X1)
	n		<- nrow(mod$data$X1)
	cat("GEE-based multivariate matched proportions model\n")
	cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))

	## odds ratio ##
	rho		<- mod$rho
	sdv		<- mod$sdv
	wald	<- mod$wald
	pval	<- mod$pval
	tab1	<- cbind(rho, sdv, wald, pval)
	colnames(tab1)	<- c('Estimate', 'Std.err', 'Wald', 'Pr(>|W|)')
	rownames(tab1)	<- 1:K
	cat("\nDifference in marginal probability\n\n")
	print(round(tab1, 3))
	cat(paste('Notes: Estimates based on GEE with family: ', mod$model$family$family, ', link: ', mod$model$family$link, '.\nStd.err estimated using robust estimator.\n', sep = ''))
	
	cat("\nTest of simultaneous marginal homogeneity\n")
	cat(paste("Score statistic:", round(mod$smh$Tgs, 3), "on", K, "df, p-value =", round(mod$smh$pTgs, 3)))

}

coef.gbmp		<- function(mod){
	dmpt	<- c(mod$rho)
		
	return(dmpt)
}


### Westfall, Troendle, and Pennello (2010) ###

bootmmp	<- function(X1, X2, B, int.hyp = FALSE){
	X		<- cbind(X1, X2)
	K		<- ncol(X1)
	n		<- nrow(X1)
	
	L		<- cbind(diag(K), -diag(K))
	
	D		<- t(L%*%t(X))
	Dbar	<- apply(D, 2, mean)
	S		<- diag(apply((t(t(D) - Dbar))^2, 2, mean))
	if(length(which(diag(S) == 0)) > 0){
	  Z <- NaN
	} else{
  	Z		<- sqrt(n)*solve(sqrt(S))%*%Dbar
	}
	
	Dboot	<- matrix(0, nrow = B, ncol = K)
	if(int.hyp){
		Zboot	<- matrix(0, nrow = B, ncol = K)
		Zmax	<- vector('numeric', length = B)
	}

	for(b in 1:B){
		idb			<- sample(1:n, n, replace = TRUE)
		Ds			<- D[idb,]
		Dbb			<- Dboot[b,] <- apply(Ds, 2, mean)
		if(int.hyp){
			Sb			<- sqrt(apply((t(t(Ds) - Dbb))^2, 2, mean))
			Zboot[b,]	<- c(sqrt(n)*ifelse(Sb > 0, Dbb/Sb, 0))
			Zmax[b]		<- max(abs(Zboot[b,]))
		}	
	}
	
	resMat	<- cbind(t(apply(Dboot, 2, quantile, probs = c(0.5, 0.025, 0.975))), apply(Dboot > 0, 2, mean))

	if(int.hyp){
		l	<- list(resMat = resMat, Dboot = Dboot, Zboot = Zboot, Zmax = Zmax, Z = Z, S = S, Dbar = Dbar, data = list(X1 = X1, X2 = X2))	
	} else {
		l	<- list(resMat = resMat, Dboot = Dboot, data = list(X1 = X1, X2 = X2))
	}
	
	class(l)	<- 'bootmmp'
	
	return(l)
}

print.bootmmp	<- function(mod){
	K		<- ncol(mod$data$X1)
	n		<- nrow(mod$data$X1)
	cat("Bootstraped multivariate matched proportions model\n")
	cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))

	## difference in marginal probabilities ##
	dmp				<- matrix(c(mod$resMat[,1]), ncol = K, nrow = 1)
	colnames(dmp)	<- 1:K
	rownames(dmp)	<- ""
	cat("\nDifference in marginal probability\n")
	print(round(dmp, digits = 3))
	
}

summary.bootmmp	<- function(mod, ci.alpha = 0.05){
	K		<- ncol(mod$data$X1)
	n		<- nrow(mod$data$X1)
	cat("Bootstraped multivariate matched proportions model\n")
	cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))

	## difference in marginal probabilities ##
	low				<- ci.alpha/2
	high			<- 1-(ci.alpha/2)
	tab1			<- mod$resMat
	colnames(tab1)	<- c('Median', paste(100*low,"%", sep = ""), paste(100*high,"%", sep = ""), "P(d > 0)")
	rownames(tab1)	<- 1:K
	cat("\nDifference in marginal probability\n\n")
	print(round(tab1, 3))

}

coef.bootmmp	<- function(mod){
	dmpt	<- c(mod$resMat[,1])
		
	return(dmpt)
}
