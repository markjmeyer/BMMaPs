# Bayesian Multivariate Matched Proportions (BMMaPs)
This repository contains code for Running Bayesian Multivariate Matched Proportions (BMMaPs) as described by Meyer and Knutson (2020), currently under review, and CSP 2020 in Sacramento.

# File details
The file bmmaps.R is the source file for all functions related to running the models.

# Brief Function Description
All functions two matrices, X1 and X2, as arguments. X1 is an n x K matrix of the first set of binary multivariate outcomes. X2 is an n x K matrix of the second or paired set of binary multivariate outcomes. Each function as additional arguments including prior specification and link.


Bayesian models:

bim(X1, X2, B, burnin, ...) - Bayesian Independent Multinomial with Dirichlet prior model

bspam(X1, X2, B, burnin, ...) - Bayesian SParsity Adjusted Matched-proportions model

bmir(X1, X2, link = 'logit', ...) - Bayesian Multivariate Item Response model


Non-Bayesian models:

cmhlc(X1, X2, ...) - CMH, Liu and Chang (2016) based model

gbmp(X1, X2, ...) - GEE, Klingenberg and Agresti (2006) based model

bootmmp(X1, X2, B, ...) - Bootstrap, Westfall, Troendle, and Pennello (2010) based model
