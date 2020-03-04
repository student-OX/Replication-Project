############ DATA FRAME FOR IMPUTED DATA
markov.analysis5 <- analysis5


############ OLB026S
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB026S), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat) 

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat) 


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB026S
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB026S))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB026S
markov.analysis5$OLB026S[ind] <- imp.dat$ymis



############ OLB018A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB018A), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB018A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB018A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB018A
markov.analysis5$OLB018A[ind] <- imp.dat$ymis



############ OLB001I
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB001I), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB001I
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB001I))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB001I
markov.analysis5$OLB001I[ind] <- imp.dat$ymis



############ OLB033A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB033A), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB033A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB033A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB033A
markov.analysis5$OLB033A[ind] <- imp.dat$ymis




############ OLB033B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB033B), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB033B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB033B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB033B
markov.analysis5$OLB033B[ind] <- imp.dat$ymis





############ OLB033C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB033C), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB033C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB033C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB033C
markov.analysis5$OLB033C[ind] <- imp.dat$ymis




############ OLB033D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB033D), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB033D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB033D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB033D
markov.analysis5$OLB033D[ind] <- imp.dat$ymis




############ OLB033F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB033F), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB033F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB033F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB033F
markov.analysis5$OLB033F[ind] <- imp.dat$ymis



############ OLB033G
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB033G), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB033G
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB033G))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB033G
markov.analysis5$OLB033G[ind] <- imp.dat$ymis




############ OLB031L
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB031L), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB031L
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB031L))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB031L
markov.analysis5$OLB031L[ind] <- imp.dat$ymis



############ OLB026R
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB026R), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB026R
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB026R))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB026R
markov.analysis5$OLB026R[ind] <- imp.dat$ymis



############ OLB018H
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB018H), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB018H
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB018H))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB018H
markov.analysis5$OLB018H[ind] <- imp.dat$ymis



############ OLB018B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB018B), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB018B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB018B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB018B
markov.analysis5$OLB018B[ind] <- imp.dat$ymis



OLB026J	OLB018C	OLB026K	OLB018D	OLB026X	OLB018E	OLB002C	OLB033E	OLB018F	OLB001N	OLB012B


############ OLB001R
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB001R), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB001R
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB001R))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB001R
markov.analysis5$OLB001R[ind] <- imp.dat$ymis


############ OLB026J
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB026J), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB026J
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB026J))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB026J
markov.analysis5$OLB026J[ind] <- imp.dat$ymis


############ OLB018C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB018C), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB018C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB018C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB018C
markov.analysis5$OLB018C[ind] <- imp.dat$ymis



############ OLB026K
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB026K), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB026K
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB026K))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB026K
markov.analysis5$OLB026K[ind] <- imp.dat$ymis



############ OLB018D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB018D), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB018D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB018D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB018D
markov.analysis5$OLB018D[ind] <- imp.dat$ymis



############ OLB026X
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB026X), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB026X
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB026X))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB026X
markov.analysis5$OLB026X[ind] <- imp.dat$ymis




############ OLB018E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB018E), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB018E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB018E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB018E
markov.analysis5$OLB018E[ind] <- imp.dat$ymis



############ OLB002C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB002C), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB002C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB002C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB002C
markov.analysis5$OLB002C[ind] <- imp.dat$ymis



############ OLB033E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB033E), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB033E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB033E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB033E
markov.analysis5$OLB033E[ind] <- imp.dat$ymis




############ OLB018F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB018F), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB018F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB018F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB018F
markov.analysis5$OLB018F[ind] <- imp.dat$ymis


############ OLB001N
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB001N), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB001N
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB001N))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB001N
markov.analysis5$OLB001N[ind] <- imp.dat$ymis



############ OLB012B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis5.mat <- matrix(c(analysis5$sex, analysis5$OLB012B), ncol=2, byrow=FALSE)

N <- nrow(analysis5.mat)

nvar <- ncol(analysis5.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis5.mat)


#  View missing values and missing data patterns
prelim.mi$nmis
prelim.mi$r


# Use EM algorithm to estimate matrix of incomplete data and obtain the form of the starting values
theta.init <- em.norm(prelim.mi, maxits = 10)


# Extract the parameters to view the format
theta.init <- getparam.norm(prelim.mi,theta.init,corr=TRUE)

theta.init$mu    #  the mean vector
theta.init$sdv   #  the standard deviations
theta.init$r     #  the correlation matrix



#  Convert parameters to the form used by the function em.norm
theta.init <- makeparam.norm(prelim.mi,theta.init)


### Create empty imputed datasets and save each dataset as a large matrix
imputed.data <- NULL


## Set parameters to be saved as lists
mu.list <- vector("list",M)
ses.list <- vector("list",M)

## Set asymptomatic covariance matrices for use with multivariate mi inference function
covs.list <- vector("list",M)
covparms.list <- vector("list",M)


#  Set the seed for the MCMC

rngseed(1234)

for (m in 1:M){
  
  #  Impute the missing values to create a single data set using MCMC.
  #  Use the da.norm (for "data augmentation") function in the norm package
  #  to impute the missing values. Here, steps is apparently the number
  #  of iterations in the chain between draws - the documentation does not give
  #  details, and it is not clear if there is a burn-in period
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  #  misobs is a vector containing the imputed missing values in the order
  #  they are missing in the original data set, so insert them into a copy
  #  of the original data set.
  
  misobs <- imp.dat$ymis
  this.imp <- analysis5.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set in case you want to do other stuff with it
  
  imputed.data <- rbind(imputed.data,this.imp)
  
  #  Reconfigure the imputed data set to 1 record per observation 
  #  for use with gls()
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured
  #  imputed data set
  
  id <- factor(this.imp.alt[,1])
  ind <- factor(this.imp.alt[,2])
  y <- this.imp.alt[,3]
  time <- rep(seq(1,2,1),N)
  
  gls.fit <- gls(y ~ -1 + ind,correlation=corSymm(form = ~time | id),
                 weights=varIdent(form= ~1 | ind),method="ML")
  
  #  Get the estimates of the mean parameters and their standard errors
  #  and also the full asymptotic covariance matrices
  
  this.mu <- gls.fit$coef
  this.ses <- sqrt(diag(gls.fit$varBeta))
  this.cov <- gls.fit$varBeta
  
  #  Save these in the lists, as the mi.inference function
  #  requires that the parameter estimates from each imputed data set and
  #  their standard errors be in lists
  
  mu.list[[m]] <- this.mu
  ses.list[[m]] <- this.ses
  covs.list[[m]] <- this.cov
  
  #  Save the estimates of the variance and covariance parameters   
  this.covmat <- getVarCov(gls.fit)
  this.covparms <- c(this.covmat[1,1],this.covmat[1,2],this.covmat[2,2])
  
  covparms.list[[m]] <- this.covparms
  
  #  Unfortunately, getting their standard errors is hard -- they are computed
  #  but are available in a different (and not clear) parameterization via 
  
  #  gls.fit$apVar
  
  #  We won't bother trying to extract them; usually these aren't of
  #  interest anyway 
  
}  

#  Now call the function mi.inference() to obtain the multiple
#  imputation estimate for mu and the ses via Rubin's formula.  This
#  function requires the parameters and their standard errors to be in
#  lists, as noted aboe

mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#  Results of running this as above

#> mi.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
#ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
#ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  If we are willing to do a little of our own programming, we
#  can write a function to do the full multivariate analysis.  This
#  one takes as input list of the parameters and asymptotic covariance
#  matrices from each imputation and returns the imputation estimator,
#  the entire covariance matrix of this estimator, the within and among
#  components thereof, and the associated standard errors/confidence limits. 

mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
  #  Get the mean of estimates and mean of (within) covariance matrices 
  #  (presumably these lists are of the same length; don't bother
  #  checking this)
  
  m <- length(est)
  qmat <- simplify2array(est)
  qbar <- apply(qmat,1,mean)
  wcov <- Reduce('+',covmat)/m
  
  #  Get among-imputation covariance matrix    
  
  bcov <- sweep(data.matrix(qmat),1,qbar)%*%t(sweep(data.matrix(qmat),1,qbar))/(m-1)
  
  #   Rubin covariance matrix and diagonal elements of each component
  
  qcovmat <- wcov + (1+1/m)*bcov
  bm <- apply(simplify2array(est),1,var)  #  should = diag(bcov)
  ubar <- diag(wcov)  
  
  #   This code is from mi.inference - CIs, DFs, etc for each component
  
  tm <- ubar + ((1 + (1/m)) * bm)  #  should - diag(qcovmat)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  
  #   First 3 elements are the mean of estimates, their SEs,   
  #   entire covariance matrix using Rubin's formula; last 2
  #   are the within and among components 
  
  result <- list(est = qbar, std.err = sqrt(tm), cov.mat = qcovmat,
                 df = nu, signif = pval,lower = low, upper = up, r =
                   rem, fminf = fminf, within = wcov, between = bcov)
  result
}

#  Call the multivariate function 

mi.mv.mu <- mi.mv.inference(mu.list,covs.list,confidence=0.95)

#  Display the results -- should be identical to mi.inference

mi.mv.results <- cbind(mi.mv.mu$est,mi.mv.mu$std.err,mi.mv.mu$lower,mi.mv.mu$upper,
                       mi.mv.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.mv.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

#> mi.mv.results
#          Est     StdErr    Lower    Upper       DF      Min      Max
# ind1 4.986591 0.03711372 4.913196 5.059986 135.8841 4.953036 5.018024
# ind2 7.923550 0.03669786 7.851316 7.995784 283.9978 7.901227 7.950898

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to OLB012B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis5$OLB012B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to OLB012B
markov.analysis5$OLB012B[ind] <- imp.dat$ymis




