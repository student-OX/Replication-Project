############ DATA FRAME FOR IMPUTED DATA
markov.analysis2 <- analysis2


############ LLB001a
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB001a), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  this.imp <- analysis2.mat
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


############# Add imputed values to LLB001a
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB001a))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB001a
markov.analysis2$LLB001a[ind] <- imp.dat$ymis


############ LLB041A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB041A), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB041A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB041A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB001a
markov.analysis2$LLB041A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB041A)


############ LLB019A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB019A), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB019A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB019A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB001a
markov.analysis2$LLB019A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB019A)



############ LLB027I
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027I), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027I
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027I))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB001a
markov.analysis2$LLB027I[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027I)


############ LLB019F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB019F), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB019F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB019F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB019F
markov.analysis2$LLB019F[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB019F)


############ LLB027A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027A), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB027A
markov.analysis2$LLB027A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027A)





############ LLB001a
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB001a), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB001a
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB001a))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB001a
markov.analysis2$LLB001a[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB001a)



############ LLB035A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB035A), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB035A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB035A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB035A
markov.analysis2$LLB035A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB035A)


############ LLB035B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB035B), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB035B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB035B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB035B
markov.analysis2$LLB035B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB035B)



############ LLB035C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB035C), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB035C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB035C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB035C
markov.analysis2$LLB035C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB035C)




############ LLB035D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB035D), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB035D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB035D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB035D
markov.analysis2$LLB035D[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB035D)




############ LLB035F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB035F), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB035F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB035F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB035F
markov.analysis2$LLB035F[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB035F)


############ LLB035G
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB035G), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB035G
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB035G))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB035G
markov.analysis2$LLB035G[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB035G)



############ LLB041B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB041B), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB041B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB041B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB041B
markov.analysis2$LLB041B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB041B)




############ LLB027N
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027N), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027N
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027N))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB027N
markov.analysis2$LLB027N[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027N)


############ LLB019B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB019B), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB019B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB019B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB019B
markov.analysis2$LLB019B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB019B)




############ LLB027J
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027J), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027J
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027J))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB027J
markov.analysis2$LLB027J[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027J)


############ LLB019G
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB019G), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB019G
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB019G))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB019G
markov.analysis2$LLB019G[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB019G)




############ LLB027B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027B), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB027B
markov.analysis2$LLB027B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027B)





############ LLB001b
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB001b), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB001b
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB001b))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB001b
markov.analysis2$LLB001b[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB001b)



############ LLB027J
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027J), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027J
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027J))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB027J
markov.analysis2$LLB027J[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027J)




############ LLB041C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB041C), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB041C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB041C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB041C
markov.analysis2$LLB041C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB041C)




############ LLB019C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB019C), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB019C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB019C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB019C
markov.analysis2$LLB019C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB019C)


############ LLB027K
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027K), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027K
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027K))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB027K
markov.analysis2$LLB027K[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027K)


############ LLB019H
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB019H), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB019H
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB019H))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB019H
markov.analysis2$LLB019H[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB019H)




############ LLB027C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027C), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB027C
markov.analysis2$LLB027C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027C)



############ LLB001c
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB001c), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB001c
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB001c))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB001c
markov.analysis2$LLB001c[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB001c)



############ LLB041D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB041D), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB041D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB041D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB041D
markov.analysis2$LLB041D[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB041D)



############ LLB019D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB019D), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB019D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB019D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB019D
markov.analysis2$LLB019D[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB019D)



############ LLB027L
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027L), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027L
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027L))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB027L
markov.analysis2$LLB027L[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027L)



############ LLB019I
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB019I), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB019I
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB019I))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB019I
markov.analysis2$LLB019I[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB019I)



############ LLB027D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027D), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB027D
markov.analysis2$LLB027D[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027D)


############ LLB041E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB041E), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB041E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB041E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB041E
markov.analysis2$LLB041E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB041E)



############ LLB019E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB019E), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB019E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB019E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB019E
markov.analysis2$LLB019E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB019E)


############ LLB027M
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027M), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027M
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027M))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB027M
markov.analysis2$LLB027M[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027M)



############# Add imputed values to LLB019E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB019E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB019E
markov.analysis2$LLB019E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB019E)


############ LLB027M
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027M), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027M
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027M))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB027M
markov.analysis2$LLB027M[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027M)




############ LLB019J
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB019J), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB019J
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB019J))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB019J
markov.analysis2$LLB019J[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB019J)



############ LLB027E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027E), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB027E
markov.analysis2$LLB027E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027E)



############ LLB001e
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB001e), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB001e
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB001e))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB001e
markov.analysis2$LLB001e[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB001e)





############ LLB027M
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027M), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027M
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027M))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB027M
markov.analysis2$LLB027M[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027M)




############ LLB035E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB035E), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB035E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB035E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB035E
markov.analysis2$LLB035E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB035E)




############ LLB019K
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB019K), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB019K
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB019K))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB019K
markov.analysis2$LLB019K[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB019K)





############ LLB027F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027F), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB027F
markov.analysis2$LLB027F[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027F)



############ LLB001f
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB001f), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB001f
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB001f))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB001f
markov.analysis2$LLB001f[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB001f)




############ LLB001g
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB001g), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB001g
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB001g))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB001g
markov.analysis2$LLB001g[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB001g)


############ LLB001A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB001A), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB001A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB001A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB001A
markov.analysis2$LLB001A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB001A)



############ LLB001B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB001B), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB001B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB001B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB001B
markov.analysis2$LLB001B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB001B)



############ LLB001C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB001C), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB001C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB001C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB001C
markov.analysis2$LLB001C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB001C)



############ LLB001F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB001F), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB001F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB001F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB001F
markov.analysis2$LLB001F[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB001F)




############ LLB001G
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB001G), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB001G
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB001G))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB001G
markov.analysis2$LLB001G[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB001G)

 


############ LLB027J.1
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027J.1), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027J.1
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027J.1))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB027J.1
markov.analysis2$LLB027J.1[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027J.1)



############ LLB027M.1
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis2$sex, analysis2$LLB027M.1), ncol=2, byrow=FALSE)

N <- nrow(analysis2.mat)

nvar <- ncol(analysis2.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis2.mat)


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
  
  imp.dat <- da.norm(prelim.mi,theta.init,steps=1,return.ymis=TRUE)   
  
  misobs <- imp.dat$ymis
  this.imp <- analysis2.mat
  this.imp[is.na(this.imp)] <- misobs
  imputed.data <- rbind(imputed.data,this.imp)
  
  this.imp.alt <- NULL
  for (i in 1:N){
    this <-cbind(rep(i,nvar),seq(1,nvar,1),this.imp[i,1:2,drop=TRUE])
    this.imp.alt <- rbind(this.imp.alt,this)
  }
  
  #  Call gls() to fit the multivariate normal model to the reconfigured imputed data set
  
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

#> mi.results
mi.mv.inference <- function (est, covmat, confidence = 0.95){
  
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
  
  #   First 3 elements are the mean of estimates, their SEs entire covariance matrix using Rubin's formula; last 2
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

#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat

## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to LLB027M.1
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis2$LLB027M.1))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to LLB027M.1
markov.analysis2$LLB027M.1[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis2$LLB027M.1)

