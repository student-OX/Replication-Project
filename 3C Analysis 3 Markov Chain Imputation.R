############ DATA FRAME FOR IMPUTED DATA
markov.analysis3 <- analysis3


############ MLB001A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB001A), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB001A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB001A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB001A
markov.analysis3$MLB001A[ind] <- imp.dat$ymis


############ MLB041A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB041A), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB041A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB041A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB001a
markov.analysis3$MLB041A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB041A)


############ MLB019A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB019A), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB019A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB019A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB001a
markov.analysis3$MLB019A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB019A)



############ MLB027I
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB027I), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB027I
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027I))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB001a
markov.analysis3$MLB027I[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027I)


############ MLB019F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB019F), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB019F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB019F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB019F
markov.analysis3$MLB019F[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB019F)


############ MLB027A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB027A), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB027A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB027A
markov.analysis3$MLB027A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027A)





############ MLB001a
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB001a), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB001a
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB001a))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB001a
markov.analysis3$MLB001a[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB001a)



############ MLB035A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB035A), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB035A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB035A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB035A
markov.analysis3$MLB035A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB035A)


############ MLB035B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB035B), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB035B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB035B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB035B
markov.analysis3$MLB035B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB035B)



############ MLB035C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB035C), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB035C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB035C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB035C
markov.analysis3$MLB035C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB035C)




############ MLB035D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB035D), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB035D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB035D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB035D
markov.analysis3$MLB035D[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB035D)




############ MLB035F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB035F), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB035F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB035F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB035F
markov.analysis3$MLB035F[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB035F)


############ MLB035G
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB035G), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB035G
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB035G))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB035G
markov.analysis3$MLB035G[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB035G)



############ MLB041B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB041B), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB041B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB041B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB041B
markov.analysis3$MLB041B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB041B)




############ MLB027N
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB027N), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB027N
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027N))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB027N
markov.analysis3$MLB027N[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027N)


############ MLB019B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB019B), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB019B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB019B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB019B
markov.analysis3$MLB019B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB019B)




############ MLB027J
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB027J), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB027J
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027J))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB027J
markov.analysis3$MLB027J[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027J)


############ MLB019G
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB019G), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB019G
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB019G))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB019G
markov.analysis3$MLB019G[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB019G)




############ MLB027B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB027B), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB027B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB027B
markov.analysis3$MLB027B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027B)





############ MLB001b
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB001b), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB001b
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB001b))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB001b
markov.analysis3$MLB001b[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB001b)



############ MLB027J
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB027J), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB027J
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027J))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB027J
markov.analysis3$MLB027J[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027J)




############ MLB041C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB041C), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB041C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB041C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB041C
markov.analysis3$MLB041C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB041C)




############ MLB019C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB019C), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB019C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB019C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB019C
markov.analysis3$MLB019C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB019C)


############ MLB027K
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB027K), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB027K
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027K))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB027K
markov.analysis3$MLB027K[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027K)


############ MLB019H
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB019H), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB019H
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB019H))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB019H
markov.analysis3$MLB019H[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB019H)




############ MLB027C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB027C), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB027C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB027C
markov.analysis3$MLB027C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027C)



############ MLB001C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB001C), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB001C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB001C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB001C
markov.analysis3$MLB001C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB001C)



############ MLB041D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB041D), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB041D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB041D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB041D
markov.analysis3$MLB041D[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB041D)



############ MLB019D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB019D), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB019D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB019D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB019D
markov.analysis3$MLB019D[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB019D)



############ MLB027L
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB027L), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB027L
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027L))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB027L
markov.analysis3$MLB027L[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027L)



############ MLB019I
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB019I), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB019I
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB019I))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB019I
markov.analysis3$MLB019I[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB019I)



############ MLB027D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB027D), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB027D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB027D
markov.analysis3$MLB027D[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027D)


############ MLB041E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB041E), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB041E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB041E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB041E
markov.analysis3$MLB041E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB041E)



############ MLB019E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB019E), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB019E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB019E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB019E
markov.analysis3$MLB019E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB019E)


############ MLB027M
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB027M), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB027M
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027M))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB027M
markov.analysis3$MLB027M[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027M)



############# Add imputed values to MLB019E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB019E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB019E
markov.analysis3$MLB019E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB019E)


############ MLB027M
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB027M), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB027M
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027M))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB027M
markov.analysis3$MLB027M[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027M)




############ MLB019J
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB019J), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB019J
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB019J))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB019J
markov.analysis3$MLB019J[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB019J)



############ MLB027E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB027E), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB027E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB027E
markov.analysis3$MLB027E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027E)



############ MLB001E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB001E), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB001E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB001E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB001E
markov.analysis3$MLB001E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB001E)





############ MLB027M
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB027M), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB027M
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027M))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB027M
markov.analysis3$MLB027M[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027M)




############ MLB035E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB035E), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB035E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB035E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB035E
markov.analysis3$MLB035E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB035E)




############ MLB019K
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB019K), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB019K
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB019K))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB019K
markov.analysis3$MLB019K[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB019K)





############ MLB027F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB027F), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB027F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB027F
markov.analysis3$MLB027F[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027F)



############ MLB001F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB001F), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB001F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB001F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB001F
markov.analysis3$MLB001F[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB001F)




############ MLB001G
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis3.mat <- matrix(c(analysis3$sex, analysis3$MLB001G), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB001G
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB001G))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB001G
markov.analysis3$MLB001G[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB001G)



############ MLB027M.1
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis3$sex, analysis3$MLB027M.1), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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
  this.imp <- analysis3.mat
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


############# Add imputed values to MLB027M.1
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027M.1))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB027M.1
markov.analysis3$MLB027M.1[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027M.1)



############ MLB001B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis3$sex, analysis3$MLB001B), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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


############# Add imputed values to MLB001B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB001B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB001B
markov.analysis3$MLB001B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB001B)




############ MLB027J.1
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis2.mat <- matrix(c(analysis3$sex, analysis3$MLB027J.1), ncol=2, byrow=FALSE)

N <- nrow(analysis3.mat)

nvar <- ncol(analysis3.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis3.mat)


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


############# Add imputed values to MLB027J.1
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis3$MLB027J.1))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to MLB027J.1
markov.analysis3$MLB027J.1[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis3$MLB027J.1)
