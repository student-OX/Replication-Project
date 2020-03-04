############ DATA FRAME FOR IMPUTED DATA
markov.analysis4 <- analysis4

############ NLB001A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB001A), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB001A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB001A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB001A
markov.analysis4$NLB001A[ind] <- imp.dat$ymis

# View imputed data
summary(markov.analysis4$NLB001A)

############ NLB041A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB041A), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB041A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB041A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB001a
markov.analysis4$NLB041A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB041A)


############ NLB019A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB019A), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB019A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB019A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB001a
markov.analysis4$NLB019A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB019A)



############ NLB027I
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027I), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027I
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027I))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB001a
markov.analysis4$NLB027I[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027I)


############ NLB019F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB019F), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB019F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB019F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB019F
markov.analysis4$NLB019F[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB019F)


############ NLB027A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027A), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB027A
markov.analysis4$NLB027A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027A)





############ NLB001a
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB001a), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB001a
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB001a))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB001a
markov.analysis4$NLB001a[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB001a)



############ NLB035A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB035A), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB035A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB035A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB035A
markov.analysis4$NLB035A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB035A)


############ NLB035B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB035B), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB035B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB035B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB035B
markov.analysis4$NLB035B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB035B)



############ NLB035C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB035C), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB035C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB035C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB035C
markov.analysis4$NLB035C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB035C)




############ NLB035D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB035D), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB035D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB035D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB035D
markov.analysis4$NLB035D[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB035D)




############ NLB035F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB035F), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB035F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB035F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB035F
markov.analysis4$NLB035F[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB035F)


############ NLB035G
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB035G), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB035G
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB035G))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB035G
markov.analysis4$NLB035G[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB035G)



############ NLB041B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB041B), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB041B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB041B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB041B
markov.analysis4$NLB041B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB041B)




############ NLB027N
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027N), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027N
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027N))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB027N
markov.analysis4$NLB027N[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027N)


############ NLB019B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB019B), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB019B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB019B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB019B
markov.analysis4$NLB019B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB019B)




############ NLB027J
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027J), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027J
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027J))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB027J
markov.analysis4$NLB027J[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027J)


############ NLB019G
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB019G), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB019G
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB019G))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB019G
markov.analysis4$NLB019G[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB019G)




############ NLB027B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027B), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB027B
markov.analysis4$NLB027B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027B)





############ NLB001b
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB001b), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB001b
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB001b))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB001b
markov.analysis4$NLB001b[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB001b)



############ NLB027J
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027J), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027J
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027J))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB027J
markov.analysis4$NLB027J[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027J)




############ NLB041C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB041C), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB041C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB041C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB041C
markov.analysis4$NLB041C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB041C)




############ NLB019C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB019C), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB019C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB019C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB019C
markov.analysis4$NLB019C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB019C)


############ NLB027K
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027K), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027K
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027K))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB027K
markov.analysis4$NLB027K[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027K)


############ NLB019H
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB019H), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB019H
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB019H))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB019H
markov.analysis4$NLB019H[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB019H)




############ NLB027C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027C), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB027C
markov.analysis4$NLB027C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027C)



############ NLB001C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB001C), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB001C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB001C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB001C
markov.analysis4$NLB001C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB001C)



############ NLB041D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB041D), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB041D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB041D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB041D
markov.analysis4$NLB041D[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB041D)



############ NLB019D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB019D), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB019D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB019D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB019D
markov.analysis4$NLB019D[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB019D)



############ NLB027L
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027L), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027L
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027L))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB027L
markov.analysis4$NLB027L[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027L)



############ NLB019I
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB019I), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB019I
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB019I))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB019I
markov.analysis4$NLB019I[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB019I)



############ NLB027D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027D), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB027D
markov.analysis4$NLB027D[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027D)


############ NLB041E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB041E), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB041E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB041E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB041E
markov.analysis4$NLB041E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB041E)



############ NLB019E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB019E), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB019E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB019E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB019E
markov.analysis4$NLB019E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB019E)


############ NLB027M
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027M), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027M
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027M))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB027M
markov.analysis4$NLB027M[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027M)



############# Add imputed values to NLB019E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB019E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB019E
markov.analysis4$NLB019E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB019E)


############ NLB027M
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027M), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027M
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027M))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB027M
markov.analysis4$NLB027M[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027M)




############ NLB019J
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB019J), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB019J
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB019J))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB019J
markov.analysis4$NLB019J[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB019J)



############ NLB027E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027E), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB027E
markov.analysis4$NLB027E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027E)



############ NLB001E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB001E), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB001E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB001E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB001E
markov.analysis4$NLB001E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB001E)





############ NLB027M
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027M), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027M
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027M))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB027M
markov.analysis4$NLB027M[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027M)




############ NLB035E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB035E), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB035E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB035E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB035E
markov.analysis4$NLB035E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB035E)




############ NLB019K
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB019K), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB019K
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB019K))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB019K
markov.analysis4$NLB019K[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB019K)





############ NLB027F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027F), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB027F
markov.analysis4$NLB027F[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027F)



############ NLB001F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB001F), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB001F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB001F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB001F
markov.analysis4$NLB001F[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB001F)




############ NLB001G
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB001G), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB001G
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB001G))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB001G
markov.analysis4$NLB001G[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB001G)


############ NLB001B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB001B), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB001B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB001B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB001B
markov.analysis4$NLB001B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB001B)


############ NLB027J.1
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027J.1), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027J.1
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027J.1))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB027J.1
markov.analysis4$NLB027J.1[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027J.1)




############ NLB027M.1
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis4.mat <- matrix(c(analysis4$sex, analysis4$NLB027M.1), ncol=2, byrow=FALSE)

N <- nrow(analysis4.mat)

nvar <- ncol(analysis4.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis4.mat)


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
  this.imp <- analysis4.mat
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


############# Add imputed values to NLB027M.1
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis4$NLB027M.1))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to NLB027M.1
markov.analysis4$NLB027M.1[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis4$NLB027M.1)

