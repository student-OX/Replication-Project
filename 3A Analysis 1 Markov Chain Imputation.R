############## MARKOV CHAIN MULTIVARIATE IMPUTATION OF MISSING PSYCOLOGICAL VARIABLES
##### Load Packages
library(norm)
library(mice)
library(nlme)


############ DATA FRAME FOR IMPUTED DATA
markov.analysis1 <- analysis1

############ KLB001a
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB001a), ncol=2, byrow=FALSE)
                        
N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
  this.imp[is.na(this.imp)] <- misobs
  
  #  Save this imputed data set 
  
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
  
}  


mi.mu <- mi.inference(mu.list,ses.list,confidence=0.95)    #  confidence=0.95 is the default

#  Display the results in the same format as SAS proc mianalyze

mu.imp <- simplify2array(mu.list)

mi.results <- cbind(mi.mu$est,mi.mu$std.err,mi.mu$lower,mi.mu$upper,mi.mu$df,apply(mu.imp,1,min),apply(mu.imp,1,max))
colnames(mi.results) <- c("Est","StdErr","Lower","Upper","DF","Min","Max")

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


#  get the components of the full Rubin covariance matrix

within.cov <- mi.mv.mu$within
between.cov <- mi.mv.mu$between
Rubin.cov <- mi.mv.mu$cov.mat


## View imputed values and length
summary(imp.dat$ymis)
length(imp.dat$ymis)


############# Add imputed values to KLB001a
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB001a))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB001a
markov.analysis1$KLB001a[ind] <- imp.dat$ymis


############ KLB041A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB041A), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB041A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB041A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB001a
markov.analysis1$KLB041A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB041A)


############ KLB019A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB019A), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB019A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB019A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB001a
markov.analysis1$KLB019A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB019A)



############ KLB027I
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027I), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027I
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027I))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB001a
markov.analysis1$KLB027I[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027I)


############ KLB019F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB019F), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB019F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB019F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB019F
markov.analysis1$KLB019F[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB019F)


############ KLB027A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027A), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB027A
markov.analysis1$KLB027A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027A)





############ KLB001a
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB001a), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB001a
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB001a))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB001a
markov.analysis1$KLB001a[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB001a)



############ KLB035A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB035A), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB035A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB035A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB035A
markov.analysis1$KLB035A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB035A)


############ KLB035B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB035B), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB035B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB035B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB035B
markov.analysis1$KLB035B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB035B)



############ KLB035C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB035C), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB035C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB035C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB035C
markov.analysis1$KLB035C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB035C)




############ KLB035D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB035D), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB035D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB035D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB035D
markov.analysis1$KLB035D[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB035D)




############ KLB035F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB035F), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB035F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB035F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB035F
markov.analysis1$KLB035F[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB035F)


############ KLB035G
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB035G), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB035G
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB035G))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB035G
markov.analysis1$KLB035G[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB035G)



############ KLB041B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB041B), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB041B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB041B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB041B
markov.analysis1$KLB041B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB041B)




############ KLB027N
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027N), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027N
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027N))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB027N
markov.analysis1$KLB027N[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027N)


############ KLB019B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB019B), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB019B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB019B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB019B
markov.analysis1$KLB019B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB019B)




############ KLB027J
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027J), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027J
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027J))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB027J
markov.analysis1$KLB027J[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027J)


############ KLB019G
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB019G), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB019G
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB019G))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB019G
markov.analysis1$KLB019G[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB019G)




############ KLB027B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027B), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB027B
markov.analysis1$KLB027B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027B)





############ KLB001b
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB001b), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB001b
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB001b))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB001b
markov.analysis1$KLB001b[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB001b)



############ KLB027J
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027J), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027J
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027J))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB027J
markov.analysis1$KLB027J[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027J)




############ KLB041C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB041C), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB041C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB041C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB041C
markov.analysis1$KLB041C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB041C)




############ KLB019C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB019C), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB019C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB019C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB019C
markov.analysis1$KLB019C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB019C)


############ KLB027K
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027K), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027K
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027K))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB027K
markov.analysis1$KLB027K[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027K)


############ KLB019H
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB019H), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB019H
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB019H))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB019H
markov.analysis1$KLB019H[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB019H)




############ KLB027C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027C), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB027C
markov.analysis1$KLB027C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027C)



############ KLB001c
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB001c), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB001c
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB001c))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB001c
markov.analysis1$KLB001c[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB001c)



############ KLB041D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB041D), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB041D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB041D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB041D
markov.analysis1$KLB041D[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB041D)



############ KLB019D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB019D), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB019D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB019D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB019D
markov.analysis1$KLB019D[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB019D)



############ KLB027L
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027L), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027L
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027L))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB027L
markov.analysis1$KLB027L[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027L)



############ KLB019I
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB019I), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB019I
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB019I))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB019I
markov.analysis1$KLB019I[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB019I)



############ KLB027D
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027D), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027D
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027D))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB027D
markov.analysis1$KLB027D[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027D)


############ KLB041E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB041E), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB041E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB041E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB041E
markov.analysis1$KLB041E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB041E)



############ KLB019E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB019E), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB019E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB019E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB019E
markov.analysis1$KLB019E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB019E)


############ KLB027M
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027M), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027M
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027M))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB027M
markov.analysis1$KLB027M[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027M)



############# Add imputed values to KLB019E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB019E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB019E
markov.analysis1$KLB019E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB019E)


############ KLB027M
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027M), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027M
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027M))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB027M
markov.analysis1$KLB027M[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027M)




############ KLB019J
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB019J), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB019J
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB019J))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB019J
markov.analysis1$KLB019J[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB019J)



############ KLB027E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027E), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB027E
markov.analysis1$KLB027E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027E)



############ KLB001e
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB001e), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB001e
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB001e))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB001e
markov.analysis1$KLB001e[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB001e)





############ KLB027M
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027M), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027M
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027M))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB027M
markov.analysis1$KLB027M[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027M)




############ KLB035E
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB035E), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB035E
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB035E))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB035E
markov.analysis1$KLB035E[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB035E)




############ KLB019K
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB019K), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB019K
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB019K))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB019K
markov.analysis1$KLB019K[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB019K)





############ KLB027F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027F), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB027F
markov.analysis1$KLB027F[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027F)



############ KLB001f
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB001f), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB001f
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB001f))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB001f
markov.analysis1$KLB001f[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB001f)




############ KLB001g
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB001g), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB001g
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB001g))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB001g
markov.analysis1$KLB001g[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB001g)


############ KLB027J.1
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027J.1), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027J.1
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027J.1))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB027J.1
markov.analysis1$KLB027J.1[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027J.1)



############ KLB027M.1
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB027M.1), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB027M.1
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB027M.1))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB027M.1
markov.analysis1$KLB027M.1[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB027M.1)


############ KLB001h
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB001h), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB001h
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB001h))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB001h
markov.analysis1$KLB001h[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB001h)


############ KLB001A
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB001A), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB001A
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB001A))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB001A
markov.analysis1$KLB001A[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB001A)



############ KLB001B
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB001B), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB001B
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB001B))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB001B
markov.analysis1$KLB001B[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB001B)



############ KLB001C
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB001C), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB001C
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB001C))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB001C
markov.analysis1$KLB001C[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB001C)



############ KLB001F
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB001F), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB001F
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB001F))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB001F
markov.analysis1$KLB001F[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB001F)




############ KLB001G
#### Number of imputed data sets
M <- 5

#  Import harmonised dataset as a matrix
analysis1.mat <- matrix(c(analysis1$sex, analysis1$KLB001G), ncol=2, byrow=FALSE)

N <- nrow(analysis1.mat)

nvar <- ncol(analysis1.mat)

#  Run prelim.norm() to set up the data for the algorithm
prelim.mi <- prelim.norm(analysis1.mat)


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
  this.imp <- analysis1.mat
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


############# Add imputed values to KLB001G
# Create vector for output values
vec <- rep(TRUE, 6529)

# Record index location of missing values
ind <- which(is.na(markov.analysis1$KLB001G))

# Expand vector to match the length of analysis rows
vec.ind <- rep(NA, length(vec))
vec.ind[ind] <- imp.dat$ymis

# Add imputed values to KLB001G
markov.analysis1$KLB001G[ind] <- imp.dat$ymis

# Check to see if all NAs have been replaced
summary(markov.analysis1$KLB001G)


