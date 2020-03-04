############ DATA FRAME FOR IMPUTED DATA
bootstrap.analysis5 <- analysis5

############ PREDICTOR VECTORS FOR IMPUTED VALUES (SEX)
x <- matrix(c(bootstrap.analysis5$sex), ncol=1, byrow=FALSE)


############ OLB026S
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB026S))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB026S))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB026S, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB026S[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB026S)



############ OLB018A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB018A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB018A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB018A, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB018A[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB018A)



############ OLB001I 
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB001I))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB001I))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB001I, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB001I[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB001I)


############ OLB033A  
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB033A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB033A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB033A, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB033A[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB033A)


############ OLB033B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB033B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB033B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB033B, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB033B[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB033B)


############ OLB033C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB033C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB033C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB033C, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB033C[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB033C)


############ OLB033D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB033D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB033D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB033D, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB033D[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB033D)


############ OLB033F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB033F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB033F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB033F, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB033F[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB033F)


############ OLB033G
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB033G))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB033G))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB033G, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB033G[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB033G)


############ OLB031L
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB031L))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB031L))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB031L, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB031L[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB031L)


############ OLB026R
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB026R))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB026R))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB026R, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB026R[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB026R)


############ OLB018H
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB018H))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB018H))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB018H, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB018H[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB018H)


############ OLB018B 
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB018B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB018B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB018B, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB018B[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB018B)


############ OLB001R
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB001R))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB001R))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB001R, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB001R[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB001R)


############ OLB026J
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB026J))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB026J))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB026J, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB026J[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB026J)



############ OLB018C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB018C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB018C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB018C, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB018C[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB018C)



############ OLB026K
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB026K))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB026K))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB026K, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB026K[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB026K)



############ OLB018D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB018D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB018D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB018D, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB018D[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB018D)



############ OLB026X
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB026X))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB026X))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB026X, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB026X[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB026X)



############ OLB018E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB018E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB018E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB018E, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB018E[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB018E)


############ OLB002C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB002C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB002C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB002C, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB002C[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB002C)


############ OLB033E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB033E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB033E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB033E, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB033E[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB033E)


############ OLB018F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB018F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB018F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB018F, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB018F[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB018F)


############ OLB001N
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB001N))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB001N))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB001N, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB001N[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB001N)


############ OLB012B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis5$OLB012B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis5$OLB012B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis5$OLB012B, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis5$OLB012B[ind] <- imp

## View imputation results
summary(bootstrap.analysis5$OLB012B)


