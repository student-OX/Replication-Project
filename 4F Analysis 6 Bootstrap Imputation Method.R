############ DATA FRAME FOR IMPUTED DATA
bootstrap.analysis6 <- analysis6

############ PREDICTOR VECTORS FOR IMPUTED VALUES (SEX)
x <- matrix(c(bootstrap.analysis6$sex), ncol=1, byrow=FALSE)


############ PLB018A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB018A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB018A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB018A, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB018A[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB018A)



############ PLB001I 
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB001I))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB001I))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB001I, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB001I[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB001I)


############ PLB033A  
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB033A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB033A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB033A, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB033A[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB033A)



############ PLB033C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB033C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB033C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB033C, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB033C[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB033C)


############ PLB033D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB033D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB033D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB033D, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB033D[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB033D)


############ PLB033F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB033F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB033F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB033F, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB033F[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB033F)


############ PLB033G
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB033G))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB033G))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB033G, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB033G[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB033G)


############ PLB031L
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB031L))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB031L))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB031L, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB031L[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB031L)


############ PLB026R
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB026R))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB026R))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB026R, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB026R[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB026R)


############ PLB018H
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB018H))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB018H))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB018H, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB018H[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB018H)


############ PLB018B 
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB018B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB018B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB018B, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB018B[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB018B)


############ PLB001R
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB001R))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB001R))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB001R, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB001R[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB001R)


############ PLB026J
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB026J))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB026J))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB026J, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB026J[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB026J)



############ PLB018C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB018C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB018C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB018C, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB018C[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB018C)



############ PLB026K
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB026K))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB026K))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB026K, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB026K[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB026K)



############ PLB018D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB018D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB018D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB018D, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB018D[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB018D)



############ PLB026X
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB026X))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB026X))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB026X, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB026X[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB026X)



############ PLB018E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB018E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB018E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB018E, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB018E[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB018E)


############ PLB002C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB002C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB002C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB002C, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB002C[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB002C)


############ PLB033E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB033E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB033E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB033E, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB033E[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB033E)


############ PLB018F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB018F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB018F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB018F, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB018F[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB018F)


############ PLB001N
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB001N))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB001N))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB001N, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB001N[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB001N)


############ PLB012B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis6$PLB012B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis6$PLB012B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis6$PLB012B, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis6$PLB012B[ind] <- imp

## View imputation results
summary(bootstrap.analysis6$PLB012B)



######### REMOVE IRRELEVANT VARIABLES
bootstrap.analysis6 <- bootstrap.analysis6[-c(18)]


