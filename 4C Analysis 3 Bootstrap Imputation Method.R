### LOAD PACKAGES
library(mice)
library(nlme)


############ DATA FRAME FOR IMPUTED DATA
bootstrap.analysis3 <- analysis3

############ PREDICTOR VECTORS FOR IMPUTED VALUES (SEX)
x <- matrix(c(bootstrap.analysis3$sex), ncol=1, byrow=FALSE)


############ MLB041A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB041A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB041A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB041A, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB041A[ind] <- imp


############ MLB019A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB019A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB019A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB019A, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB019A[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB019A)


############ MLB027I
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB027I))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB027I))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB027I, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB027I[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB027I)



############ MLB019F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB019F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB019F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB019F, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB019F[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB019F)



############ MLB027A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB027A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB027A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB027A, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB027A[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB027A)



############ MLB001A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB001A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB001A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB001A, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB001A[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB001A)



############ MLB035A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB035A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB035A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB035A, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB035A[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB035A)



############ MLB035B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB035B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB035B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB035B, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB035B[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB035B)



############ MLB035C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB035C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB035C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB035C, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB035C[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB035C)


############ MLB035D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB035D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB035D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB035D, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB035D[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB035D)



############ MLB035F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB035F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB035F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB035F, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB035F[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB035F)


############ MLB035G
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB035G))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB035G))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB035G, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB035G[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB035G)



############ MLB041B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB041B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB041B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB041B, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB041B[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB041B)



############ MLB027N
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB027N))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB027N))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB027N, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB027N[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB027N)



############ MLB019B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB019B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB019B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB019B, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB019B[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB019B)


############ MLB027J
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB027J))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB027J))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB027J, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB027J[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB027J)



############ MLB019G
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB019G))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB019G))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB019G, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB019G[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB019G)


############ MLB027B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB027B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB027B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB027B, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB027B[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB027B)



############ MLB001B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB001B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB001B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB001B, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB001B[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB001B)



############ MLB027J
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB027J))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB027J))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB027J, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB027J[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB027J)


############ MLB041C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB041C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB041C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB041C, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB041C[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB041C)


############ MLB019C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB019C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB019C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB019C, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB019C[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB019C)


############ MLB027K
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB027K))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB027K))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB027K, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB027K[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB027K)


############ MLB019H
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB019H))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB019H))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB019H, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB019H[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB019H)


############ MLB027C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB027C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB027C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB027C, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB027C[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB027C)


############ MLB001C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB001C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB001C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB001C, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB001C[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB001C)


############ MLB041D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB041D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB041D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB041D, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB041D[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB041D)


############ MLB019D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB019D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB019D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB019D, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB019D[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB019D)


############ MLB027L
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB027L))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB027L))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB027L, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB027L[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB027L)


############ MLB019I
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB019I))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB019I))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB019I, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB019I[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB019I)


############ MLB027D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB027D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB027D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB027D, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB027D[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB027D)


############ MLB041E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB041E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB041E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB041E, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB041E[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB041E)


############ MLB019E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB019E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB019E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB019E, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB019E[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB019E)


############ MLB027M
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB027M))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB027M))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB027M, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB027M[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB027M)



############ MLB019J
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB019J))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB019J))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB019J, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB019J[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB019J)


############ MLB027E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB027E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB027E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB027E, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB027E[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB027E)


############ MLB001e
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB001e))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB001e))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB001e, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB001e[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB001e)


############ MLB027M
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB027M))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB027M))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB027M, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB027M[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB027M)


############ MLB035E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB035E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB035E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB035E, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB035E[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB035E)


############ MLB019K
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB019K))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB019K))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB019K, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB019K[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB019K)


############ MLB027F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB027F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB027F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB027F, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB027F[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB027F)


############ MLB001F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB001F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB001F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB001F, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB001F[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB001F)


############ MLB001G
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB001G))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB001G))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB001G, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB001G[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB001G)




############ MLB027J.1 
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB027J.1))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB027J.1))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB027J.1, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB027J.1[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB027J.1)


############ MLB001E 
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB001E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB001E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB001E, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB001E[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB001E)


############ MLB027M.1 
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis3$MLB027M.1))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis3$MLB027M.1))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis3$MLB027M.1, ry, x)


## Add imputed values to bootstrap.analysis3 data frame
bootstrap.analysis3$MLB027M.1[ind] <- imp


# View imputed data
summary(bootstrap.analysis3$MLB027M.1)



