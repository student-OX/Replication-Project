### LOAD PACKAGES
library(mice)
library(nlme)


############ DATA FRAME FOR IMPUTED DATA
bootstrap.analysis2 <- analysis2

############ PREDICTOR VECTORS FOR IMPUTED VALUES (SEX)
x <- matrix(c(bootstrap.analysis2$sex), ncol=1, byrow=FALSE)


############ LLB041A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB041A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB041A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB041A, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB041A[ind] <- imp


############ LLB019A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB019A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB019A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB019A, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB019A[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB019A)


############ LLB027I
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB027I))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB027I))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB027I, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB027I[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB027I)



############ LLB019F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB019F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB019F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB019F, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB019F[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB019F)



############ LLB027A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB027A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB027A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB027A, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB027A[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB027A)



############ LLB001A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB001A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB001A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB001A, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB001A[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB001A)



############ LLB035A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB035A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB035A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB035A, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB035A[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB035A)



############ LLB035B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB035B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB035B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB035B, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB035B[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB035B)



############ LLB035C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB035C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB035C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB035C, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB035C[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB035C)


############ LLB035D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB035D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB035D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB035D, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB035D[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB035D)



############ LLB035F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB035F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB035F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB035F, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB035F[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB035F)


############ LLB035G
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB035G))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB035G))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB035G, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB035G[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB035G)



############ LLB041B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB041B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB041B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB041B, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB041B[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB041B)



############ LLB027N
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB027N))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB027N))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB027N, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB027N[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB027N)



############ LLB019B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB019B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB019B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB019B, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB019B[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB019B)


############ LLB027J
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB027J))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB027J))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB027J, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB027J[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB027J)



############ LLB019G
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB019G))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB019G))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB019G, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB019G[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB019G)


############ LLB027B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB027B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB027B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB027B, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB027B[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB027B)



############ LLB001B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB001B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB001B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB001B, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB001B[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB001B)



############ LLB027J
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB027J))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB027J))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB027J, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB027J[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB027J)


############ LLB041C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB041C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB041C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB041C, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB041C[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB041C)


############ LLB019C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB019C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB019C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB019C, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB019C[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB019C)


############ LLB027K
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB027K))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB027K))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB027K, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB027K[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB027K)


############ LLB019H
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB019H))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB019H))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB019H, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB019H[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB019H)


############ LLB027C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB027C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB027C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB027C, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB027C[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB027C)


############ LLB001C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB001C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB001C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB001C, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB001C[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB001C)


############ LLB041D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB041D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB041D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB041D, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB041D[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB041D)


############ LLB019D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB019D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB019D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB019D, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB019D[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB019D)


############ LLB027L
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB027L))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB027L))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB027L, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB027L[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB027L)


############ LLB019I
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB019I))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB019I))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB019I, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB019I[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB019I)


############ LLB027D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB027D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB027D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB027D, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB027D[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB027D)


############ LLB041E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB041E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB041E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB041E, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB041E[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB041E)


############ LLB019E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB019E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB019E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB019E, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB019E[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB019E)


############ LLB027M
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB027M))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB027M))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB027M, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB027M[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB027M)



############ LLB019J
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB019J))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB019J))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB019J, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB019J[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB019J)


############ LLB027E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB027E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB027E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB027E, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB027E[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB027E)


############ LLB001e
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB001e))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB001e))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB001e, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB001e[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB001e)


############ LLB027M
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB027M))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB027M))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB027M, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB027M[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB027M)


############ LLB035E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB035E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB035E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB035E, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB035E[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB035E)


############ LLB019K
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB019K))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB019K))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB019K, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB019K[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB019K)


############ LLB027F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB027F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB027F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB027F, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB027F[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB027F)


############ LLB001F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB001F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB001F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB001F, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB001F[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB001F)


############ LLB001G
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB001G))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB001G))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB001G, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB001G[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB001G)



############ LLB001a
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB001a))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB001a))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB001a, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB001a[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB001a)


############### LLB001b 
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB001b))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB001b))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB001b, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB001b[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB001b)


############### LLB027J.1  
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB027J.1))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB027J.1))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB027J.1, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB027J.1[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB027J.1)



############### LLB001c   
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB001c))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB001c))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB001c, ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB001c[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB001c)



############### LLB027M.1   
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB027M.1))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB027M.1))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB027M.1 , ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB027M.1[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB027M.1)


############### LLB001f   
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB001f))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB001f))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB001f , ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB001f[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB001f)



############### LLB001g   
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis2$LLB001g))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis2$LLB001g))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis2$LLB001g , ry, x)


## Add imputed values to bootstrap.analysis2 data frame
bootstrap.analysis2$LLB001g[ind] <- imp


# View imputed data
summary(bootstrap.analysis2$LLB001g)

