### LOAD PACKAGES
library(mice)
library(nlme)

###Warning message: In bootstrap.analysis1$KLB001F[ind] <- imp :number of items to replace is not a multiple of replacement length  KLB035F KLB035G KLB027J.1 KLB019C KLB027M.1 KLB027F KLB001h KLB001B KLB001C KLB001F KLB001G  

############ DATA FRAME FOR IMPUTED DATA
bootstrap.analysis1 <- analysis1

############ PREDICTOR VECTORS FOR IMPUTED VALUES (SEX)
x <- matrix(c(bootstrap.analysis1$sex), ncol=1, byrow=FALSE)


############ KLB041A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB041A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB041A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB041A, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB041A[ind] <- imp


############ KLB019A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB019A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB019A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB019A, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB019A[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB019A)


############ KLB027I
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB027I))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB027I))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB027I, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB027I[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB027I)



############ KLB019F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB019F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB019F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB019F, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB019F[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB019F)



############ KLB027A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB027A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB027A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB027A, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB027A[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB027A)



############ KLB001a
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB001a))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB001a))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB001a, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB001a[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB001a)



############ KLB035A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB035A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB035A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB035A, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB035A[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB035A)



############ KLB035B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB035B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB035B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB035B, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB035B[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB035B)



############ KLB035C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB035C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB035C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB035C, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB035C[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB035C)


############ KLB035D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB035D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB035D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB035D, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB035D[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB035D)



############ KLB035F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB035F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB035F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB035F, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB035F[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB035F)


############ KLB035G
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB035G))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB035G))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB035G, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB035G[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB035G)



############ KLB041B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB041B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB041B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB041B, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB041B[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB041B)



############ KLB027N
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB027N))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB027N))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB027N, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB027N[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB027N)



############ KLB019B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB019B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB019B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB019B, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB019B[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB019B)


############ KLB027J.1
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB027J.1))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB027J.1))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB027J.1, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB027J.1[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB027J.1)



############ KLB019G
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB019G))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB019G))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB019G, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB019G[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB019G)


############ KLB027B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB027B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB027B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB027B, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB027B[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB027B)



############ KLB001b
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB001b))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB001b))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB001b, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB001b[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB001b)



############ KLB027J
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB027J))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB027J))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB027J, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB027J[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB027J)


############ KLB041C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB041C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB041C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB041C, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB041C[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB041C)


############ KLB019C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB019C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB019C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB019C, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB019C[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB019C)


############ KLB027K
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB027K))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB027K))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB027K, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB027K[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB027K)


############ KLB019H
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB019H))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB019H))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB019H, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB019H[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB019H)


############ KLB027C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB027C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB027C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB027C, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB027C[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB027C)


############ KLB001c
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB001c))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB001c))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB001c, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB001c[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB001c)


############ KLB041D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB041D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB041D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB041D, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB041D[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB041D)


############ KLB019D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB019D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB019D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB019D, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB019D[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB019D)


############ KLB027L
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB027L))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB027L))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB027L, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB027L[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB027L)


############ KLB019I
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB019I))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB019I))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB019I, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB019I[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB019I)


############ KLB027D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB027D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB027D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB027D, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB027D[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB027D)


############ KLB041E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB041E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB041E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB041E, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB041E[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB041E)


############ KLB019E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB019E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB019E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB019E, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB019E[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB019E)


############ KLB027M.1
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB027M.1))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB027M.1))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB027M.1, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB027M.1[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB027M.1)



############ KLB019J
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB019J))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB019J))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB019J, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB019J[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB019J)


############ KLB027E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB027E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB027E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB027E, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB027E[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB027E)


############ KLB001e
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB001e))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB001e))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB001e, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB001e[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB001e)


############ KLB027M
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB027M))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB027M))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB027M, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB027M[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB027M)


############ KLB035E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB035E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB035E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB035E, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB035E[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB035E)


############ KLB019K
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB019K))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB019K))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB019K, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB019K[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB019K)


############ KLB027F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB027F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB027F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB027F, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB027F[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB027F)


############ KLB001F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB001F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB001F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB001F, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB001F[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB001F)


############ KLB001G
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB001G))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB001G))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB001G, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB001G[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB001G)



############ KLB001h 
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB001h ))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB001h))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB001h, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB001h [ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB001h)



############ KLB001B 
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis1$KLB001B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis1$KLB001B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis1$KLB001B, ry, x)


## Add imputed values to bootstrap.analysis1 data frame
bootstrap.analysis1$KLB001B[ind] <- imp


# View imputed data
summary(bootstrap.analysis1$KLB001B)


