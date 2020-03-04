############ DATA FRAME FOR IMPUTED DATA
bootstrap.analysis4 <- analysis4

############ PREDICTOR VECTORS FOR IMPUTED VALUES (SEX)
x <- matrix(c(bootstrap.analysis4$sex), ncol=1, byrow=FALSE)


############ NLB041A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB041A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB041A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB041A, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB041A[ind] <- imp

## View imputation results
summary(bootstrap.analysis4$NLB041A)


############ NLB019A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB019A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB019A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB019A, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB019A[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB019A)


############ NLB027I
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB027I))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB027I))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB027I, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB027I[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB027I)



############ NLB019F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB019F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB019F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB019F, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB019F[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB019F)



############ NLB027A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB027A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB027A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB027A, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB027A[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB027A)



############ NLB001A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB001A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB001A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB001A, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB001A[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB001A)



############ NLB035A
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB035A))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB035A))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB035A, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB035A[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB035A)



############ NLB035B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB035B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB035B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB035B, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB035B[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB035B)



############ NLB035C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB035C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB035C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB035C, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB035C[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB035C)


############ NLB035D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB035D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB035D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB035D, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB035D[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB035D)



############ NLB035F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB035F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB035F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB035F, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB035F[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB035F)


############ NLB035G
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB035G))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB035G))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB035G, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB035G[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB035G)



############ NLB041B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB041B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB041B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB041B, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB041B[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB041B)



############ NLB027N
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB027N))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB027N))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB027N, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB027N[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB027N)



############ NLB019B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB019B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB019B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB019B, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB019B[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB019B)


############ NLB027J
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB027J))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB027J))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB027J, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB027J[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB027J)



############ NLB019G
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB019G))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB019G))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB019G, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB019G[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB019G)


############ NLB027B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB027B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB027B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB027B, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB027B[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB027B)



############ NLB001B
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB001B))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB001B))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB001B, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB001B[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB001B)



############ NLB027J
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB027J))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB027J))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB027J, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB027J[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB027J)


############ NLB041C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB041C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB041C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB041C, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB041C[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB041C)


############ NLB019C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB019C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB019C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB019C, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB019C[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB019C)


############ NLB027K
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB027K))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB027K))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB027K, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB027K[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB027K)


############ NLB019H
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB019H))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB019H))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB019H, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB019H[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB019H)


############ NLB027C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB027C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB027C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB027C, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB027C[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB027C)


############ NLB001C
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB001C))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB001C))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB001C, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB001C[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB001C)


############ NLB041D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB041D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB041D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB041D, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB041D[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB041D)


############ NLB019D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB019D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB019D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB019D, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB019D[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB019D)


############ NLB027L
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB027L))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB027L))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB027L, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB027L[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB027L)


############ NLB019I
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB019I))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB019I))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB019I, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB019I[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB019I)


############ NLB027D
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB027D))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB027D))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB027D, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB027D[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB027D)


############ NLB041E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB041E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB041E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB041E, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB041E[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB041E)


############ NLB019E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB019E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB019E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB019E, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB019E[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB019E)


############ NLB027M
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB027M))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB027M))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB027M, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB027M[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB027M)



############ NLB019J
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB019J))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB019J))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB019J, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB019J[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB019J)


############ NLB027E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB027E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB027E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB027E, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB027E[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB027E)


############ NLB001e
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB001e))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB001e))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB001e, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB001e[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB001e)


############ NLB027M
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB027M))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB027M))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB027M, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB027M[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB027M)


############ NLB035E
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB035E))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB035E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB035E, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB035E[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB035E)


############ NLB019K
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB019K))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB019K))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB019K, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB019K[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB019K)


############ NLB027F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB027F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB027F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB027F, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB027F[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB027F)


############ NLB001F
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB001F))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB001F))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB001F, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB001F[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB001F)


############ NLB001G
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB001G))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB001G))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB001G, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB001G[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB001G)



############ NLB027J.1 
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB027J.1))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB027J.1))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB027J.1, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB027J.1[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB027J.1)



############ NLB001E  
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB001E ))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB001E))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB001E, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB001E[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB001E)


############ NLB027M.1  
## Logical vector of length length(y) indicating the the subset y[ry] of elements in y to which the imputation model is fitted
ry <- rep(TRUE, length(bootstrap.analysis4$NLB027M.1))

## Index location of missing values
ind <- which(is.na(bootstrap.analysis4$NLB027M.1))

## Mark missing values as FALSE in logical vector
ry[ind] <- FALSE


## Impute missing data
imp <- mice.impute.norm.boot(bootstrap.analysis4$NLB027M.1, ry, x)


## Add imputed values to bootstrap.analysis4 data frame
bootstrap.analysis4$NLB027M.1[ind] <- imp


# View imputed data
summary(bootstrap.analysis4$NLB027M.1)

