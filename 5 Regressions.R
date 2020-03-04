########### CORRELATION MATRICES FOR PSYCHO-SOCIAL VARIABLES
##### Non-Imputed data
psyc.nonImp <- data.frame(analysis1[11:54], analysis1[97:101])

# Set NAs to 0
psyc.nonImp[is.na(psyc.nonImp)] <- 0

# Correlation matrix
corNonImp <- cor(psyc.nonImp)

##### Markov Chain Imputed Datasets
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
corMarkov.analysis1 <- cor(markov.analysis1[11:54])


#### 2008 baseline covariates, 2010 and 2012 stroke incidence
corMarkov.analysis2 <- data.frame(markov.analysis2[13:55], markov.analysis2[98:102])

corMarkov.analysis2 <- cor(corMarkov.analysis2)



#### 2010 baseline covariates, 2012 and 2014 stroke incidence
corMarkov.analysis3 <- cor(markov.analysis3[12:54])


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
corMarkov.analysis4 <- cor(markov.analysis4[11:53])


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
corMarkov.analysis5 <- cor(markov.analysis5[5:29])


#### 2016 baseline covariates, 2018 stroke incidence
corMarkov.analysis6 <- cor(markov.analysis6[5:17], markov.analysis6[19:28])
corMarkov.analysis6 <- cor(corMarkov.analysis6)




##### Bootstrap Imputed Datasets
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
corBoot.analysis1 <- cor(bootstrap.analysis1[11:54])


#### 2008 baseline covariates, 2010 and 2012 stroke incidence
corBoot.analysis2 <- data.frame(bootstrap.analysis2[13:55], bootstrap.analysis2[98:102])

corBoot.analysis2 <- cor(corBoot.analysis2)


#### 2010 baseline covariates, 2012 and 2014 stroke incidence
corBoot.analysis3 <- cor(bootstrap.analysis3[12:54])


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
corBoot.analysis4 <- cor(bootstrap.analysis4[11:53])


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
corBoot.analysis5 <- cor(bootstrap.analysis5[5:29])


#### 2016 baseline covariates, 2018 stroke incidence
corBoot.analysis6 <- cor(bootstrap.analysis6[5:27])



######### PERFORM KAISER-MEYER-OLKIN (KMO) TEST
library(psych)
library(GPArotation)

KMO(r=corMarkov.analysis1)

KMO(r=corMarkov.analysis2)

KMO(r=corMarkov.analysis3)

KMO(r=corMarkov.analysis4)

KMO(r=corMarkov.analysis5)

KMO(r=corMarkov.analysis6)



######### DETERMINE NUMBER OF FACTORS

## 2006 baseline covariates, 2008 and 2010 stroke incidence: Five factors
Fig6a <- fa.parallel(corMarkov.analysis1, fm = "minres", fa = "fa")



######### FACTOR ANALYSIS OF PSYCOLOGICAL VARIABLES
## 2006 baseline covariates, 2008 and 2010 stroke incidence: 2 factors
## Non imputed data
fa.analysis1 <-  fa(r=corNonImp, nfactors = 5, rotate="none",fm="pa")

fa.diagram(fa.analysis1, rsize=5)


# Markov imputed data
fa.markov.analysis1 <-  fa(r=corMarkov.analysis1, nfactors = 5, rotate="none",fm="pa")

print(fa.markov.analysis1)

fa.diagram(fa.markov.analysis1)


# Bootsrap imputed data
fa.boot.analysis1 <-  fa(r=corBoot.analysis1, nfactors = 5, rotate="none",fm="pa")

print(fa.boot.analysis1)

fa.diagram(fa.boot.analysis1)



#### 2008 baseline covariates, 2010 and 2012 stroke incidence: 5 factors
fa.markov.analysis2 <-  fa(r=corMarkov.analysis2, nfactors = 5, rotate="none",fm="pa")

print(fa.markov.analysis2)

fa.diagram(fa.markov.analysis2)



#### 2010 baseline covariates, 2012 and 2014 stroke incidence: 3 factors
fa.markov.analysis3 <-  fa(r=corMarkov.analysis3, nfactors = 3, rotate="none",fm="pa")

print(fa.markov.analysis3)

fa.diagram(fa.markov.analysis3)



#### 2012 baseline covariates, 2014 and 2016 stroke incidence: 3 factors
fa.markov.analysis4 <-  fa(r=corMarkov.analysis4, nfactors = 3, rotate="none",fm="pa")

print(fa.markov.analysis4)

fa.diagram(fa.markov.analysis4)


#### 2014 baseline covariates, 2016 and 2018 stroke incidence: 5 factors
fa.markov.analysis5 <-  fa(r=corMarkov.analysis5, nfactors = 5, rotate="none",fm="pa")

print(fa.markov.analysis5)

fa.diagram(fa.markov.analysis5)


#### 2016 baseline covariates, 2018 stroke incidence: 2 factors
fa.markov.analysis6 <-  fa(r=corMarkov.analysis6, nfactors = 2, rotate="none",fm="pa")

print(fa.markov.analysis6)

fa.diagram(fa.markov.analysis6)




####### LOGISTIC REGRESSIONS
library(logistf)
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
## 1.1 Without imputation, choosing more appropriate variables and addressing variable transformation issues
mod1.1 <- logistf(stroke ~ KIMpurposeinlife + age  + sex + marstat + race + education + wealth + func014 + func016 + func021 + func023 + func025 + numCigs + mildexercise + modexercise + highexercise + alcohol  + diabetes + bmi + heartdis + KLB041A	+ KLB019A	+ KLB027I	+ KLB019F	+ KLB027A	+ KLB001a	+ KLB041B	+ KLB027N	+ KLB019B	+ KLB027J	+ KLB019G	+ KLB027B	+ KLB001b	+ KLB027J	+ KLB041C	+ KLB019C	+ KLB027K	+ KLB019H	+ KLB027C	+ KLB001c	+ KLB041D	+ KLB019D	+ KLB027L	+ KLB019I	+ KLB027D	+ KLB041E	+ KLB019E	+ KLB027M	+ KLB019J	+ KLB027E	+ KLB001e	+ KLB027M	+ KLB019K	+ KLB027F	+ KLB001f	+ KLB001g, data = analysis1)



## 1.2 Without imputation, choosing more appropriate variables and addressing variable transformation issues and factor analysed 
# Get factor scores
scores.analysis1 <-factor.scores(corNonImp, fa.analysis1)

# View factor weights
scores.analysis1$weights


# Create interaction terms
analysis1$KLB041A  <-  analysis1$KLB041A*0.022099254
analysis1$KLB019A  <-  analysis1$KLB019A*0.020915051
analysis1$KLB027I  <-  analysis1$KLB027I*0.060887889
analysis1$KLB019F  <-  analysis1$KLB019F*0.02177315
analysis1$KLB027A  <-  analysis1$KLB027A*0.062960155
analysis1$KLB001a  <-  analysis1$KLB001a*0.124778728
analysis1$KLB035A  <-  analysis1$KLB035A*0.037686744
analysis1$KLB035B  <-  analysis1$KLB035B*0.024947812
analysis1$KLB035C  <-  analysis1$KLB035C*0.037159444
analysis1$KLB035D  <-  analysis1$KLB035D*0.023397105
analysis1$KLB035F  <-  analysis1$KLB035F*0.003075364
analysis1$KLB035G  <-  analysis1$KLB035G*0.040021234
analysis1$KLB041B  <-  analysis1$KLB041B*0.027918603
analysis1$KLB027N  <-  analysis1$KLB027N*0.044883651
analysis1$KLB019B  <-  analysis1$KLB019B*0.02670577
analysis1$KLB027J  <-  analysis1$KLB027J*0.053061674
analysis1$KLB019G  <-  analysis1$KLB019G*0.029066359
analysis1$KLB027B  <-  analysis1$KLB027B*0.065114858
analysis1$KLB001b  <-  analysis1$KLB001b*-0.025838499
analysis1$KLB027J.1  <-  analysis1$KLB027J.1*0.05697613
analysis1$KLB041C  <-  analysis1$KLB041C*0.018300887
analysis1$KLB019C  <-  analysis1$KLB019C*0.012413776
analysis1$KLB027K  <-  analysis1$KLB027K*0.033007901
analysis1$KLB019H  <-  analysis1$KLB019H*0.032345725
analysis1$KLB027C  <-  analysis1$KLB027C*0.024193204
analysis1$KLB001c  <-  analysis1$KLB001c*-0.066826599
analysis1$KLB041D  <-  analysis1$KLB041D*0.01675385
analysis1$KLB019D  <-  analysis1$KLB019D*0.022026544
analysis1$KLB027L  <-  analysis1$KLB027L*0.009430357
analysis1$KLB019I  <-  analysis1$KLB019I*0.046398451
analysis1$KLB027D  <-  analysis1$KLB027D*0.045208671
analysis1$KLB041E  <-  analysis1$KLB041E*0.021110281
analysis1$KLB019E  <-  analysis1$KLB019E*0.013505573
analysis1$KLB027M  <-  analysis1$KLB027M*0.062841035
analysis1$KLB019J  <-  analysis1$KLB019J*0.026591422
analysis1$KLB027E  <-  analysis1$KLB027E*0.054380244
analysis1$KLB001e  <-  analysis1$KLB001e*0.007402667
analysis1$KLB027M.1  <-  analysis1$KLB027M.1*0.053952592
analysis1$KLB035E  <-  analysis1$KLB035E*0.023826164
analysis1$KLB019K  <-  analysis1$KLB019K*0.028273106
analysis1$KLB027F  <-  analysis1$KLB027F*0.039643467
analysis1$KLB001f  <-  analysis1$KLB001f*0.032006463
analysis1$KLB001g  <-  analysis1$KLB001g*-0.008175871
analysis1$KLB001h  <-  analysis1$KLB001h*0.006872656
analysis1$KLB001A  <-  analysis1$KLB001A*0.142883882
analysis1$KLB001B  <-  analysis1$KLB001B*-0.013202316
analysis1$KLB001C  <-  analysis1$KLB001C*-0.058010098
analysis1$KLB001F  <-  analysis1$KLB001F*0.044236622
analysis1$KLB001G  <-  analysis1$KLB001G*0.007892469


## Regression
mod1.2 <- logistf(stroke ~ KIMpurposeinlife + age  + sex + marstat + race + education + wealth + func014 + func016 + func021 + func023 + func025 + numCigs + mildexercise + modexercise + highexercise + alcohol  + diabetes + bmi + heartdis + KLB041A	+ KLB019A	+ KLB027I	+ KLB019F	+ KLB027A	+ KLB001a	+ KLB041B	+ KLB027N	+ KLB019B	+ KLB027J	+ KLB019G	+ KLB027B	+ KLB001b	+ KLB027J	+ KLB041C	+ KLB019C	+ KLB027K	+ KLB019H	+ KLB027C	+ KLB001c	+ KLB041D	+ KLB019D	+ KLB027L	+ KLB019I	+ KLB027D	+ KLB041E	+ KLB019E	+ KLB027M	+ KLB019J	+ KLB027E	+ KLB001e	+ KLB027M	+ KLB019K	+ KLB027F	+ KLB001f	+ KLB001g, data = analysis1,  family = "binomial")


## 1.3 Markov chain imputation, choosing more appropriate variables and addressing variable transformation issues
mod1.3 <- logistf(stroke ~ KIMpurposeinlife + age  + sex + marstat + race + education + wealth + func014 + func016 + func021 + func023 + func025 + numCigs + mildexercise + modexercise + highexercise + alcohol  + diabetes + bmi + heartdis + KLB041A	+ KLB019A	+ KLB027I	+ KLB019F	+ KLB027A	+ KLB001a	+ KLB041B	+ KLB027N	+ KLB019B	+ KLB027J	+ KLB019G	+ KLB027B	+ KLB001b	+ KLB027J	+ KLB041C	+ KLB019C	+ KLB027K	+ KLB019H	+ KLB027C	+ KLB001c	+ KLB041D	+ KLB019D	+ KLB027L	+ KLB019I	+ KLB027D	+ KLB041E	+ KLB019E	+ KLB027M	+ KLB019J	+ KLB027E	+ KLB001e	+ KLB027M	+ KLB019K	+ KLB027F	+ KLB001f	+ KLB001g, data = markov.analysis1,  family = binomial(logit))


## 1.4 Bootstrap imputation, choosing more appropriate variables and addressing variable transformation issues
mod1.4 <- logistf(stroke ~ KIMpurposeinlife + age  + sex + marstat + race + education + wealth + func014 + func016 + func021 + func023 + func025 + numCigs + mildexercise + modexercise + highexercise + alcohol  + diabetes + bmi + heartdis + KLB041A	+ KLB019A	+ KLB027I	+ KLB019F	+ KLB027A	+ KLB001a	+ KLB041B	+ KLB027N	+ KLB019B	+ KLB027J	+ KLB019G	+ KLB027B	+ KLB001b	+ KLB027J	+ KLB041C	+ KLB019C	+ KLB027K	+ KLB019H	+ KLB027C	+ KLB001c	+ KLB041D	+ KLB019D	+ KLB027L	+ KLB019I	+ KLB027D	+ KLB041E	+ KLB019E	+ KLB027M	+ KLB019J	+ KLB027E	+ KLB001e	+ KLB027M	+ KLB019K	+ KLB027F	+ KLB001f	+ KLB001g, data = bootstrap.analysis1,  family = "binomial")


## 1.5 Markov chain imputation, interaction term for stress using factor analysed variables, more appropriate variables and addressing transformation issues
# Get factor scores
scores.corMarkov.analysis1 <-factor.scores(corMarkov.analysis1, fa.markov.analysis1)

# View factor weights
scores.corMarkov.analysis1$weights

# Create interaction terms
markov.analysis1$KLB041A  <-  markov.analysis1$KLB041A*-0.044349633
markov.analysis1$KLB019A  <-  markov.analysis1$KLB019A*-0.022925923
markov.analysis1$KLB027I  <-  markov.analysis1$KLB027I*0.054913451
markov.analysis1$KLB019F  <-  markov.analysis1$KLB019F*0.048452391
markov.analysis1$KLB027A  <-  markov.analysis1$KLB027A*0.066553387
markov.analysis1$KLB001a  <-  markov.analysis1$KLB001a*-0.00718568
markov.analysis1$KLB035A  <-  markov.analysis1$KLB035A*0.037873051
markov.analysis1$KLB035B  <-  markov.analysis1$KLB035B*0.034422908
markov.analysis1$KLB035C  <-  markov.analysis1$KLB035C*0.03223822
markov.analysis1$KLB035D  <-  markov.analysis1$KLB035D*0.035070877
markov.analysis1$KLB035F  <-  markov.analysis1$KLB035F*-0.025024443
markov.analysis1$KLB035G  <-  markov.analysis1$KLB035G*0.038122571
markov.analysis1$KLB041B  <-  markov.analysis1$KLB041B*-0.066842288
markov.analysis1$KLB027N  <-  markov.analysis1$KLB027N*0.05504158
markov.analysis1$KLB019B  <-  markov.analysis1$KLB019B*-0.02905499
markov.analysis1$KLB027J  <-  markov.analysis1$KLB027J*0.08393913
markov.analysis1$KLB019G  <-  markov.analysis1$KLB019G*0.033189216
markov.analysis1$KLB027B  <-  markov.analysis1$KLB027B*0.075566421
markov.analysis1$KLB001b  <-  markov.analysis1$KLB001b*-0.016806558
markov.analysis1$KLB027J.1  <-  markov.analysis1$KLB027J.1*0.082351395
markov.analysis1$KLB041C  <-  markov.analysis1$KLB041C*-0.049862844
markov.analysis1$KLB019C  <-  markov.analysis1$KLB019C*-0.034510593
markov.analysis1$KLB027K  <-  markov.analysis1$KLB027K*0.043855118
markov.analysis1$KLB019H  <-  markov.analysis1$KLB019H*0.032053185
markov.analysis1$KLB027C  <-  markov.analysis1$KLB027C*0.060425727
markov.analysis1$KLB001c  <-  markov.analysis1$KLB001c*-0.026064935
markov.analysis1$KLB041D  <-  markov.analysis1$KLB041D*-0.038773094
markov.analysis1$KLB019D  <-  markov.analysis1$KLB019D*-0.03129089
markov.analysis1$KLB027L  <-  markov.analysis1$KLB027L*-0.052977168
markov.analysis1$KLB019I  <-  markov.analysis1$KLB019I*0.038122261
markov.analysis1$KLB027D  <-  markov.analysis1$KLB027D*0.05815278
markov.analysis1$KLB041E  <-  markov.analysis1$KLB041E*-0.038291807
markov.analysis1$KLB019E  <-  markov.analysis1$KLB019E*-0.042733709
markov.analysis1$KLB027M  <-  markov.analysis1$KLB027M*0.071262516
markov.analysis1$KLB019J  <-  markov.analysis1$KLB019J*0.051592894
markov.analysis1$KLB027E  <-  markov.analysis1$KLB027E*0.069975427
markov.analysis1$KLB001e  <-  markov.analysis1$KLB001e*-0.020959781
markov.analysis1$KLB027M.1  <-  markov.analysis1$KLB027M.1*0.076229901
markov.analysis1$KLB035E  <-  markov.analysis1$KLB035E*0.039003778
markov.analysis1$KLB019K  <-  markov.analysis1$KLB019K*0.045240567
markov.analysis1$KLB027F  <-  markov.analysis1$KLB027F*0.066423495
markov.analysis1$KLB001f  <-  markov.analysis1$KLB001f*-0.02415495
markov.analysis1$KLB001g  <-  markov.analysis1$KLB001g*-0.01547244
markov.analysis1$KLB001h  <-  markov.analysis1$KLB001h*0.008067607

# Markov Regression model with interaction term for factor loadings
mod1.5 <- logistf(stroke ~ KIMpurposeinlife + age  + sex + marstat + race + education + wealth + func014 + func016 + func021 + func023 + func025 + numCigs + mildexercise + modexercise + highexercise + alcohol  + diabetes + bmi + heartdis + KLB041A	+ KLB019A	+ KLB027I	+ KLB019F	+ KLB027A	+ KLB001a	+ KLB041B	+ KLB027N	+ KLB019B	+ KLB027J	+ KLB019G	+ KLB027B	+ KLB001b	+ KLB027J.1	+ KLB041C	+ KLB019C	+ KLB027K	+ KLB019H	+ KLB027C	+ KLB001c	+ KLB041D	+ KLB019D	+ KLB027L	+ KLB019I	+ KLB027D	+ KLB041E	+ KLB019E	+ KLB027M	+ KLB019J	+ KLB027E	+ KLB001e	+ KLB027M.1	+ KLB019K	+ KLB027F	+ KLB001f	+ KLB001g	+ KLB001h, data = markov.analysis1,  family = "binomial")



## 1.6 Bootstrap imputation, interaction term for stress using factor analysed variables, more appropriate variables and addressing transformation issues
# Get factor scores
scores.corBoot.analysis1 <-factor.scores(corBoot.analysis1, fa.boot.analysis1)

# View factor weights
scores.corBoot.analysis1$weights

# Create interaction terms
bootstrap.analysis1$KLB041A  <-  bootstrap.analysis1$KLB041A*-0.04567426
bootstrap.analysis1$KLB019A  <-  bootstrap.analysis1$KLB019A*-0.021420989
bootstrap.analysis1$KLB027I  <-  bootstrap.analysis1$KLB027I*0.055084555
bootstrap.analysis1$KLB019F  <-  bootstrap.analysis1$KLB019F*0.049133548
bootstrap.analysis1$KLB027A  <-  bootstrap.analysis1$KLB027A*0.069190299
bootstrap.analysis1$KLB001a  <-  bootstrap.analysis1$KLB001a*-0.007131205
bootstrap.analysis1$KLB035A  <-  bootstrap.analysis1$KLB035A*0.036572692
bootstrap.analysis1$KLB035B  <-  bootstrap.analysis1$KLB035B*0.033882544
bootstrap.analysis1$KLB035C  <-  bootstrap.analysis1$KLB035C*0.032178522
bootstrap.analysis1$KLB035D  <-  bootstrap.analysis1$KLB035D*0.036146984
bootstrap.analysis1$KLB035F  <-  bootstrap.analysis1$KLB035F*-0.02574598
bootstrap.analysis1$KLB035G  <-  bootstrap.analysis1$KLB035G*0.036992265
bootstrap.analysis1$KLB041B  <-  bootstrap.analysis1$KLB041B*-0.066314002
bootstrap.analysis1$KLB027N  <-  bootstrap.analysis1$KLB027N*0.055621168
bootstrap.analysis1$KLB019B  <-  bootstrap.analysis1$KLB019B*-0.029360858
bootstrap.analysis1$KLB027J  <-  bootstrap.analysis1$KLB027J*0.084599231
bootstrap.analysis1$KLB019G  <-  bootstrap.analysis1$KLB019G*0.032349009
bootstrap.analysis1$KLB027B  <-  bootstrap.analysis1$KLB027B*0.074461521
bootstrap.analysis1$KLB001b  <-  bootstrap.analysis1$KLB001b*-0.016186133
bootstrap.analysis1$KLB027J.1  <-  bootstrap.analysis1$KLB027J.1*0.07669932
bootstrap.analysis1$KLB041C  <-  bootstrap.analysis1$KLB041C*-0.050374448
bootstrap.analysis1$KLB019C  <-  bootstrap.analysis1$KLB019C*-0.034491549
bootstrap.analysis1$KLB027K  <-  bootstrap.analysis1$KLB027K*0.044442935
bootstrap.analysis1$KLB019H  <-  bootstrap.analysis1$KLB019H*0.032289226
bootstrap.analysis1$KLB027C  <-  bootstrap.analysis1$KLB027C*0.059453092
bootstrap.analysis1$KLB001c  <-  bootstrap.analysis1$KLB001c*-0.025782601
bootstrap.analysis1$KLB041D  <-  bootstrap.analysis1$KLB041D*-0.038306564
bootstrap.analysis1$KLB019D  <-  bootstrap.analysis1$KLB019D*-0.031418005
bootstrap.analysis1$KLB027L  <-  bootstrap.analysis1$KLB027L*-0.052811245
bootstrap.analysis1$KLB019I  <-  bootstrap.analysis1$KLB019I*0.038544076
bootstrap.analysis1$KLB027D  <-  bootstrap.analysis1$KLB027D*0.05880768
bootstrap.analysis1$KLB041E  <-  bootstrap.analysis1$KLB041E*-0.03800501
bootstrap.analysis1$KLB019E  <-  bootstrap.analysis1$KLB019E*-0.043703248
bootstrap.analysis1$KLB027M  <-  bootstrap.analysis1$KLB027M*0.08257601
bootstrap.analysis1$KLB019J  <-  bootstrap.analysis1$KLB019J*0.051110618
bootstrap.analysis1$KLB027E  <-  bootstrap.analysis1$KLB027E*0.067783173
bootstrap.analysis1$KLB001e  <-  bootstrap.analysis1$KLB001e*-0.020767955
bootstrap.analysis1$KLB027M.1  <-  bootstrap.analysis1$KLB027M.1*0.070286009
bootstrap.analysis1$KLB035E  <-  bootstrap.analysis1$KLB035E*0.038580923
bootstrap.analysis1$KLB019K  <-  bootstrap.analysis1$KLB019K*0.044629548
bootstrap.analysis1$KLB027F  <-  bootstrap.analysis1$KLB027F*0.068694183
bootstrap.analysis1$KLB001f  <-  bootstrap.analysis1$KLB001f*-0.023858093
bootstrap.analysis1$KLB001g  <-  bootstrap.analysis1$KLB001g*-0.016324628
bootstrap.analysis1$KLB001h  <-  bootstrap.analysis1$KLB001h*0.009533624


# Regression
mod1.6 <- logistf(stroke ~ KIMpurposeinlife + age  + sex + marstat + race + education + wealth + func014 + func016 + func021 + func023 + func025 + numCigs + mildexercise + modexercise + highexercise + alcohol  + diabetes + bmi + heartdis + KLB041A	+ KLB019A	+ KLB027I	+ KLB019F	+ KLB027A	+ KLB001a	+ KLB041B	+ KLB027N	+ KLB019B	+ KLB027J	+ KLB019G	+ KLB027B	+ KLB001b	+ KLB027J.1	+ KLB041C	+ KLB019C	+ KLB027K	+ KLB019H	+ KLB027C	+ KLB001c	+ KLB041D	+ KLB019D	+ KLB027L	+ KLB019I	+ KLB027D	+ KLB041E	+ KLB019E	+ KLB027M	+ KLB019J	+ KLB027E	+ KLB001e	+ KLB027M.1	+ KLB019K	+ KLB027F	+ KLB001f	+ KLB001g	+ KLB001h, data = bootstrap.analysis1,  family = "binomial")



###### PANEL LOGIT MODEL
#### Combine all datasets into one data frame

library(questionr)

# Non-imputed data
names(analysis1)<-gsub("KLB","LB",names(analysis1))

names(analysis2)<-gsub("LLB","LB",names(analysis2))

names(analysis3)<-gsub("MLB","LB",names(analysis3))

names(analysis4)<-gsub("NLB","LB",names(analysis4))

names(analysis5)<-gsub("OLB","LB",names(analysis5))

names(analysis6)<-gsub("PLB","LB",names(analysis6))



# Markov chain imputed data
names(markov.analysis1)<-gsub("KLB","LB",names(markov.analysis1))

names(markov.analysis2)<-gsub("LLB","LB",names(markov.analysis2))

names(markov.analysis3)<-gsub("MLB","LB",names(markov.analysis3))

names(markov.analysis4)<-gsub("NLB","LB",names(markov.analysis4))

names(markov.analysis5)<-gsub("OLB","LB",names(markov.analysis5))

names(markov.analysis6)<-gsub("PLB","LB",names(markov.analysis6))



# Boostrap imputed data
names(bootstrap.analysis1)<-gsub("KLB","LB",names(bootstrap.analysis1))

names(bootstrap.analysis2)<-gsub("LLB","LB",names(bootstrap.analysis2))

names(bootstrap.analysis3)<-gsub("MLB","LB",names(bootstrap.analysis3))

names(bootstrap.analysis4)<-gsub("NLB","LB",names(bootstrap.analysis4))

names(bootstrap.analysis5)<-gsub("OLB","LB",names(bootstrap.analysis5))

names(bootstrap.analysis6)<-gsub("PLB","LB",names(bootstrap.analysis6))




##### COMBINE DATA FRAMES
library(dplyr)
library(data.table)

# Create unique IDs
## Non-imputed data
analysis1$ID <- cumsum(!duplicated(analysis1[1:2]))

analysis2$ID <- cumsum(!duplicated(analysis2[1:2]))

analysis3$ID <- cumsum(!duplicated(analysis3[1:2]))

analysis4$ID <- cumsum(!duplicated(analysis4[1:2]))

analysis5$ID <- cumsum(!duplicated(analysis5[1:2]))

analysis6$ID <- cumsum(!duplicated(analysis6[1:2]))

## Markov imputed data
markov.analysis1$ID <- cumsum(!duplicated(markov.analysis1[1:2]))

markov.analysis2$ID <- cumsum(!duplicated(markov.analysis2[1:2]))

markov.analysis3$ID <- cumsum(!duplicated(markov.analysis3[1:2]))

markov.analysis4$ID <- cumsum(!duplicated(markov.analysis4[1:2]))

markov.analysis5$ID <- cumsum(!duplicated(markov.analysis5[1:2]))

markov.analysis6$ID <- cumsum(!duplicated(markov.analysis6[1:2]))

## Boostrap imputed dataset
bootstrap.analysis1$ID <- cumsum(!duplicated(bootstrap.analysis1[1:2]))

bootstrap.analysis2$ID <- cumsum(!duplicated(bootstrap.analysis2[1:2]))

bootstrap.analysis3$ID <- cumsum(!duplicated(bootstrap.analysis3[1:2]))

bootstrap.analysis4$ID <- cumsum(!duplicated(bootstrap.analysis4[1:2]))

bootstrap.analysis5$ID <- cumsum(!duplicated(bootstrap.analysis5[1:2]))

bootstrap.analysis6$ID <- cumsum(!duplicated(bootstrap.analysis6[1:2]))



## Non-imputed data
allHRSnon.imp <- bind_rows(analysis1, analysis2, analysis3, analysis4, analysis5, analysis6)

wave.yr <- c(rep(2006, 6529), rep(2008, 5584), rep(2010, 7655), rep(2012, 10405), rep(2014, 9107), rep(2016, 10894))

allHRSnon.imp <- cbind(allHRSnon.imp, wave.yr)

allHRSnon.imp <- pdata.frame(allHRSnon.imp, index=c("wave.yr","ID"), drop.index=TRUE, row.names=TRUE)


## Markov imputed data
allHRSmarkov <- bind_rows(markov.analysis1, markov.analysis2, markov.analysis3, markov.analysis4, markov.analysis5, markov.analysis6)

allHRSmarkov <- cbind(allHRSmarkov, wave.yr)

allHRSmarkov <- pdata.frame(allHRSmarkov, index=c("wave.yr","ID"), drop.index=TRUE, row.names=TRUE)



## Bootstrap imputed data
allHRSboot <- bind_rows(bootstrap.analysis1, bootstrap.analysis2, bootstrap.analysis3, bootstrap.analysis4, bootstrap.analysis5, bootstrap.analysis6)

allHRSboot <- cbind(allHRSboot, wave.yr)

allHRSboot <- pdata.frame(allHRSboot, index=c("wave.yr","ID"), drop.index=TRUE, row.names=TRUE)



  
########### CORRELATION MATRICES FOR PANEL DATA
corMarkov <- cor(allHRSmarkov[11:50])

corBoot <- cor(allHRSboot[11:50])

########## FACTOR ANALYSE PSYCHO-SOCIAL VARIABLES
#### Non-imputed data
# Correlation matrices for markov imputed dataset
pa.psyc <- data.frame(allHRSnon.imp[11:54], allHRSnon.imp[139:163])

# Set all NA to 0
pa.psyc[is.na(pa.psyc)] <- 0

# Correlation matrix
cor.pa.psyc <- cor(pa.psyc)

## Factor analysis
fa.pa.psyc <-  fa(r=cor.pa.psyc, nfactors = 5, rotate="none",fm="pa")

print(fa.pa.psyc)

fa.diagram(fa.pa.psyc)


## Markov chain imputed dataset
fa.markov <-  fa(r=corMarkov, nfactors = 5, rotate="none",fm="pa")

print(fa.markov)

fa.diagram(fa.markov)


#### Markov chain imputed data
# Correlation matrices for markov imputed dataset
pa.markov.psyc <- data.frame(allHRSmarkov[11:54], allHRSmarkov[97:101])

# Set all NA to 0
pa.markov.psyc[is.na(pa.markov.psyc)] <- 0

# Correlation matrix
corMarkov <- cor(pa.markov.psyc)


## Markov chain imputed dataset
fa.markov <-  fa(r=corMarkov, nfactors = 5, rotate="none",fm="pa")

print(fa.markov)

fa.diagram(fa.markov)


####### Bootstap imputed data
# Correlation matrices for bootstrap imputed dataset
pa.bootstrap.psyc <- data.frame(allHRSboot[11:54], allHRSboot[139:163])

# Set all NA to 0
pa.bootstrap.psyc[is.na(pa.bootstrap.psyc)] <- 0

# Correlation matrix
corBoot <- cor(pa.bootstrap.psyc)


## Bootsrap imputed dataset
fa.boot <-  fa(r=corBoot, nfactors = 5, rotate="none",fm="pa")

print(fa.boot)

fa.diagram(fa.boot)



##### PANEL REGRESSION MODELS
library(plm)

############ Interaction terms
### Non-imputed data
# Get factor scores
scores.fa.pa.psyc <-factor.scores(cor.pa.psyc, fa.pa.psyc)

# View factor weights
scores.fa.pa.psyc$weights


# Interaction terms
allHRSnon.imp$LB041A  <-  allHRSnon.imp$LB041A*0.0190867077
allHRSnon.imp$LB019A  <-  allHRSnon.imp$LB019A*0.0170609115
allHRSnon.imp$LB027I  <-  allHRSnon.imp$LB027I*0.0495217836
allHRSnon.imp$LB019F  <-  allHRSnon.imp$LB019F*0.0091559117
allHRSnon.imp$LB027A  <-  allHRSnon.imp$LB027A*0.0313475039
allHRSnon.imp$LB001a  <-  allHRSnon.imp$LB001a*0.0074000983
allHRSnon.imp$LB035A  <-  allHRSnon.imp$LB035A*0.0252438647
allHRSnon.imp$LB035B  <-  allHRSnon.imp$LB035B*0.0119972654
allHRSnon.imp$LB035C  <-  allHRSnon.imp$LB035C*0.0285400257
allHRSnon.imp$LB035D  <-  allHRSnon.imp$LB035D*0.0110713643
allHRSnon.imp$LB035F  <-  allHRSnon.imp$LB035F*0.0046102488
allHRSnon.imp$LB035G  <-  allHRSnon.imp$LB035G*0.0294958965
allHRSnon.imp$LB041B  <-  allHRSnon.imp$LB041B*0.0232476403
allHRSnon.imp$LB027N  <-  allHRSnon.imp$LB027N*0.0289086296
allHRSnon.imp$LB019B  <-  allHRSnon.imp$LB019B*0.0200555592
allHRSnon.imp$LB027J  <-  allHRSnon.imp$LB027J*-0.0186263254
allHRSnon.imp$LB019G  <-  allHRSnon.imp$LB019G*0.0191746015
allHRSnon.imp$LB027B  <-  allHRSnon.imp$LB027B*0.0253842821
allHRSnon.imp$LB001b  <-  allHRSnon.imp$LB001b*0.0140988627
allHRSnon.imp$LB027J.1  <-  allHRSnon.imp$LB027J.1*0.0537576403
allHRSnon.imp$LB041C  <-  allHRSnon.imp$LB041C*0.0174850388
allHRSnon.imp$LB019C  <-  allHRSnon.imp$LB019C*0.0106372583
allHRSnon.imp$LB027K  <-  allHRSnon.imp$LB027K*0.0132575726
allHRSnon.imp$LB019H  <-  allHRSnon.imp$LB019H*0.0188875273
allHRSnon.imp$LB027C  <-  allHRSnon.imp$LB027C*0.014402493
allHRSnon.imp$LB001c  <-  allHRSnon.imp$LB001c*0.0138934301
allHRSnon.imp$LB041D  <-  allHRSnon.imp$LB041D*0.0174912096
allHRSnon.imp$LB019D  <-  allHRSnon.imp$LB019D*0.0193869513
allHRSnon.imp$LB027L  <-  allHRSnon.imp$LB027L*0.0118558247
allHRSnon.imp$LB019I  <-  allHRSnon.imp$LB019I*0.0348360128
allHRSnon.imp$LB027D  <-  allHRSnon.imp$LB027D*0.0236680157
allHRSnon.imp$LB041E  <-  allHRSnon.imp$LB041E*0.018081511
allHRSnon.imp$LB019E  <-  allHRSnon.imp$LB019E*0.0129767545
allHRSnon.imp$LB027M  <-  allHRSnon.imp$LB027M*0.0344407424
allHRSnon.imp$LB019J  <-  allHRSnon.imp$LB019J*0.0176543198
allHRSnon.imp$LB027E  <-  allHRSnon.imp$LB027E*0.017404097
allHRSnon.imp$LB001e  <-  allHRSnon.imp$LB001e*0.0109200454
allHRSnon.imp$LB027M.1  <-  allHRSnon.imp$LB027M.1*0.0264008537
allHRSnon.imp$LB035E  <-  allHRSnon.imp$LB035E*0.0136561853
allHRSnon.imp$LB019K  <-  allHRSnon.imp$LB019K*0.0135823397
allHRSnon.imp$LB027F  <-  allHRSnon.imp$LB027F*0.0162983131
allHRSnon.imp$LB001f  <-  allHRSnon.imp$LB001f*0.0139298589
allHRSnon.imp$LB001g  <-  allHRSnon.imp$LB001g*0.0030319251
allHRSnon.imp$LB001h  <-  allHRSnon.imp$LB001h*-0.0002525305
allHRSnon.imp$LB026S  <-  allHRSnon.imp$LB026S*-0.0004093334
allHRSnon.imp$LB018A  <-  allHRSnon.imp$LB018A*-0.0127918261
allHRSnon.imp$LB001I  <-  allHRSnon.imp$LB001I*-0.0046966956
allHRSnon.imp$LB033A  <-  allHRSnon.imp$LB033A*-0.0288703031
allHRSnon.imp$LB033B  <-  allHRSnon.imp$LB033B*-0.0131744989
allHRSnon.imp$LB033C  <-  allHRSnon.imp$LB033C*-0.0267520513
allHRSnon.imp$LB033D  <-  allHRSnon.imp$LB033D*-0.012344937
allHRSnon.imp$LB033F  <-  allHRSnon.imp$LB033F*-0.0108074306
allHRSnon.imp$LB033G  <-  allHRSnon.imp$LB033G*-0.0210819517
allHRSnon.imp$LB031L  <-  allHRSnon.imp$LB031L*-0.0226780226
allHRSnon.imp$LB026R  <-  allHRSnon.imp$LB026R*-0.0423246391
allHRSnon.imp$LB018H  <-  allHRSnon.imp$LB018H*-0.0083699343
allHRSnon.imp$LB018B  <-  allHRSnon.imp$LB018B*-0.0194300728
allHRSnon.imp$LB001R  <-  allHRSnon.imp$LB001R*-0.005498494
allHRSnon.imp$LB026J  <-  allHRSnon.imp$LB026J*-0.0214741486
allHRSnon.imp$LB018C  <-  allHRSnon.imp$LB018C*-0.018273239
allHRSnon.imp$LB026K  <-  allHRSnon.imp$LB026K*-0.0293955536
allHRSnon.imp$LB018D  <-  allHRSnon.imp$LB018D*-0.0348140529
allHRSnon.imp$LB026X  <-  allHRSnon.imp$LB026X*-0.022708698
allHRSnon.imp$LB018E  <-  allHRSnon.imp$LB018E*-0.0175204032
allHRSnon.imp$LB002C  <-  allHRSnon.imp$LB002C*-0.0189648546
allHRSnon.imp$LB033E  <-  allHRSnon.imp$LB033E*-0.0159007972
allHRSnon.imp$LB018F  <-  allHRSnon.imp$LB018F*-0.0130001371
allHRSnon.imp$LB001N  <-  allHRSnon.imp$LB001N*-0.0085322701
allHRSnon.imp$LB012B  <-  allHRSnon.imp$LB012B*-0.0177190885



### Markov imputed data
# Get factor scores
scores.fa.markov <-factor.scores(corMarkov, fa.markov)

# View factor weights
scores.fa.markov$weights


# Interaction terms
allHRSmarkov$LB041A  <-  allHRSmarkov$LB041A*0.0293847
allHRSmarkov$LB019A  <-  allHRSmarkov$LB019A*0.02520319
allHRSmarkov$LB027I  <-  allHRSmarkov$LB027I*0.06813449
allHRSmarkov$LB019F  <-  allHRSmarkov$LB019F*0.01548507
allHRSmarkov$LB027A  <-  allHRSmarkov$LB027A*0.03130092
allHRSmarkov$LB001a  <-  allHRSmarkov$LB001a*0.01021483
allHRSmarkov$LB035A  <-  allHRSmarkov$LB035A*0.03533286
allHRSmarkov$LB035B  <-  allHRSmarkov$LB035B*0.01528324
allHRSmarkov$LB035C  <-  allHRSmarkov$LB035C*0.04104408
allHRSmarkov$LB035D  <-  allHRSmarkov$LB035D*0.01413707
allHRSmarkov$LB035F  <-  allHRSmarkov$LB035F*0.005191974
allHRSmarkov$LB035G  <-  allHRSmarkov$LB035G*0.03845109
allHRSmarkov$LB041B  <-  allHRSmarkov$LB041B*0.04073584
allHRSmarkov$LB027N  <-  allHRSmarkov$LB027N*0.05302067
allHRSmarkov$LB019B  <-  allHRSmarkov$LB019B*0.0286091
allHRSmarkov$LB027J  <-  allHRSmarkov$LB027J*0.03799974
allHRSmarkov$LB019G  <-  allHRSmarkov$LB019G*0.02759492
allHRSmarkov$LB027B  <-  allHRSmarkov$LB027B*0.03491695
allHRSmarkov$LB001b  <-  allHRSmarkov$LB001b*0.02327305
allHRSmarkov$LB027J.1  <-  allHRSmarkov$LB027J.1*0.00628966
allHRSmarkov$LB041C  <-  allHRSmarkov$LB041C*0.02638563
allHRSmarkov$LB019C  <-  allHRSmarkov$LB019C*0.01556632
allHRSmarkov$LB027K  <-  allHRSmarkov$LB027K*0.01886593
allHRSmarkov$LB019H  <-  allHRSmarkov$LB019H*0.02718491
allHRSmarkov$LB027C  <-  allHRSmarkov$LB027C*0.02893426
allHRSmarkov$LB001c  <-  allHRSmarkov$LB001c*0.01826416
allHRSmarkov$LB041D  <-  allHRSmarkov$LB041D*0.02653583
allHRSmarkov$LB019D  <-  allHRSmarkov$LB019D*0.02733848
allHRSmarkov$LB027L  <-  allHRSmarkov$LB027L*0.02139964
allHRSmarkov$LB019I  <-  allHRSmarkov$LB019I*0.04310159
allHRSmarkov$LB027D  <-  allHRSmarkov$LB027D*0.03073444
allHRSmarkov$LB041E  <-  allHRSmarkov$LB041E*0.02621988
allHRSmarkov$LB019E  <-  allHRSmarkov$LB019E*0.02041907
allHRSmarkov$LB027M  <-  allHRSmarkov$LB027M*0.05946887
allHRSmarkov$LB019J  <-  allHRSmarkov$LB019J*0.02560605
allHRSmarkov$LB027E  <-  allHRSmarkov$LB027E*0.01932771
allHRSmarkov$LB001e  <-  allHRSmarkov$LB001e*0.007956669
allHRSmarkov$LB027M.1  <-  allHRSmarkov$LB027M.1*0.04387439
allHRSmarkov$LB035E  <-  allHRSmarkov$LB035E*0.01776397
allHRSmarkov$LB019K  <-  allHRSmarkov$LB019K*0.01901691
allHRSmarkov$LB027F  <-  allHRSmarkov$LB027F*0.02513294
allHRSmarkov$LB001f  <-  allHRSmarkov$LB001f*0.0363402
allHRSmarkov$LB001g  <-  allHRSmarkov$LB001g*-0.005695116
allHRSmarkov$LB001h  <-  allHRSmarkov$LB001h*-0.01162157
allHRSmarkov$LB001A  <-  allHRSmarkov$LB001A*-0.002494386
allHRSmarkov$LB001B  <-  allHRSmarkov$LB001B*0.00004556611
allHRSmarkov$LB001C  <-  allHRSmarkov$LB001C*-0.004819121
allHRSmarkov$LB001F  <-  allHRSmarkov$LB001F*-0.0003504254
allHRSmarkov$LB001G  <-  allHRSmarkov$LB001G*-0.01566898



### Bootstrap imputed data
# Get factor scores
scores.fa.boot <-factor.scores(corBoot, fa.boot)

# View factor weights
scores.fa.boot$weights


# Interaction terms
allHRSboot$LB041A  <-  allHRSboot$LB041A*0.0229937704
allHRSboot$LB019A  <-  allHRSboot$LB019A*0.0142004729
allHRSboot$LB027I  <-  allHRSboot$LB027I*0.0548408383
allHRSboot$LB019F  <-  allHRSboot$LB019F*0.0062409069
allHRSboot$LB027A  <-  allHRSboot$LB027A*0.0240639065
allHRSboot$LB001a  <-  allHRSboot$LB001a*0.0046331234
allHRSboot$LB035A  <-  allHRSboot$LB035A*0.029287254
allHRSboot$LB035B  <-  allHRSboot$LB035B*0.010016152
allHRSboot$LB035C  <-  allHRSboot$LB035C*0.0314751139
allHRSboot$LB035D  <-  allHRSboot$LB035D*0.0104077965
allHRSboot$LB035F  <-  allHRSboot$LB035F*0.0052603073
allHRSboot$LB035G  <-  allHRSboot$LB035G*0.0210907548
allHRSboot$LB041B  <-  allHRSboot$LB041B*0.0264006127
allHRSboot$LB027N  <-  allHRSboot$LB027N*0.0213275809
allHRSboot$LB019B  <-  allHRSboot$LB019B*0.0161924656
allHRSboot$LB027J  <-  allHRSboot$LB027J*0.0095841026
allHRSboot$LB019G  <-  allHRSboot$LB019G*0.0163227948
allHRSboot$LB027B  <-  allHRSboot$LB027B*0.0236455173
allHRSboot$LB001b  <-  allHRSboot$LB001b*0.0124004295
allHRSboot$LB027J.1  <-  allHRSboot$LB027J.1*0.0255908102
allHRSboot$LB041C  <-  allHRSboot$LB041C*0.0208366912
allHRSboot$LB019C  <-  allHRSboot$LB019C*0.0102555285
allHRSboot$LB027K  <-  allHRSboot$LB027K*0.0176764444
allHRSboot$LB019H  <-  allHRSboot$LB019H*0.0171456535
allHRSboot$LB027C  <-  allHRSboot$LB027C*0.0185663526
allHRSboot$LB001c  <-  allHRSboot$LB001c*0.0105416578
allHRSboot$LB041D  <-  allHRSboot$LB041D*0.003687403
allHRSboot$LB019D  <-  allHRSboot$LB019D*0.0183657486
allHRSboot$LB027L  <-  allHRSboot$LB027L*0.0155377842
allHRSboot$LB019I  <-  allHRSboot$LB019I*0.0320885641
allHRSboot$LB027D  <-  allHRSboot$LB027D*0.0219698697
allHRSboot$LB041E  <-  allHRSboot$LB041E*0.0083815078
allHRSboot$LB019E  <-  allHRSboot$LB019E*0.013506215
allHRSboot$LB027M  <-  allHRSboot$LB027M*0.0199016038
allHRSboot$LB019J  <-  allHRSboot$LB019J*0.0157360909
allHRSboot$LB027E  <-  allHRSboot$LB027E*0.0097308298
allHRSboot$LB001e  <-  allHRSboot$LB001e*0.0050110548
allHRSboot$LB027M.1  <-  allHRSboot$LB027M.1*0.0623737653
allHRSboot$LB035E  <-  allHRSboot$LB035E*0.0134102232
allHRSboot$LB019K  <-  allHRSboot$LB019K*0.0117649442
allHRSboot$LB027F  <-  allHRSboot$LB027F*0.0194766404
allHRSboot$LB001f  <-  allHRSboot$LB001f*0.0194368439
allHRSboot$LB001g  <-  allHRSboot$LB001g*0.0013207345
allHRSboot$LB001h  <-  allHRSboot$LB001h*0.0001165801
allHRSboot$LB026S  <-  allHRSboot$LB026S*0.0134628382
allHRSboot$LB018A  <-  allHRSboot$LB018A*-0.0099820938
allHRSboot$LB001I  <-  allHRSboot$LB001I*-0.0039960179
allHRSboot$LB033A  <-  allHRSboot$LB033A*-0.0264330226
allHRSboot$LB033B  <-  allHRSboot$LB033B*-0.031249251
allHRSboot$LB033C  <-  allHRSboot$LB033C-0.0260517705
allHRSboot$LB033D  <-  allHRSboot$LB033D-0.0082591583
allHRSboot$LB033F  <-  allHRSboot$LB033F-0.0083637751
allHRSboot$LB033G  <-  allHRSboot$LB033G-0.0235369891
allHRSboot$LB031L  <-  allHRSboot$LB031L-0.017974738
allHRSboot$LB026R  <-  allHRSboot$LB026R-0.038448691
allHRSboot$LB018H  <-  allHRSboot$LB018H-0.0059177688
allHRSboot$LB018B  <-  allHRSboot$LB018B-0.0199945322
allHRSboot$LB001R  <-  allHRSboot$LB001R-0.0051803433
allHRSboot$LB026J  <-  allHRSboot$LB026J-0.0169661613
allHRSboot$LB018C  <-  allHRSboot$LB018C-0.0185227102
allHRSboot$LB026K  <-  allHRSboot$LB026K-0.0339759168
allHRSboot$LB018D  <-  allHRSboot$LB018D-0.0307033135
allHRSboot$LB026X  <-  allHRSboot$LB026X-0.0225638844
allHRSboot$LB018E  <-  allHRSboot$LB018E-0.0146943275
allHRSboot$LB002C  <-  allHRSboot$LB002C-0.0186208455
allHRSboot$LB033E  <-  allHRSboot$LB033E-0.0185762369
allHRSboot$LB018F  <-  allHRSboot$LB018F-0.0125387891
allHRSboot$LB001N  <-  allHRSboot$LB001N-0.0082508083
allHRSboot$LB012B  <-  allHRSboot$LB012B-0.0318671084




############ Non-imputed data
### Fixed effects model
mod2.1 <- plm(stroke~KIMpurposeinlife + age + sex + stroke ~ KIMpurposeinlife + age  + sex + marstat + race + education + wealth + func014 + func016 + func021 + func023 + func025 + numCigs + mildexercise + modexercise + highexercise + alcohol  + diabetes + bmi + LB041A	+ LB019A + LB027I + LB019F + LB027A + LB001a + LB041B + LB027N + LB019B + LB027J + LB019G + LB027B + LB027B + LB001b + LB027J.1	+ LB041C + LB019C	+ LB027K 	+ LB019H 	+ LB027C + LB001c + LB041D + LB019D	+ LB027L	+ LB019I 	+ LB027D	+ LB041E  + LB019E	+ LB027M	+ LB019J + LB027E	+ LB001e	+ LB027M.1 	+ LB019K	+ LB027F	+ LB001f 	+ LB001g, data = allHRSnon.imp, model = "within")
              


######## Markov chain imputed data
### Fixed effects model
mod2.2<- plm(stroke~KIMpurposeinlife + age + sex + stroke ~ KIMpurposeinlife + age  + sex + marstat + race + education + wealth + func014 + func016 + func021 + func023 + func025 + numCigs + mildexercise + modexercise + highexercise + alcohol  + diabetes + bmi + LB041A	+ LB019A + LB027I + LB019F + LB027A + LB001a + LB041B + LB027N + LB019B + LB027J + LB019G + LB027B + LB027B + LB001b + LB027J.1	+ LB041C + LB019C	+ LB027K 	+ LB019H 	+ LB027C + LB001c + LB041D + LB019D	+ LB027L	+ LB019I 	+ LB027D	+ LB041E  + LB019E	+ LB027M	+ LB019J + LB027E	+ LB001e	+ LB027M.1 	+ LB019K	+ LB027F	+ LB001f, data = allHRSmarkov, model = "within")


######## Bootstrap imputed data
### Fixed effects model
mod2.3 <- plm(stroke~KIMpurposeinlife + age + sex + stroke ~ KIMpurposeinlife + age  + sex + marstat + race + education + wealth + func014 + func016 + func021 + func023 + func025 + numCigs + mildexercise + modexercise + highexercise + alcohol  + diabetes + bmi + LB041A	+ LB019A + LB027I + LB019F + LB027A + LB001a + LB041B + LB027N + LB019B + LB027J + LB019G + LB027B + LB027B + LB001b + LB027J.1	+ LB041C + LB019C	+ LB027K 	+ LB019H 	+ LB027C + LB001c + LB041D + LB019D	+ LB027L	+ LB019I 	+ LB027D	+ LB041E  + LB019E	+ LB027M	+ LB019J + LB027E	+ LB001e	+ LB027M.1 	+ LB019K	+ LB027F	+ LB001f 	+ LB001g, data = allHRSboot, model = "within")




### PLOT REGRESSION ESTIMATORS
library(jtools)
library(ggstance)
library(dotwhisker)
library(dplyr)

## Data frame of regression estimators
coef.vec <- c(mod1.1$coefficients[2], mod1.2$coefficients[2], mod1.3$coefficients[2], mod1.4$coefficients[2], mod1.5$coefficients[2], mod1.6$coefficients[2], mod2.1$coefficients[1], mod2.2$coefficients[1], mod2.3$coefficients[1])

se.vec <- c(9.635697, 9.635663, 6.451053, 6.667106, 6.565470, 6.782787, 8.5990, 5.2675, 5.4604)

var.names <- c("1.1 Core", "1.2 Core Factor", "1.3 Markov Core", "1.4 Bootstrap Core", "1.5 Markov Factor", "1.6 Bootstrap Factor", "2.1 Core Fixed Effects Factor", "2.2 Markov Fixed Effects Factor", "2.3 Bootstrap Fixed Effects Factor")




# Format data as tidy dataframe
results_df <- data.frame(term=var.names, estimate=coef.vec,
                         std.error=se.vec)

# Draw dot-and-whisker plot
results_df %>% dwplot + theme_bw() + theme(legend.position="none") +
  ggtitle("Model Specifications of the Relationship Between Purpose in Life and Stroke") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  annotate("text", x = 1.05, y = 8, size = 4, hjust = 0,
           label = "n = 10,894") + theme(plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5))



