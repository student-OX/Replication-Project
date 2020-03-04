######## Create analysis datasets for each cross-section
#### Stroke incidence data frames
stroke08 <- data.frame("HHID"= HRS2008$HHID, "PN"= HRS2008$PN, "LC053"= HRS2008$LC053, "VC053"= HRS2008$VC053)

stroke10 <- data.frame("HHID"= HRS2010$HHID, "PN"= HRS2010$PN, "MC053"= HRS2010$MC053, "WC053"= HRS2010$WC053)

stroke12 <- data.frame("HHID"= HRS2012$HHID, "PN"= HRS2012$PN, "NC053"= HRS2012$NC053, "XC053"= HRS2012$XC053)

stroke14 <- data.frame("HHID"= HRS2014$HHID, "PN"= HRS2014$PN, "OC053"= HRS2014$OC053, "YC053"= HRS2014$YC053)

stroke16 <- data.frame("HHID"= HRS2016$HHID, "PN"= HRS2016$PN, "PC053"= HRS2016$PC053, "ZC053"= HRS2016$ZC053)

stroke18 <- data.frame("HHID"= HRS2018$HHID, "PN"= HRS2018$PN, "QC053"= HRS2018$QC053)


#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1 <- merge(HRS2006, stroke08, by= intersect(names(HRS2006), names(stroke08)))

analysis1 <- merge(analysis1, stroke10, by= intersect(names(analysis1), names(stroke10)))


#### 2008 baseline covariates, 2010 and 2012 stroke incidence

analysis2 <- merge(HRS2008, stroke10, by= intersect(names(HRS2008), names(stroke10)))

analysis2 <- merge(analysis2, stroke12, by= intersect(names(analysis2), names(stroke12)))


#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3 <- merge(HRS2010, stroke12, by= intersect(names(HRS2010), names(stroke12)))

analysis3 <- merge(analysis3, stroke14, by= intersect(names(analysis3), names(stroke14)))



#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4 <- merge(HRS2012, stroke14, by= intersect(names(HRS2012), names(stroke14)))

analysis4 <- merge(analysis4, stroke16, by= intersect(names(analysis4), names(stroke16)))




#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5 <- merge(HRS2014, stroke16, by= intersect(names(HRS2014), names(stroke16)))

analysis5 <- merge(analysis5, stroke18, by= intersect(names(analysis5), names(stroke18)))


#### 2016 baseline covariates and 2018 stroke incidence
analysis6 <- merge(HRS2016, stroke18, by= intersect(names(HRS2016), names(stroke18)))



###### Restrict sample to respondents who only took part in face-to-face enhanced self-report questionnaire
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1 <- subset(analysis1, KLBRTYPE ==4)



#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2 <- subset(analysis2, LLBRTYPE ==4)



#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3 <- subset(analysis3, MLBRTYPE ==4)


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4 <- subset(analysis4, NLBRTYPE ==4)


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5 <- subset(analysis5, OLBRTYPE ==4)


#### 2016 baseline covariates and 2018 stroke incidence
analysis6 <- subset(analysis6, PLBRTYPE ==4)




###### Restrict sample to respondents with self-reported history of being stroke-free at baseline 
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1 <- subset(analysis1, KC053 == 4 | KC053 == 5)

## Frequency table
w = table(analysis1$KC053)

as.data.frame(w)


#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2 <- subset(analysis2, LC053 == 4 | LC053 == 5)

## Frequency table
w = table(analysis2$LC053)

as.data.frame(w)



#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3 <- subset(analysis3, MC053 == 4 | MC053 == 5)

## Frequency table
w = table(analysis3$MC053)

as.data.frame(w)


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4 <- subset(analysis4, NC053 == 4 | NC053 == 5)

## Frequency table
w = table(analysis4$NC053)

as.data.frame(w)


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5 <- subset(analysis5, OC053 == 4 | OC053 == 5)

## Frequency table
w = table(analysis5$OC053)

as.data.frame(w)




#### 2016 baseline covariates and 2018 stroke incidence
analysis6 <- subset(analysis6, PC053 == 4 | PC053 == 5)

## Frequency table
w = table(analysis6$PC053)

as.data.frame(w)




#### Count final sample size
#### (6529) 2006 baseline covariates, 2008 and 2010 stroke incidence

#### (5584) 2008 baseline covariates, 2010 and 2012 stroke incidence


#### (7655) 2010 baseline covariates, 2012 and 2014 stroke incidence


#### (10405) 2012 baseline covariates, 2014 and 2016 stroke incidence


#### (9107) 2014 baseline covariates, 2016 and 2018 stroke incidence


#### (10894) 2016 baseline covariates and 2018 stroke incidence



################### RECODE AND RENAME VARIABLES

############# STROKE

#### 2006 baseline covariates, 2008 and 2010 stroke incidence
# Create dummy variable for stroke incidence using reported stroke information in 2008 and 2010 (LC053, VC053, MC053, WC053)
analysis1$stroke <- analysis1$LC053

analysis1$stroke[analysis1$stroke==2] <- 1

analysis1$stroke[analysis1$stroke==3] <- 1

analysis1$stroke[analysis1$stroke==4] <- 0

analysis1$stroke[analysis1$stroke==5] <- 0

analysis1$stroke[analysis1$stroke==5] <- 0

analysis1$stroke[analysis1$stroke==8] <- NA

analysis1$stroke[analysis1$stroke==9] <- NA

# VC053
analysis1$stroke[analysis1$VC053==1] <- 1

analysis1$stroke[analysis1$VC053==2] <- 1

analysis1$stroke[analysis1$VC053==3] <- 1

# MC053
analysis1$stroke[analysis1$MC053==1] <- 1

analysis1$stroke[analysis1$MC053==2] <- 1

analysis1$stroke[analysis1$MC053==3] <- 1


# WC053
analysis1$stroke[analysis1$WC053==1] <- 1

analysis1$stroke[analysis1$WC053==2] <- 1

analysis1$stroke[analysis1$WC053==3] <- 1




#### 2008 baseline covariates, 2010 and 2012 stroke incidence
# Create dummy variable for stroke incidence using reported stroke information in 2010 and 2012 (MC053, WC053, NC053, XC053)
analysis2$stroke <- analysis2$MC053

analysis2$stroke[analysis2$stroke==2] <- 1

analysis2$stroke[analysis2$stroke==3] <- 1

analysis2$stroke[analysis2$stroke==4] <- 0

analysis2$stroke[analysis2$stroke==5] <- 0

analysis2$stroke[analysis2$stroke==5] <- 0

analysis2$stroke[analysis2$stroke==8] <- NA

analysis2$stroke[analysis2$stroke==9] <- NA

# WC035
analysis2$stroke[analysis2$WC053==1] <- 1

analysis2$stroke[analysis2$WC053==2] <- 1

analysis2$stroke[analysis2$WC053==3] <- 1

# NC053
analysis2$stroke[analysis2$NC053==1] <- 1

analysis2$stroke[analysis2$NC053==2] <- 1

analysis2$stroke[analysis2$NC053==3] <- 1


# XC053
analysis2$stroke[analysis2$XC053==1] <- 1

analysis2$stroke[analysis2$XC053==2] <- 1

analysis2$stroke[analysis2$XC053==3] <- 1



#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3$stroke <- analysis3$NC053

analysis3$stroke[analysis3$stroke==2] <- 1

analysis3$stroke[analysis3$stroke==3] <- 1

analysis3$stroke[analysis3$stroke==4] <- 0

analysis3$stroke[analysis3$stroke==5] <- 0

analysis3$stroke[analysis3$stroke==5] <- 0

analysis3$stroke[analysis3$stroke==8] <- NA

analysis3$stroke[analysis3$stroke==9] <- NA

# XC053
analysis3$stroke[analysis3$XC053==1] <- 1

analysis3$stroke[analysis3$XC053==2] <- 1

analysis3$stroke[analysis3$XC053==3] <- 1


# OC053
analysis3$stroke[analysis3$OC053==1] <- 1

analysis3$stroke[analysis3$OC053==2] <- 1

analysis3$stroke[analysis3$OC053==3] <- 1


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4$stroke <- analysis4$OC053

analysis4$stroke[analysis4$stroke==2] <- 1

analysis4$stroke[analysis4$stroke==3] <- 1

analysis4$stroke[analysis4$stroke==4] <- 0

analysis4$stroke[analysis4$stroke==5] <- 0

analysis4$stroke[analysis4$stroke==5] <- 0

analysis4$stroke[analysis4$stroke==8] <- NA

analysis4$stroke[analysis4$stroke==9] <- NA

# YC053
analysis4$stroke[analysis4$YC053==1] <- 1

analysis4$stroke[analysis4$YC053==2] <- 1

analysis4$stroke[analysis4$YC053==3] <- 1

# PC053
analysis4$stroke[analysis4$PC053==1] <- 1

analysis4$stroke[analysis4$PC053==2] <- 1

analysis4$stroke[analysis4$PC053==3] <- 1

# ZC053
analysis4$stroke[analysis4$ZC053==1] <- 1

analysis4$stroke[analysis4$ZC053==2] <- 1

analysis4$stroke[analysis4$ZC053==3] <- 1


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5$stroke <- analysis5$PC053

analysis5$stroke[analysis5$stroke==2] <- 1

analysis5$stroke[analysis5$stroke==3] <- 1

analysis5$stroke[analysis5$stroke==4] <- 0

analysis5$stroke[analysis5$stroke==5] <- 0

analysis5$stroke[analysis5$stroke==5] <- 0

analysis5$stroke[analysis5$stroke==8] <- NA

analysis5$stroke[analysis5$stroke==9] <- NA

# ZC053
analysis5$stroke[analysis5$ZC053==1] <- 1

analysis5$stroke[analysis5$ZC053==2] <- 1

analysis5$stroke[analysis5$ZC053==3] <- 1

# QC053
analysis5$stroke[analysis5$QC053==1] <- 1

analysis5$stroke[analysis5$QC053==2] <- 1

analysis5$stroke[analysis5$QC053==3] <- 1



#### 2016 baseline covariates and 2018 stroke incidence
analysis6$stroke <- analysis6$QC053

analysis6$stroke[analysis6$stroke==2] <- 1

analysis6$stroke[analysis6$stroke==3] <- 1

analysis6$stroke[analysis6$stroke==4] <- 0

analysis6$stroke[analysis6$stroke==5] <- 0

analysis6$stroke[analysis6$stroke==5] <- 0

analysis6$stroke[analysis6$stroke==8] <- NA

analysis6$stroke[analysis6$stroke==9] <- NA




############# AGE AND SEX

#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1$age <- analysis1$KA019

analysis1$sex <- analysis1$KX060_R



#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2$age <- analysis2$LA019

analysis2$sex <- analysis2$LX060_R



#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3$age <- analysis3$MA019

analysis3$sex <- analysis3$MX060_R


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4$age <- analysis4$NA019

analysis4$sex <- analysis4$NX060_R




#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5$age <- analysis5$OA019

analysis5$sex <- analysis5$OX060_R


#### 2016 baseline covariates and 2018 stroke incidence
analysis6$age <- analysis6$PA019

analysis6$sex <- analysis6$PX060_R


############# MARITAL STATUS (1- MARRIED, 2-NOT MARRIED)
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1$marstat <- analysis1$KX065_R

analysis1$marstat[analysis1$marstat==2] <- 1

analysis1$marstat[analysis1$marstat==2] <- 1

analysis1$marstat[analysis1$marstat==3] <- 1

analysis1$marstat[analysis1$marstat==6] <- 0

analysis1$marstat[analysis1$KX065_R=="."] <- NA


## Frequency table
w = table(analysis1$marstat)

as.data.frame(w)


#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2$marstat <- analysis2$LX065_R

analysis2$marstat[analysis2$marstat==2] <- 1

analysis2$marstat[analysis2$marstat==2] <- 1

analysis2$marstat[analysis2$marstat==3] <- 1

analysis2$marstat[analysis2$marstat==6] <- 0

analysis2$marstat[analysis2$LX065_R=="."] <- NA

## Frequency table
w = table(analysis2$marstat)

as.data.frame(w)


#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3$marstat <- analysis3$MX065_R

analysis3$marstat[analysis3$marstat==2] <- 1

analysis3$marstat[analysis3$marstat==2] <- 1

analysis3$marstat[analysis3$marstat==3] <- 1

analysis3$marstat[analysis3$marstat==4] <- 1

analysis3$marstat[analysis3$marstat==6] <- 0

analysis3$marstat[analysis3$MX065_R=="."] <- NA

## Frequency table
w = table(analysis3$marstat)

as.data.frame(w)


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4$marstat <- analysis4$NX065_R

analysis4$marstat[analysis4$marstat==2] <- 1

analysis4$marstat[analysis4$marstat==2] <- 1

analysis4$marstat[analysis4$marstat==3] <- 1

analysis4$marstat[analysis4$marstat==4] <- 1

analysis4$marstat[analysis4$marstat==6] <- 0

analysis4$marstat[analysis4$NX065_R=="."] <- NA

## Frequency table
w = table(analysis4$marstat)

as.data.frame(w)


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5$marstat <- analysis5$OX065_R

analysis5$marstat[analysis5$marstat==2] <- 1

analysis5$marstat[analysis5$marstat==2] <- 1

analysis5$marstat[analysis5$marstat==3] <- 1

analysis5$marstat[analysis5$marstat==4] <- 1

analysis5$marstat[analysis5$marstat==6] <- 0

analysis5$marstat[analysis5$OX065_R=="."] <- NA

## Frequency table
w = table(analysis5$marstat)

as.data.frame(w)


#### 2016 baseline covariates and 2018 stroke incidence
analysis6$marstat <- analysis6$PX065_R

analysis6$marstat[analysis6$marstat==2] <- 1

analysis6$marstat[analysis6$marstat==3] <- 1

analysis6$marstat[analysis6$marstat==6] <- 0

analysis6$marstat[analysis6$PX065_R=="."] <- NA


## Frequency table
w = table(analysis6$marstat)

as.data.frame(w)


############# RACE/ETHNICITY- KIM ET AL, 2013 (1- CAUCASIAN, 2- AFRICAN AMERICAN, 3- HISPANIC)

#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1$KIMrace <- analysis1$KB089M1M

analysis1$KIMrace[analysis1$KIMrace==97] <- 4

analysis1$KIMrace[analysis1$KB028==1] <- 3

analysis1$KIMrace[analysis1$KIMrace==98] <- NA

analysis1$KIMrace[analysis1$KIMrace==99] <- NA

## Frequency table
w = table(analysis1$KIMrace)

as.data.frame(w)



#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2$KIMrace <- analysis2$LB089M1M

analysis2$KIMrace[analysis2$KIMrace==97] <- 4

analysis2$KIMrace[analysis2$LB028==1] <- 3

## Frequency table
w = table(analysis2$KIMrace)

as.data.frame(w)


#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3$KIMrace <- analysis3$MB089M1M

analysis3$KIMrace[analysis3$KIMrace==97] <- 4

analysis3$KIMrace[analysis3$KIMrace==98] <- NA

analysis3$KIMrace[analysis3$KIMrace==99] <- NA

analysis3$KIMrace[analysis3$MB028==1] <- 3

## Frequency table
w = table(analysis3$KIMrace)

as.data.frame(w)


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4$KIMrace <- analysis4$NB089M1M

analysis4$KIMrace[analysis4$KIMrace==97] <- 4

analysis4$KIMrace[analysis4$KIMrace==98] <- NA

analysis4$KIMrace[analysis4$NB028==1] <- 3


## Frequency table
w = table(analysis4$KIMrace)

as.data.frame(w)

#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5$KIMrace <- analysis5$OB089M1M

analysis5$KIMrace[analysis5$KIMrace==97] <- 4

analysis5$KIMrace[analysis5$OB028==1] <- 3


## Frequency table
w = table(analysis5$KIMrace)

as.data.frame(w)

#### 2016 baseline covariates and 2018 stroke incidence
analysis6$KIMrace <- analysis6$PB089M1M

analysis6$KIMrace[analysis6$KIMrace==97] <- 4

analysis6$KIMrace[analysis6$KIMrace==98] <- NA

analysis6$KIMrace[analysis6$KIMrace==99] <- NA

analysis6$KIMrace[analysis6$PB028==1] <- 3


## Frequency table
w = table(analysis6$KIMrace)

as.data.frame(w)



############# RACE/ETHNICITY (1- CAUCASIAN, 2- AFRICAN AMERICAN, 3- HISPANIC, 4- OTHER, 5- MIXED)

#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1$race <- analysis1$KB089M1M

analysis1$race[analysis1$race==97] <- 4

analysis1$race[analysis1$KB028==1] <- 3

analysis1$race[analysis1$KB089M1M==1 & analysis1$KB089M2M==2] <- 5

analysis1$race[analysis1$KB089M1M==2 & analysis1$KB089M2M==1] <- 5

analysis1$race[analysis1$KB089M1M==1 & analysis1$KB089M2M==97] <- 5

analysis1$race[analysis1$KB089M1M==97 & analysis1$KB089M2M==1] <- 5

analysis1$race[analysis1$KB089M1M==1 & analysis1$KB089M3M==2] <- 5

analysis1$race[analysis1$KB089M1M==2 & analysis1$KB089M3M==1] <- 5

analysis1$race[analysis1$race==98] <- NA

analysis1$race[analysis1$race==99] <- NA

## Frequency table
w = table(analysis1$race)

as.data.frame(w)


#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2$race <- analysis2$LB089M1M

analysis2$race[analysis2$race==97] <- 4

analysis2$race[analysis2$LB028==1] <- 3

analysis2$race[analysis2$LB089M1M==1 & analysis2$LB089M2M==2] <- 5

analysis2$race[analysis2$LB089M1M==2 & analysis2$LB089M2M==1] <- 5

analysis2$race[analysis2$LB089M1M==1 & analysis2$LB089M2M==97] <- 5

analysis2$race[analysis2$LB089M1M==97 & analysis2$LB089M2M==1] <- 5

analysis2$race[analysis2$LB089M1M==1 & analysis2$LB089M3M==2] <- 5

analysis2$race[analysis2$LB089M1M==2 & analysis2$LB089M3M==1] <- 5


## Frequency table
w = table(analysis2$race)

as.data.frame(w)

#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3$race <- analysis3$MB089M1M

analysis3$race[analysis3$race==97] <- 4

analysis3$race[analysis3$race==98] <- NA

analysis3$race[analysis3$race==99] <- NA

analysis3$race[analysis3$MB028==1] <- 3

analysis3$race[analysis3$MB089M1M==1 & analysis3$MB089M2M==2] <- 5

analysis3$race[analysis3$MB089M1M==2 & analysis3$MB089M2M==1] <- 5

analysis3$race[analysis3$MB089M1M==1 & analysis3$MB089M2M==97] <- 5

analysis3$race[analysis3$MB089M1M==97 & analysis3$MB089M2M==1] <- 5

analysis3$race[analysis3$MB089M1M==1 & analysis3$MB089M3M==2] <- 5

analysis3$race[analysis3$MB089M1M==2 & analysis3$MB089M3M==1] <- 5


## Frequency table
w = table(analysis3$race)

as.data.frame(w)


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4$race <- analysis4$NB089M1M

analysis4$race[analysis4$race==97] <- 4

analysis4$race[analysis4$race==98] <- NA

analysis4$race[analysis4$race==99] <- NA

analysis4$race[analysis4$NB028==1] <- 3

analysis4$race[analysis4$NB089M1M==1 & analysis4$NB089M2M==2] <- 5

analysis4$race[analysis4$NB089M1M==2 & analysis4$NB089M2M==1] <- 5

analysis4$race[analysis4$NB089M1M==1 & analysis4$NB089M2M==97] <- 5

analysis4$race[analysis4$NB089M1M==97 & analysis4$NB089M2M==1] <- 5

analysis4$race[analysis4$NB089M1M==1 & analysis4$NB089M3M==2] <- 5

analysis4$race[analysis4$NB089M1M==2 & analysis4$NB089M3M==1] <- 5


## Frequency table
w = table(analysis4$race)

as.data.frame(w)


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5$race <- analysis5$OB089M1M

analysis5$race[analysis5$race==97] <- 4

analysis5$race[analysis5$race==98] <- NA

analysis5$race[analysis5$race==99] <- NA

analysis5$race[analysis5$OB028==1] <- 3

analysis5$race[analysis5$OB089M1M==1 & analysis5$OB089M2M==2] <- 5

analysis5$race[analysis5$OB089M1M==2 & analysis5$OB089M2M==1] <- 5

analysis5$race[analysis5$OB089M1M==1 & analysis5$OB089M2M==97] <- 5

analysis5$race[analysis5$OB089M1M==97 & analysis5$OB089M2M==1] <- 5

analysis5$race[analysis5$OB089M1M==1 & analysis5$OB089M3M==2] <- 5

analysis5$race[analysis5$OB089M1M==2 & analysis5$OB089M3M==1] <- 5


## Frequency table
w = table(analysis5$race)

as.data.frame(w)


#### 2016 baseline covariates and 2018 stroke incidence
analysis6$race <- analysis6$PB089M1M

analysis6$race[analysis6$race==97] <- 4

analysis6$race[analysis6$race==98] <- NA

analysis6$race[analysis6$race==99] <- NA

analysis6$race[analysis6$PB028==1] <- 3

analysis6$race[analysis6$PB089M1M==1 & analysis6$PB089M2M==2] <- 5

analysis6$race[analysis6$PB089M1M==2 & analysis6$PB089M2M==1] <- 5

analysis6$race[analysis6$PB089M1M==1 & analysis6$PB089M2M==97] <- 5

analysis6$race[analysis6$PB089M1M==97 & analysis6$PB089M2M==1] <- 5

analysis6$race[analysis6$PB089M1M==1 & analysis6$PB089M3M==2] <- 5

analysis6$race[analysis6$PB089M1M==2 & analysis6$PB089M3M==1] <- 5


## Frequency table
w = table(analysis6$race)

as.data.frame(w)



############# EDUCATION (1- LESS THAN HIGH SCHOOL, 2- HIGH SCHOOL, 3- COLLEGE AND ABOVE)

#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1$KIMeducation <- analysis1$KZ216

analysis1$KIMeducation[analysis1$KIMeducation <12 ] <- 1

analysis1$KIMeducation[analysis1$KIMeducation >11 & analysis1$KIMeducation < 15] <- 2

analysis1$KIMeducation[analysis1$KIMeducation >14 ] <- 3

## Frequency table
w = table(analysis1$KIMeducation)

as.data.frame(w)


#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2$KIMeducation <- analysis2$LZ216

analysis2$KIMeducation[analysis2$KIMeducation <12 ] <- 1

analysis2$KIMeducation[analysis2$KIMeducation >11 & analysis2$KIMeducation < 15] <- 2

analysis2$KIMeducation[analysis2$KIMeducation >14 ] <- 3


## Frequency table
w = table(analysis2$KIMeducation)

as.data.frame(w)


#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3$KIMeducation <- analysis3$MZ216

analysis3$KIMeducation[analysis3$KIMeducation == 99] <- NA

analysis3$KIMeducation[analysis3$KIMeducation <12 ] <- 1

analysis3$KIMeducation[analysis3$KIMeducation >11 & analysis3$KIMeducation < 15] <- 2

analysis3$KIMeducation[analysis3$KIMeducation >14 ] <- 3


## Frequency table
w = table(analysis3$KIMeducation)

as.data.frame(w)


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4$KIMeducation <- analysis4$NZ216

analysis4$KIMeducation[analysis4$KIMeducation <12 ] <- 1

analysis4$KIMeducation[analysis4$KIMeducation >11 & analysis4$KIMeducation < 15] <- 2

analysis4$KIMeducation[analysis4$KIMeducation >14 ] <- 3


## Frequency table
w = table(analysis4$KIMeducation)

as.data.frame(w)


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5$KIMeducation <- analysis5$OZ216

analysis5$KIMeducation[analysis5$KIMeducation <12 ] <- 1

analysis5$KIMeducation[analysis5$KIMeducation >11 & analysis5$KIMeducation < 15] <- 2

analysis5$KIMeducation[analysis5$KIMeducation >14 ] <- 3


## Frequency table
w = table(analysis5$KIMeducation)

as.data.frame(w)


#### 2016 baseline covariates and 2018 stroke incidence
analysis6$KIMeducation <- analysis6$PZ216

analysis6$KIMeducation[analysis6$KIMeducation <12 ] <- 1

analysis6$KIMeducation[analysis6$KIMeducation >11 & analysis6$KIMeducation < 15] <- 2

analysis6$KIMeducation[analysis6$KIMeducation >14 ] <- 3


## Frequency table
w = table(analysis6$KIMeducation)

as.data.frame(w)


############# EDUCATION- YEARS OF EDUCATION
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1$education <- analysis1$KZ216

analysis1$education[analysis1$education==97] <- NA

## Frequency table
w = table(analysis1$education)

as.data.frame(w)


#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2$education <- analysis2$LZ216

## Frequency table
w = table(analysis2$education)

as.data.frame(w)


#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3$education <- analysis3$MZ216

analysis3$education[analysis3$education==99] <- NA

## Frequency table
w = table(analysis3$education)

as.data.frame(w)


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4$education <- analysis4$NZ216

## Frequency table
w = table(analysis4$education)

as.data.frame(w)


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5$education <- analysis5$OZ216

## Frequency table
w = table(analysis5$education)

as.data.frame(w)


#### 2016 baseline covariates and 2018 stroke incidence
analysis6$education <- analysis6$PZ216

## Frequency table
w = table(analysis6$education)

as.data.frame(w)



############# TOTAL WEALTH- CONTINUOUS
library(raster)
library(dplyr)
library(ryouready)


## Coerce from integer to numeric
## test[47:144] <- lapply(test[47:144], as.numeric)

##### Missing wealth values
wealth.na <- "99999998= NA; 99999999=NA; 999998=NA; 999999=NA; 9999998=NA; 9999999=NA; 9998=NA; 9999=NA; 998=NA; 999=NA"


### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1 <- recode2(analysis1, vars=46:143, recodes=wealth.na)

analysis1$income <- rowSums(analysis1 [, c(46:72, 75:76, 78:79, 82, 84:116, 119:135, 138:139, 142:143)], na.rm=TRUE)

analysis1$expenses <- rowSums(analysis1[, c(73, 74, 77, 80, 81, 83, 117, 118, 136, 137, 140)], na.rm=TRUE)

analysis1$wealth <- analysis1$income - analysis1$expenses


#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2 <- recode2(analysis2, vars=47:144, recodes=wealth.na)

analysis2$income <- rowSums(analysis2 [, c(47:73, 76,77, 79, 80, 83, 85:117, 120:136, 139:140, 142:144)], na.rm=TRUE)

analysis2$expenses <- rowSums(analysis2[, c(74, 75, 78, 82, 84, 118, 119, 137, 138, 141)], na.rm=TRUE)

analysis2$wealth <- analysis2$income - analysis2$expenses


#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3 <- recode2(analysis3, vars=50:146, recodes=wealth.na)

analysis3$income <- rowSums(analysis3 [, c(50:76, 79, 80, 82, 83, 86, 88:114)], na.rm=TRUE)

analysis3$expenses <- rowSums(analysis3[, c(77, 78, 81, 84, 85, 87, 121, 122, 139, 140, 143)], na.rm=TRUE)

analysis3$wealth <- analysis3$income - analysis3$expenses





#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4 <- recode2(analysis4, vars=45:141, recodes=wealth.na)

analysis4$income <- rowSums(analysis4 [, c(45:71, 74, 77, 78, 81, 83:109)], na.rm=TRUE)

analysis4$expenses <- rowSums(analysis4[, c(72, 73, 76, 79, 80, 82, 116, 117, 134, 135, 138)], na.rm=TRUE)

analysis4$wealth <- analysis4$income - analysis4$expenses


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5 <- recode2(analysis5, vars=45:137, recodes=wealth.na)

analysis5.names <- names(analysis5)

grepl(pattern= "OQ478", analysis5.names)

analysis5$income <- rowSums(analysis5 [, c(45:71, 74:75, 77:78, 81, 83:108)], na.rm=TRUE)

analysis5$expenses <- rowSums(analysis5[, c(72, 73, 76, 79, 80, 82, 115, 116, 131, 134)], na.rm=TRUE)

analysis5$wealth <- analysis5$income - analysis5$expenses


#### 2016 baseline covariates and 2018 stroke incidence
analysis6 <- recode2(analysis6, vars=45:137, recodes=wealth.na)

analysis6.names <- names(analysis6)

grepl(pattern= "PQ478", analysis6.names)

analysis6$income <- rowSums(analysis6 [, c(45:71, 74:75, 77:78, 81, 83:108)], na.rm=TRUE)

analysis6$expenses <- rowSums(analysis6[, c(72, 73, 76, 79, 80, 82, 115, 116, 131, 134)], na.rm=TRUE)

analysis6$wealth <- analysis6$income - analysis6$expenses



############# TOTAL WEALTH- QUINTILES

#### 2006 baseline covariates, 2008 and 2010 stroke incidence
quantile(analysis1$wealth)

analysis1$KIMwealth[analysis1$wealth <9495] <- 1

analysis1$KIMwealth[analysis1$wealth >9494 & analysis1$wealth <72325 ] <- 2

analysis1$KIMwealth[analysis1$wealth >72324 & analysis1$wealth <275056 ] <- 3

analysis1$KIMwealth[analysis1$wealth >275055 & analysis1$wealth <82122531 ] <- 4


#### 2008 baseline covariates, 2010 and 2012 stroke incidence
quantile(analysis2$wealth)

analysis2$KIMwealth[analysis2$wealth <154100] <- 1

analysis2$KIMwealth[analysis2$wealth >154099 & analysis2$wealth <20002007] <- 2

analysis2$KIMwealth[analysis2$wealth >20002006 & analysis2$wealth <120009607] <- 3

analysis2$KIMwealth[analysis2$wealth >120009606 & analysis2$wealth <1410002168] <- 4

  
#### 2010 baseline covariates, 2012 and 2014 stroke incidence
quantile(analysis3$wealth)

analysis3$KIMwealth[analysis3$wealth <2800] <- 1

analysis3$KIMwealth[analysis3$wealth >2799 & analysis3$wealth <50825] <- 2

analysis3$KIMwealth[analysis3$wealth >50824 & analysis3$wealth <172975] <- 3

analysis3$KIMwealth[analysis3$wealth >172974 & analysis3$wealth <13423733] <- 4



## Frequency table
w = table(analysis3$KIMwealth)

as.data.frame(w)


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
quantile(analysis4$wealth)

analysis4$KIMwealth[analysis4$wealth <1539 ] <- 1

analysis4$KIMwealth[analysis4$wealth >1538 & analysis4$wealth <45540] <- 2

analysis4$KIMwealth[analysis4$wealth >45539 & analysis4$wealth <189070] <- 3

analysis4$KIMwealth[analysis4$wealth >189069 & analysis4$wealth <51995357] <- 4



## Frequency table
w = table(analysis4$KIMwealth)

as.data.frame(w)


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
quantile(analysis5$wealth)

analysis5$KIMwealth[analysis5$wealth <2475 ] <- 1

analysis5$KIMwealth[analysis5$wealth >2474 & analysis5$wealth <66360] <- 2

analysis5$KIMwealth[analysis5$wealth >66359 & analysis5$wealth <342840] <- 3

analysis5$KIMwealth[analysis5$wealth >342839 & analysis5$wealth <1004103159] <- 4



## Frequency table
w = table(analysis5$KIMwealth)

as.data.frame(w)


#### 2016 baseline covariates and 2018 stroke incidence
quantile(analysis6$wealth)

analysis6$KIMwealth[analysis6$wealth <2301.5 ] <- 1

analysis6$KIMwealth[analysis6$wealth >2301.4 & analysis6$wealth <70402.5] <- 2

analysis6$KIMwealth[analysis6$wealth >70402.4 & analysis6$wealth < 232769] <- 3

analysis6$KIMwealth[analysis6$wealth >232768 & analysis6$wealth <30000501] <- 4



## Frequency table
w = table(analysis6$KIMwealth)

as.data.frame(w)




############# FUNCTIONAL STATUS (0- NO DIFFICULTY, 1- HAS DIFFICULTY, 1- CAN'T DO)

#### 2006 baseline covariates, 2008 and 2010 stroke incidence
# KG014
analysis1$func014[analysis1$KG014==1] <- 1

analysis1$func014[analysis1$KG014==5] <- 0

analysis1$func014[analysis1$KG014==6] <- 1

# KG016
analysis1$func016[analysis1$KG016==1] <- 1

analysis1$func016[analysis1$KG016==5] <- 0

analysis1$func016[analysis1$KG016==6] <- 1


# KG021
analysis1$func021[analysis1$KG021==1] <- 1

analysis1$func021[analysis1$KG021==5] <- 0

analysis1$func021[analysis1$KG021==6] <- 1


# KG023
analysis1$func023[analysis1$KG023==1] <- 1

analysis1$func023[analysis1$KG023==5] <- 0

analysis1$func023[analysis1$KG023==6] <- 1


# KG025
analysis1$func025[analysis1$KG025==1] <- 1

analysis1$func025[analysis1$KG025==5] <- 0

analysis1$func025[analysis1$KG025==6] <- 1



## Frequency table
w = table(analysis1$func014)

as.data.frame(w)



#### 2008 baseline covariates, 2010 and 2012 stroke incidence
# LG014
analysis2$func014[analysis2$LG014==1] <- 1

analysis2$func014[analysis2$LG014==5] <- 0

analysis2$func014[analysis2$LG014==6] <- 1

# LG016
analysis2$func016[analysis2$LG016==1] <- 1

analysis2$func016[analysis2$LG016==5] <- 0

analysis2$func016[analysis2$LG016==6] <- 1


# LG021
analysis2$func021[analysis2$LG021==1] <- 1

analysis2$func021[analysis2$LG021==5] <- 0

analysis2$func021[analysis2$LG021==6] <- 1


# LG023
analysis2$func023[analysis2$LG023==1] <- 1

analysis2$func023[analysis2$LG023==5] <- 0

analysis2$func023[analysis2$LG023==6] <- 1


# LG025
analysis2$func025[analysis2$LG025==1] <- 1

analysis2$func025[analysis2$LG025==5] <- 0

analysis2$func025[analysis2$LG025==6] <- 1



## Frequency table
w = table(analysis2$func025)

as.data.frame(w)




#### 2010 baseline covariates, 2012 and 2014 stroke incidence
# MG014
analysis3$func014[analysis3$MG014==1] <- 1

analysis3$func014[analysis3$MG014==5] <- 0

analysis3$func014[analysis3$MG014==6] <- 1

# MG016
analysis3$func016[analysis3$MG016==1] <- 1

analysis3$func016[analysis3$MG016==5] <- 0

analysis3$func016[analysis3$MG016==6] <- 1


# MG021
analysis3$func021[analysis3$MG021==1] <- 1

analysis3$func021[analysis3$MG021==5] <- 0

analysis3$func021[analysis3$MG021==6] <- 1


# MG023
analysis3$func023[analysis3$MG023==1] <- 1

analysis3$func023[analysis3$MG023==5] <- 0

analysis3$func023[analysis3$MG023==6] <- 1


# MG025
analysis3$func025[analysis3$MG025==1] <- 1

analysis3$func025[analysis3$MG025==5] <- 0

analysis3$func025[analysis3$MG025==6] <- 1



## Frequency table
w = table(analysis3$func025)

as.data.frame(w)


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
# NG014
analysis4$func014[analysis4$NG014==1] <- 1

analysis4$func014[analysis4$NG014==5] <- 0

analysis4$func014[analysis4$NG014==6] <- 1

## Frequency table
w = table(analysis4$func014)

as.data.frame(w)


# NG016
analysis4$func016[analysis4$NG016==1] <- 1

analysis4$func016[analysis4$NG016==5] <- 0

analysis4$func016[analysis4$NG016==6] <- 1

## Frequency table
w = table(analysis4$func016)

as.data.frame(w)


# NG021
analysis4$func021[analysis4$NG021==1] <- 1

analysis4$func021[analysis4$NG021==5] <- 0

analysis4$func021[analysis4$NG021==6] <- 1


## Frequency table
w = table(analysis4$func021)

as.data.frame(w)


# NG023
analysis4$func023[analysis4$NG023==1] <- 1

analysis4$func023[analysis4$NG023==5] <- 0

analysis4$func023[analysis4$NG023==6] <- 1

## Frequency table
w = table(analysis4$func023)

as.data.frame(w)


# NG025
analysis4$func025[analysis4$NG025==1] <- 1

analysis4$func025[analysis4$NG025==5] <- 0

analysis4$func025[analysis4$NG025==6] <- 1


## Frequency table
w = table(analysis4$func025)

as.data.frame(w)


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
# OG014
analysis5$func014[analysis5$OG014==1] <- 1

analysis5$func014[analysis5$OG014==5] <- 0

analysis5$func014[analysis5$OG014==6] <- 1

## Frequency table
w = table(analysis5$func014)

as.data.frame(w)


# OG016
analysis5$func016[analysis5$OG016==1] <- 1

analysis5$func016[analysis5$OG016==5] <- 0

analysis5$func016[analysis5$OG016==6] <- 1

## Frequency table
w = table(analysis5$func016)

as.data.frame(w)


# OG021
analysis5$func021[analysis5$OG021==1] <- 1

analysis5$func021[analysis5$OG021==5] <- 0

analysis5$func021[analysis5$OG021==6] <- 1


## Frequency table
w = table(analysis5$func021)

as.data.frame(w)


# OG023
analysis5$func023[analysis5$OG023==1] <- 1

analysis5$func023[analysis5$OG023==5] <- 0

analysis5$func023[analysis5$OG023==6] <- 1

## Frequency table
w = table(analysis5$OG023)

as.data.frame(w)


# OG025
analysis5$func025[analysis5$OG025==1] <- 1

analysis5$func025[analysis5$OG025==5] <- 0

analysis5$func025[analysis5$OG025==6] <- 1


## Frequency table
w = table(analysis5$func025)

as.data.frame(w)


#### 2016 baseline covariates and 2018 stroke incidence
# PG014
analysis6$func014[analysis6$PG014==1] <- 1

analysis6$func014[analysis6$PG014==5] <- 0

analysis6$func014[analysis6$PG014==6] <- 1

## Frequency table
w = table(analysis6$func014)

as.data.frame(w)


# PG016
analysis6$func016[analysis6$PG016==1] <- 1

analysis6$func016[analysis6$PG016==5] <- 0

analysis6$func016[analysis6$PG016==6] <- 1

## Frequency table
w = table(analysis6$func016)

as.data.frame(w)


# PG021
analysis6$func021[analysis6$PG021==1] <- 1

analysis6$func021[analysis6$PG021==5] <- 0

analysis6$func021[analysis6$PG021==6] <- 1


## Frequency table
w = table(analysis6$func021)

as.data.frame(w)


# PG023
analysis6$func023[analysis6$PG023==1] <- 1

analysis6$func023[analysis6$PG023==5] <- 0

analysis6$func023[analysis6$PG023==6] <- 1


## Frequency table
w = table(analysis6$func023)

as.data.frame(w)


# PG025
analysis6$func025[analysis6$PG025==1] <- 1

analysis6$func025[analysis6$PG025==5] <- 0

analysis6$func025[analysis6$PG025==6] <- 1



## Frequency table
w = table(analysis6$func025)

as.data.frame(w)


#### MEAN FUNCTIONAL STATUS
analysis1$func <- rowMeans(as.data.frame(analysis1$func014 + analysis1$func016 + analysis1$func021 + analysis1$func023 + analysis1$func025), na.rm=FALSE)

analysis2$func <- rowMeans(as.data.frame(analysis2$func014 + analysis2$func016 + analysis2$func021 + analysis2$func023 + analysis2$func025), na.rm=FALSE)

analysis3$func <- rowMeans(as.data.frame(analysis3$func014 + analysis3$func016 + analysis3$func021 + analysis3$func023 + analysis3$func025), na.rm=FALSE)

analysis4$func <- rowMeans(as.data.frame(analysis4$func014 + analysis4$func016 + analysis4$func021 + analysis4$func023 + analysis4$func025), na.rm=FALSE)

analysis5$func <- rowMeans(as.data.frame(analysis5$func014 + analysis5$func016 + analysis5$func021 + analysis5$func023 + analysis5$func025), na.rm=FALSE)

analysis6$func <- rowMeans(as.data.frame(analysis6$func014 + analysis6$func016 + analysis6$func021 + analysis6$func023 + analysis6$func025), na.rm=FALSE)


## Frequency table
w = table(analysis6$func)

as.data.frame(w)



############# SMOKING STATUS (1- CURRENTLY SMOKES, 0- DOES NOT SMOKE)

#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1$KIMsmoke[analysis1$KC117==1] <- 1

analysis1$KIMsmoke[analysis1$KC117==5] <- 0

#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2$KIMsmoke[analysis2$LC117==1] <- 1

analysis2$KIMsmoke[analysis2$LC117==5] <- 0

#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3$KIMsmoke[analysis3$MC117==1] <- 1

analysis3$KIMsmoke[analysis3$MC117==5] <- 0


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4$KIMsmoke[analysis4$NC117==1] <- 1

analysis4$KIMsmoke[analysis4$NC117==5] <- 0



#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5$KIMsmoke[analysis5$OC117==1] <- 1

analysis5$KIMsmoke[analysis5$OC117==5] <- 0


#### 2016 baseline covariates and 2018 stroke incidence
analysis6$KIMsmoke[analysis6$PC117==1] <- 1

analysis6$KIMsmoke[analysis6$PC117==5] <- 0



############# NUMBER OF CIGARETTES SMOKED
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1$numCigs <- analysis1$KC118

analysis1$numCigs[analysis1$numCigs==998] <- NA

analysis1$numCigs[analysis1$numCigs==999] <- NA


## Frequency table
w = table(analysis1$numCigs)

as.data.frame(w)



#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2$numCigs <- analysis2$LC118

analysis2$numCigs[analysis2$numCigs==998] <- NA

analysis2$numCigs[analysis2$numCigs==999] <- NA


#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3$numCigs <- analysis3$MC118

analysis3$numCigs[analysis3$numCigs==998] <- NA

analysis3$numCigs[analysis3$numCigs==999] <- NA


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4$numCigs <- analysis4$NC118

analysis4$numCigs[analysis4$numCigs==998] <- NA

analysis4$numCigs[analysis4$numCigs==999] <- NA


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5$numCigs <- analysis5$OC118

analysis5$numCigs[analysis5$numCigs==998] <- NA

analysis5$numCigs[analysis5$numCigs==999] <- NA


#### 2016 baseline covariates and 2018 stroke incidence
analysis6$numCigs <- analysis6$PC118

analysis6$numCigs[analysis6$numCigs==998] <- NA

analysis6$numCigs[analysis6$numCigs==999] <- NA




############# EXERCISE- KIM ET AL

#### 2006 baseline covariates, 2008 and 2010 stroke incidence
# Never exercise
analysis1$KIMneverexercise[analysis1$KC223==4] <- 1

analysis1$KIMneverexercise[analysis1$KC223==1] <- 0

analysis1$KIMneverexercise[analysis1$KC223==2] <- 0

analysis1$KIMneverexercise[analysis1$KC223==3] <- 0

analysis1$KIMneverexercise[analysis1$KC223==7] <- 0

analysis1$KIMneverexercise[analysis1$KC223==8] <- NA

analysis1$KIMneverexercise[analysis1$KC223==9] <- NA

## Frequency table
w = table(analysis1$KIMneverexercise)

as.data.frame(w)

# Low exercise
analysis1$KIMlowexercise[analysis1$KC225==3] <- 1

analysis1$KIMlowexercise[analysis1$KC225==1] <- 0

analysis1$KIMlowexercise[analysis1$KC225==2] <- 0

analysis1$KIMlowexercise[analysis1$KC225==4] <- 0

analysis1$KIMlowexercise[analysis1$KC225==7] <- 0


## Frequency table
w = table(analysis1$KIMlowexercise)

as.data.frame(w)

# Moderate exercise
analysis1$KIMmodexercise[analysis1$KC224==1] <- 1

analysis1$KIMmodexercise[analysis1$KC224==2] <- 0

analysis1$KIMmodexercise[analysis1$KC224==3] <- 0

analysis1$KIMmodexercise[analysis1$KC224==4] <- 0

analysis1$KIMmodexercise[analysis1$KC224==7] <- 0



## Frequency table
w = table(analysis1$KIMmodexercise)

as.data.frame(w)

# High exercise
analysis1$KIMhighexercise[analysis1$KC223==7] <- 1

analysis1$KIMhighexercise[analysis1$KC223==1] <- 0

analysis1$KIMhighexercise[analysis1$KC223==2] <- 0

analysis1$KIMhighexercise[analysis1$KC223==3] <- 0

analysis1$KIMhighexercise[analysis1$KC223==4] <- 0



## Frequency table
w = table(analysis1$KIMhighexercise)

as.data.frame(w)




#### 2008 baseline covariates, 2010 and 2012 stroke incidence
# Never exercise
analysis2$KIMneverexercise[analysis2$LC223==4] <- 1

analysis2$KIMneverexercise[analysis2$LC223==1] <- 0

analysis2$KIMneverexercise[analysis2$LC223==2] <- 0

analysis2$KIMneverexercise[analysis2$LC223==3] <- 0

analysis2$KIMneverexercise[analysis2$LC223==7] <- 0

analysis2$KIMneverexercise[analysis2$LC223==8] <- NA

analysis2$KIMneverexercise[analysis2$LC223==9] <- NA

## Frequency table
w = table(analysis2$KIMneverexercise)

as.data.frame(w)

# Low exercise
analysis2$KIMlowexercise[analysis2$LC225==3] <- 1

analysis2$KIMlowexercise[analysis2$LC225==1] <- 0

analysis2$KIMlowexercise[analysis2$LC225==2] <- 0

analysis2$KIMlowexercise[analysis2$LC225==4] <- 0

analysis2$KIMlowexercise[analysis2$LC225==7] <- 0


## Frequency table
w = table(analysis2$KIMlowexercise)

as.data.frame(w)

# Moderate exercise
analysis2$KIMmodexercise[analysis2$LC224==1] <- 1

analysis2$KIMmodexercise[analysis2$LC224==2] <- 0

analysis2$KIMmodexercise[analysis2$LC224==3] <- 0

analysis2$KIMmodexercise[analysis2$LC224==4] <- 0

analysis2$KIMmodexercise[analysis2$LC224==7] <- 0



## Frequency table
w = table(analysis2$KIMmodexercise)

as.data.frame(w)

# High exercise
analysis2$KIMhighexercise[analysis2$LC223==7] <- 1

analysis2$KIMhighexercise[analysis2$LC223==1] <- 0

analysis2$KIMhighexercise[analysis2$LC223==2] <- 0

analysis2$KIMhighexercise[analysis2$LC223==3] <- 0

analysis2$KIMhighexercise[analysis2$LC223==4] <- 0



## Frequency table
w = table(analysis2$KIMhighexercise)

as.data.frame(w)



#### 2010 baseline covariates, 2012 and 2014 stroke incidence
# Never exercise
analysis3$KIMneverexercise[analysis3$MC223==4] <- 1

analysis3$KIMneverexercise[analysis3$MC223==1] <- 0

analysis3$KIMneverexercise[analysis3$MC223==2] <- 0

analysis3$KIMneverexercise[analysis3$MC223==3] <- 0

analysis3$KIMneverexercise[analysis3$MC223==7] <- 0

analysis3$KIMneverexercise[analysis3$MC223==8] <- NA

analysis3$KIMneverexercise[analysis3$MC223==9] <- NA

## Frequency table
w = table(analysis3$KIMneverexercise)

as.data.frame(w)

# Low exercise
analysis3$KIMlowexercise[analysis3$MC225==3] <- 1

analysis3$KIMlowexercise[analysis3$MC225==1] <- 0

analysis3$KIMlowexercise[analysis3$MC225==2] <- 0

analysis3$KIMlowexercise[analysis3$MC225==4] <- 0

analysis3$KIMlowexercise[analysis3$MC225==7] <- 0


## Frequency table
w = table(analysis3$KIMlowexercise)

as.data.frame(w)

# Moderate exercise
analysis3$KIMmodexercise[analysis3$MC224==1] <- 1

analysis3$KIMmodexercise[analysis3$MC224==2] <- 0

analysis3$KIMmodexercise[analysis3$MC224==3] <- 0

analysis3$KIMmodexercise[analysis3$MC224==4] <- 0

analysis3$KIMmodexercise[analysis3$MC224==7] <- 0



## Frequency table
w = table(analysis3$KIMmodexercise)

as.data.frame(w)

# High exercise
analysis3$KIMhighexercise[analysis3$MC223==7] <- 1

analysis3$KIMhighexercise[analysis3$MC223==1] <- 0

analysis3$KIMhighexercise[analysis3$MC223==2] <- 0

analysis3$KIMhighexercise[analysis3$MC223==3] <- 0

analysis3$KIMhighexercise[analysis3$MC223==4] <- 0



## Frequency table
w = table(analysis3$KIMhighexercise)

as.data.frame(w)


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
# Never exercise
analysis4$KIMneverexercise[analysis4$NC223==4] <- 1

analysis4$KIMneverexercise[analysis4$NC223==1] <- 0

analysis4$KIMneverexercise[analysis4$NC223==2] <- 0

analysis4$KIMneverexercise[analysis4$NC223==3] <- 0

analysis4$KIMneverexercise[analysis4$NC223==7] <- 0

analysis4$KIMneverexercise[analysis4$NC223==8] <- NA

analysis4$KIMneverexercise[analysis4$NC223==9] <- NA

## Frequency table
w = table(analysis4$KIMneverexercise)

as.data.frame(w)

# Low exercise
analysis4$KIMlowexercise[analysis4$NC225==3] <- 1

analysis4$KIMlowexercise[analysis4$NC225==1] <- 0

analysis4$KIMlowexercise[analysis4$NC225==2] <- 0

analysis4$KIMlowexercise[analysis4$NC225==4] <- 0

analysis4$KIMlowexercise[analysis4$NC225==7] <- 0


## Frequency table
w = table(analysis4$KIMlowexercise)

as.data.frame(w)

# Moderate exercise
analysis4$KIMmodexercise[analysis4$NC224==1] <- 1

analysis4$KIMmodexercise[analysis4$NC224==2] <- 0

analysis4$KIMmodexercise[analysis4$NC224==3] <- 0

analysis4$KIMmodexercise[analysis4$NC224==4] <- 0

analysis4$KIMmodexercise[analysis4$NC224==7] <- 0



## Frequency table
w = table(analysis4$KIMmodexercise)

as.data.frame(w)

# High exercise
analysis4$KIMhighexercise[analysis4$NC223==7] <- 1

analysis4$KIMhighexercise[analysis4$NC223==1] <- 0

analysis4$KIMhighexercise[analysis4$NC223==2] <- 0

analysis4$KIMhighexercise[analysis4$NC223==3] <- 0

analysis4$KIMhighexercise[analysis4$NC223==4] <- 0



## Frequency table
w = table(analysis4$KIMhighexercise)

as.data.frame(w)


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
# Never exercise
analysis5$KIMneverexercise[analysis5$OC223==4] <- 1

analysis5$KIMneverexercise[analysis5$OC223==1] <- 0

analysis5$KIMneverexercise[analysis5$OC223==2] <- 0

analysis5$KIMneverexercise[analysis5$OC223==3] <- 0

analysis5$KIMneverexercise[analysis5$OC223==7] <- 0

analysis5$KIMneverexercise[analysis5$OC223==8] <- NA

analysis5$KIMneverexercise[analysis5$OC223==9] <- NA

## Frequency table
w = table(analysis5$KIMneverexercise)

as.data.frame(w)

# Low exercise
analysis5$KIMlowexercise[analysis5$OC225==3] <- 1

analysis5$KIMlowexercise[analysis5$OC225==1] <- 0

analysis5$KIMlowexercise[analysis5$OC225==2] <- 0

analysis5$KIMlowexercise[analysis5$OC225==4] <- 0

analysis5$KIMlowexercise[analysis5$OC225==7] <- 0


## Frequency table
w = table(analysis5$KIMlowexercise)

as.data.frame(w)

# Moderate exercise
analysis5$KIMmodexercise[analysis5$OC224==1] <- 1

analysis5$KIMmodexercise[analysis5$OC224==2] <- 0

analysis5$KIMmodexercise[analysis5$OC224==3] <- 0

analysis5$KIMmodexercise[analysis5$OC224==4] <- 0

analysis5$KIMmodexercise[analysis5$OC224==7] <- 0



## Frequency table
w = table(analysis5$KIMmodexercise)

as.data.frame(w)


# High exercise
analysis5$KIMhighexercise[analysis5$OC223==7] <- 1

analysis5$KIMhighexercise[analysis5$OC223==1] <- 0

analysis5$KIMhighexercise[analysis5$OC223==2] <- 0

analysis5$KIMhighexercise[analysis5$OC223==3] <- 0

analysis5$KIMhighexercise[analysis5$OC223==4] <- 0



## Frequency table
w = table(analysis5$KIMhighexercise)

as.data.frame(w)



#### 2016 baseline covariates and 2018 stroke incidence
# Never exercise
analysis6$KIMneverexercise[analysis6$PC223==4] <- 1

analysis6$KIMneverexercise[analysis6$PC223==1] <- 0

analysis6$KIMneverexercise[analysis6$PC223==2] <- 0

analysis6$KIMneverexercise[analysis6$PC223==3] <- 0

analysis6$KIMneverexercise[analysis6$PC223==7] <- 0

analysis6$KIMneverexercise[analysis6$PC223==8] <- NA

analysis6$KIMneverexercise[analysis6$PC223==9] <- NA

## Frequency table
w = table(analysis6$KIMneverexercise)

as.data.frame(w)


# Low exercise
analysis6$KIMlowexercise[analysis6$PC225==3] <- 1

analysis6$KIMlowexercise[analysis6$PC225==1] <- 0

analysis6$KIMlowexercise[analysis6$PC225==2] <- 0

analysis6$KIMlowexercise[analysis6$PC225==4] <- 0

analysis6$KIMlowexercise[analysis6$PC225==7] <- 0


## Frequency table
w = table(analysis6$KIMlowexercise)

as.data.frame(w)

# Moderate exercise
analysis6$KIMmodexercise[analysis6$PC224==1] <- 1

analysis6$KIMmodexercise[analysis6$PC224==2] <- 0

analysis6$KIMmodexercise[analysis6$PC224==3] <- 0

analysis6$KIMmodexercise[analysis6$PC224==4] <- 0

analysis6$KIMmodexercise[analysis6$PC224==7] <- 0



## Frequency table
w = table(analysis6$KIMmodexercise)

as.data.frame(w)


# High exercise
analysis6$KIMhighexercise[analysis6$PC223==7] <- 1

analysis6$KIMhighexercise[analysis6$PC223==1] <- 0

analysis6$KIMhighexercise[analysis6$PC223==2] <- 0

analysis6$KIMhighexercise[analysis6$PC223==3] <- 0

analysis6$KIMhighexercise[analysis6$PC223==4] <- 0



## Frequency table
w = table(analysis6$KIMhighexercise)

as.data.frame(w)


############# EXERCISE- OMITTED VALUES INCLUDED (1- HARDLY EVER OR NEVER, 2- ONE TO THREE TIMES A MONTH, 3- ONCE A WEEK, 4- MORE THAN ONCE A WEEK, 5- EVERYDAY)
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
# Mild exercise
analysis1$mildexercise[analysis1$KC225== 1] <- 4

analysis1$mildexercise[analysis1$KC225== 2] <- 3

analysis1$mildexercise[analysis1$KC225== 3] <- 2

analysis1$mildexercise[analysis1$KC225== 4] <- 1

analysis1$mildexercise[analysis1$KC225== 7] <- 5

## Frequency table
w = table(analysis1$mildexercise)

as.data.frame(w)


# Moderate activity
analysis1$modexercise[analysis1$KC224== 1] <- 4

analysis1$modexercise[analysis1$KC224== 2] <- 3

analysis1$modexercise[analysis1$KC224== 3] <- 2

analysis1$modexercise[analysis1$KC224== 4] <- 1

analysis1$modexercise[analysis1$KC224== 7] <- 5

## Frequency table
w = table(analysis1$modexercise)

as.data.frame(w)

# High activity
analysis1$highexercise[analysis1$KC223== 1] <- 4

analysis1$highexercise[analysis1$KC223== 2] <- 3

analysis1$highexercise[analysis1$KC223== 3] <- 2

analysis1$highexercise[analysis1$KC223== 4] <- 1

analysis1$highexercise[analysis1$KC223== 7] <- 5

## Frequency table
w = table(analysis1$highexercise)

as.data.frame(w)



#### 2008 baseline covariates, 2010 and 2012 stroke incidence
# Mild exercise
analysis2$mildexercise[analysis2$LC225== 1] <- 4

analysis2$mildexercise[analysis2$LC225== 2] <- 3

analysis2$mildexercise[analysis2$LC225== 3] <- 2

analysis2$mildexercise[analysis2$LC225== 4] <- 1

analysis2$mildexercise[analysis2$LC225== 7] <- 5

## Frequency table
w = table(analysis2$mildexercise)

as.data.frame(w)


# Moderate activity
analysis2$modexercise[analysis2$LC224== 1] <- 4

analysis2$modexercise[analysis2$LC224== 2] <- 3

analysis2$modexercise[analysis2$LC224== 3] <- 2

analysis2$modexercise[analysis2$LC224== 4] <- 1

analysis2$modexercise[analysis2$LC224== 7] <- 5

## Frequency table
w = table(analysis2$modexercise)

as.data.frame(w)

# High activity
analysis2$highexercise[analysis2$LC223== 1] <- 4

analysis2$highexercise[analysis2$LC223== 2] <- 3

analysis2$highexercise[analysis2$LC223== 3] <- 2

analysis2$highexercise[analysis2$LC223== 4] <- 1

analysis2$highexercise[analysis2$LC223== 7] <- 5

## Frequency table
w = table(analysis2$highexercise)

as.data.frame(w)



#### 2010 baseline covariates, 2012 and 2014 stroke incidence
# Mild exercise
analysis3$mildexercise[analysis3$MC225== 1] <- 4

analysis3$mildexercise[analysis3$MC225== 2] <- 3

analysis3$mildexercise[analysis3$MC225== 3] <- 2

analysis3$mildexercise[analysis3$MC225== 4] <- 1

analysis3$mildexercise[analysis3$MC225== 7] <- 5

## Frequency table
w = table(analysis3$mildexercise)

as.data.frame(w)


# Moderate activity
analysis3$modexercise[analysis3$MC224== 1] <- 4

analysis3$modexercise[analysis3$MC224== 2] <- 3

analysis3$modexercise[analysis3$MC224== 3] <- 2

analysis3$modexercise[analysis3$MC224== 4] <- 1

analysis3$modexercise[analysis3$MC224== 7] <- 5

## Frequency table
w = table(analysis3$modexercise)

as.data.frame(w)

# High activity
analysis3$highexercise[analysis3$MC223== 1] <- 4

analysis3$highexercise[analysis3$MC223== 2] <- 3

analysis3$highexercise[analysis3$MC223== 3] <- 2

analysis3$highexercise[analysis3$MC223== 4] <- 1

analysis3$highexercise[analysis3$MC223== 7] <- 5

## Frequency table
w = table(analysis3$highexercise)

as.data.frame(w)


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
# Mild exercise
analysis4$mildexercise[analysis4$NC225== 1] <- 4

analysis4$mildexercise[analysis4$NC225== 2] <- 3

analysis4$mildexercise[analysis4$NC225== 3] <- 2

analysis4$mildexercise[analysis4$NC225== 4] <- 1

analysis4$mildexercise[analysis4$NC225== 7] <- 5

## Frequency table
w = table(analysis4$mildexercise)

as.data.frame(w)


# Moderate activity
analysis4$modexercise[analysis4$NC224== 1] <- 4

analysis4$modexercise[analysis4$NC224== 2] <- 3

analysis4$modexercise[analysis4$NC224== 3] <- 2

analysis4$modexercise[analysis4$NC224== 4] <- 1

analysis4$modexercise[analysis4$NC224== 7] <- 5

## Frequency table
w = table(analysis4$modexercise)

as.data.frame(w)

# High activity
analysis4$highexercise[analysis4$NC223== 1] <- 4

analysis4$highexercise[analysis4$NC223== 2] <- 3

analysis4$highexercise[analysis4$NC223== 3] <- 2

analysis4$highexercise[analysis4$NC223== 4] <- 1

analysis4$highexercise[analysis4$NC223== 7] <- 5

## Frequency table
w = table(analysis4$highexercise)

as.data.frame(w)


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
# Mild exercise
analysis5$mildexercise[analysis5$OC225== 1] <- 4

analysis5$mildexercise[analysis5$OC225== 2] <- 3

analysis5$mildexercise[analysis5$OC225== 3] <- 2

analysis5$mildexercise[analysis5$OC225== 4] <- 1

analysis5$mildexercise[analysis5$OC225== 7] <- 5

## Frequency table
w = table(analysis5$mildexercise)

as.data.frame(w)


# Moderate activity
analysis5$modexercise[analysis5$OC224== 1] <- 4

analysis5$modexercise[analysis5$OC224== 2] <- 3

analysis5$modexercise[analysis5$OC224== 3] <- 2

analysis5$modexercise[analysis5$OC224== 4] <- 1

analysis5$modexercise[analysis5$OC224== 7] <- 5

## Frequency table
w = table(analysis5$modexercise)

as.data.frame(w)

# High activity
analysis5$highexercise[analysis5$OC223== 1] <- 4

analysis5$highexercise[analysis5$OC223== 2] <- 3

analysis5$highexercise[analysis5$OC223== 3] <- 2

analysis5$highexercise[analysis5$OC223== 4] <- 1

analysis5$highexercise[analysis5$OC223== 7] <- 5

## Frequency table
w = table(analysis5$highexercise)

as.data.frame(w)


#### 2016 baseline covariates and 2018 stroke incidence
# Mild exercise
analysis6$mildexercise[analysis6$PC225== 1] <- 4

analysis6$mildexercise[analysis6$PC225== 2] <- 3

analysis6$mildexercise[analysis6$PC225== 3] <- 2

analysis6$mildexercise[analysis6$PC225== 4] <- 1

analysis6$mildexercise[analysis6$PC225== 7] <- 5

## Frequency table
w = table(analysis6$mildexercise)

as.data.frame(w)


# Moderate activity
analysis6$modexercise[analysis6$PC224== 1] <- 4

analysis6$modexercise[analysis6$PC224== 2] <- 3

analysis6$modexercise[analysis6$PC224== 3] <- 2

analysis6$modexercise[analysis6$PC224== 4] <- 1

analysis6$modexercise[analysis6$PC224== 7] <- 5

## Frequency table
w = table(analysis6$modexercise)

as.data.frame(w)

# High activity
analysis6$highexercise[analysis6$PC223== 1] <- 4

analysis6$highexercise[analysis6$PC223== 2] <- 3

analysis6$highexercise[analysis6$PC223== 3] <- 2

analysis6$highexercise[analysis6$PC223== 4] <- 1

analysis6$highexercise[analysis6$PC223== 7] <- 5

## Frequency table
w = table(analysis6$highexercise)

as.data.frame(w)




############# ALCOHOL- KIM ET AL

#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1$KIMalcohol[analysis1$KC128==1]<- 1

analysis1$KIMalcohol[analysis1$KC128==5]<- 0

## Frequency table
w = table(analysis1$KIMalcohol)

as.data.frame(w)



#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2$KIMalcohol[analysis2$LC128==1]<- 1

analysis2$KIMalcohol[analysis2$LC128==5]<- 0

## Frequency table
w = table(analysis2$KIMalcohol)

as.data.frame(w)


#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3$KIMalcohol[analysis3$MC128==1]<- 1

analysis3$KIMalcohol[analysis3$MC128==5]<- 0


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4$KIMalcohol[analysis4$NC128==1]<- 1

analysis4$KIMalcohol[analysis4$NC128==5]<- 0

#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5$KIMalcohol[analysis5$OC128==1]<- 1

analysis5$KIMalcohol[analysis5$OC128==5]<- 0

#### 2016 baseline covariates and 2018 stroke incidence
analysis6$KIMalcohol[analysis6$PC128==1]<- 1

analysis6$KIMalcohol[analysis6$PC128==5]<- 0


############# ALCOHOL- NUMBER OF DAYS PER WEEK ALCOHOL CONSUMED

#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1$alcohol <- analysis1$KC129 

analysis1$alcohol[analysis1$alcohol==8] <- NA

analysis1$alcohol[analysis1$alcohol==9] <- NA

## Frequency table
w = table(analysis1$alcohol)

as.data.frame(w)

#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2$alcohol <- analysis2$LC129 

analysis2$alcohol[analysis2$alcohol==8] <- NA

analysis2$alcohol[analysis2$alcohol==9] <- NA

## Frequency table
w = table(analysis2$alcohol)

as.data.frame(w)

#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3$alcohol <- analysis3$MC129 

analysis3$alcohol[analysis3$alcohol==8] <- NA

analysis3$alcohol[analysis3$alcohol==9] <- NA


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4$alcohol <- analysis4$NC129 

analysis4$alcohol[analysis4$alcohol==8] <- NA

analysis4$alcohol[analysis4$alcohol==9] <- NA


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5$alcohol <- analysis5$OC129 

analysis5$alcohol[analysis5$alcohol==8] <- NA

analysis5$alcohol[analysis5$alcohol==9] <- NA


#### 2016 baseline covariates and 2018 stroke incidence
analysis6$alcohol <- analysis6$PC129 

analysis6$alcohol[analysis6$alcohol==8] <- NA

analysis6$alcohol[analysis6$alcohol==9] <- NA


############# HYPERTENSION
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1$hypertension[analysis1$KC005==1] <- 1

analysis1$hypertension[analysis1$KC005==3] <- 1

analysis1$hypertension[analysis1$KC005==5] <- 0

analysis1$hypertension[analysis1$KC005==4] <- 0

## Frequency table
w = table(analysis1$hypertension)

as.data.frame(w)


#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2$hypertension[analysis2$LC005==1] <- 1

analysis2$hypertension[analysis2$LC005==3] <- 1

analysis2$hypertension[analysis2$LC005==5] <- 0

analysis2$hypertension[analysis2$LC005==4] <- 0

## Frequency table
w = table(analysis2$hypertension)

as.data.frame(w)


#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3$hypertension[analysis3$MC005==1] <- 1

analysis3$hypertension[analysis3$MC005==3] <- 1

analysis3$hypertension[analysis3$MC005==5] <- 0

analysis3$hypertension[analysis3$MC005==4] <- 0


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4$hypertension[analysis4$NC005==1] <- 1

analysis4$hypertension[analysis4$NC005==3] <- 1

analysis4$hypertension[analysis4$NC005==5] <- 0

analysis4$hypertension[analysis4$NC005==4] <- 0



#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5$hypertension[analysis5$OC005==1] <- 1

analysis5$hypertension[analysis5$OC005==3] <- 1

analysis5$hypertension[analysis5$OC005==5] <- 0

analysis5$hypertension[analysis5$OC005==4] <- 0


#### 2016 baseline covariates and 2018 stroke incidence
analysis6$hypertension[analysis6$PC005==1] <- 1

analysis6$hypertension[analysis6$PC005==3] <- 1

analysis6$hypertension[analysis6$PC005==5] <- 0

analysis6$hypertension[analysis6$PC005==4] <- 0




############# DIABETES
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1$diabetes[analysis1$KC010==1] <- 1

analysis1$diabetes[analysis1$KC010==3] <- 1

analysis1$diabetes[analysis1$KC010==5] <- 0

analysis1$diabetes[analysis1$KC010==4] <- 0


## Frequency table
w = table(analysis1$diabetes)

as.data.frame(w)

#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2$diabetes[analysis2$LC010==1] <- 1

analysis2$diabetes[analysis2$LC010==3] <- 1

analysis2$diabetes[analysis2$LC010==5] <- 0

analysis2$diabetes[analysis2$LC010==4] <- 0


## Frequency table
w = table(analysis2$diabetes)

as.data.frame(w)

#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3$diabetes[analysis3$MC010==1] <- 1

analysis3$diabetes[analysis3$MC010==3] <- 1

analysis3$diabetes[analysis3$MC010==5] <- 0

analysis3$diabetes[analysis3$MC010==4] <- 0


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4$diabetes[analysis4$NC010==1] <- 1

analysis4$diabetes[analysis4$NC010==3] <- 1

analysis4$diabetes[analysis4$NC010==5] <- 0

analysis4$diabetes[analysis4$NC010==4] <- 0


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5$diabetes[analysis5$OC010==1] <- 1

analysis5$diabetes[analysis5$OC010==3] <- 1

analysis5$diabetes[analysis5$OC010==5] <- 0

analysis5$diabetes[analysis5$OC010==4] <- 0


#### 2016 baseline covariates and 2018 stroke incidence
analysis6$diabetes[analysis6$PC010==1] <- 1

analysis6$diabetes[analysis6$PC010==3] <- 1

analysis6$diabetes[analysis6$PC010==5] <- 0

analysis6$diabetes[analysis6$PC010==4] <- 0




############# SYSTOLIC BLOOD PRESSURE, KIM ET AL
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1$KIMsystolic <- analysis1$KI859

analysis1$KIMsystolic[analysis1$KIMsystolic==993] <- NA

analysis1$KIMsystolic[analysis1$KIMsystolic==996] <- NA

analysis1$KIMsystolic[analysis1$KIMsystolic==997] <- NA

analysis1$KIMsystolic[analysis1$KIMsystolic==998] <- NA

analysis1$KIMsystolic[analysis1$KIMsystolic==999] <- NA

## Frequency table
w = table(analysis1$KIMsystolic)

as.data.frame(w)


#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2$KIMsystolic <- analysis2$LI859

analysis2$KIMsystolic[analysis2$KIMsystolic==993] <- NA

analysis2$KIMsystolic[analysis2$KIMsystolic==996] <- NA

analysis2$KIMsystolic[analysis2$KIMsystolic==997] <- NA

analysis2$KIMsystolic[analysis2$KIMsystolic==998] <- NA

analysis2$KIMsystolic[analysis2$KIMsystolic==999] <- NA

## Frequency table
w = table(analysis2$KIMsystolic)

as.data.frame(w)


#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3$KIMsystolic <- analysis3$MI859

analysis3$KIMsystolic[analysis3$KIMsystolic==993] <- NA

analysis3$KIMsystolic[analysis3$KIMsystolic==996] <- NA

analysis3$KIMsystolic[analysis3$KIMsystolic==997] <- NA

analysis3$KIMsystolic[analysis3$KIMsystolic==998] <- NA

analysis3$KIMsystolic[analysis3$KIMsystolic==999] <- NA

# Frequency table
w = table(analysis3$KIMsystolic)

as.data.frame(w)


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4$KIMsystolic <- analysis4$NI859

analysis4$KIMsystolic[analysis4$KIMsystolic==993] <- NA

analysis4$KIMsystolic[analysis4$KIMsystolic==996] <- NA

analysis4$KIMsystolic[analysis4$KIMsystolic==997] <- NA

analysis4$KIMsystolic[analysis4$KIMsystolic==998] <- NA

analysis4$KIMsystolic[analysis4$KIMsystolic==999] <- NA

# Frequency table
w = table(analysis4$KIMsystolic)

as.data.frame(w)


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5$KIMsystolic <- analysis5$OI859

analysis5$KIMsystolic[analysis5$KIMsystolic==993] <- NA

analysis5$KIMsystolic[analysis5$KIMsystolic==996] <- NA

analysis5$KIMsystolic[analysis5$KIMsystolic==997] <- NA

analysis5$KIMsystolic[analysis5$KIMsystolic==998] <- NA

analysis5$KIMsystolic[analysis5$KIMsystolic==999] <- NA

# Frequency table
w = table(analysis5$KIMsystolic)

as.data.frame(w)



#### 2016 baseline covariates and 2018 stroke incidence
analysis6$KIMsystolic <- analysis6$PI859

analysis6$KIMsystolic[analysis6$KIMsystolic==993] <- NA

analysis6$KIMsystolic[analysis6$KIMsystolic==996] <- NA

analysis6$KIMsystolic[analysis6$KIMsystolic==997] <- NA

analysis6$KIMsystolic[analysis6$KIMsystolic==998] <- NA

analysis6$KIMsystolic[analysis6$KIMsystolic==999] <- NA


# Frequency table
w = table(analysis6$KIMsystolic)

as.data.frame(w)



############# SYSTOLIC BLOOD PRESSURE- ALL THREE MEASUREMENTS
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
# Systolic 1
analysis1$systolic1 <- analysis1$KI859

analysis1$systolic1[analysis1$systolic1==993] <- NA

analysis1$systolic1[analysis1$systolic1==996] <- NA

analysis1$systolic1[analysis1$systolic1==997] <- NA

analysis1$systolic1[analysis1$systolic1==998] <- NA

analysis1$systolic1[analysis1$systolic1==999] <- NA

# Systolic 2
analysis1$systolic2 <- analysis1$KI864

analysis1$systolic2[analysis1$systolic2==993] <- NA

analysis1$systolic2[analysis1$systolic2==996] <- NA

analysis1$systolic2[analysis1$systolic2==997] <- NA

analysis1$systolic2[analysis1$systolic2==998] <- NA

analysis1$systolic2[analysis1$systolic2==999] <- NA


# Systolic 3
analysis1$systolic3 <- analysis1$KI869

analysis1$systolic3[analysis1$systolic3==993] <- NA

analysis1$systolic3[analysis1$systolic3==996] <- NA

analysis1$systolic3[analysis1$systolic3==997] <- NA

analysis1$systolic3[analysis1$systolic3==998] <- NA

analysis1$systolic3[analysis1$systolic3==999] <- NA

# MEAN SYSTOLIC BLOOD PRESSURE
analysis1$systolic <- (analysis1$systolic1 + analysis1$systolic2 + analysis1$systolic3)/3



#### 2008 baseline covariates, 2010 and 2012 stroke incidence
# Systolic 1
analysis2$systolic1 <- analysis2$LI859

analysis2$systolic1[analysis2$systolic1==993] <- NA

analysis2$systolic1[analysis2$systolic1==996] <- NA

analysis2$systolic1[analysis2$systolic1==997] <- NA

analysis2$systolic1[analysis2$systolic1==998] <- NA

analysis2$systolic1[analysis2$systolic1==999] <- NA

# Systolic 2
analysis2$systolic2 <- analysis2$LI864

analysis2$systolic2[analysis2$systolic2==993] <- NA

analysis2$systolic2[analysis2$systolic2==996] <- NA

analysis2$systolic2[analysis2$systolic2==997] <- NA

analysis2$systolic2[analysis2$systolic2==998] <- NA

analysis2$systolic2[analysis2$systolic2==999] <- NA


# Systolic 3
analysis2$systolic3 <- analysis2$LI869

analysis2$systolic3[analysis2$systolic3==993] <- NA

analysis2$systolic3[analysis2$systolic3==996] <- NA

analysis2$systolic3[analysis2$systolic3==997] <- NA

analysis2$systolic3[analysis2$systolic3==998] <- NA

analysis2$systolic3[analysis2$systolic3==999] <- NA

# MEAN SYSTOLIC BLOOD PRESSURE
analysis2$systolic <- (analysis2$systolic1 + analysis2$systolic2 + analysis2$systolic3)/3



#### 2010 baseline covariates, 2012 and 2014 stroke incidence
# Systolic 1
analysis3$systolic1 <- analysis3$MI859

analysis3$systolic1[analysis3$systolic1==993] <- NA

analysis3$systolic1[analysis3$systolic1==996] <- NA

analysis3$systolic1[analysis3$systolic1==997] <- NA

analysis3$systolic1[analysis3$systolic1==998] <- NA

analysis3$systolic1[analysis3$systolic1==999] <- NA

# Systolic 2
analysis3$systolic2 <- analysis3$MI864

analysis3$systolic2[analysis3$systolic2==993] <- NA

analysis3$systolic2[analysis3$systolic2==996] <- NA

analysis3$systolic2[analysis3$systolic2==997] <- NA

analysis3$systolic2[analysis3$systolic2==998] <- NA

analysis3$systolic2[analysis3$systolic2==999] <- NA


# Systolic 3
analysis3$systolic3 <- analysis3$MI869

analysis3$systolic3[analysis3$systolic3==993] <- NA

analysis3$systolic3[analysis3$systolic3==996] <- NA

analysis3$systolic3[analysis3$systolic3==997] <- NA

analysis3$systolic3[analysis3$systolic3==998] <- NA

analysis3$systolic3[analysis3$systolic3==999] <- NA

# MEAN SYSTOLIC BLOOD PRESSURE
analysis3$systolic <- (analysis3$systolic1 + analysis3$systolic2 + analysis3$systolic3)/3


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
# Systolic 1
analysis4$systolic1 <- analysis4$NI859

analysis4$systolic1[analysis4$systolic1==993] <- NA

analysis4$systolic1[analysis4$systolic1==996] <- NA

analysis4$systolic1[analysis4$systolic1==997] <- NA

analysis4$systolic1[analysis4$systolic1==998] <- NA

analysis4$systolic1[analysis4$systolic1==999] <- NA

# Systolic 2
analysis4$systolic2 <- analysis4$NI864

analysis4$systolic2[analysis4$systolic2==993] <- NA

analysis4$systolic2[analysis4$systolic2==996] <- NA

analysis4$systolic2[analysis4$systolic2==997] <- NA

analysis4$systolic2[analysis4$systolic2==998] <- NA

analysis4$systolic2[analysis4$systolic2==999] <- NA


# Systolic 3
analysis4$systolic3 <- analysis4$NI869

analysis4$systolic3[analysis4$systolic3==993] <- NA

analysis4$systolic3[analysis4$systolic3==996] <- NA

analysis4$systolic3[analysis4$systolic3==997] <- NA

analysis4$systolic3[analysis4$systolic3==998] <- NA

analysis4$systolic3[analysis4$systolic3==999] <- NA

# MEAN SYSTOLIC BLOOD PRESSURE
analysis4$systolic <- (analysis4$systolic1 + analysis4$systolic2 + analysis4$systolic3)/3


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
# Systolic 1
analysis5$systolic1 <- analysis5$OI859

analysis5$systolic1[analysis5$systolic1==993] <- NA

analysis5$systolic1[analysis5$systolic1==996] <- NA

analysis5$systolic1[analysis5$systolic1==997] <- NA

analysis5$systolic1[analysis5$systolic1==998] <- NA

analysis5$systolic1[analysis5$systolic1==999] <- NA

# Systolic 2
analysis5$systolic2 <- analysis5$OI864

analysis5$systolic2[analysis5$systolic2==993] <- NA

analysis5$systolic2[analysis5$systolic2==996] <- NA

analysis5$systolic2[analysis5$systolic2==997] <- NA

analysis5$systolic2[analysis5$systolic2==998] <- NA

analysis5$systolic2[analysis5$systolic2==999] <- NA


# Systolic 3
analysis5$systolic3 <- analysis5$OI869

analysis5$systolic3[analysis5$systolic3==993] <- NA

analysis5$systolic3[analysis5$systolic3==996] <- NA

analysis5$systolic3[analysis5$systolic3==997] <- NA

analysis5$systolic3[analysis5$systolic3==998] <- NA

analysis5$systolic3[analysis5$systolic3==999] <- NA

# MEAN SYSTOLIC BLOOD PRESSURE
analysis5$systolic <- (analysis5$systolic1 + analysis5$systolic2 + analysis5$systolic3)/3


#### 2016 baseline covariates and 2018 stroke incidence
# Systolic 1
analysis6$systolic1 <- analysis6$PI859

analysis6$systolic1[analysis6$systolic1==993] <- NA

analysis6$systolic1[analysis6$systolic1==996] <- NA

analysis6$systolic1[analysis6$systolic1==997] <- NA

analysis6$systolic1[analysis6$systolic1==998] <- NA

analysis6$systolic1[analysis6$systolic1==999] <- NA

# Systolic 2
analysis6$systolic2 <- analysis6$PI864

analysis6$systolic2[analysis6$systolic2==993] <- NA

analysis6$systolic2[analysis6$systolic2==996] <- NA

analysis6$systolic2[analysis6$systolic2==997] <- NA

analysis6$systolic2[analysis6$systolic2==998] <- NA

analysis6$systolic2[analysis6$systolic2==999] <- NA


# Systolic 3
analysis6$systolic3 <- analysis6$PI869

analysis6$systolic3[analysis6$systolic3==993] <- NA

analysis6$systolic3[analysis6$systolic3==996] <- NA

analysis6$systolic3[analysis6$systolic3==997] <- NA

analysis6$systolic3[analysis6$systolic3==998] <- NA

analysis6$systolic3[analysis6$systolic3==999] <- NA


# MEAN SYSTOLIC BLOOD PRESSURE
analysis6$systolic <- (analysis6$systolic1 + analysis6$systolic2 + analysis6$systolic3)/3




############# BODY MASS INDEX
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
# Height
analysis1$heightin <- analysis1$KI834

analysis1$heightin[analysis1$heightin==996] <- NA

analysis1$heightm <- analysis1$heightin*0.0254

# Frequency table
w <- table(analysis1$heightin)

as.data.frame(w)

# Weight
analysis1$weightlbs <- analysis1$KC139

analysis1$weightlbs[analysis1$weightlbs==999] <- NA

analysis1$weightlbs[analysis1$weightlbs==998] <- NA

analysis1$weightkg <- analysis1$weightlbs*0.453592

# Frequency table
w <- table(analysis1$weightlbs)

as.data.frame(w)


# BMI
analysis1$bmi <- analysis1$weightkg/(analysis1$heightm^2)



#### 2008 baseline covariates, 2010 and 2012 stroke incidence
# Height
analysis2$heightin <- analysis2$LI834

analysis2$heightin[analysis2$heightin==996] <- NA

analysis2$heightm <- analysis2$heightin*0.0254

# Frequency table
w <- table(analysis2$heightin)

as.data.frame(w)

# Weight
analysis2$weightlbs <- analysis2$LC139

analysis2$weightlbs[analysis2$weightlbs==999] <- NA

analysis2$weightlbs[analysis2$weightlbs==998] <- NA

analysis2$weightkg <- analysis2$weightlbs*0.453592

# Frequency table
w <- table(analysis2$weightlbs)

as.data.frame(w)


# BMI
analysis2$bmi <- analysis2$weightkg/(analysis2$heightm^2)


#### 2010 baseline covariates, 2012 and 2014 stroke incidence
# Height
analysis3$heightin <- analysis3$MI834

analysis3$heightin[analysis3$heightin==996] <- NA

analysis3$heightm <- analysis3$heightin*0.0254

# Weight
analysis3$weightlbs <- analysis3$MC139

analysis3$weightlbs[analysis3$weightlbs==999] <- NA
  
analysis3$weightlbs[analysis3$weightlbs==998] <- NA

analysis3$weightkg <- analysis3$weightlbs*0.453592

# BMI
analysis3$bmi <- analysis3$weightkg/(analysis3$heightm^2)


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
# Height (no data)
analysis4$heightin <- analysis4$NI834

analysis4$heightin[analysis4$heightin==996] <- NA

analysis4$heightm <- analysis4$heightin*0.0254

# Weight
analysis4$weightlbs <- analysis4$NC139

analysis4$weightlbs[analysis4$weightlbs==999] <- NA

analysis4$weightlbs[analysis4$weightlbs==998] <- NA

analysis4$weightkg <- analysis4$weightlbs*0.453592

# BMI
analysis4$bmi <- analysis4$weightkg/(analysis4$heightm^2)


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
# Height
analysis5$heightin <- analysis5$OI834

analysis5$heightin[analysis5$heightin==996] <- NA

analysis5$heightm <- analysis5$heightin*0.0254

# Weight
analysis5$weightlbs <- analysis5$OC139

analysis5$weightlbs[analysis5$weightlbs==999] <- NA

analysis5$weightlbs[analysis5$weightlbs==998] <- NA

analysis5$weightkg <- analysis5$weightlbs*0.453592

# BMI 
analysis5$bmi <- analysis5$weightkg/(analysis5$heightm^2)


#### 2016 baseline covariates and 2018 stroke incidence
# Height 
analysis6$heightin <- analysis6$PI834

analysis6$heightin[analysis6$heightin==996] <- NA

analysis6$heightm <- analysis6$heightin*0.0254

# Weight
analysis6$weightlbs <- analysis6$PC139

analysis6$weightlbs[analysis6$weightlbs==999] <- NA

analysis6$weightlbs[analysis6$weightlbs==998] <- NA

analysis6$weightkg <- analysis6$weightlbs*0.453592

# BMI
analysis6$bmi <- analysis6$weightkg/(analysis6$heightm^2)


############# BMI CATEGORICAL
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1$KIMbmi[analysis1$bmi < 25] <- 1

analysis1$KIMbmi[analysis1$bmi > 24.9999 & analysis1$bmi < 30] <- 2

analysis1$KIMbmi[analysis1$bmi > 29.9999 ] <- 3


#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2$KIMbmi[analysis2$bmi < 25] <- 1

analysis2$KIMbmi[analysis2$bmi > 24.9999 & analysis2$bmi < 30] <- 2

analysis2$KIMbmi[analysis2$bmi > 29.9999 ] <- 3



#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3$KIMbmi[analysis3$bmi < 25] <- 1

analysis3$KIMbmi[analysis3$bmi > 24.9999 & analysis3$bmi < 30] <- 2

analysis3$KIMbmi[analysis3$bmi > 29.9999 ] <- 3


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4$KIMbmi[analysis4$bmi < 25] <- 1

analysis4$KIMbmi[analysis4$bmi > 24.9999 & analysis4$bmi < 30] <- 2

analysis4$KIMbmi[analysis4$bmi > 29.9999 ] <- 3


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5$KIMbmi[analysis5$bmi < 25] <- 1

analysis5$KIMbmi[analysis5$bmi > 24.9999 & analysis5$bmi < 30] <- 2

analysis5$KIMbmi[analysis5$bmi > 29.9999 ] <- 3

#### 2016 baseline covariates and 2018 stroke incidence
analysis6$KIMbmi[analysis6$bmi < 25] <- 1

analysis6$KIMbmi[analysis6$bmi > 24.9999 & analysis6$bmi < 30] <- 2

analysis6$KIMbmi[analysis6$bmi > 29.9999 ] <- 3




############# HEART DISEASE
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1$heartdis[analysis1$KC036== 1] <- 1

analysis1$heartdis[analysis1$KC036== 3] <- 1

analysis1$heartdis[analysis1$KC036== 5] <- 0

analysis1$heartdis[analysis1$KC036== 4] <- 0

# Frequency table
w <- table(analysis1$heartdis)
as.data.frame(w)

#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2$heartdis[analysis2$LC036== 1] <- 1

analysis2$heartdis[analysis2$LC036== 3] <- 1

analysis2$heartdis[analysis2$LC036== 5] <- 0

analysis2$heartdis[analysis2$LC036== 4] <- 0

# Frequency table
w <- table(analysis2$heartdis)
as.data.frame(w)


#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3$heartdis[analysis3$MC036== 1] <- 1

analysis3$heartdis[analysis3$MC036== 3] <- 1

analysis3$heartdis[analysis3$MC036== 5] <- 0

analysis3$heartdis[analysis3$MC036== 4] <- 0


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4$heartdis[analysis4$NC036== 1] <- 1

analysis4$heartdis[analysis4$NC036== 3] <- 1

analysis4$heartdis[analysis4$NC036== 5] <- 0

analysis4$heartdis[analysis4$NC036== 4] <- 0


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5$heartdis[analysis5$OC036== 1] <- 1

analysis5$heartdis[analysis5$OC036== 3] <- 1

analysis5$heartdis[analysis5$OC036== 5] <- 0

analysis5$heartdis[analysis5$OC036== 4] <- 0

#### 2016 baseline covariates and 2018 stroke incidence
analysis6$heartdis[analysis6$PC036== 1] <- 1

analysis6$heartdis[analysis6$PC036== 3] <- 1

analysis6$heartdis[analysis6$PC036== 5] <- 0

analysis6$heartdis[analysis6$PC036== 4] <- 0

#### Frequency table
w = table(analysis3$heartdis)

as.data.frame(w)



############# PSYCHOLOGICAL VARIABLES REVERSE CODED
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
### Depression
analysis1$KLB027I <- 6 - analysis1$KLB027I

analysis1$KLB027J <- 6 - analysis1$KLB027J

analysis1$KLB027K <- 6 - analysis1$KLB027K

analysis1$KLB027L <- 6 - analysis1$KLB027L

analysis1$KLB027M <- 6 - analysis1$KLB027M

analysis1$KLB027N <- 6 - analysis1$KLB027N


### Optimism
analysis1$KLB019F <- 6 - analysis1$KLB019F

analysis1$KLB019J <- 6 - analysis1$KLB019J

analysis1$KLB019K <- 6 - analysis1$KLB019K


### Positive Affect
analysis1$KLB027A <- 6 - analysis1$KLB027A

analysis1$KLB027B <- 6 - analysis1$KLB027B

analysis1$KLB027C <- 6 - analysis1$KLB027C

analysis1$KLB027D <- 6 - analysis1$KLB027D

analysis1$KLB027E <- 6 - analysis1$KLB027E

analysis1$KLB027F <- 6 - analysis1$KLB027F

### Negative Affect
analysis1$KLB027I <- 6 - analysis1$KLB027I

analysis1$KLB027J <- 6 - analysis1$KLB027J

analysis1$KLB027K <- 6 - analysis1$KLB027K

analysis1$KLB027M <- 6 - analysis1$KLB027M

analysis1$KLB027N <- 6 - analysis1$KLB027N



### Social participation
analysis1$KLB001A <- 6 - analysis1$KLB001a

analysis1$KLB001B <- 6 - analysis1$KLB001b

analysis1$KLB001C <- 6 - analysis1$KLB001c

analysis1$KLB001F <- 6 - analysis1$KLB001f

analysis1$KLB001G <- 6 - analysis1$KLB001g


### Purpose in life
analysis1$KLB035B <- 6 - analysis1$KLB035B

analysis1$KLB035D <- 6 - analysis1$KLB035D

analysis1$KLB035E <- 6 - analysis1$KLB035E




#### 2008 baseline covariates, 2010 and 2012 stroke incidence
### Depression
analysis2$LLB027I <- 6 - analysis2$LLB027I

analysis2$LLB027J <- 6 - analysis2$LLB027J

analysis2$LLB027K <- 6 - analysis2$LLB027K

analysis2$LLB027L <- 6 - analysis2$LLB027L

analysis2$LLB027M <- 6 - analysis2$LLB027M

analysis2$LLB027N <- 6 - analysis2$LLB027N


### Optimism
analysis2$LLB019F <- 6 - analysis2$LLB019F

analysis2$LLB019J <- 6 - analysis2$LLB019J

analysis2$LLB019K <- 6 - analysis2$LLB019K


### Positive Affect
analysis2$LLB027A <- 6 - analysis2$LLB027A

analysis2$LLB027B <- 6 - analysis2$LLB027B

analysis2$LLB027C <- 6 - analysis2$LLB027C

analysis2$LLB027D <- 6 - analysis2$LLB027D

analysis2$LLB027E <- 6 - analysis2$LLB027E

analysis2$LLB027F <- 6 - analysis2$LLB027F

### Negative Affect
analysis2$LLB027I <- 6 - analysis2$LLB027I

analysis2$LLB027J <- 6 - analysis2$LLB027J

analysis2$LLB027K <- 6 - analysis2$LLB027K

analysis2$LLB027M <- 6 - analysis2$LLB027M

analysis2$LLB027N <- 6 - analysis2$LLB027N



### Social participation
analysis2$LLB001A <- 6 - analysis2$LLB001a

analysis2$LLB001B <- 6 - analysis2$LLB001b

analysis2$LLB001C <- 6 - analysis2$LLB001c

analysis2$LLB001F <- 6 - analysis2$LLB001f

analysis2$LLB001G <- 6 - analysis2$LLB001g


### Purpose in life
analysis2$LLB035B <- 6 - analysis2$LLB035B

analysis2$LLB035D <- 6 - analysis2$LLB035D

analysis2$LLB035E <- 6 - analysis2$LLB035E




#### 2010 baseline covariates, 2012 and 2014 stroke incidence
### Depression
analysis3$MLB027I <- 6 - analysis3$MLB027I

analysis3$MLB027J <- 6 - analysis3$MLB027J

analysis3$MLB027K <- 6 - analysis3$MLB027K

analysis3$MLB027L <- 6 - analysis3$MLB027L

analysis3$MLB027M <- 6 - analysis3$MLB027M

analysis3$MLB027N <- 6 - analysis3$MLB027N


### Optimism
analysis3$MLB019F <- 6 - analysis3$MLB019F

analysis3$MLB019J <- 6 - analysis3$MLB019J

analysis3$MLB019K <- 6 - analysis3$MLB019K


### Positive Affect
analysis3$MLB027A <- 6 - analysis3$MLB027A

analysis3$MLB027B <- 6 - analysis3$MLB027B

analysis3$MLB027C <- 6 - analysis3$MLB027C

analysis3$MLB027D <- 6 - analysis3$MLB027D

analysis3$MLB027E <- 6 - analysis3$MLB027E

analysis3$MLB027F <- 6 - analysis3$MLB027F

### Negative Affect
analysis3$MLB027I <- 6 - analysis3$MLB027I

analysis3$MLB027J <- 6 - analysis3$MLB027J

analysis3$MLB027K <- 6 - analysis3$MLB027K

analysis3$MLB027M <- 6 - analysis3$MLB027M

analysis3$MLB027N <- 6 - analysis3$MLB027N



### Social participation
analysis3$MLB001A <- 6 - analysis3$MLB001A

analysis3$MLB001B <- 6 - analysis3$MLB001B

analysis3$MLB001C <- 6 - analysis3$MLB001C

analysis3$MLB001F <- 6 - analysis3$MLB001F

analysis3$MLB001G <- 6 - analysis3$MLB001G


### Purpose in life
analysis3$MLB035B <- 6 - analysis3$MLB035B

analysis3$MLB035D <- 6 - analysis3$MLB035D

analysis3$MLB035E <- 6 - analysis3$MLB035E



#### 2012 baseline covariates, 2014 and 2016 stroke incidence
### Depression
analysis4$NLB027I <- 6 - analysis4$NLB027I

analysis4$NLB027J <- 6 - analysis4$NLB027J

analysis4$NLB027K <- 6 - analysis4$NLB027K

analysis4$NLB027L <- 6 - analysis4$NLB027L

analysis4$NLB027M <- 6 - analysis4$NLB027M

analysis4$NLB027N <- 6 - analysis4$NLB027N


### Optimism
analysis4$NLB019F <- 6 - analysis4$NLB019F

analysis4$NLB019J <- 6 - analysis4$NLB019J

analysis4$NLB019K <- 6 - analysis4$NLB019K


### Positive Affect
analysis4$NLB027A <- 6 - analysis4$NLB027A

analysis4$NLB027B <- 6 - analysis4$NLB027B

analysis4$NLB027C <- 6 - analysis4$NLB027C

analysis4$NLB027D <- 6 - analysis4$NLB027D

analysis4$NLB027E <- 6 - analysis4$NLB027E

analysis4$NLB027F <- 6 - analysis4$NLB027F

### Negative Affect
analysis4$NLB027I <- 6 - analysis4$NLB027I

analysis4$NLB027J <- 6 - analysis4$NLB027J

analysis4$NLB027K <- 6 - analysis4$NLB027K

analysis4$NLB027M <- 6 - analysis4$NLB027M

analysis4$NLB027N <- 6 - analysis4$NLB027N



### Social participation
analysis4$NLB001A <- 6 - analysis4$NLB001A

analysis4$NLB001B <- 6 - analysis4$NLB001B

analysis4$NLB001C <- 6 - analysis4$NLB001C

analysis4$NLB001F <- 6 - analysis4$NLB001F

analysis4$NLB001G <- 6 - analysis4$NLB001G


### Purpose in life
analysis4$NLB035B <- 6 - analysis4$NLB035B

analysis4$NLB035D <- 6 - analysis4$NLB035D

analysis4$NLB035E <- 6 - analysis4$NLB035E


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
### Anxiety
analysis5$OLB031L <- 5 - analysis5$OLB031L

### Social participation
analysis5$OLB001R <- 8 - analysis5$OLB001R


### Frustration
analysis5$OLB026J <- 6 - analysis5$OLB026J


## Positive affect
analysis5$OLB026K <- 6 - analysis5$OLB026K

analysis5$OLB026X <- 6 - analysis5$OLB026X 


## Optimism
analysis5$OLB018E <- 7 - analysis5$OLB018E

analysis5$OLB018F <- 7 - analysis5$OLB018F


## Purpose in life
analysis5$OLB033E <- 7 - analysis5$OLB033E 

## Social pariticipation
analysis5$OLB001N <- 8 - analysis5$OLB001N

analysis5$OLB012B <- 8 - analysis5$OLB012B



#### 2016 baseline covariates and 2018 stroke incidence
### Anxiety
analysis6$PLB031L <- 5 - analysis6$PLB031L

### Social participation
analysis6$PLB001R <- 8 - analysis6$PLB001R


### Frustration
analysis6$PLB026J <- 6 - analysis6$PLB026J


## Positive affect
analysis6$PLB026K <- 6 - analysis6$PLB026K

analysis6$PLB026X <- 6 - analysis6$PLB026X 


## Optimism
analysis6$PLB018E <- 7 - analysis6$PLB018E

analysis6$PLB018F <- 7 - analysis6$PLB018F


## Purpose in life
analysis6$PLB033E <- 7 - analysis6$PLB033E 

## Social pariticipation
analysis6$PLB001N <- 8 - analysis6$PLB001N

analysis6$PLB012B <- 8 - analysis6$PLB012B




############# PSYCHOLOGICAL VARIABLES, KIM ET AL
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
### Anxiety
analysis1$KIManxiety <- (analysis1$KLB041A + analysis1$KLB041B + analysis1$KLB041C + analysis1$KLB041D + analysis1$KLB041E)/5

### Cynical hostility
analysis1$KIMcynical <- (analysis1$KLB019A + analysis1$KLB019B + analysis1$KLB019C + analysis1$KLB019D + analysis1$KLB019E)/5

### Depression
analysis1$KIMdepression <- (analysis1$KLB027I + analysis1$KLB027J + analysis1$KLB027K + analysis1$KLB027L + analysis1$KLB027M + analysis1$KLB027N)/6

### Optimism
analysis1$KIMoptimism <- (analysis1$KLB019F + analysis1$KLB019G + analysis1$KLB019H + analysis1$KLB019I + analysis1$KLB019J + analysis1$KLB019K)/6

### Positive Affect
analysis1$KIMposaf <- (analysis1$KLB027A + analysis1$KLB027B + analysis1$KLB027C + analysis1$KLB027D + analysis1$KLB027E + analysis1$KLB027F)/6

### Negative Affect
analysis1$KIMnegaf <- (analysis1$KLB027I + analysis1$KLB027J + analysis1$KLB027K + analysis1$KLB027L + analysis1$KLB027M + analysis1$KLB027N)/6

### Social participation
analysis1$KIMparticipation <- (analysis1$KLB001a + analysis1$KLB001b + analysis1$KLB001c + analysis1$KLB001e + analysis1$KLB001f + analysis1$KLB001g)/8

### Purpose in life
analysis1$KIMpurposeinlife <- (analysis1$KLB035A + analysis1$KLB035B + analysis1$KLB035C + analysis1$KLB035D + analysis1$KLB035E + analysis1$KLB035F + analysis1$KLB035G)/7




#### 2008 baseline covariates, 2010 and 2012 stroke incidence
### Anxiety
analysis2$KIManxiety <- (analysis2$LLB041A + analysis2$LLB041B + analysis2$LLB041C + analysis2$LLB041D + analysis2$LLB041E)/5

### Cynical hostility
analysis2$KIMcynical <- (analysis2$LLB019A + analysis2$LLB019B + analysis2$LLB019C + analysis2$LLB019D + analysis2$LLB019E)/5

### Depression
analysis2$KIMdepression <- (analysis2$LLB027I + analysis2$LLB027J + analysis2$LLB027K + analysis2$LLB027L + analysis2$LLB027M + analysis2$LLB027N)/6

### Optimism
analysis2$KIMoptimism <- (analysis2$LLB019F + analysis2$LLB019G + analysis2$LLB019H + analysis2$LLB019I + analysis2$LLB019J + analysis2$LLB019K)/6

### Positive Affect
analysis2$KIMposaf <- (analysis2$LLB027A + analysis2$LLB027B + analysis2$LLB027C + analysis2$LLB027D + analysis2$LLB027E + analysis2$LLB027F)/6

### Negative Affect
analysis2$KIMnegaf <- (analysis2$LLB027I + analysis2$LLB027J + analysis2$LLB027K + analysis2$LLB027L + analysis2$LLB027M + analysis2$LLB027N)/6

### Social participation
analysis2$KIMparticipation <- (analysis2$LLB001a + analysis2$LLB001b + analysis2$LLB001c + analysis2$LLB001e + analysis2$LLB001f + analysis2$LLB001g)/8

### Purpose in life
analysis2$KIMpurposeinlife <- (analysis2$LLB035A + analysis2$LLB035B + analysis2$LLB035C + analysis2$LLB035D + analysis2$LLB035E + analysis2$LLB035F + analysis2$LLB035G)/7



#### 2010 baseline covariates, 2012 and 2014 stroke incidence
### Anxiety
analysis3$KIManxiety <- (analysis3$MLB041A + analysis3$MLB041B + analysis3$MLB041C + analysis3$MLB041D + analysis3$MLB041E)/5

### Cynical hostility
analysis3$KIMcynical <- (analysis3$MLB019A + analysis3$MLB019B + analysis3$MLB019C + analysis3$MLB019D + analysis3$MLB019E)/5

### Depression
analysis3$KIMdepression <- (analysis3$MLB027I + analysis3$MLB027J + analysis3$MLB027K + analysis3$MLB027L + analysis3$MLB027M + analysis3$MLB027N)/6

### Optimism
analysis3$KIMoptimism <- (analysis3$MLB019F + analysis3$MLB019G + analysis3$MLB019H + analysis3$MLB019I + analysis3$MLB019J + analysis3$MLB019K)/6

### Positive Affect
analysis3$KIMposaf <- (analysis3$MLB027A + analysis3$MLB027B + analysis3$MLB027C + analysis3$MLB027D + analysis3$MLB027E + analysis3$MLB027F)/6

### Negative Affect
analysis3$KIMnegaf <- (analysis3$MLB027I + analysis3$MLB027J + analysis3$MLB027K + analysis3$MLB027L + analysis3$MLB027M + analysis3$MLB027N)/6

### Social participation
analysis3$KIMparticipation <- (analysis3$MLB001A + analysis3$MLB001B + analysis3$MLB001C + analysis3$MLB001E + analysis3$MLB001F + analysis3$MLB001G)/8

### Purpose in life
analysis3$KIMpurposeinlife <- (analysis3$MLB035A + analysis3$MLB035B + analysis3$MLB035C + analysis3$MLB035D + analysis3$MLB035E + analysis3$MLB035F + analysis3$MLB035G)/7


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
### Anxiety
analysis4$KIManxiety <- (analysis4$NLB041A + analysis4$NLB041B + analysis4$NLB041C + analysis4$NLB041D + analysis4$NLB041E)/5

### Cynical hostility
analysis4$KIMcynical <- (analysis4$NLB019A + analysis4$NLB019B + analysis4$NLB019C + analysis4$NLB019D + analysis4$NLB019E)/5

### Depression
analysis4$KIMdepression <- (analysis4$NLB027I + analysis4$NLB027J + analysis4$NLB027K + analysis4$NLB027L + analysis4$NLB027M + analysis4$NLB027N)/6

### Optimism
analysis4$KIMoptimism <- (analysis4$NLB019F + analysis4$NLB019G + analysis4$NLB019H + analysis4$NLB019I + analysis4$NLB019J + analysis4$NLB019K)/6

### Positive Affect
analysis4$KIMposaf <- (analysis4$NLB027A + analysis4$NLB027B + analysis4$NLB027C + analysis4$NLB027D + analysis4$NLB027E + analysis4$NLB027F)/6

### Negative Affect
analysis4$KIMnegaf <- (analysis4$NLB027I + analysis4$NLB027J + analysis4$NLB027K + analysis4$NLB027L + analysis4$NLB027M + analysis4$NLB027N)/6

### Social participation
analysis4$KIMparticipation <- (analysis4$NLB001A + analysis4$NLB001B + analysis4$NLB001C + analysis4$NLB001E + analysis4$NLB001F + analysis4$NLB001G)/8

### Purpose in life
analysis4$KIMpurposeinlife <- (analysis4$NLB035A + analysis4$NLB035B + analysis4$NLB035C + analysis4$NLB035D + analysis4$NLB035E + analysis4$NLB035F + analysis4$NLB035G)/7



#### 2014 baseline covariates and 2016, 2018 stroke incidence
## Depression
analysis5$depression <- (analysis5$OLB026S	+ analysis5$OLB018H + analysis5$OLB026K)/3

## Optimism
analysis5$optimism <- (analysis5$OLB018A + analysis5$OLB018B + analysis5$OLB018C + analysis5$OLB018D + analysis5$OLB018E + analysis5$OLB018F)/6

## Social participation
analysis5$participation <- (analysis5$OLB001I + analysis5$OLB001R + analysis5$OLB001N + analysis5$OLB012B)/4

## Purpose in life
analysis5$purposeinlife <- (analysis5$OLB033A  + analysis5$OLB033B + analysis5$OLB033C + analysis5$OLB033D + analysis5$OLB033F + analysis5$OLB033G + analysis5$OLB033E)/7

## Anxiety
analysis5$anxiety <- (analysis5$OLB031L + analysis5$OLB026R)/2

## Frustration
analysis5$frustration <- (analysis5$OLB026J + analysis5$OLB026X)/2



#### 2016 baseline covariates and 2018 stroke incidence
## Depression
analysis6$depression <- (analysis6$PLB018H + analysis6$PLB026K)/2

## Optimism
analysis6$optimism <- (analysis6$PLB018A + analysis6$PLB018B + analysis6$PLB018C + analysis6$PLB018D + analysis6$PLB018E + analysis6$PLB018F)/6

## Social participation
analysis6$participation <- (analysis6$PLB001I + analysis6$PLB001R + analysis6$PLB001N + analysis6$PLB012B)/4

## Purpose in life
analysis6$purposeinlife <- (analysis6$PLB033A  + analysis6$PLB033C + analysis6$PLB033D + analysis6$PLB033F + analysis6$PLB033G + analysis6$PLB033E)/7

## Anxiety
analysis6$anxiety <- (analysis6$PLB031L + analysis6$PLB026R)/2

## Frustration
analysis6$frustration <- (analysis6$PLB026J + analysis6$PLB026X)/2




############### SAVE ANALYSIS VARIABLES
#### Varibles for analysis1 to analysis4
analysis1to4.vars <- c("HHID", "PN", "MSUBHH", "age", "sex", "stroke", "marstat", "KIMrace", "race", "KIMeducation", "education", "expenses", "income", "wealth", "KIMwealth", "func", "func014", "func016",  "func021",  "func023", "func025", "KIMsmoke", "numCigs", "KIMneverexercise", "KIMlowexercise", "KIMmodexercise", "KIMhighexercise", "mildexercise", "modexercise", "highexercise", "KIMalcohol", "alcohol", "hypertension", "diabetes", "KIMsystolic", "systolic", "bmi", "KIMbmi", "heartdis", "KIManxiety", "KIMcynical", "KIMdepression", "KIMoptimism", "KIMposaf", "KIMnegaf", "KIMparticipation", "KIMpurposeinlife" )

analysis5to6.vars <- c("HHID", "PN", "MSUBHH", "age", "sex", "stroke", "marstat", "KIMrace", "race", "KIMeducation", "education", "expenses", "income", "wealth", "KIMwealth", "func", "func014", "func016",  "func021",  "func023", "func025", "KIMsmoke", "numCigs", "KIMneverexercise", "KIMlowexercise", "KIMmodexercise", "KIMhighexercise", "mildexercise", "modexercise", "highexercise", "KIMalcohol", "alcohol", "hypertension", "diabetes", "KIMsystolic", "systolic", "bmi", "KIMbmi", "heartdis",  "depression", "optimism", "participation", "purposeinlife",  "anxiety", "frustration")

psychology1to4.vars <- c("NLB041A",	"NLB019A",	"NLB027I",	"NLB019F",	"NLB027A",	"NLB001A",	"NLB035A ",	"NLB035B",	"NLB035C",	"NLB035D",	"NLB035F",	"NLB035G",	"NLB041B",	"NLB027N",	"NLB019B",	"NLB027J",	"NLB019G",	"NLB027B",	"NLB001B",	"NLB027J",	"NLB041C",	"NLB019C",	"NLB027K",	"NLB019H",	"NLB027C",	"NLB001C",	"NLB041D",	"NLB019D",	"NLB027L",	"NLB019I",	"NLB027D",	"NLB041E",	"NLB019E",	"NLB027M",	"NLB019J",	"NLB027E",	"NLB001E",	"NLB027M",	"NLB035E ",	"NLB019K",	"NLB027F",	"NLB001F",	"NLB001G")

psychology5.vars <- c("OLB001I",	"OLB026S",	"OLB018H",	"OLB018B",	"OLB018E",	"OLB018F",	"OLB033C",	"OLB033D",	"OLB033F",	"OLB033G",	"OLB026R",	"OLB026K",	"OLB018C",	"OLB033E",	"OLB026J",	"OLB018A",	"OLB018D",	"OLB031L",	"OLB026X",	"OLB001R",	"OLB001N",	"OLB012B",	"OLB033A",	"OLB033B")

psychology6.vars <- c("OLB018H",	"OLB018B",	"OLB018A",	"OLB026K",	"OLB033F",	"OLB033G", "OLB033D",	"OLB018C",	"OLB018D",	"OLB018E",	"OLB018F",	"OLB001I",	"OLB001R",	"OLB001N",	"OLB012B",	"OLB033A",	"OLB033C",	"OLB033E",	"OLB031L",	"OLB026R",	"OLB026J",	"OLB026X")


######## CLEAR WORKSPACE
rm(a, dam, HRS2006, HRS2006to2018, HRS2008, HRS2009, HRS2010, HRS2012, HRS2014, HRS2016, HRS2018, mydata, stroke10, stroke12, stroke14, stroke16, stroke18, test, test1, analysis.vars, i, rec, tar, val, w, wealth.na)


######## Remove variables not needed
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
analysis1 <- analysis1[-c(4:38, 46:144, 189:191, 193:197, 206:207)]

#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2 <- analysis2[-c(4:38, 47:144, 189:191, 193:197, 206:207)]


#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3 <- analysis3[-c(4:42, 50:146, 191:193, 195:199, 208, 209)]



#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4 <- analysis4[-c(4:38, 45:141, 186:194)]


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5 <- analysis5[-c(4:137, 164:171)]


#### 2016 baseline covariates and 2018 stroke incidence
analysis6 <- analysis6[-c(4:137, 164:169)]



