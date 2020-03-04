########### CORRELATION MATRICES FOR PSYCHO-SOCIAL VARIABLES
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


########### FIG 1- HEATMAP OF COLLINEARITY IN PSYCHO-SOCIAL VARIABLES
##### Load packages
library(corrplot)

##### Markov Chain Imputed Datasets
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
corrplot(corMarkov.analysis1, method = "circle", tl.cex=0.5, tl.col="black")


######### FIG 2- ILLUSTRATION OF RESIDUALS FROM MODELS WITH AND WITHOUT APPROPRIATE CONTROL VARIABLES

#### 2006 baseline covariates, 2008 and 2010 stroke incidence
## Transform sex, marital status and race into nominal variables
# Non imputed dataset
analysis1$sex <- factor(analysis1$sex)

analysis1$marstat <-factor(analysis1$marstat)

analysis1$race <- factor(analysis1$race)


# Markov chain imputed dataset
markov.analysis1$sex <- factor(markov.analysis1$sex)

markov.analysis1$marstat <-factor(markov.analysis1$marstat)

markov.analysis1$race <- factor(markov.analysis1$race)


# Bootstrap imputed dataset
bootstrap.analysis1$sex <- factor(markov.analysis1$sex)

bootstrap.analysis1$marstat <-factor(markov.analysis1$marstat)

bootstrap.analysis1$race <- factor(markov.analysis1$race)



## Logistic regression
# Kim et al original model, no imputation
KIM.model <- glm(stroke ~ KIMpurposeinlife + age  + sex + race + marstat + KIMrace + KIMeducation + KIMwealth + func + KIMsmoke + KIMneverexercise + KIMlowexercise + KIMmodexercise + KIMhighexercise + KIMalcohol + hypertension + KIMsystolic+ diabetes + KIMbmi + heartdis + KIManxiety + KIMcynical + KIMdepression + KIMoptimism + KIMposaf + KIMnegaf + KIMparticipation, data = markov.analysis1,  family = "binomial")

# Non-imputed data
non.fig2 <- glm(stroke ~ KIMpurposeinlife + age  + sex + marstat + race + KIMeducation + KIMwealth + func + numCigs + KIMneverexercise + KIMlowexercise + KIMmodexercise + KIMhighexercise + alcohol  + diabetes + KIMbmi + heartdis + KIManxiety + KIMcynical + KIMdepression + KIMoptimism + KIMposaf + KIMnegaf + KIMparticipation, data = analysis1,  family = "binomial")

# Markov imputation
markov.fig2 <- glm(stroke ~ KIMpurposeinlife + age  + sex + marstat + race + KIMeducation + KIMwealth + func + numCigs + KIMneverexercise + KIMlowexercise + KIMmodexercise + KIMhighexercise + alcohol  + diabetes + KIMbmi + heartdis + KIManxiety + KIMcynical + KIMdepression + KIMoptimism + KIMposaf + KIMnegaf + KIMparticipation, data = markov.analysis1,  family = "binomial")



## Visualise residuals
library(ggplot2)
library(ggeffects)
library(devtools)

# Kim et al original model, no imputation
resid.KIM.fig2 <- ggplot(KIM.model, aes(x= predict(KIM.model), y=residuals(KIM.model))) + geom_point(aes(col=residuals(KIM.model), alpha=0.0)) + geom_smooth(method='lm',formula=y~x)

resid.KIM.fig2 + xlab("Predicted Values") + ylab("Residuals")


# Non-imputed data
resid.non.fig2 <- ggplot(non.fig2, aes(x= predict(non.fig2), y=residuals(non.fig2))) + geom_point(aes(col=residuals(non.fig2), alpha=0.0)) + geom_smooth(method='lm',formula=y~x)

resid.non.fig2 + xlab("Predicted Values") + ylab("Residuals")


# Markov imputation
resid.markov.fig2 <- ggplot(markov.fig2, aes(x= predict(markov.fig2), y=residuals(markov.fig2))) + geom_point(aes(col=residuals(markov.fig2), alpha=0.0)) + geom_smooth(method='lm',formula=y~x)

resid.markov.fig2 + xlab("Predicted Values") + ylab("Residuals")



########## TABLE 2- VARIANCE INFLATION FACTOR- PSYCHO-SOCIAL VARIABLES
## Load packages
library(regclass)

# Kim et al original model, Markov imputation
KIM.model.vif <- glm(stroke ~ KIMpurposeinlife + age  + sex + race + marstat + KIMrace + KIMeducation + KIMwealth + func + KIMsmoke + KIMneverexercise + KIMlowexercise + KIMmodexercise + KIMhighexercise + KIMalcohol + hypertension + KIMsystolic+ diabetes + KIMbmi + heartdis +  KIManxiety + KIMcynical + KIMdepression + KIMoptimism + KIMposaf + KIMparticipation, data = markov.analysis1,  family = "binomial")


VIF(KIM.model.vif)


########## FIG 3-  ILLUSTRATION OF NON-RESPONSE RATES
## Load packages
library(visdat)
library(naniar)
library(ggpubr)
library(ggplot2)

## View all psych-social missing data
all.miss <- vis_miss(analysis1[11:54],large_data_size =400000)


## g_miss_fct- missingness in purpose in life by sex, race and marital status
# Subset psychosocial variables
analysis1.psycho.mis <- data.frame(analysis1[17:22], analysis1[49], analysis1[57:64])

# Plot missing data- sex
gg_miss_fct(x = analysis1.psycho.mis, fct = sex)

# Plot missing data- race
gg_miss_fct(x = analysis1.psycho.mis, fct = KIMrace)

# Plot missing data- marital status
gg_miss_fct(x = analysis1.psycho.mis, fct = marstat)

# Plot missing data- education
gg_miss_fct(x = analysis1.psycho.mis, fct = education)


## Scatter plot- missingness in purpose in life by age and income
# Age
AgemisPlot035A <- ggplot(analysis1.psycho.mis, aes(x = Making.plans, y = age)) + geom_miss_point() 

AgemisPlot035B <- ggplot(analysis1.psycho.mis, aes(x = Activities.trivial, y = age)) + geom_miss_point() 

AgemisPlot035C <- ggplot(analysis1.psycho.mis, aes(x = Active.planning, y = age)) + geom_miss_point() 

AgemisPlot035D <- ggplot(analysis1.psycho.mis, aes(x = Accomplish, y = age)) + geom_miss_point() 

AgemisPlot035E <- ggplot(analysis1.psycho.mis, aes(x = Done.all.to.do.in.life, y = age)) + geom_miss_point() 

AgemisPlot035G <- ggplot(analysis1.psycho.mis, aes(x = Sense.of.direction, y = age)) + geom_miss_point() 

ggarrange(AgemisPlot035A, AgemisPlot035B, AgemisPlot035C, AgemisPlot035D, AgemisPlot035E, AgemisPlot035G)


# Income
EducmisPlot035A <- ggplot(analysis1.psycho.mis, aes(x = Making.plans, y = wealth)) + geom_miss_point() 

EducmisPlot035B <- ggplot(analysis1.psycho.mis, aes(x = Activities.trivial, y = wealth)) + geom_miss_point() 

EducmisPlot035C <- ggplot(analysis1.psycho.mis, aes(x = Active.planning, y = wealth)) + geom_miss_point() 

EducmisPlot035D <- ggplot(analysis1.psycho.mis, aes(x = Accomplish, y = wealth)) + geom_miss_point() 

EducmisPlot035E <- ggplot(analysis1.psycho.mis, aes(x = Done.all.to.do.in.life, y = wealth)) + geom_miss_point() 

EducmisPlot035F <- ggplot(analysis1.psycho.mis, aes(x = One.day.a.time, y = wealth)) + geom_miss_point() 

EducmisPlot035G <- ggplot(analysis1.psycho.mis, aes(x = Sense.of.direction, y = wealth)) + geom_miss_point() 

ggarrange(EducmisPlot035A, EducmisPlot035B, EducmisPlot035C, EducmisPlot035D, EducmisPlot035E, EducmisPlot035F, EducmisPlot035G)




######### FIG 4- RESIDUALS FROM MODELS WITH DICHOTOMISED AND CONTINUOUS VARIABLES
## Logistic regression
# Non-imputed data
non.fig4 <- glm(stroke ~ KIMpurposeinlife + age  + sex + marstat + KIMrace + education + wealth + func014 + func016 + func021 + func023 + func025 + KIMsmoke + mildexercise + modexercise + highexercise + KIMalcohol  + diabetes + bmi + heartdis + KIMsystolic + hypertension + KIManxiety + KIMcynical + KIMdepression + KIMoptimism + KIMposaf + KIMnegaf + KIMparticipation, data = analysis1,  family = "binomial")

# Markov imputation
markov.fig4 <- glm(stroke ~ KIMpurposeinlife + age  + sex + marstat + KIMrace + education + wealth + func014 + func016 + func021 + func023 + func025 + KIMsmoke + mildexercise + modexercise + highexercise + KIMalcohol  + diabetes + bmi + heartdis + KIMsystolic + hypertension + KIManxiety + KIMcynical + KIMdepression + KIMoptimism + KIMposaf + KIMnegaf + KIMparticipation, data = markov.analysis1,  family = "binomial")


## Plot residuals
# Kim et al original model, no imputation
resid.KIM.fig2 <- ggplot(KIM.model, aes(x= predict(KIM.model), y=residuals(KIM.model))) + geom_point(aes(col=residuals(KIM.model), alpha=0.0)) + geom_smooth(method='lm',formula=y~x)

resid.KIM.fig2 + xlab("Predicted Values") + ylab("Residuals")


# Non-imputed data
resid.non.fig4 <- ggplot(non.fig4, aes(x= predict(non.fig4), y=residuals(non.fig4))) + geom_point(aes(col=residuals(non.fig4), alpha=0.0)) + geom_smooth(method='lm',formula=y~x)

resid.non.fig4 + xlab("Predicted Values") + ylab("Residuals")


# Markov imputation
resid.markov.fig4 <- ggplot(markov.fig4, aes(x= predict(markov.fig4), y=residuals(markov.fig4))) + geom_point(aes(col=residuals(markov.fig4), alpha=0.0)) + geom_smooth(method='lm',formula=y~x)

resid.markov.fig4 + xlab("Predicted Values") + ylab("Residuals")


######### FIG 5- FREQUENCY DISTRIBUTIONS OF PURPOSE IN LIFE USING MARKOV CHAIN IMPUTATION, SINGLE IMPUTATION AND BOOSTRAPPING
sub <- data.frame(analysis1[109], markov.analysis1[109], bootstrap.analysis1[109])

#### 2006 baseline covariates, 2008 and 2010 stroke incidence
fig5a <- ggplot(sub, aes(x=KIMpurposeinlife)) + geom_histogram(binwidth=0.2, fill="steelblue4", color="black") + geom_freqpoly(aes(x=sub$KIMpurposeinlife.1), binwidth=0.2, size=0.8) + geom_freqpoly(aes(x=sub$KIMpurposeinlife.2), binwidth=0.2, size=0.2, color="orange1")

fig5a <- fig5a + xlab("Purpose in Life") + ylab("Frequency")


#### 2008 baseline covariates, 2010 and 2012 stroke incidence
fig5b <- ggplot(analysis2, aes(x=KIMpurposeinlife)) + geom_histogram(binwidth=0.2, fill="steelblue4", color="black") + geom_freqpoly(aes(x=bootstrap.analysis2$KIMpurposeinlife), binwidth=0.2, size=0.8) + geom_freqpoly(aes(x=markov.analysis2$KIMpurposeinlife), binwidth=0.2, size=0.2, color="orange1")

fig5b <- fig5b + xlab("Purpose in Life") + ylab("Frequency")



#### 2010 baseline covariates, 2012 and 2014 stroke incidence
fig5c <- ggplot(analysis3, aes(x=KIMpurposeinlife)) + geom_histogram(binwidth=0.2, fill="steelblue4", color="black") + geom_freqpoly(aes(x=bootstrap.analysis3$KIMpurposeinlife), binwidth=0.2, size=0.8) + geom_freqpoly(aes(x=markov.analysis3$KIMpurposeinlife), binwidth=0.2, size=0.2, color="orange1")

fig5c <- fig5c + xlab("Purpose in Life") + ylab("Frequency")



#### 2012 baseline covariates, 2014 and 2016 stroke incidence
fig5d <- ggplot(analysis4, aes(x=KIMpurposeinlife)) + geom_histogram(binwidth=0.2, fill="steelblue4", color="black") + geom_freqpoly(aes(x=bootstrap.analysis4$KIMpurposeinlife), binwidth=0.2, size=0.8) + geom_freqpoly(aes(x=markov.analysis4$KIMpurposeinlife), binwidth=0.2, size=0.2, color="orange1")

fig5d <- fig5d + xlab("Purpose in Life") + ylab("Frequency")


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
fig5e <- ggplot(analysis5, aes(x=purposeinlife)) + geom_histogram(binwidth=0.2, fill="steelblue4", color="black") + geom_freqpoly(aes(x=bootstrap.analysis5$purposeinlife), binwidth=0.2, size=0.8) + geom_freqpoly(aes(x=markov.analysis5$purposeinlife), binwidth=0.2, size=0.2, color="orange1")

fig5e <- fig5e + xlab("Purpose in Life") + ylab("Frequency")


#### 2016 baseline covariates, 2018 stroke incidence
fig5f <- ggplot(analysis6, aes(x=purposeinlife)) + geom_histogram(binwidth=0.2, fill="steelblue4", color="black") + geom_freqpoly(aes(x=bootstrap.analysis6$purposeinlife), binwidth=0.2, size=0.8) + geom_freqpoly(aes(x=markov.analysis6$purposeinlife), binwidth=0.2, size=0.2, color="orange1")

fig5f <- fig5f + xlab("Purpose in Life") + ylab("Frequency")


## Combine outputs
ggarrange(fig5a, fig5b, fig5c, fig5d, fig5e, fig5f)


########### FIG 6 PA

################# APPENDIX
######### DETERMINE NUMBER OF FACTORS
###### Parallel Analysis screen plots

#### 2008 baseline covariates, 2010 and 2012 stroke incidence: 2 factors
Fig6b <- fa.parallel(corMarkov.analysis2, fm = "minres", fa = "fa")



#### 2010 baseline covariates, 2012 and 2014 stroke incidence: 3 factors
Fig6c <- fa.parallel(corMarkov.analysis3, fm = "minres", fa = "fa")



#### 2012 baseline covariates, 2014 and 2016 stroke incidence: 3 factors
Fig6d <- fa.parallel(corMarkov.analysis4, fm = "minres", fa = "fa")


#### 2014 baseline covariates, 2016 and 2018 stroke incidence: 5 factors
Fig6e <- fa.parallel(corMarkov.analysis5, fm = "minres", fa = "fa")


#### 2016 baseline covariates, 2018 stroke incidence: 2 factors
Fig6f <- fa.parallel(corMarkov.analysis6, fm = "minres", fa = "fa")


####### FIG 1 HEATMAP OF COLINEARITY IN PSYCHO-SOCIAL VARIABLES
#### 2008 baseline covariates, 2010 and 2012 stroke incidence
corrplot(corMarkov.analysis2, method = "circle", tl.cex=0.5, tl.col="black")


#### 2010 baseline covariates, 2012 and 2014 stroke incidence
corrplot(corMarkov.analysis3, method = "circle", tl.cex=0.5, tl.col="black")



#### 2012 baseline covariates, 2014 and 2016 stroke incidence
corrplot(corMarkov.analysis4, method = "circle", tl.cex=0.5, tl.col="black")


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
corrplot(corMarkov.analysis5, method = "circle", tl.cex=0.5, tl.col="black")


#### 2016 baseline covariates, 2018 stroke incidence
corrplot(corMarkov.analysis6, method = "circle", tl.cex=0.5, tl.col="black")


