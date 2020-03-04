#### CLEAR WORKSPACE
rm(list = ls(all.names = TRUE))


## Load Workspace of imputed and non-imputed datasets
load("C:/Users/Student/OneDrive - Nexus365/Replication Project/Practical/Specification Curve/Replication Project.RData")


#### CREATE KIM ET AL PSYCHO-SOCIAL VARIABLES ON IMPUTED DATASETS
#### 2006 baseline covariates, 2008 and 2010 stroke incidence

## Markov Chain Imputed Dataset
### Anxiety
markov.analysis1$KIManxiety <- (markov.analysis1$KLB041A + markov.analysis1$KLB041B + markov.analysis1$KLB041C + markov.analysis1$KLB041D + markov.analysis1$KLB041E)/5

### Cynical hostility
markov.analysis1$KIMcynical <- (markov.analysis1$KLB019A + markov.analysis1$KLB019B + markov.analysis1$KLB019C + markov.analysis1$KLB019D + markov.analysis1$KLB019E)/5

### Depression
markov.analysis1$KIMdepression <- (markov.analysis1$KLB027I + markov.analysis1$KLB027J + markov.analysis1$KLB027K + markov.analysis1$KLB027L + markov.analysis1$KLB027M + markov.analysis1$KLB027N)/6

### Optimism
markov.analysis1$KIMoptimism <- (markov.analysis1$KLB019F + markov.analysis1$KLB019G + markov.analysis1$KLB019H + markov.analysis1$KLB019I + markov.analysis1$KLB019J + markov.analysis1$KLB019K)/6

### Positive Affect
markov.analysis1$KIMposaf <- (markov.analysis1$KLB027A + markov.analysis1$KLB027B + markov.analysis1$KLB027C + markov.analysis1$KLB027D + markov.analysis1$KLB027E + markov.analysis1$KLB027F)/6

### Negative Affect
markov.analysis1$KIMnegaf <- (markov.analysis1$KLB027I + markov.analysis1$KLB027J + markov.analysis1$KLB027K + markov.analysis1$KLB027L + markov.analysis1$KLB027M + markov.analysis1$KLB027N)/6

### Social participation
markov.analysis1$KIMparticipation <- (markov.analysis1$KLB001a + markov.analysis1$KLB001b + markov.analysis1$KLB001c + markov.analysis1$KLB001e + markov.analysis1$KLB001f + markov.analysis1$KLB001g)/8

### Purpose in life
markov.analysis1$KIMpurposeinlife <- (markov.analysis1$KLB035A + markov.analysis1$KLB035B + markov.analysis1$KLB035C + markov.analysis1$KLB035D + markov.analysis1$KLB035E + markov.analysis1$KLB035F + markov.analysis1$KLB035G)/7


## Bootstrap Imputed Dataset
### Anxiety
bootstrap.analysis1$KIManxiety <- (bootstrap.analysis1$KLB041A + bootstrap.analysis1$KLB041B + bootstrap.analysis1$KLB041C + bootstrap.analysis1$KLB041D + bootstrap.analysis1$KLB041E)/5

### Cynical hostility
bootstrap.analysis1$KIMcynical <- (bootstrap.analysis1$KLB019A + bootstrap.analysis1$KLB019B + bootstrap.analysis1$KLB019C + bootstrap.analysis1$KLB019D + bootstrap.analysis1$KLB019E)/5

### Depression
bootstrap.analysis1$KIMdepression <- (bootstrap.analysis1$KLB027I + bootstrap.analysis1$KLB027J + bootstrap.analysis1$KLB027K + bootstrap.analysis1$KLB027L + bootstrap.analysis1$KLB027M + bootstrap.analysis1$KLB027N)/6

### Optimism
bootstrap.analysis1$KIMoptimism <- (bootstrap.analysis1$KLB019F + bootstrap.analysis1$KLB019G + bootstrap.analysis1$KLB019H + bootstrap.analysis1$KLB019I + bootstrap.analysis1$KLB019J + bootstrap.analysis1$KLB019K)/6

### Positive Affect
bootstrap.analysis1$KIMposaf <- (bootstrap.analysis1$KLB027A + bootstrap.analysis1$KLB027B + bootstrap.analysis1$KLB027C + bootstrap.analysis1$KLB027D + bootstrap.analysis1$KLB027E + bootstrap.analysis1$KLB027F)/6

### Negative Affect
bootstrap.analysis1$KIMnegaf <- (bootstrap.analysis1$KLB027I + bootstrap.analysis1$KLB027J + bootstrap.analysis1$KLB027K + bootstrap.analysis1$KLB027L + bootstrap.analysis1$KLB027M + bootstrap.analysis1$KLB027N)/6

### Social participation
bootstrap.analysis1$KIMparticipation <- (bootstrap.analysis1$KLB001a + bootstrap.analysis1$KLB001b + bootstrap.analysis1$KLB001c + bootstrap.analysis1$KLB001e + bootstrap.analysis1$KLB001f + bootstrap.analysis1$KLB001g)/8

### Purpose in life
bootstrap.analysis1$KIMpurposeinlife <- (bootstrap.analysis1$KLB035A + bootstrap.analysis1$KLB035B + bootstrap.analysis1$KLB035C + bootstrap.analysis1$KLB035D + bootstrap.analysis1$KLB035E + bootstrap.analysis1$KLB035F + bootstrap.analysis1$KLB035G)/7




#### 2008 baseline covariates, 2010 and 2012 stroke incidence
## Markov Chain Imputed Dataset
### Anxiety
markov.analysis2$KIManxiety <- (markov.analysis2$LLB041A + markov.analysis2$LLB041B + markov.analysis2$LLB041C + markov.analysis2$LLB041D + markov.analysis2$LLB041E)/5

### Cynical hostility
markov.analysis2$KIMcynical <- (markov.analysis2$LLB019A + markov.analysis2$LLB019B + markov.analysis2$LLB019C + markov.analysis2$LLB019D + markov.analysis2$LLB019E)/5

### Depression
markov.analysis2$KIMdepression <- (markov.analysis2$LLB027I + markov.analysis2$LLB027J + markov.analysis2$LLB027K + markov.analysis2$LLB027L + markov.analysis2$LLB027M + markov.analysis2$LLB027N)/6

### Optimism
markov.analysis2$KIMoptimism <- (markov.analysis2$LLB019F + markov.analysis2$LLB019G + markov.analysis2$LLB019H + markov.analysis2$LLB019I + markov.analysis2$LLB019J + markov.analysis2$LLB019K)/6

### Positive Affect
markov.analysis2$KIMposaf <- (markov.analysis2$LLB027A + markov.analysis2$LLB027B + markov.analysis2$LLB027C + markov.analysis2$LLB027D + markov.analysis2$LLB027E + markov.analysis2$LLB027F)/6

### Negative Affect
markov.analysis2$KIMnegaf <- (markov.analysis2$LLB027I + markov.analysis2$LLB027J + markov.analysis2$LLB027K + markov.analysis2$LLB027L + markov.analysis2$LLB027M + markov.analysis2$LLB027N)/6

### Social participation
markov.analysis2$KIMparticipation <- (markov.analysis2$LLB001a + markov.analysis2$LLB001b + markov.analysis2$LLB001c + markov.analysis2$LLB001e + markov.analysis2$LLB001f + markov.analysis2$LLB001g)/8

### Purpose in life
markov.analysis2$KIMpurposeinlife <- (markov.analysis2$LLB035A + markov.analysis2$LLB035B + markov.analysis2$LLB035C + markov.analysis2$LLB035D + markov.analysis2$LLB035E + markov.analysis2$LLB035F + markov.analysis2$LLB035G)/7



## Bootstrap Imputed Dataset
### Anxiety
bootstrap.analysis2$KIManxiety <- (bootstrap.analysis2$LLB041A + bootstrap.analysis2$LLB041B + bootstrap.analysis2$LLB041C + bootstrap.analysis2$LLB041D + bootstrap.analysis2$LLB041E)/5

### Cynical hostility
bootstrap.analysis2$KIMcynical <- (bootstrap.analysis2$LLB019A + bootstrap.analysis2$LLB019B + bootstrap.analysis2$LLB019C + bootstrap.analysis2$LLB019D + bootstrap.analysis2$LLB019E)/5

### Depression
bootstrap.analysis2$KIMdepression <- (bootstrap.analysis2$LLB027I + bootstrap.analysis2$LLB027J + bootstrap.analysis2$LLB027K + bootstrap.analysis2$LLB027L + bootstrap.analysis2$LLB027M + bootstrap.analysis2$LLB027N)/6

### Optimism
bootstrap.analysis2$KIMoptimism <- (bootstrap.analysis2$LLB019F + bootstrap.analysis2$LLB019G + bootstrap.analysis2$LLB019H + bootstrap.analysis2$LLB019I + bootstrap.analysis2$LLB019J + bootstrap.analysis2$LLB019K)/6

### Positive Affect
bootstrap.analysis2$KIMposaf <- (bootstrap.analysis2$LLB027A + bootstrap.analysis2$LLB027B + bootstrap.analysis2$LLB027C + bootstrap.analysis2$LLB027D + bootstrap.analysis2$LLB027E + bootstrap.analysis2$LLB027F)/6

### Negative Affect
bootstrap.analysis2$KIMnegaf <- (bootstrap.analysis2$LLB027I + bootstrap.analysis2$LLB027J + bootstrap.analysis2$LLB027K + bootstrap.analysis2$LLB027L + bootstrap.analysis2$LLB027M + bootstrap.analysis2$LLB027N)/6

### Social participation
bootstrap.analysis2$KIMparticipation <- (bootstrap.analysis2$LLB001a + bootstrap.analysis2$LLB001b + bootstrap.analysis2$LLB001c + bootstrap.analysis2$LLB001e + bootstrap.analysis2$LLB001f + bootstrap.analysis2$LLB001g)/8

### Purpose in life
bootstrap.analysis2$KIMpurposeinlife <- (bootstrap.analysis2$LLB035A + bootstrap.analysis2$LLB035B + bootstrap.analysis2$LLB035C + bootstrap.analysis2$LLB035D + bootstrap.analysis2$LLB035E + bootstrap.analysis2$LLB035F + bootstrap.analysis2$LLB035G)/7



#### 2010 baseline covariates, 2012 and 2014 stroke incidence
## Markov Chain Imputed Dataset
### Anxiety
markov.analysis3$KIManxiety <- (markov.analysis3$MLB041A + markov.analysis3$MLB041B + markov.analysis3$MLB041C + markov.analysis3$MLB041D + markov.analysis3$MLB041E)/5

### Cynical hostility
markov.analysis3$KIMcynical <- (markov.analysis3$MLB019A + markov.analysis3$MLB019B + markov.analysis3$MLB019C + markov.analysis3$MLB019D + markov.analysis3$MLB019E)/5

### Depression
markov.analysis3$KIMdepression <- (markov.analysis3$MLB027I + markov.analysis3$MLB027J + markov.analysis3$MLB027K + markov.analysis3$MLB027L + markov.analysis3$MLB027M + markov.analysis3$MLB027N)/6

### Optimism
markov.analysis3$KIMoptimism <- (markov.analysis3$MLB019F + markov.analysis3$MLB019G + markov.analysis3$MLB019H + markov.analysis3$MLB019I + markov.analysis3$MLB019J + markov.analysis3$MLB019K)/6

### Positive Affect
markov.analysis3$KIMposaf <- (markov.analysis3$MLB027A + markov.analysis3$MLB027B + markov.analysis3$MLB027C + markov.analysis3$MLB027D + markov.analysis3$MLB027E + markov.analysis3$MLB027F)/6

### Negative Affect
markov.analysis3$KIMnegaf <- (markov.analysis3$MLB027I + markov.analysis3$MLB027J + markov.analysis3$MLB027K + markov.analysis3$MLB027L + markov.analysis3$MLB027M + markov.analysis3$MLB027N)/6

### Social participation
markov.analysis3$KIMparticipation <- (markov.analysis3$MLB001A + markov.analysis3$MLB001B + markov.analysis3$MLB001C + markov.analysis3$MLB001E + markov.analysis3$MLB001F + markov.analysis3$MLB001G)/8

### Purpose in life
markov.analysis3$KIMpurposeinlife <- (markov.analysis3$MLB035A + markov.analysis3$MLB035B + markov.analysis3$MLB035C + markov.analysis3$MLB035D + markov.analysis3$MLB035E + markov.analysis3$MLB035F + markov.analysis3$MLB035G)/7


## Bootstrap Imputed Dataset
### Anxiety
bootstrap.analysis3$KIManxiety <- (bootstrap.analysis3$MLB041A + bootstrap.analysis3$MLB041B + bootstrap.analysis3$MLB041C + bootstrap.analysis3$MLB041D + bootstrap.analysis3$MLB041E)/5

### Cynical hostility
bootstrap.analysis3$KIMcynical <- (bootstrap.analysis3$MLB019A + bootstrap.analysis3$MLB019B + bootstrap.analysis3$MLB019C + bootstrap.analysis3$MLB019D + bootstrap.analysis3$MLB019E)/5

### Depression
bootstrap.analysis3$KIMdepression <- (bootstrap.analysis3$MLB027I + bootstrap.analysis3$MLB027J + bootstrap.analysis3$MLB027K + bootstrap.analysis3$MLB027L + bootstrap.analysis3$MLB027M + bootstrap.analysis3$MLB027N)/6

### Optimism
bootstrap.analysis3$KIMoptimism <- (bootstrap.analysis3$MLB019F + bootstrap.analysis3$MLB019G + bootstrap.analysis3$MLB019H + bootstrap.analysis3$MLB019I + bootstrap.analysis3$MLB019J + bootstrap.analysis3$MLB019K)/6

### Positive Affect
bootstrap.analysis3$KIMposaf <- (bootstrap.analysis3$MLB027A + bootstrap.analysis3$MLB027B + bootstrap.analysis3$MLB027C + bootstrap.analysis3$MLB027D + bootstrap.analysis3$MLB027E + bootstrap.analysis3$MLB027F)/6

### Negative Affect
bootstrap.analysis3$KIMnegaf <- (bootstrap.analysis3$MLB027I + bootstrap.analysis3$MLB027J + bootstrap.analysis3$MLB027K + bootstrap.analysis3$MLB027L + bootstrap.analysis3$MLB027M + bootstrap.analysis3$MLB027N)/6

### Social participation
bootstrap.analysis3$KIMparticipation <- (bootstrap.analysis3$MLB001A + bootstrap.analysis3$MLB001B + bootstrap.analysis3$MLB001C + bootstrap.analysis3$MLB001E + bootstrap.analysis3$MLB001F + bootstrap.analysis3$MLB001G)/8

### Purpose in life
bootstrap.analysis3$KIMpurposeinlife <- (bootstrap.analysis3$MLB035A + bootstrap.analysis3$MLB035B + bootstrap.analysis3$MLB035C + bootstrap.analysis3$MLB035D + bootstrap.analysis3$MLB035E + bootstrap.analysis3$MLB035F + bootstrap.analysis3$MLB035G)/7



#### 2012 baseline covariates, 2014 and 2016 stroke incidence
## Markov Chain Imputed Dataset
markov.analysis4$KIManxiety <- (markov.analysis4$NLB041A + markov.analysis4$NLB041B + markov.analysis4$NLB041C + markov.analysis4$NLB041D + markov.analysis4$NLB041E)/5

### Cynical hostility
markov.analysis4$KIMcynical <- (markov.analysis4$NLB019A + markov.analysis4$NLB019B + markov.analysis4$NLB019C + markov.analysis4$NLB019D + markov.analysis4$NLB019E)/5

### Depression
markov.analysis4$KIMdepression <- (markov.analysis4$NLB027I + markov.analysis4$NLB027J + markov.analysis4$NLB027K + markov.analysis4$NLB027L + markov.analysis4$NLB027M + markov.analysis4$NLB027N)/6

### Optimism
markov.analysis4$KIMoptimism <- (markov.analysis4$NLB019F + markov.analysis4$NLB019G + markov.analysis4$NLB019H + markov.analysis4$NLB019I + markov.analysis4$NLB019J + markov.analysis4$NLB019K)/6

### Positive Affect
markov.analysis4$KIMposaf <- (markov.analysis4$NLB027A + markov.analysis4$NLB027B + markov.analysis4$NLB027C + markov.analysis4$NLB027D + markov.analysis4$NLB027E + markov.analysis4$NLB027F)/6

### Negative Affect
markov.analysis4$KIMnegaf <- (markov.analysis4$NLB027I + markov.analysis4$NLB027J + markov.analysis4$NLB027K + markov.analysis4$NLB027L + markov.analysis4$NLB027M + markov.analysis4$NLB027N)/6

### Social participation
markov.analysis4$KIMparticipation <- (markov.analysis4$NLB001A + markov.analysis4$NLB001B + markov.analysis4$NLB001C + markov.analysis4$NLB001E + markov.analysis4$NLB001F + markov.analysis4$NLB001G)/8

### Purpose in life
markov.analysis4$KIMpurposeinlife <- (markov.analysis4$NLB035A + markov.analysis4$NLB035B + markov.analysis4$NLB035C + markov.analysis4$NLB035D + markov.analysis4$NLB035E + markov.analysis4$NLB035F + markov.analysis4$NLB035G)/7



## Bootstrap Imputed Dataset
## Anxiety
bootstrap.analysis4$KIManxiety <- (bootstrap.analysis4$NLB041A + bootstrap.analysis4$NLB041B + bootstrap.analysis4$NLB041C + bootstrap.analysis4$NLB041D + bootstrap.analysis4$NLB041E)/5

### Cynical hostility
bootstrap.analysis4$KIMcynical <- (bootstrap.analysis4$NLB019A + bootstrap.analysis4$NLB019B + bootstrap.analysis4$NLB019C + bootstrap.analysis4$NLB019D + bootstrap.analysis4$NLB019E)/5

### Depression
bootstrap.analysis4$KIMdepression <- (bootstrap.analysis4$NLB027I + bootstrap.analysis4$NLB027J + bootstrap.analysis4$NLB027K + bootstrap.analysis4$NLB027L + bootstrap.analysis4$NLB027M + bootstrap.analysis4$NLB027N)/6

### Optimism
bootstrap.analysis4$KIMoptimism <- (bootstrap.analysis4$NLB019F + bootstrap.analysis4$NLB019G + bootstrap.analysis4$NLB019H + bootstrap.analysis4$NLB019I + bootstrap.analysis4$NLB019J + bootstrap.analysis4$NLB019K)/6

### Positive Affect
bootstrap.analysis4$KIMposaf <- (bootstrap.analysis4$NLB027A + bootstrap.analysis4$NLB027B + bootstrap.analysis4$NLB027C + bootstrap.analysis4$NLB027D + bootstrap.analysis4$NLB027E + bootstrap.analysis4$NLB027F)/6

### Negative Affect
bootstrap.analysis4$KIMnegaf <- (bootstrap.analysis4$NLB027I + bootstrap.analysis4$NLB027J + bootstrap.analysis4$NLB027K + bootstrap.analysis4$NLB027L + bootstrap.analysis4$NLB027M + bootstrap.analysis4$NLB027N)/6

### Social participation
bootstrap.analysis4$KIMparticipation <- (bootstrap.analysis4$NLB001A + bootstrap.analysis4$NLB001B + bootstrap.analysis4$NLB001C + bootstrap.analysis4$NLB001E + bootstrap.analysis4$NLB001F + bootstrap.analysis4$NLB001G)/8

### Purpose in life
bootstrap.analysis4$KIMpurposeinlife <- (bootstrap.analysis4$NLB035A + bootstrap.analysis4$NLB035B + bootstrap.analysis4$NLB035C + bootstrap.analysis4$NLB035D + bootstrap.analysis4$NLB035E + bootstrap.analysis4$NLB035F + bootstrap.analysis4$NLB035G)/7


#### 2014 baseline covariates, 2016 and 2018 stroke incidence
## Markov Chain Imputed Dataset
## Depression
markov.analysis5$depression <- (markov.analysis5$OLB026S	+ markov.analysis5$OLB018H + markov.analysis5$OLB026K)/3

## Optimism
markov.analysis5$optimism <- (markov.analysis5$OLB018A + markov.analysis5$OLB018B + markov.analysis5$OLB018C + markov.analysis5$OLB018D + markov.analysis5$OLB018E + markov.analysis5$OLB018F)/6

## Social participation
markov.analysis5$participation <- (markov.analysis5$OLB001I + markov.analysis5$OLB001R + markov.analysis5$OLB001N + markov.analysis5$OLB012B)/4

## Purpose in life
markov.analysis5$purposeinlife <- (markov.analysis5$OLB033A  + markov.analysis5$OLB033B + markov.analysis5$OLB033C + markov.analysis5$OLB033D + markov.analysis5$OLB033F + markov.analysis5$OLB033G + markov.analysis5$OLB033E)/7

## Anxiety
markov.analysis5$anxiety <- (markov.analysis5$OLB031L + markov.analysis5$OLB026R)/2

## Frustration
markov.analysis5$frustration <- (markov.analysis5$OLB026J + markov.analysis5$OLB026X)/2



## Bootstrap Imputed Dataset
## Depression
bootstrap.analysis5$depression <- (bootstrap.analysis5$OLB026S	+ bootstrap.analysis5$OLB018H + bootstrap.analysis5$OLB026K)/3

## Optimism
bootstrap.analysis5$optimism <- (bootstrap.analysis5$OLB018A + bootstrap.analysis5$OLB018B + bootstrap.analysis5$OLB018C + bootstrap.analysis5$OLB018D + bootstrap.analysis5$OLB018E + bootstrap.analysis5$OLB018F)/6

## Social participation
bootstrap.analysis5$participation <- (bootstrap.analysis5$OLB001I + bootstrap.analysis5$OLB001R + bootstrap.analysis5$OLB001N + bootstrap.analysis5$OLB012B)/4

## Purpose in life
bootstrap.analysis5$purposeinlife <- (bootstrap.analysis5$OLB033A  + bootstrap.analysis5$OLB033B + bootstrap.analysis5$OLB033C + bootstrap.analysis5$OLB033D + bootstrap.analysis5$OLB033F + bootstrap.analysis5$OLB033G + bootstrap.analysis5$OLB033E)/7

## Anxiety
bootstrap.analysis5$anxiety <- (bootstrap.analysis5$OLB031L + bootstrap.analysis5$OLB026R)/2

## Frustration
bootstrap.analysis5$frustration <- (bootstrap.analysis5$OLB026J + bootstrap.analysis5$OLB026X)/2



#### 2016 baseline covariates, 2018 stroke incidence
## Markov Chain Imputed Dataset
## Depression
markov.analysis6$depression <- (markov.analysis6$PLB018H + markov.analysis6$PLB026K)/2

## Optimism
markov.analysis6$optimism <- (markov.analysis6$PLB018A + markov.analysis6$PLB018B + markov.analysis6$PLB018C + markov.analysis6$PLB018D + markov.analysis6$PLB018E + markov.analysis6$PLB018F)/6

## Social participation
markov.analysis6$participation <- (markov.analysis6$PLB001I + markov.analysis6$PLB001R + markov.analysis6$PLB001N + markov.analysis6$PLB012B)/4

## Purpose in life
markov.analysis6$purposeinlife <- (markov.analysis6$PLB033A  + markov.analysis6$PLB033C + markov.analysis6$PLB033D + markov.analysis6$PLB033F + markov.analysis6$PLB033G + markov.analysis6$PLB033E)/7

## Anxiety
markov.analysis6$anxiety <- (markov.analysis6$PLB031L + markov.analysis6$PLB026R)/2

## Frustration
markov.analysis6$frustration <- (markov.analysis6$PLB026J + markov.analysis6$PLB026X)/2



## Bootstrap Imputed Dataset
## Depression
bootstrap.analysis6$depression <- (bootstrap.analysis6$PLB018H + bootstrap.analysis6$PLB026K)/2

## Optimism
bootstrap.analysis6$optimism <- (bootstrap.analysis6$PLB018A + bootstrap.analysis6$PLB018B + bootstrap.analysis6$PLB018C + bootstrap.analysis6$PLB018D + bootstrap.analysis6$PLB018E + bootstrap.analysis6$PLB018F)/6

## Social participation
bootstrap.analysis6$participation <- (bootstrap.analysis6$PLB001I + bootstrap.analysis6$PLB001R + bootstrap.analysis6$PLB001N + bootstrap.analysis6$PLB012B)/4

## Purpose in life
bootstrap.analysis6$purposeinlife <- (bootstrap.analysis6$PLB033A  + bootstrap.analysis6$PLB033C + bootstrap.analysis6$PLB033D + bootstrap.analysis6$PLB033F + bootstrap.analysis6$PLB033G + bootstrap.analysis6$PLB033E)/7

## Anxiety
bootstrap.analysis6$anxiety <- (bootstrap.analysis6$PLB031L + bootstrap.analysis6$PLB026R)/2

## Frustration
bootstrap.analysis6$frustration <- (bootstrap.analysis6$PLB026J + bootstrap.analysis6$PLB026X)/2




#### RENAME PSYCHO-SOCIAL VARIABLES
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
colnames(analysis1)[11] <- "Fear Worst"

analysis1 <- rename.variable(analysis1, "KLB019A", "Dislike helping")
analysis1  <-  rename.variable(analysis1,  "KLB027I", "Nothing cheer")
analysis1  <-  rename.variable(analysis1,  "KLB019F", "Things go wrong")
analysis1  <-  rename.variable(analysis1,  "KLB027A", "Cheerful")
analysis1  <-  rename.variable(analysis1,  "KLB001a", "Newspaper")
analysis1  <-  rename.variable(analysis1,  "KLB035A", "Making plans")
analysis1  <-  rename.variable(analysis1,  "KLB035B", "Activities trivial")
analysis1  <-  rename.variable(analysis1,  "KLB035C", "Active planning")
analysis1  <-  rename.variable(analysis1,  "KLB035D", "Accomplish")
analysis1  <-  rename.variable(analysis1,  "KLB035F", "One day a time")
analysis1  <-  rename.variable(analysis1,  "KLB035G", "Sense of direction")
analysis1  <-  rename.variable(analysis1,  "KLB041B", "Nervous")
analysis1  <-  rename.variable(analysis1,  "KLB027N", "Nervous")
analysis1  <-  rename.variable(analysis1,  "KLB019B", "Unfair profit")
analysis1  <-  rename.variable(analysis1,  "KLB027J", "Hopeless")
analysis1  <-  rename.variable(analysis1,  "KLB019G", "Optimistic about future")
analysis1  <-  rename.variable(analysis1,  "KLB027B", "Good spririts")
analysis1  <-  rename.variable(analysis1,  "KLB001b", "Has hobby")
analysis1  <-  rename.variable(analysis1,  "KLB027J", "Frustrated")
analysis1  <-  rename.variable(analysis1,  "KLB041C", "Trembling hands")
analysis1  <-  rename.variable(analysis1,  "KLB019C", "People care")
analysis1  <-  rename.variable(analysis1,  "KLB027K", "Restless")
analysis1  <-  rename.variable(analysis1,  "KLB019H", "Expect best")
analysis1  <-  rename.variable(analysis1,  "KLB027C", "Happy")
analysis1  <-  rename.variable(analysis1,  "KLB001c", "Vacation")
analysis1  <-  rename.variable(analysis1,  "KLB041D", "Fear dying")
analysis1  <-  rename.variable(analysis1,  "KLB019D", "Lie to get ahead")
analysis1  <-  rename.variable(analysis1,  "KLB027L", "Everything an effort")
analysis1  <-  rename.variable(analysis1,  "KLB019I", "Expect good things")
analysis1  <-  rename.variable(analysis1,  "KLB027D", "Calm")
analysis1  <-  rename.variable(analysis1,  "KLB041E", "Faint")
analysis1  <-  rename.variable(analysis1,  "KLB019E", "Sceptical people nice")
analysis1  <-  rename.variable(analysis1,  "KLB027M", "Worthless")
analysis1  <-  rename.variable(analysis1,  "KLB019J", "Things won't go well")
analysis1  <-  rename.variable(analysis1,  "KLB027E", "Satisfied")
analysis1  <-  rename.variable(analysis1,  "KLB001e", "Daytrip")
analysis1  <-  rename.variable(analysis1,  "KLB027M", "Worthless")
analysis1  <-  rename.variable(analysis1,  "KLB035E", "Done all to do in life")
analysis1  <-  rename.variable(analysis1,  "KLB019K", "Can't count on good")
analysis1  <-  rename.variable(analysis1,  "KLB027F", "Full of life")
analysis1  <-  rename.variable(analysis1,  "KLB001f", "Use Internet")
analysis1  <-  rename.variable(analysis1,  "KLB001g", "Use cell phone")


analysis1  <-  rename.variable(analysis1,  "KLB001A", "Newspaper")
analysis1  <-  rename.variable(analysis1,  "KLB001B", "Has hobby")
analysis1  <-  rename.variable(analysis1,  "KLB001C", "Vacation")
analysis1  <-  rename.variable(analysis1,  "KLB001F", "Use Internet")
analysis1  <-  rename.variable(analysis1,  "KLB001G", "Use cell phone")
analysis1  <-  rename.variable(analysis1,  "KLB027J.1", "Hopeless")
analysis1  <-  rename.variable(analysis1,  "KLB027M.1", "Worthless")
analysis1  <-  rename.variable(analysis1,  "KLB001h", "Hopeless")


#### 2008 baseline covariates, 2010 and 2012 stroke incidence
analysis2  <-  rename.variable(analysis2,  "LLB041A", "Fear worst")
analysis2  <-  rename.variable(analysis2,  "LLB019A", "Dislike helping")
analysis2  <-  rename.variable(analysis2,  "LLB027I", "Nothing cheer")
analysis2  <-  rename.variable(analysis2,  "LLB019F", "Things go wrong")
analysis2  <-  rename.variable(analysis2,  "LLB027A", "Cheerful")
analysis2  <-  rename.variable(analysis2,  "LLB001a", "Newspaper")
analysis2  <-  rename.variable(analysis2,  "LLB035A", "Making plans")
analysis2  <-  rename.variable(analysis2,  "LLB035B", "Activities trivial")
analysis2  <-  rename.variable(analysis2,  "LLB035C", "Active planning")
analysis2  <-  rename.variable(analysis2,  "LLB035D", "Accomplish")
analysis2  <-  rename.variable(analysis2,  "LLB035F", "One day a time")
analysis2  <-  rename.variable(analysis2,  "LLB035G", "Sense of direction")

analysis2  <-  rename.variable(analysis2,  "LLB041B", "Nervous")
analysis2  <-  rename.variable(analysis2,  "LLB027N", "Nervous")
analysis2  <-  rename.variable(analysis2,  "LLB019B", "Unfair profit")
analysis2  <-  rename.variable(analysis2,  "LLB027J", "Hopeless")
analysis2  <-  rename.variable(analysis2,  "LLB019G", "Optimistic about future")
analysis2  <-  rename.variable(analysis2,  "LLB027B", "Good spririts")
analysis2  <-  rename.variable(analysis2,  "LLB001b", "Has hobby")
analysis2  <-  rename.variable(analysis2,  "LLB027J", "Frustrated")

analysis2  <-  rename.variable(analysis2,  "LLB041C", "Trembling hands")
analysis2  <-  rename.variable(analysis2,  "LLB019C", "People care")
analysis2  <-  rename.variable(analysis2,  "LLB027K", "Restless")
analysis2  <-  rename.variable(analysis2,  "LLB019H", "Expect best")
analysis2  <-  rename.variable(analysis2,  "LLB027C", "Happy")
analysis2  <-  rename.variable(analysis2,  "LLB001c", "Vacation")
analysis2  <-  rename.variable(analysis2,  "LLB041D", "Fear dying")
analysis2  <-  rename.variable(analysis2,  "LLB019D", "Lie to get ahead")
analysis2  <-  rename.variable(analysis2,  "LLB027L", "Everything an effort")
analysis2  <-  rename.variable(analysis2,  "LLB019I", "Expect good things")
analysis2  <-  rename.variable(analysis2,  "LLB027D", "Calm")
analysis2  <-  rename.variable(analysis2,  "LLB041E", "Faint")
analysis2  <-  rename.variable(analysis2,  "LLB019E", "Sceptical people nice")
analysis2  <-  rename.variable(analysis2,  "LLB027M", "Worthless")
analysis2  <-  rename.variable(analysis2,  "LLB019J", "Things won't go well")
analysis2  <-  rename.variable(analysis2,  "LLB027E", "Satisfied")
analysis2  <-  rename.variable(analysis2,  "LLB001e", "Daytrip")
analysis2  <-  rename.variable(analysis2,  "LLB027M", "Worthless")
analysis2  <-  rename.variable(analysis2,  "LLB035E", "Done all to do in life")
analysis2  <-  rename.variable(analysis2,  "LLB019K", "Can't count on good")
analysis2  <-  rename.variable(analysis2,  "LLB027F", "Full of life")
analysis2  <-  rename.variable(analysis2,  "LLB001f", "Use Internet")
analysis2  <-  rename.variable(analysis2,  "LLB001g", "Use cell phone")

analysis2  <-  rename.variable(analysis2,  "LLB027J.1", "Hopeless")
analysis2  <-  rename.variable(analysis2,  "LLB027M.1", "Worthless")
analysis2  <-  rename.variable(analysis2,  "LLB001A", "Newspaper")
analysis2  <-  rename.variable(analysis2,  "LLB001B", "Has hobby")
analysis2  <-  rename.variable(analysis2,  "LLB001C", "Vacation")
analysis2  <-  rename.variable(analysis2,  "LLB001F", "Use Internet")
analysis2  <-  rename.variable(analysis2,  "LLB001G", "Use cell phone")




#### 2010 baseline covariates, 2012 and 2014 stroke incidence
analysis3  <-  rename.variable(analysis3,  "MLB041A", "Fear worst")
analysis3  <-  rename.variable(analysis3,  "MLB019A", "Dislike helping")
analysis3  <-  rename.variable(analysis3,  "MLB027I", "Nothing cheer")
analysis3  <-  rename.variable(analysis3,  "MLB019F", "Things go wrong")
analysis3  <-  rename.variable(analysis3,  "MLB027A", "Cheerful")
analysis3  <-  rename.variable(analysis3,  "MLB001a", "Newspaper")
analysis3  <-  rename.variable(analysis3,  "MLB035A", "Making plans")
analysis3  <-  rename.variable(analysis3,  "MLB035B", "Activities trivial")
analysis3  <-  rename.variable(analysis3,  "MLB035C", "Active planning")
analysis3  <-  rename.variable(analysis3,  "MLB035D", "Accomplish")
analysis3  <-  rename.variable(analysis3,  "MLB035F", "One day a time")
analysis3  <-  rename.variable(analysis3,  "MLB035G", "Sense of direction")

analysis3  <-  rename.variable(analysis3,  "MLB041B", "Nervous")
analysis3  <-  rename.variable(analysis3,  "MLB027N", "Nervous")
analysis3  <-  rename.variable(analysis3,  "MLB019B", "Unfair profit")
analysis3  <-  rename.variable(analysis3,  "MLB027J", "Hopeless")
analysis3  <-  rename.variable(analysis3,  "MLB019G", "Optimistic about future")
analysis3  <-  rename.variable(analysis3,  "MLB027B", "Good spririts")
analysis3  <-  rename.variable(analysis3,  "MLB001b", "Has hobby")
analysis3  <-  rename.variable(analysis3,  "MLB027J", "Frustrated")

analysis3  <-  rename.variable(analysis3,  "MLB041C", "Trembling hands")
analysis3  <-  rename.variable(analysis3,  "MLB019C", "People care")
analysis3  <-  rename.variable(analysis3,  "MLB027K", "Restless")
analysis3  <-  rename.variable(analysis3,  "MLB019H", "Expect best")
analysis3  <-  rename.variable(analysis3,  "MLB027C", "Happy")
analysis3  <-  rename.variable(analysis3,  "MLB001c", "Vacation")
analysis3  <-  rename.variable(analysis3,  "MLB041D", "Fear dying")
analysis3  <-  rename.variable(analysis3,  "MLB019D", "Lie to get ahead")
analysis3  <-  rename.variable(analysis3,  "MLB027L", "Everything an effort")
analysis3  <-  rename.variable(analysis3,  "MLB019I", "Expect good things")
analysis3  <-  rename.variable(analysis3,  "MLB027D", "Calm")
analysis3  <-  rename.variable(analysis3,  "MLB041E", "Faint")
analysis3  <-  rename.variable(analysis3,  "MLB019E", "Sceptical people nice")
analysis3  <-  rename.variable(analysis3,  "MLB027M", "Worthless")
analysis3  <-  rename.variable(analysis3,  "MLB019J", "Things won't go well")
analysis3  <-  rename.variable(analysis3,  "MLB027E", "Satisfied")
analysis3  <-  rename.variable(analysis3,  "MLB001e", "Daytrip")
analysis3  <-  rename.variable(analysis3,  "MLB027M", "Worthless")
analysis3  <-  rename.variable(analysis3,  "MLB035E", "Done all to do in life")
analysis3  <-  rename.variable(analysis3,  "MLB019K", "Can't count on good")
analysis3  <-  rename.variable(analysis3,  "MLB027F", "Full of life")
analysis3  <-  rename.variable(analysis3,  "MLB001f", "Use Internet")
analysis3  <-  rename.variable(analysis3,  "MLB001g", "Use cell phone")


analysis3  <-  rename.variable(analysis3,  "MLB001A", "Newspaper")
analysis3  <-  rename.variable(analysis3,  "MLB001B", "Has hobby")
analysis3  <-  rename.variable(analysis3,  "MLB027J.1", "Hopeless")
analysis3  <-  rename.variable(analysis3,  "MLB001E", "Daytrip")
analysis3  <-  rename.variable(analysis3,  "MLB027M.1", "Worthless")
analysis3  <-  rename.variable(analysis3,  "MLB001F", "Use Internet")
analysis3  <-  rename.variable(analysis3,  "MLB001G", "Use cell phone")



#### 2012 baseline covariates, 2014 and 2016 stroke incidence
analysis4  <-  rename.variable(analysis4,  "NLB041A", "Fear worst")
analysis4  <-  rename.variable(analysis4,  "NLB019A", "Dislike helping")
analysis4  <-  rename.variable(analysis4,  "NLB027I", "Nothing cheer")
analysis4  <-  rename.variable(analysis4,  "NLB019F", "Things go wrong")
analysis4  <-  rename.variable(analysis4,  "NLB027A", "Cheerful")
analysis4  <-  rename.variable(analysis4,  "NLB001a", "Newspaper")
analysis4  <-  rename.variable(analysis4,  "NLB035A", "Making plans")
analysis4  <-  rename.variable(analysis4,  "NLB035B", "Activities trivial")
analysis4  <-  rename.variable(analysis4,  "NLB035C", "Active planning")
analysis4  <-  rename.variable(analysis4,  "NLB035D", "Accomplish")
analysis4  <-  rename.variable(analysis4,  "NLB035F", "One day a time")
analysis4  <-  rename.variable(analysis4,  "NLB035G", "Sense of direction")

analysis4  <-  rename.variable(analysis4,  "NLB041B", "Nervous")
analysis4  <-  rename.variable(analysis4,  "NLB027N", "Nervous")
analysis4  <-  rename.variable(analysis4,  "NLB019B", "Unfair profit")
analysis4  <-  rename.variable(analysis4,  "NLB027J", "Hopeless")
analysis4  <-  rename.variable(analysis4,  "NLB019G", "Optimistic about future")
analysis4  <-  rename.variable(analysis4,  "NLB027B", "Good spririts")
analysis4  <-  rename.variable(analysis4,  "NLB001b", "Has hobby")
analysis4  <-  rename.variable(analysis4,  "NLB027J", "Frustrated")

analysis4  <-  rename.variable(analysis4,  "NLB041C", "Trembling hands")
analysis4  <-  rename.variable(analysis4,  "NLB019C", "People care")
analysis4  <-  rename.variable(analysis4,  "NLB027K", "Restless")
analysis4  <-  rename.variable(analysis4,  "NLB019H", "Expect best")
analysis4  <-  rename.variable(analysis4,  "NLB027C", "Happy")
analysis4  <-  rename.variable(analysis4,  "NLB001c", "Vacation")
analysis4  <-  rename.variable(analysis4,  "NLB041D", "Fear dying")
analysis4  <-  rename.variable(analysis4,  "NLB019D", "Lie to get ahead")
analysis4  <-  rename.variable(analysis4,  "NLB027L", "Everything an effort")
analysis4  <-  rename.variable(analysis4,  "NLB019I", "Expect good things")
analysis4  <-  rename.variable(analysis4,  "NLB027D", "Calm")
analysis4  <-  rename.variable(analysis4,  "NLB041E", "Faint")
analysis4  <-  rename.variable(analysis4,  "NLB019E", "Sceptical people nice")
analysis4  <-  rename.variable(analysis4,  "NLB027M", "Worthless")
analysis4  <-  rename.variable(analysis4,  "NLB019J", "Things won't go well")
analysis4  <-  rename.variable(analysis4,  "NLB027E", "Satisfied")
analysis4  <-  rename.variable(analysis4,  "NLB001e", "Daytrip")
analysis4  <-  rename.variable(analysis4,  "NLB027M", "Worthless")
analysis4  <-  rename.variable(analysis4,  "NLB035E", "Done all to do in life")
analysis4  <-  rename.variable(analysis4,  "NLB019K", "Can't count on good")
analysis4  <-  rename.variable(analysis4,  "NLB027F", "Full of life")
analysis4  <-  rename.variable(analysis4,  "NLB001f", "Use Internet")
analysis4  <-  rename.variable(analysis4,  "NLB001g", "Use cell phone")


analysis4  <-  rename.variable(analysis4,  "NLB001A", "Newspaper")
analysis4  <-  rename.variable(analysis4,  "NLB001B", "Has hobby")
analysis4  <-  rename.variable(analysis4,  "NLB027J.1", "Hopeless")
analysis4  <-  rename.variable(analysis4,  "NLB001E", "Daytrip")
analysis4  <-  rename.variable(analysis4,  "NLB027M.1", "Worthless")
analysis4  <-  rename.variable(analysis4,  "NLB001F", "Use Internet")
analysis4  <-  rename.variable(analysis4,  "NLB001G", "Use cell phone")



#### 2014 baseline covariates, 2016 and 2018 stroke incidence
analysis5  <-  rename.variable(analysis5,  "OLB026S",  "Nothing cheer")
analysis5  <-  rename.variable(analysis5,  "OLB018A",  "Things go wrong")
analysis5  <-  rename.variable(analysis5,  "OLB001I ",  "Newspaper")
analysis5  <-  rename.variable(analysis5,  "OLB033A  ",  "Making plans")
analysis5  <-  rename.variable(analysis5,  "OLB033B ",  "Activities trivial")
analysis5  <-  rename.variable(analysis5,  "OLB033C",  "Active planning")
analysis5  <-  rename.variable(analysis5,  "OLB033D",  "Accomplish")
analysis5  <-  rename.variable(analysis5,  "OLB033F",  "One day a time")
analysis5  <-  rename.variable(analysis5,  "OLB031L ",  "Nervous")
analysis5  <-  rename.variable(analysis5,  "OLB026R",  "Nervous")
analysis5  <-  rename.variable(analysis5,  "OLB018H",  "Hopeless")
analysis5  <-  rename.variable(analysis5,  "OLB018B ",  "Optimistic about future")
analysis5  <-  rename.variable(analysis5,  "OLB026J",  "Frustrated")
analysis5  <-  rename.variable(analysis5,  "OLB018C",  "Expect best")
analysis5  <-  rename.variable(analysis5,  "OLB026K",  "Happy")
analysis5  <-  rename.variable(analysis5,  "OLB018D",  "Expect good things")
analysis5  <-  rename.variable(analysis5,  "OLB026X",  "Calm")
analysis5  <-  rename.variable(analysis5,  "OLB018E",  "Things won't go well")
analysis5  <-  rename.variable(analysis5,  "OLB002C",  "Satisfied")
analysis5  <-  rename.variable(analysis5,  "OLB033E",  "Done all to do in life")
analysis5  <-  rename.variable(analysis5,  "OLB018F",  "Can't count on good")
analysis5  <-  rename.variable(analysis5,  "OLB001N",  "Use Internet")
analysis5  <-  rename.variable(analysis5,  "OLB012B",  "Use cell phone")

analysis5  <-  rename.variable(analysis5,  "OLB001I",  "Newspaper")
analysis5  <-  rename.variable(analysis5,  "OLB033A",  "Making plans")
analysis5  <-  rename.variable(analysis5,  "OLB033B",  "Activities trivial")
analysis5  <-  rename.variable(analysis5,  "OLB033G",  "Sense of direction")
analysis5  <-  rename.variable(analysis5,  "OLB001R",  "Has hobby")
analysis5  <-  rename.variable(analysis5,  "OLB031L",  "Nervous")
analysis5  <-  rename.variable(analysis5,  "OLB018B",  "Optimistic about future")


#### 2016 baseline covariates, 2018 stroke incidence
analysis6  <-  rename.variable(analysis6,  "PLB018A",  "Things go wrong")
analysis6  <-  rename.variable(analysis6,  "PLB001I ",  "Newspaper")
analysis6  <-  rename.variable(analysis6,  "PLB033A  ",  "Making plans")
analysis6  <-  rename.variable(analysis6,  "PLB033C",  "Active planning")
analysis6  <-  rename.variable(analysis6,  "PLB033D",  "Accomplish")
analysis6  <-  rename.variable(analysis6,  "PLB033F",  "One day a time")
analysis6  <-  rename.variable(analysis6,  "PLB031L ",  "Nervous")
analysis6  <-  rename.variable(analysis6,  "PLB026R",  "Nervous")
analysis6  <-  rename.variable(analysis6,  "PLB018H",  "Hopeless")
analysis6  <-  rename.variable(analysis6,  "PLB018B ",  "Optimistic about future")
analysis6  <-  rename.variable(analysis6,  "PLB026J",  "Frustrated")
analysis6  <-  rename.variable(analysis6,  "PLB018C",  "Expect best")
analysis6  <-  rename.variable(analysis6,  "PLB026K",  "Happy")
analysis6  <-  rename.variable(analysis6,  "PLB018D",  "Expect good things")
analysis6  <-  rename.variable(analysis6,  "PLB026X",  "Calm")
analysis6  <-  rename.variable(analysis6,  "PLB018E",  "Things won't go well")
analysis6  <-  rename.variable(analysis6,  "PLB002C",  "Satisfied")
analysis6  <-  rename.variable(analysis6,  "PLB033E",  "Done all to do in life")
analysis6  <-  rename.variable(analysis6,  "PLB018F",  "Can't count on good")
analysis6  <-  rename.variable(analysis6,  "PLB001N",  "Use Internet")
analysis6  <-  rename.variable(analysis6,  "PLB012B",  "Use cell phone")

analysis6  <-  rename.variable(analysis6,  "PLB001I",  "Newspaper")
analysis6  <-  rename.variable(analysis6,  "PLB033A",  "Making plans")
analysis6  <-  rename.variable(analysis6,  "PLB033B",  "Activities trivial")
analysis6  <-  rename.variable(analysis6,  "PLB033G",  "Sense of direction")
analysis6  <-  rename.variable(analysis6,  "PLB001R",  "Has hobby")
analysis6  <-  rename.variable(analysis6,  "PLB031L",  "Nervous")
analysis6  <-  rename.variable(analysis6,  "PLB018B",  "Optimistic about future")

analysis6  <-  rename.variable(analysis6,  "PLB042A",  "Frustration")


##### MARKOV CHAIN IMPUTED DATASETS
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
colnames(markov.analysis1)[11] <- "Fear Worst"

markov.analysis1 <- rename.variable(markov.analysis1, "KLB019A", "Dislike helping")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB027I", "Nothing cheer")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB019F", "Things go wrong")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB027A", "Cheerful")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB001a", "Newspaper")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB035A", "Making plans")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB035B", "Activities trivial")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB035C", "Active planning")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB035D", "Accomplish")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB035F", "One day a time")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB035G", "Sense of direction")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB041B", "Nervous")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB027N", "Nervous")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB019B", "Unfair profit")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB027J", "Hopeless")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB019G", "Optimistic about future")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB027B", "Good spririts")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB001b", "Has hobby")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB027J", "Frustrated")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB041C", "Trembling hands")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB019C", "People care")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB027K", "Restless")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB019H", "Expect best")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB027C", "Happy")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB001c", "Vacation")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB041D", "Fear dying")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB019D", "Lie to get ahead")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB027L", "Everything an effort")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB019I", "Expect good things")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB027D", "Calm")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB041E", "Faint")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB019E", "Sceptical people nice")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB027M", "Worthless")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB019J", "Things won't go well")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB027E", "Satisfied")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB001e", "Daytrip")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB027M", "Worthless")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB035E", "Done all to do in life")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB019K", "Can't count on good")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB027F", "Full of life")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB001f", "Use Internet")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB001g", "Use cell phone")


markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB001A", "Newspaper")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB001B", "Has hobby")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB001C", "Vacation")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB001F", "Use Internet")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB001G", "Use cell phone")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB027J.1", "Hopeless")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB027M.1", "Worthless")
markov.analysis1  <-  rename.variable(markov.analysis1,  "KLB001h", "Hopeless")


#### 2008 baseline covariates, 2010 and 2012 stroke incidence
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB041A", "Fear worst")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB019A", "Dislike helping")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB027I", "Nothing cheer")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB019F", "Things go wrong")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB027A", "Cheerful")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB001a", "Newspaper")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB035A", "Making plans")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB035B", "Activities trivial")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB035C", "Active planning")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB035D", "Accomplish")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB035F", "One day a time")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB035G", "Sense of direction")

markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB041B", "Nervous")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB027N", "Nervous")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB019B", "Unfair profit")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB027J", "Hopeless")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB019G", "Optimistic about future")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB027B", "Good spririts")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB001b", "Has hobby")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB027J", "Frustrated")

markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB041C", "Trembling hands")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB019C", "People care")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB027K", "Restless")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB019H", "Expect best")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB027C", "Happy")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB001c", "Vacation")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB041D", "Fear dying")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB019D", "Lie to get ahead")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB027L", "Everything an effort")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB019I", "Expect good things")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB027D", "Calm")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB041E", "Faint")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB019E", "Sceptical people nice")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB027M", "Worthless")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB019J", "Things won't go well")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB027E", "Satisfied")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB001e", "Daytrip")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB027M", "Worthless")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB035E", "Done all to do in life")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB019K", "Can't count on good")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB027F", "Full of life")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB001f", "Use Internet")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB001g", "Use cell phone")

markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB027J.1", "Hopeless")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB027M.1", "Worthless")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB001A", "Newspaper")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB001B", "Has hobby")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB001C", "Vacation")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB001F", "Use Internet")
markov.analysis2  <-  rename.variable(markov.analysis2,  "LLB001G", "Use cell phone")




#### 2010 baseline covariates, 2012 and 2014 stroke incidence
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB041A", "Fear worst")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB019A", "Dislike helping")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB027I", "Nothing cheer")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB019F", "Things go wrong")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB027A", "Cheerful")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB001a", "Newspaper")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB035A", "Making plans")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB035B", "Activities trivial")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB035C", "Active planning")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB035D", "Accomplish")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB035F", "One day a time")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB035G", "Sense of direction")

markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB041B", "Nervous")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB027N", "Nervous")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB019B", "Unfair profit")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB027J", "Hopeless")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB019G", "Optimistic about future")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB027B", "Good spririts")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB001b", "Has hobby")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB027J", "Frustrated")

markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB041C", "Trembling hands")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB019C", "People care")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB027K", "Restless")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB019H", "Expect best")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB027C", "Happy")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB001c", "Vacation")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB041D", "Fear dying")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB019D", "Lie to get ahead")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB027L", "Everything an effort")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB019I", "Expect good things")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB027D", "Calm")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB041E", "Faint")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB019E", "Sceptical people nice")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB027M", "Worthless")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB019J", "Things won't go well")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB027E", "Satisfied")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB001e", "Daytrip")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB027M", "Worthless")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB035E", "Done all to do in life")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB019K", "Can't count on good")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB027F", "Full of life")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB001f", "Use Internet")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB001g", "Use cell phone")


markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB001A", "Newspaper")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB001B", "Has hobby")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB027J.1", "Hopeless")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB001E", "Daytrip")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB027M.1", "Worthless")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB001F", "Use Internet")
markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB001G", "Use cell phone")


markov.analysis3  <-  rename.variable(markov.analysis3,  "MLB001C", "Vacation")


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB041A", "Fear worst")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB019A", "Dislike helping")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB027I", "Nothing cheer")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB019F", "Things go wrong")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB027A", "Cheerful")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB001a", "Newspaper")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB035A", "Making plans")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB035B", "Activities trivial")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB035C", "Active planning")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB035D", "Accomplish")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB035F", "One day a time")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB035G", "Sense of direction")

markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB041B", "Nervous")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB027N", "Nervous")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB019B", "Unfair profit")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB027J", "Hopeless")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB019G", "Optimistic about future")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB027B", "Good spririts")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB001b", "Has hobby")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB027J", "Frustrated")

markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB041C", "Trembling hands")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB019C", "People care")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB027K", "Restless")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB019H", "Expect best")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB027C", "Happy")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB001c", "Vacation")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB041D", "Fear dying")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB019D", "Lie to get ahead")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB027L", "Everything an effort")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB019I", "Expect good things")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB027D", "Calm")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB041E", "Faint")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB019E", "Sceptical people nice")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB027M", "Worthless")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB019J", "Things won't go well")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB027E", "Satisfied")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB001e", "Daytrip")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB027M", "Worthless")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB035E", "Done all to do in life")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB019K", "Can't count on good")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB027F", "Full of life")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB001f", "Use Internet")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB001g", "Use cell phone")


markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB001A", "Newspaper")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB001B", "Has hobby")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB027J.1", "Hopeless")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB001E", "Daytrip")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB027M.1", "Worthless")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB001F", "Use Internet")
markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB001G", "Use cell phone")

markov.analysis4  <-  rename.variable(markov.analysis4,  "NLB001C", "Vacation")



#### 2014 baseline covariates, 2016 and 2018 stroke incidence
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB026S",  "Nothing cheer")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB018A",  "Things go wrong")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB001I ",  "Newspaper")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB033A  ",  "Making plans")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB033B ",  "Activities trivial")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB033C",  "Active planning")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB033D",  "Accomplish")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB033F",  "One day a time")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB031L ",  "Nervous")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB026R",  "Nervous")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB018H",  "Hopeless")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB018B ",  "Optimistic about future")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB026J",  "Frustrated")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB018C",  "Expect best")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB026K",  "Happy")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB018D",  "Expect good things")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB026X",  "Calm")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB018E",  "Things won't go well")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB002C",  "Satisfied")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB033E",  "Done all to do in life")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB018F",  "Can't count on good")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB001N",  "Use Internet")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB012B",  "Use cell phone")

markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB001I",  "Newspaper")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB033A",  "Making plans")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB033B",  "Activities trivial")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB033G",  "Sense of direction")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB001R",  "Has hobby")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB031L",  "Nervous")
markov.analysis5  <-  rename.variable(markov.analysis5,  "OLB018B",  "Optimistic about future")



#### 2016 baseline covariates, 2018 stroke incidence
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB018A",  "Things go wrong")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB001I ",  "Newspaper")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB033A  ",  "Making plans")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB033C",  "Active planning")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB033D",  "Accomplish")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB033F",  "One day a time")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB031L ",  "Nervous")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB026R",  "Nervous")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB018H",  "Hopeless")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB018B ",  "Optimistic about future")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB026J",  "Frustrated")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB018C",  "Expect best")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB026K",  "Happy")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB018D",  "Expect good things")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB026X",  "Calm")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB018E",  "Things won't go well")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB002C",  "Satisfied")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB033E",  "Done all to do in life")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB018F",  "Can't count on good")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB001N",  "Use Internet")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB012B",  "Use cell phone")

markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB001I",  "Newspaper")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB033A",  "Making plans")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB033B",  "Activities trivial")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB033G",  "Sense of direction")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB001R",  "Has hobby")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB031L",  "Nervous")
markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB018B",  "Optimistic about future")

markov.analysis6  <-  rename.variable(markov.analysis6,  "PLB042A",  "Frustration")



##### BOOTSTRAP IMPUTED DATASETS
#### 2006 baseline covariates, 2008 and 2010 stroke incidence
colnames(bootstrap.analysis1)[11] <- "Fear Worst"

bootstrap.analysis1 <- rename.variable(bootstrap.analysis1, "KLB019A", "Dislike helping")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB027I", "Nothing cheer")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB019F", "Things go wrong")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB027A", "Cheerful")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB001a", "Newspaper")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB035A", "Making plans")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB035B", "Activities trivial")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB035C", "Active planning")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB035D", "Accomplish")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB035F", "One day a time")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB035G", "Sense of direction")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB041B", "Nervous")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB027N", "Nervous")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB019B", "Unfair profit")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB027J", "Hopeless")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB019G", "Optimistic about future")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB027B", "Good spririts")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB001b", "Has hobby")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB027J", "Frustrated")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB041C", "Trembling hands")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB019C", "People care")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB027K", "Restless")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB019H", "Expect best")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB027C", "Happy")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB001c", "Vacation")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB041D", "Fear dying")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB019D", "Lie to get ahead")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB027L", "Everything an effort")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB019I", "Expect good things")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB027D", "Calm")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB041E", "Faint")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB019E", "Sceptical people nice")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB027M", "Worthless")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB019J", "Things won't go well")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB027E", "Satisfied")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB001e", "Daytrip")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB027M", "Worthless")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB035E", "Done all to do in life")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB019K", "Can't count on good")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB027F", "Full of life")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB001f", "Use Internet")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB001g", "Use cell phone")


bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB001A", "Newspaper")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB001B", "Has hobby")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB001C", "Vacation")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB001F", "Use Internet")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB001G", "Use cell phone")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB027J.1", "Hopeless")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB027M.1", "Worthless")
bootstrap.analysis1  <-  rename.variable(bootstrap.analysis1,  "KLB001h", "Hopeless")


#### 2008 baseline covariates, 2010 and 2012 stroke incidence
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB041A", "Fear worst")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB019A", "Dislike helping")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB027I", "Nothing cheer")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB019F", "Things go wrong")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB027A", "Cheerful")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB001a", "Newspaper")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB035A", "Making plans")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB035B", "Activities trivial")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB035C", "Active planning")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB035D", "Accomplish")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB035F", "One day a time")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB035G", "Sense of direction")

bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB041B", "Nervous")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB027N", "Nervous")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB019B", "Unfair profit")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB027J", "Hopeless")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB019G", "Optimistic about future")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB027B", "Good spririts")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB001b", "Has hobby")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB027J", "Frustrated")

bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB041C", "Trembling hands")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB019C", "People care")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB027K", "Restless")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB019H", "Expect best")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB027C", "Happy")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB001c", "Vacation")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB041D", "Fear dying")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB019D", "Lie to get ahead")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB027L", "Everything an effort")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB019I", "Expect good things")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB027D", "Calm")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB041E", "Faint")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB019E", "Sceptical people nice")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB027M", "Worthless")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB019J", "Things won't go well")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB027E", "Satisfied")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB001e", "Daytrip")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB027M", "Worthless")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB035E", "Done all to do in life")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB019K", "Can't count on good")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB027F", "Full of life")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB001f", "Use Internet")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB001g", "Use cell phone")

bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB027J.1", "Hopeless")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB027M.1", "Worthless")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB001A", "Newspaper")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB001B", "Has hobby")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB001C", "Vacation")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB001F", "Use Internet")
bootstrap.analysis2  <-  rename.variable(bootstrap.analysis2,  "LLB001G", "Use cell phone")




#### 2010 baseline covariates, 2012 and 2014 stroke incidence
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB041A", "Fear worst")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB019A", "Dislike helping")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB027I", "Nothing cheer")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB019F", "Things go wrong")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB027A", "Cheerful")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB001a", "Newspaper")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB035A", "Making plans")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB035B", "Activities trivial")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB035C", "Active planning")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB035D", "Accomplish")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB035F", "One day a time")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB035G", "Sense of direction")

bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB041B", "Nervous")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB027N", "Nervous")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB019B", "Unfair profit")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB027J", "Hopeless")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB019G", "Optimistic about future")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB027B", "Good spririts")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB001b", "Has hobby")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB027J", "Frustrated")

bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB041C", "Trembling hands")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB019C", "People care")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB027K", "Restless")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB019H", "Expect best")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB027C", "Happy")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB001c", "Vacation")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB041D", "Fear dying")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB019D", "Lie to get ahead")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB027L", "Everything an effort")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB019I", "Expect good things")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB027D", "Calm")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB041E", "Faint")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB019E", "Sceptical people nice")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB027M", "Worthless")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB019J", "Things won't go well")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB027E", "Satisfied")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB001e", "Daytrip")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB027M", "Worthless")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB035E", "Done all to do in life")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB019K", "Can't count on good")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB027F", "Full of life")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB001f", "Use Internet")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB001g", "Use cell phone")


bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB001A", "Newspaper")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB001B", "Has hobby")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB027J.1", "Hopeless")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB001E", "Daytrip")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB027M.1", "Worthless")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB001F", "Use Internet")
bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB001G", "Use cell phone")


bootstrap.analysis3  <-  rename.variable(bootstrap.analysis3,  "MLB001C", "Vacation")


#### 2012 baseline covariates, 2014 and 2016 stroke incidence
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB041A", "Fear worst")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB019A", "Dislike helping")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB027I", "Nothing cheer")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB019F", "Things go wrong")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB027A", "Cheerful")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB001a", "Newspaper")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB035A", "Making plans")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB035B", "Activities trivial")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB035C", "Active planning")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB035D", "Accomplish")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB035F", "One day a time")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB035G", "Sense of direction")

bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB041B", "Nervous")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB027N", "Nervous")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB019B", "Unfair profit")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB027J", "Hopeless")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB019G", "Optimistic about future")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB027B", "Good spririts")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB001b", "Has hobby")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB027J", "Frustrated")

bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB041C", "Trembling hands")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB019C", "People care")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB027K", "Restless")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB019H", "Expect best")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB027C", "Happy")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB001c", "Vacation")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB041D", "Fear dying")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB019D", "Lie to get ahead")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB027L", "Everything an effort")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB019I", "Expect good things")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB027D", "Calm")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB041E", "Faint")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB019E", "Sceptical people nice")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB027M", "Worthless")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB019J", "Things won't go well")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB027E", "Satisfied")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB001e", "Daytrip")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB027M", "Worthless")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB035E", "Done all to do in life")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB019K", "Can't count on good")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB027F", "Full of life")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB001f", "Use Internet")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB001g", "Use cell phone")


bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB001A", "Newspaper")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB001B", "Has hobby")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB027J.1", "Hopeless")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB001E", "Daytrip")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB027M.1", "Worthless")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB001F", "Use Internet")
bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB001G", "Use cell phone")

bootstrap.analysis4  <-  rename.variable(bootstrap.analysis4,  "NLB001C", "Vacation")



#### 2014 baseline covariates, 2016 and 2018 stroke incidence
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB026S",  "Nothing cheer")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB018A",  "Things go wrong")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB001I ",  "Newspaper")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB033A  ",  "Making plans")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB033B ",  "Activities trivial")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB033C",  "Active planning")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB033D",  "Accomplish")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB033F",  "One day a time")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB031L ",  "Nervous")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB026R",  "Nervous")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB018H",  "Hopeless")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB018B ",  "Optimistic about future")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB026J",  "Frustrated")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB018C",  "Expect best")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB026K",  "Happy")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB018D",  "Expect good things")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB026X",  "Calm")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB018E",  "Things won't go well")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB002C",  "Satisfied")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB033E",  "Done all to do in life")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB018F",  "Can't count on good")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB001N",  "Use Internet")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB012B",  "Use cell phone")

bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB001I",  "Newspaper")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB033A",  "Making plans")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB033B",  "Activities trivial")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB033G",  "Sense of direction")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB001R",  "Has hobby")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB031L",  "Nervous")
bootstrap.analysis5  <-  rename.variable(bootstrap.analysis5,  "OLB018B",  "Optimistic about future")



#### 2016 baseline covariates, 2018 stroke incidence
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB018A",  "Things go wrong")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB001I ",  "Newspaper")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB033A  ",  "Making plans")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB033C",  "Active planning")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB033D",  "Accomplish")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB033F",  "One day a time")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB031L ",  "Nervous")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB026R",  "Nervous")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB018H",  "Hopeless")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB018B ",  "Optimistic about future")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB026J",  "Frustrated")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB018C",  "Expect best")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB026K",  "Happy")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB018D",  "Expect good things")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB026X",  "Calm")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB018E",  "Things won't go well")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB002C",  "Satisfied")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB033E",  "Done all to do in life")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB018F",  "Can't count on good")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB001N",  "Use Internet")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB012B",  "Use cell phone")

bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB001I",  "Newspaper")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB033A",  "Making plans")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB033B",  "Activities trivial")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB033G",  "Sense of direction")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB001R",  "Has hobby")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB031L",  "Nervous")
bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB018B",  "Optimistic about future")

bootstrap.analysis6  <-  rename.variable(bootstrap.analysis6,  "PLB042A",  "Frustration")




#### RECODE INFINITE VALUES
## Function to recode infinite values
Inf.recode <- function(data)
{
  for (i in 1:ncol(data)){
    data[,i][is.infinite(data[,i])] = 0
  }
  return(data)
}

## Markov chain imputed datasets
markov.analysis1 <- Inf.recode(markov.analysis1)

markov.analysis2 <- Inf.recode(markov.analysis2)

markov.analysis3 <- Inf.recode(markov.analysis3)

markov.analysis4 <- Inf.recode(markov.analysis4)

markov.analysis5 <- Inf.recode(markov.analysis5)

markov.analysis6 <- Inf.recode(markov.analysis6)



## Bootstrap imputed datasets
bootstrap.analysis1 <- Inf.recode(bootstrap.analysis1)

bootstrap.analysis2 <- Inf.recode(bootstrap.analysis2)

bootstrap.analysis3 <- Inf.recode(bootstrap.analysis3)

bootstrap.analysis4 <- Inf.recode(bootstrap.analysis4)

bootstrap.analysis5 <- Inf.recode(bootstrap.analysis5)

bootstrap.analysis6 <- Inf.recode(bootstrap.analysis6)
