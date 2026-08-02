#==============================================================================
#Project:    Treatment effects using R
#Author:     Shaun Scholes
#Date:       August 2026
#Purpose:    Estimate treatment effects in Stata
#Input:      cattaneo2.dta
#Output:     Estimates
#Location:   "C:/Causal_Inference/Treatment_effects/code/Treatment_effects_using R"
#==============================================================================

# Remove trailing whitespace from all R files in the project
for (f in list.files(path = "..", pattern = "\\.R$",
                     recursive=TRUE, full.names = TRUE)) {
  lines <- readLines(f, warn = FALSE)
  cleaned <- sub("\\s+$", "", lines)
  if (!identical(lines, cleaned)) writeLines(cleaned, f)
}
getwd()
.libPaths()
print(R.version.string)

drive<-"C:/Causal_Inference/Treatment_effects/"
folder<-"data"
setwd(paste0(drive, folder))

library(haven)
library(ipw)
library(tidyverse)
library(survey)
library(MatchIt)
library(tableone)
library(sandwich)
library(Matching)
library(glue)
df <-haven::read_dta('./cattaneo2.dta')

#regression for control.
lm.treat0 <- lm(bweight~1 + mage,data=df[df$mbsmoke==0,])
lm.treat0
df$yhat0 <- (3108.28) + (11.36*df$mage)
mean(df$yhat0)

#regression for treatment.
lm.treat1 <- lm(bweight~1 + mage,data=df[df$mbsmoke==1,])
lm.treat1
df$yhat1 <- (3237.091) + (-3.951*df$mage)
mean(df$yhat1)

df$diff<-df$yhat1  - df$yhat0

#ate
summary(df$diff)

#new tibble.
df %>% dplyr::summarize(
  ate = mean(diff)
)

#atet: (=ate for the observed Tx: mbsmoke=1)

df %>%
  dplyr::group_by(mbsmoke) %>%
  dplyr::summarize(atet = mean(diff))

##########################
#Regression adjustment
##########################

df <-haven::read_dta('./cattaneo2.dta')

#regression for control.
lm.treat0 <- lm(bweight~1 + fage + mage + mmarried + prenatal1
                ,data=df[df$mbsmoke==0,])
lm.treat0
df$yhat0 <- (3120.4836) + (-0.3222 * df$fage) +
  (4.9787*df$mage) + (162.2292*df$mmarried) + (55.6804*df$prenatal1)
mean(df$yhat0)

#regression for treatment.
lm.treat1 <- lm(bweight~1 + fage + mage + mmarried + prenatal1
                ,data=df[df$mbsmoke==1,])
lm.treat1
df$yhat1 <- (3267.10555) + (-0.03353 * df$fage) + (-8.51209*df$mage)
+ (131.72129*df$mmarried)
+ (33.70559*df$prenatal1)
mean(df$yhat1)

df$diff<-df$yhat1  - df$yhat0

#ate
summary(df$diff)

#new tibble.
df %>% dplyr::summarize(
  ate = mean(diff)
)

#atet: (=ate for the observed Tx: mbsmoke=1)

df %>%
  dplyr::group_by(mbsmoke) %>%
  dplyr::summarize(atet = mean(diff))

print(paste("Regression adjustment finished at", format(Sys.time(),
                                                        "%Y-%m-%d %H:%M:%S")))

##################
###IPW
##################

library(haven)
library(ipw)
library(tidyverse)
library(survey)
library(MatchIt)
library(tableone)
library(sandwich)
library(Matching)
library(glue)
library(WeightIt)
library(estimatr)

#############
#IPW.
#############
#Example 1.

df <-haven::read_dta('./cattaneo2.dta')

# Rename 'mbsmoke' to 'treat'
df <- df %>% rename(treat = mbsmoke)
# Propensity score model
ps_model <- glm(treat ~ fage + mage + mmarried + fbaby, data = df, family = binomial())
# Predict propensity scores
df$pscore <- predict(ps_model, type = "response")

# IP weights for POMs and ATE
df$ipw <- ifelse(df$treat == 1, 1 / df$pscore, 1 / (1 - df$pscore))
# Define survey design
design_ate <- svydesign(ids = ~1, weights = ~ipw, data = df)

# Potential means by treatment group
svyby(~bweight, ~treat, design_ate, svymean)

# Estimate ATE as the difference in means
ate <- svyglm(bweight ~ treat, design = design_ate)
summary(ate)

# IP weights for ATET
df$ipw2 <- ifelse(df$treat == 1, 1, df$pscore / (1 - df$pscore))
# Design for ATET
design_atet <- svydesign(ids = ~1, weights = ~ipw2, data = df)

# Estimate group means
svyby(~bweight, ~treat, design_atet, svymean)

# ATET from regression
atet <- svyglm(bweight ~ treat, design = design_atet)
summary(atet)

#alternative.
library(WeightIt)
library(estimatr)

# Weighting for ATE

#library("cobalt")
#data("lalonde", package = "cobalt")
#head(lalonde)

#library("WeightIt")
#W.out <- weightit(treat ~ age + educ + race + married + nodegree +
#                    re74 + re75,
#                  data = lalonde,
#                  estimand = "ATT",
#                  method = "glm")
#W.out #print the output
#summary(W.out)

#required treat as integer not double.
df$treat <- as.integer(df$treat)
class(df$treat)

# Drop columns by name
df <- df %>% dplyr::select(-pscore, -ipw, -ipw2)

w <- weightit(treat ~ 1, data = df,
              method = "glm", estimand = "ATE")
df$w <- w$weights

w_ate <- WeightIt::weightit(treat ~ fage + mage + mmarried + fbaby,
                            data = df,
                            method = "ps", estimand = "ATE")
df$w_ate <- w_ate$weights


lm_ate <- lm_robust(bweight ~ treat, weights = w_ate, data = df)
summary(lm_ate)


# Weighting for ATET
w_atet <- weightit(treat ~ fage + mage + mmarried + fbaby,
                   data = df, method = "ps", estimand = "ATT")
df$w_atet <- w_atet$weights

lm_atet <- lm_robust(bweight ~ treat, weights = w_atet, data = df)
summary(lm_atet)


############
#Example 2.
############

set.seed(6795)
# Set number of observations
n <- 1000
# Create dataset
id <- 1:n
treat <- runif(n) < 0.2
male <- runif(n) < 0.5
outcome <- sample(1:10, n, replace = TRUE)
x <- sample(1:10, n, replace = TRUE)
# Combine into a data frame
df <- data.frame(id, treat, male, outcome, x)
head(df)

#model for the Tx in the denominator.
tx <- ipwpoint(
  exposure = treat,
  family = "binomial",
  link = "logit",
  numerator = ~ 1,
  denominator = ~ x + male,
  data = df)
summary(tx$ipw.weights)
df$ipw <- tx$ipw.weights

#ATE:
msm <- (svyglm(outcome ~ treat,
               design = svydesign(~ 1,
                                  weights = ~ ipw,
                                  data = df)))
print(coef(msm))
print(confint(msm))


##example.
#Using WeightIt to Estimate Balancing Weights

library("cobalt")
data("lalonde", package = "cobalt")
head(lalonde)
bal.tab(treat ~ age + educ + race + married + nodegree +
          re74 + re75,
        data = lalonde,
        estimand = "ATT",
        thresholds = c(m = .05))

library("WeightIt")
W.out <- weightit(treat ~ age + educ + race + married + nodegree +
                    re74 + re75,
                  data = lalonde,
                  estimand = "ATT",
                  method = "glm")
W.out #print the output
summary(W.out)
bal.tab(W.out, stats = c("m", "v"),
        thresholds = c(m = .05))

W.out <- weightit(treat ~ age + educ + race + married + nodegree +
                    re74 + re75,
                  data = lalonde,
                  estimand = "ATT",
                  method = "ebal")
summary(W.out)
bal.tab(W.out, stats = c("m", "v"),
        thresholds = c(m = .05))
# Fit outcome model
fit <- lm_weightit(re78 ~ treat * (age + educ + race + married +
                                     nodegree + re74 + re75),
                   data = lalonde, weightit = W.out)

# G-computation for the treatment effect
library("marginaleffects")

avg_comparisons(fit, variables = "treat",
                newdata = subset(treat == 1))


print(paste("IPW estimator finished at", format(Sys.time(),
                                                "%Y-%m-%d %H:%M:%S")))

#================
#IPWRA.
#================


library(haven)
library(ipw)
library(tidyverse)
library(survey)
library(MatchIt)
library(tableone)
library(sandwich)
library(Matching)
library(twang)
library(glue)

##################
#Example 1
##################

data <-haven::read_dta('./cattaneo2.dta')
head(data)

# Rename variables
data <- data %>%
  rename(treat = mbsmoke)

# Propensity score model for the treatment variable
ps_model <- glm(treat ~ mmarried + mage + fbaby + medu,
                family = binomial(), data = data)

# Calculate propensity scores
data$pscore <- predict(ps_model, type = "response")

# Compute the inverse probability weights
data$ipw <- ifelse(data$treat == 0, 1 / (1 - data$pscore), 1 / data$pscore)

# Fit the regression for the treated group and the control group
model_treat <- lm(bweight ~ mage + prenatal1 + mmarried + fbaby,
                  data = data[data$treat == 1, ], weights = ipw)
model_control <- lm(bweight ~ mage + prenatal1 + mmarried + fbaby,
                    data = data[data$treat == 0, ], weights = ipw)
model_treat

data$pom_t<-(3201.664) + (-6.452*data$mage) + (26.590*data$prenatal1) +
  (136.652*data$mmarried) + (50.282*data$fbaby)
length(data$pom_t)

model_control


data$pom_c<-(3187.79) + (3.15*data$mage) + (66.36*data$prenatal1) +
  (156.70*data$mmarried) + (-71.50*data$fbaby)
length(data$pom_c)

# Summary statistics for treated and control group predictions
summary(data$pom_c)
summary(data$pom_t)

# Calculate the treatment effect (difference in predictions)
data$ate <- data$pom_t - data$pom_c
mean_ate <- mean(data$ate, na.rm = TRUE)
# Print the mean of the 'ate' variable
print(mean_ate)

##################
#Example 2
##################

data <-haven::read_dta('./cattaneo2.dta')
head(data)

# Rename variables
data <- data %>%
  rename(treat = mbsmoke)

# Install the twang package if not already installed
# install.packages("twang")

# Propensity score model for the treatment variable
ps_model <- glm(treat ~ fage + mage + fbaby + mmarried,
                family = binomial(), data = data)

# Calculate propensity scores
data$pscore <- predict(ps_model, type = "response")

# Compute the inverse probability weights
data$ipw <- ifelse(data$treat == 0, 1 / (1 - data$pscore), 1 / data$pscore)

# Fit the regression for the treated group and the control group
model_treat <- lm(bweight ~ fage + mage + mmarried + prenatal1,
                  data = data[data$treat == 1, ], weights = ipw)
model_control <- lm(bweight ~ fage + mage + mmarried + prenatal1 ,
                    data = data[data$treat == 0, ], weights = ipw)
model_treat

data$pom_t<-(3241.334) + (1.788*data$fage) + (-8.818*data$mage)
+ (122.586*data$mmarried) + (28.622*data$prenatal1)
length(data$pom_t)

model_control
data$pom_c<-(3114.5369) + (-0.5052*data$fage) + (5.4908*data$mage)
+ (161.6814*data$mmarried) + (54.3625*data$prenatal1)
length(data$pom_c)

# Summary statistics for treated and control group predictions
summary(data$pom_c)
summary(data$pom_t)

# Calculate the treatment effect (difference in predictions)
data$ate <- data$pom_t - data$pom_c
mean_ate <- mean(data$ate, na.rm = TRUE)
# Print the mean of the 'ate' variable
print(mean_ate)

print(paste("IPWRA finished at", format(Sys.time(),
                                        "%Y-%m-%d %H:%M:%S")))
#===========
#AIPW
#===========

# Load necessary libraries
library(haven)    # For reading Stata files
library(dplyr)    # For data manipulation
library(stats)    # For logistic regression and linear regression
library(glue)

data <-haven::read_dta('./cattaneo2.dta')
head(data)

data <- data %>%
  rename(treat = mbsmoke)  # Rename 'mbsmoke' to 'treat'

# Step 2: Estimate the propensity score using logistic regression
ps_model <- glm(treat ~ mmarried + mage + fbaby + medu,
                family = binomial(), data = data)

# Calculate propensity scores
data$pscore <- predict(ps_model, type = "response")

# Step 3: Compute the inverse probability weights (IPW) for treated
# and control groups
# IPW for control group
data$ipw0 <- ifelse(data$treat == 0, 1 / (1 - data$pscore), 0)
# IPW for treated group
data$ipw1 <- ifelse(data$treat == 1, 1 / data$pscore, 0)

# Step 4: Estimate the Potential Outcomes Model (POM) for treated group
# No weights in the outcome model
treated_data <- subset(data, treat == 1)

# Fit regression model for the treated group
model_treat <- lm(bweight ~ mage + prenatal1 + mmarried + fbaby,
                  data = treated_data)

# Generate potential outcome for treated group (POM1)
data$pom1 <- 3227.169 + (41.43991 * data$fbaby)
+ (133.6617 * data$mmarried) +
  (25.11133 * data$prenatal1) + (-7.370881 * data$mage)

# Adjust the prediction with the inverse probability weight (IPW)
data$pom1 <- data$pom1 + data$ipw1 * (data$bweight - data$pom1)

# Step 5: Estimate the Potential Outcomes Model (POM) for control group
control_data <- subset(data, treat == 0)

# Fit regression model for the control group
model_control <- lm(bweight ~ mage + prenatal1 + mmarried + fbaby,
                    data = control_data)

# Generate potential outcome for control group (POM0)
data$pom0 <- 3202.746 + (-71.3286 * data$fbaby)
+ (160.9513 * data$mmarried) +
  (64.40859 * data$prenatal1) + (2.546828 * data$mage)

# Adjust the prediction with the inverse probability weight (IPW)
data$pom0 <- data$pom0 + data$ipw0 * (data$bweight - data$pom0)

# Step 6: Calculate ATE using the potential outcomes
pom_c <- mean(data$pom0, na.rm = TRUE)  # Mean potential outcome for control
pom_t <- mean(data$pom1, na.rm = TRUE)  # Mean potential outcome for treated

# Calculate the ATE
ATE <- pom_t - pom_c

# Print the ATE
cat("AIPW =", ATE, "\n")

print(paste("AIPW finished at", format(Sys.time(),
                                       "%Y-%m-%d %H:%M:%S")))

#================
#matching
#================

library(PSW)
library(optmatch)
library(glue)

data <-haven::read_dta('./cattaneo2.dta')
head(data)

m.out <- matchit(mbsmoke ~ mmarried + mage + medu + fbaby,
                 data = data,
                 distance = "glm",
                 estimand = "ATT",
                 method = "nearest")
m.out
summary(m.out)   # 864 matched pairs
md<-match.data(m.out) # matched dataset for estimation
#View(md)
fit <- lm(bweight ~ mbsmoke, data = md, weights = weights)
fit

m.out <- matchit(mbsmoke ~ mage, data = data,
                 distance = "glm",
                 estimand = "ATT",
                 method = "nearest")
m.out
summary(m.out)
md<-match.data(m.out)
md
fit <- lm(bweight ~ mbsmoke, data = md, weights = weights)
fit

####################################################
#same result by user-supplied propensity scores
####################################################

pscore <- fitted(glm(mbsmoke ~ mage, data = data,family = binomial))
summary(pscore)
#match on the pscore
m.out <- matchit(mbsmoke ~ mage, data = data, distance = pscore)
summary(m.out)
md<-match.data(m.out)
md
fit1 <- lm(bweight ~ mbsmoke, data = md, weights = weights)
print(fit1)

print(paste("Matching estimators finished at", format(Sys.time(),
                                                      "%Y-%m-%d %H:%M:%S")))


#=============
#NNM matching
#=============

library(MatchIt)
library(glue)
library(haven)

data <-haven::read_dta('./cattaneo2.dta')
head(data)

#One-to-one matching
m.out <- MatchIt::matchit(mbsmoke ~ mmarried + mage,
                          distance = "mahalanobis",
                          method = "nearest",
                          estimand = "ATT",
                          ratio = 1, replace=TRUE,
                          data = data)
m.out
summary(m.out)
md<-match.data(m.out)
fit1 <- lm(bweight ~ mbsmoke, data = md,weights = weights)
print(fit1)

#Example 2.
data <-haven::read_dta('./cattaneo2.dta')
head(data)

#One-to-one matching
m.out <- MatchIt::matchit(mbsmoke ~ mmarried + mage + fage + medu + prenatal1,
                          exact ~ mmarried + prenatal1,
                          distance = "mahalanobis",
                          method = "nearest",
                          estimand = "ATT",
                          ratio = 1,
                          data = data)
m.out
summary(m.out)
md<-match.data(m.out)
fit1 <- lm(bweight ~ mbsmoke, data = md,weights = weights)
print(fit1)

print(paste("NNM finished at", format(Sys.time(),
                                         "%Y-%m-%d %H:%M:%S")))

#==========================
#Propensity score matching
#===========================

df <-haven::read_dta('./cattaneo2.dta')
head(df)


##################################################
#Example 1:
#matching variables:mmarried; mage; medu; fbaby
#distance= mahalanobis
#(One-to-one).
###################################################


#The match.data() output is preferred when pair membership is not
#directly included in the analysis;

#The get_matches() output is preferred when
#pair membership is to be included.

#NN Mahalanobis distance matching w/ replacement.
m.out <- matchit(mbsmoke ~ mmarried + mage + medu + fbaby,
                 data = df,
                 distance = "mahalanobis", ties=TRUE,
                 replace = TRUE)
m.out
summary(m.out)

#1151 obs matched.
#864 Tx matched with 287 control (meaning 864 pairs).
#864*2 = 1728.
#EACH TX=1 has been matched with at least 1 control.

#match data.

md <- match.data(m.out, data = df,
                 distance = "prop_score")
dim(md) #one row per matched unit [N=1151].
head(md, 10)

#get matches
gm <- get_matches(m.out, data = df,
                  distance = "prop_score")
dim(gm) #multiple rows per matched unit
head(gm, 10)
#864 pairs (1728 rows)

#Number of control units in each match stratum
table(table(gm$subclass[gm$A == 0]))
table(table(gm$subclass[gm$A == 1]))

#match.data() output
fit1md <- lm(bweight ~ mbsmoke, data = md, weights = weights)
fit1md


##################################################
#Example 2:
#matching variables:mmarried; mage; medu; fbaby
#distance= mahalanobis
#ratio=2 (2 controls for each Tx)
###################################################


# NN Mahalanobis distance matching w/ replacement.
m.out <- matchit(mbsmoke ~ mmarried + mage + medu + fbaby,
                 data = df,
                 distance = "mahalanobis", ties=TRUE,ratio=2,
                 replace = TRUE)
m.out
summary(m.out)

#1409 obs matched.
#864 Tx matched with 545 control (meaning 864 pairs).
#864*2 = 1728.
#EACH TX=1 has been matched with at least 1 control.

#match data.
md <- match.data(m.out, data = df,
                 distance = "prop_score")
dim(md) #one row per matched unit [N=1409].
head(md, 10)

#get matches

gm <- get_matches(m.out, data = df,
                  distance = "prop_score")
dim(gm) #multiple rows per matched unit
head(gm, 10)

#864 pairs (2592 rows)
#View(gm)
#(1 Tx + 1 control + 1 control) * 864.
#864*3 = 2592.

#Number of control units in each match stratum
table(table(gm$subclass[gm$mbsmoke == 0])) # 2 control for each Tx.

#match.data() output
fit1md <- lm(bweight ~ mbsmoke, data = md, weights = weights)
fit1md

##################################################
#Example 3:
#matching variables:mmarried; mage; medu; fbaby
#distance= NNM
#ratio=1
###################################################

#Nearest neighbour
m.out <- matchit(mbsmoke ~ mmarried + mage + medu + fbaby,
                 data = df,
                 method = "nearest", ties=TRUE,ratio=1,
                 replace = TRUE)
m.out
summary(m.out)

#1155 obs matched.
#864 Tx matched with 291 control (meaning 864 pairs).
#864*2 = 1728.
#EACH TX=1 has been matched with at least 1 control.

#match data.
md <- match.data(m.out, data = df,
                 distance = "prop_score")
dim(md) #one row per matched unit [N=1409].
head(md, 10)

#get matches
gm <- get_matches(m.out, data = df,
                  distance = "prop_score")
dim(gm) #multiple rows per matched unit
head(gm, 10)
#864 pairs (1728 rows)

#Number of control units in each match stratum
table(table(gm$subclass[gm$mbsmoke == 0])) # 1 control for each Tx.

#match.data() output
fit1md <- lm(bweight ~ mbsmoke, data = md, weights = weights)
fit1md


##################################################
#Example 4:
#matching variables:mmarried; mage; medu; fbaby
#distance= NNM
#variable number of matches
###################################################


m.out <- matchit(mbsmoke ~ mmarried + mage + medu + fbaby,
                 data = df,
                 ties=TRUE,
                 replace = TRUE,
                 method = "nearest",
                 ratio = 2,
                 min.controls = 1, max.controls = 80)
m.out
summary(m.out)

#1462 obs matched.
#864 Tx matched with 598 control (meaning 864 pairs).
#864*2 = 1728.
#EACH TX=1 has been matched with at least 1 control.

#match data.
md <- match.data(m.out, data = df,
                 distance = "prop_score")
dim(md) #one row per matched unit [N=1151].
head(md, 10)

#get matches
gm <- get_matches(m.out, data = df,
                  distance = "prop_score")
dim(gm) #multiple rows per matched unit
head(gm, 10)
#864 pairs (1728 rows)

#Number of control units in each match stratum
table(table(gm$subclass[gm$mbsmoke == 0]))

#match.data() output
fit1md <- lm(bweight ~ mbsmoke, data = md, weights = weights)
fit1md

print(paste("PSM matching finished at", format(Sys.time(),
                                                  "%Y-%m-%d %H:%M:%S")))

##########.
#FINISHED.
##########.



































