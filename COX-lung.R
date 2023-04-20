##### Univariate Cox regression on lung dataset

#install.packages(c("survival", "survminer"))
library("survival")
library("survminer")
View(lung)
head(lung)
?lung

# Calculate Cox PH - HR between male and female
res.cox <- coxph(Surv(time, status) ~ sex, data = lung)
res.cox
summary(res.cox)

#The Cox model gives the hazard ratio (HR) for the second group relative to the first group, that is, female versus male.
#The exponentiated coefficients (exp(coef) = exp(-0.53) = 0.59), also known as hazard ratios, give the effect size of covariates. 
#For example, being female (sex=2) reduces the hazard by a factor of 0.59, or 41%.

-------------
##### Multivariate Cox regression analysis
  
#A Cox regression of time to death on the time-constant covariates is specified as follow:

# Calculate Cox PH - HR between male and female and their ages
res.coxM <- coxph(Surv(time, status) ~ age + sex + ph.ecog, data =  lung)
summary(res.coxM)

##Visualizing the estimated distribution of survival times
ggsurvplot(survfit(Surv(time, status) ~ sex, data=lung))
ggsurvplot(survfit(Surv(time, status) ~ age, data = lung))
ggsurvplot(survfit(Surv(time, status) ~ ph.ecog, data=lung))


