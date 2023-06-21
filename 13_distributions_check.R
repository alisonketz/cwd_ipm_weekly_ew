####################################################################
###
### Data and paramters needed to calculate likelihoods
### Fake parameters or loaded from runs w/o integrated models
###
####################################################################

beta0_survival_sus <- -8.75
beta0_survival_inf <- -8
beta_male <- .5

age_effect_survival_test <- 2*exp(-.09*seq(1:nT_age_surv))
age_effect_survival_test[600:nT_age_surv] <- .00001*exp(.009*seq(600:nT_age_surv))
age_effect_survival_test <- age_effect_survival_test- mean(age_effect_survival_test)

period_effect_survival_test <- 1 * sin(2/52 * pi * (1:nT_period_overall_ext)) + rnorm(nT_period_overall_ext,0,.1)
period_effect_survival_test <- period_effect_survival_test - mean(period_effect_survival_test)

##########################################################################
### FOI Parameters
### loading age and period effects from Transmission v3 w/o Survival
###########################################################################

load("~/Documents/Transmission/Transmission_v5/01_urw2_ew_timeRW1/fit_sum.Rdata")

m_age_foi <- fit_sum[grep("m_age",rownames(fit_sum))[1:7],1]
f_age_foi <- fit_sum[grep("f_age",rownames(fit_sum))[1:7],1]
m_period_foi <- fit_sum[grep("m_time",rownames(fit_sum)),1]
f_period_foi <- fit_sum[grep("f_time",rownames(fit_sum)),1]
x <- 2002:2022
lmfperiod <- lm(f_period_foi~x)
lmmperiod <- lm(m_period_foi~x)
pred_foiperiod_f <- predict(lmfperiod,newdata=data.frame(x=1994:2001))
pred_foiperiod_m <- predict(lmmperiod,newdata=data.frame(x=1994:2001))
f_period_foi <- c(pred_foiperiod_f,f_period_foi)
m_period_foi <- c(pred_foiperiod_m,m_period_foi)
f_period_foi <- f_period_foi[-length(f_period_foi)]
m_period_foi <- m_period_foi[-length(m_period_foi)]


