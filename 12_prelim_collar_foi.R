


######################################################################
###
### Preliminary constants for running in the model
###
######################################################################

######################################################################
###
### setting up nested indexing for the harvest surveillance data, 
### for the individuals that are born prior to the start of the study
###
######################################################################


# lookup_pe_surv <- c(rep(1, intvl_step_yr_weekly * 2), 1:nT_period_overall)
# length(lookup_pe_surv)
# nT_period_overall_ext



##########################################
###
### age class indexing and period indexing
###
##########################################
n_age_lookup <- length(age_lookup_f)
age_lookup_col <- c(age_lookup_f,rep(n_ageclassf,max(d_surv$right_age_s,na.rm=TRUE) - n_age_lookup))
age_lookup_col_f <- age_lookup_col
age_lookup_col_m <- age_lookup_col
age_lookup_col_m[age_lookup_col == 7] <- 6

n_age_lookup_col <- length(age_lookup_col)
# period_lookup_col <- c(rep(24:28, each=52), rep(31,as.numeric(ceiling(difftime(death_end,"2022-01-01",units="weeks")))))#January 2017 - May 2022
# n_period_lookup_col <- length(period_lookup_col)
# 
##########################################
###
### age class indexing and period indexing
### for infected at capture
###
##########################################

# age_lookup_col_inf_f  <- age_lookup_col_inf_m <- age_lookup_col_inf <- age_lookup_col
# age_lookup_col_inf_m[age_lookup_col_inf_m == 7] <- 6
# n_age_lookup_col_inf <- length(age_lookup_col_inf)

# period_lookup_col_inf <- c(rep(1,ceiling(difftime("2002-01-01","2001-05-15",units="weeks"))),
#                            rep(1:30, each = 52), rep(31,22))
# n_period_lookup_col_inf <- length(period_lookup_col_inf)
