####################################################################################
###
### calculating age_lookup for AAH population model, 
### which is the ageclass given the number of weeks age
###
####################################################################################

age_lookup_f <- c(rep(1:4, each = intvl_step_yr_weekly),
                       rep(5, 2 * intvl_step_yr_weekly),
                       rep(6, 3 * intvl_step_yr_weekly))
age_lookup_f <- c(age_lookup_f,
                  rep(7, nT_age_surv - length(age_lookup_f)))

age_lookup_m <- c(rep(1:4, each = intvl_step_yr_weekly),
                       rep(5, 2 * intvl_step_yr_weekly),
                       rep(6, 3 * intvl_step_yr_weekly))
age_lookup_m <- c(age_lookup_m,
                  rep(6, nT_age_surv - length(age_lookup_m)))

######################################################
### age to date conversion for FOI age/period effects
######################################################
study_start_ext <- "1985-05-15"
death_end <- "2022-05-14"

cwd_df$birthweek <- (interval(study_start_ext,
                            cwd_df$birth_date) %/% weeks(1))

###############################################################
###############################################################
### For glenn and I to discuss.... 
### I had birthweek specified with "+ 1", because in the
#### but when I incorporate the age2date vector for the surveillance data
#### i've used 
#hunt_pos_age2date = d_fit_hunt_pos$birthweek - 1,
#hunt_neg_age2date = d_fit_hunt_neg$birthweek - 1,
### because in the likelihood, we index for loops and add those indexes, 
### starting with age 1, so we default add one back in in the age2date conversion
### but the problem that I have here, is that the 
### parameter indicating the number of intervals prior to the start of the study
nT_period_prestudy_ext 
### this is 470, which is when the study starts, but when I calculate this using 
### the interval function, I get 469
interval("1985-05-15","1994-05-14") %/% weeks(1) #469
### which indicates that there are 469 intervals prior to study starting, 
### the minimum birth date, "1992-05-15" 
min(cwd_df$birth_date)
### this aligns with the interval 365
interval("1985-05-15","1992-05-15") %/% weeks(1) #365
### and if I don't add 1, then this is the correct birth week, but given 
min(cwd_df$birthweek) #365
### this appears to be okay...
###############################################################
###############################################################

cwd_df$weekkill <- interval(study_start_ext,
                            cwd_df$kill_date) %/% weeks(1)
cwd_df$yearkill <- cwd_df$kill_year - year(study_start_ext) + 1

# interval("1994-05-15","1995-01-01") %/% weeks(1)

period_lookup_foi <- c(rep(1, interval("1985-05-15","1994-05-15") %/% weeks(1)),
                       rep(1, interval("1994-05-15","1995-01-01") %/% weeks(1)),
                       rep(2:n_year, each = intvl_step_yr_weekly))

period_lookup_foi <- c(period_lookup_foi, rep(n_year,nT_period_overall_ext-
                                          length(period_lookup_foi)))

(n_period_lookup <- length(period_lookup_foi))

period_lookup_foi_study <- period_lookup_foi[(nT_period_prestudy_ext+1):nT_period_overall_ext]
num_foi_cal <- table(period_lookup_foi_study)

#############################################################################################
### Creating adjacency matrix and hyper parameter
### values for the dcar_normal implementation
### for period effects for FOI 
#############################################################################################

#create num vector
num_period <- c(1, rep(2, n_year - 2), 1)

#create adjacency vector along both years
temp <- as.matrix(bandSparse(n = n_year, k = c(1), symmetric = T))
temp2 <- matrix(0, n_year, n_year)
for (i in 1:nrow(temp2)) {
  temp2[i, ] <- temp[i, ] * c(1:n_year)
}
adj_period <- t(temp2)[which(t(temp2) != 0)]
n_adj_period <- length(adj_period)
weights_period <- rep(1, length(adj_period))

#########################################################
###
### Age and period indexing for FOI for collared deer
###
#########################################################

# age_week_indx <- c(rep(1,intvl_step_yr_weekly),#fawns
#                    rep(2,intvl_step_yr_weekly),#1
#                    rep(3,intvl_step_yr_weekly),#2
#                    rep(4,intvl_step_yr_weekly),#3
#                    rep(5, intvl_step_yr_weekly * 2),#4-5
#                    rep(6, intvl_step_yr_weekly * 3),#6-8
#                    rep(7,nT_age_surv - 
#                          length(c(rep(1, intvl_step_yr_weekly),#fawns
#                                   rep(2, intvl_step_yr_weekly),#1
#                                   rep(3, intvl_step_yr_weekly),#2
#                                   rep(4, intvl_step_yr_weekly),#3
#                                   rep(5, intvl_step_yr_weekly * 2),#4-5
#                                   rep(6, intvl_step_yr_weekly * 3)))))#6_

# period_week_indx <- c(rep(1,51),#2017
#                    rep(2,52),#2018
#                    rep(3,52),#2019
#                    rep(4,52),#2020
#                    rep(5,52),#2021
#                    rep(6,nT_period_collar - length(c(rep(1,51),#2017
#                                               rep(2,52),#2018
#                                               rep(3,52),#2019
#                                               rep(4,52),#2020
#                                               rep(5,52))))#2022
#                    )

# period_week_indx_col <- period_week_indx + n_year_precollar_ext - 1

###############################################################
###
### Setting up aggregation by study
### Separating hunter harvested
### CWD test positive deer from CWD test negative
###
### d_fit_hunt_pos
### d_fit_hunt_neg
###
##############################################################

cwd_df_agg <- cwd_df %>% 
              group_by(teststatus, ageweeks, birthweek, sex, ew) %>%
              summarise(n_cases = n(), .groups = 'drop')
d_fit_hunt_neg <- cwd_df_agg[cwd_df_agg$teststatus == 0, ]
d_fit_hunt_pos <- cwd_df_agg[cwd_df_agg$teststatus == 1, ]
