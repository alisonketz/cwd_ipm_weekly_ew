####################################################################################
###
### calculating age_lookup for AAH population model, 
### which is the ageclass given the number of weeks age
###
####################################################################################

age_lookup_f <- c(rep(1:4, each = intvl_step_yr_weekly),
                       rep(5, 2 * intvl_step_yr_weekly),
                       rep(6, 3 * intvl_step_yr_weekly),
                       rep(7, intvl_step_yr_weekly))
age_lookup_m <- c(rep(1:4, each = intvl_step_yr_weekly),
                       rep(5, 2 * intvl_step_yr_weekly),
                       rep(6, 4 * intvl_step_yr_weekly))
n_age_lookup_f <- length(age_lookup_f)
n_age_lookup_m <- length(age_lookup_m)

# This is what I used before
# age_lookup_f <- c(rep(1:4, each = intvl_step_yr_weekly),
#                        rep(5, 2 * intvl_step_yr_weekly),
#                        rep(6, 3 * intvl_step_yr_weekly),
#                        7)
# age_lookup_m <- age_lookup_f
# age_lookup_m[age_lookup_m == 7] <- 6
# n_age_lookup <- length(age_lookup_f)

###################################################################################
### age to date conversion within model
###################################################################################

# birth_start <- min(cwd_df$birth_date)
study_start_ext <- "1985-05-15"
death_end <- "2022-05-15"
cwd_df$birthweek <- (interval(study_start_ext,
                            cwd_df$birth_date) %/% weeks(1)) + 1
cwd_df$weekkill <- interval(study_start_ext,
                            cwd_df$kill_date) %/% weeks(1)
cwd_df$yearkill <- cwd_df$kill_year - year(study_start_ext) + 1

period_lookup_foi <- c(rep(1, nT_period_prestudy_ext),
                   rep(1:n_year, each = intvl_step_yr_weekly))

n_period_lookup <- length(period_lookup_foi)

#############################################################################################
###
### Creating adjacency matrix and hyper parameter values for the dcar_normal implementation
### For PERIOD 
###
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

age_week_indx <- c(rep(1,52),#fawns
                   rep(2,52),#1
                   rep(3,52),#2
                   rep(4,52),#3
                   rep(5, 52 * 2),#4-5
                   rep(6, 52 * 3),#6-8
                   rep(7,nT_age_surv-length(c(rep(1,52),#fawns
                                        rep(2,52),#1
                   rep(3,52),#2
                   rep(4,52),#3
                   rep(5, 52 * 2),#4-5
                   rep(6, 52 * 3)))))

period_week_indx <- c(rep(1,51),#2017
                   rep(2,52),#2018
                   rep(3,52),#2019
                   rep(4,52),#2020
                   rep(5,52),#2021
                   rep(6,nT_period_collar - length(c(rep(1,51),#2017
                                              rep(2,52),#2018
                                              rep(3,52),#2019
                                              rep(4,52),#2020
                                              rep(5,52))))#2022
                   )

period_week_indx_col <- period_week_indx + n_year_precollar_ext - 1

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
cwd_df_agg

d_fit_hunt_neg <- cwd_df_agg[cwd_df_agg$teststatus == 0, ]
d_fit_hunt_pos <- cwd_df_agg[cwd_df_agg$teststatus == 1, ]
