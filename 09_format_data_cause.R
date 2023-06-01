###############################################################
###
### Formatting data to fit the conditional probability that 
### individuals that die during the hunting season die from 
### hunter harvest mortality
###
##############################################################

d_fit_hh <- d_surv[d_surv$censored == 0,]
age_class_indx <- c(intvl_step_yr_weekly,#fawns
                    intvl_step_yr_weekly * 2,#1
                    intvl_step_yr_weekly * 3,#2
                    intvl_step_yr_weekly * 4,#3
                    intvl_step_yr_weekly * 6,#4-5
                    intvl_step_yr_weekly * 9#6-8
                    )#9.5+

d_fit_hh$ageclassmort <- c()
for(i in 1:nrow(d_fit_hh)) {
    if(d_fit_hh$right_age_s[i] > intvl_step_yr_weekly * 9) {
        d_fit_hh$ageclassmort[i] <- 7
    }else{
        d_fit_hh$ageclassmort[i] <- 
            min(which(d_fit_hh$right_age_s[i] < age_class_indx))
    }
}

#females are 0
#males are 1
#creating the response for the cause-specific likelihood
#mort_h == 0 means that deer are not hunter harvest
#mort_h == 1 means that deer are hunter harvested 
d_fit_hh$mort_h <- 0
d_fit_hh$mort_h[d_fit_hh$p1 == 1] <- 1

records_cause <- nrow(d_fit_hh)

#change the sex for antlerless males (i.e. male fawns) to sex == 0
d_fit_hh$sex[d_fit_hh$sex==0 & d_fit_hh$ageclassmort == 1] <- 0

###########################################
###
### Read data - hunting season dates
###
###########################################

### huntseason is for duration of collar data only
d_huntseason <- read_xlsx(paste0(filepath,"Hunting_SeasonDates.xlsx"),1)
d_huntseason <- d_huntseason %>% filter(Year > 2016)

### season is specified from beginning of 
### the pop model study (May 15, 1994)
### for the full study
d_season <- read_xlsx(paste0(filepath,"Hunting_SeasonDates.xlsx"),1)
d_season <- d_season %>% filter(Year > 1993)


startdate <- min(df_cap$date_cap)
n_year <- length(unique(d_season$Year))
n_year_collar <- 5

startng <- c()
endng <- c()
startgun <- c()
endgun <- c()
for (i in 1:n_year_collar) {
    startng[i] <- interval(startdate,
        min(d_huntseason$OpenDate[d_huntseason$Year ==
                                  (i + 2016)])) %/% weeks(1) + 1
    endng[i] <- interval(startdate,
        max(d_huntseason$CloseDate[d_huntseason$Year ==
                                  (i + 2016)])) %/% weeks(1) + 1
    startgun[i] <- interval(startdate,
        d_huntseason$OpenDate[d_huntseason$Year ==
                              (i + 2016) &
        d_huntseason$SeasonType == "nineday"]) %/% weeks(1) + 1
    endgun[i] <- interval(startdate,
        d_huntseason$CloseDate[d_huntseason$Year ==
                                (i + 2016) &
        d_huntseason$SeasonType == "nineday"]) %/% weeks(1) + 1
}

study_start <- "1994-05-15"
season_ng_start <- c()
season_ng_end <- c()
for (i in 1:length(unique(d_season$Year))) {
    season_ng_start[i] <- 
        ceiling(as.duration(ymd("1994-05-15") %--%
                ymd(min(d_season$OpenDate[d_season$Year ==
                            (i + 1993)])))/dweeks(1))
    season_ng_end[i] <- 
        ceiling(as.duration(ymd("1994-05-15") %--%
            ymd(max(d_season$CloseDate[d_season$Year ==
            (i + 1993)])))/dweeks(1))
}

# pulling out approximate 9 day gun season
temp <- d_season %>% filter(SeasonType == "nineday")
nonineday <- c(1994:2021)[which(!(c(1994:2021) %in% temp$Year))]
temp2 <- d_season %>% 
            filter(Year %in% nonineday) %>% 
            filter(SeasonType != "Archery" &
                   SeasonType != "Muzzleloader" &
                   SeasonType != "Youth Gun")
temp2rm  <- c(5,6,8,11,12)
temp2 <- temp2[-temp2rm,]
temp <- rbind(temp,temp2)
temp <- temp[order(temp$Year),]

season_gun_start <- c()
season_gun_end <- c()
for (i in 1:length(unique(temp$Year))) {
    season_gun_start[i] <- 
        ceiling(as.duration(ymd("1994-05-15") %--%
            ymd(temp$OpenDate[temp$Year == (i + 1993)]))/dweeks(1))
    season_gun_end[i] <- 
        ceiling(as.duration(ymd("1994-05-15") %--%
        ymd(temp$CloseDate[temp$Year == (i + 1993)]))/dweeks(1))
}

### need 8 start/end points throughout year
### for the season indexing which corresponds to 
### start of birth pulse period
### then end non-harvest period
### then ng_harvest start
### then gun_harvest start
### gun harvest end
### ng harvest end
### postharvest start
### end of aah annual year, expressed in period effects

n_year_precollar <- 23
d_fit_season <- matrix(NA, nrow = n_year, ncol = 8)
for(t in 1:n_year) {
   d_fit_season[t, ] <- c(intvl_step_yr_weekly * (t - 1) + 1,
                          season_ng_start[t] - 1,
                          season_ng_start[t],
                          season_gun_start[t],
                          season_gun_end[t],
                          season_ng_end[t],
                          season_ng_end[t] + 1,
                          intvl_step_yr_weekly * t)
}
d_fit_season[n_year_precollar, 6] <- d_fit_season[n_year_precollar, 6] - 1

d_fit_season <- data.frame(d_fit_season)
colnames(d_fit_season) <- c("yr_start",
                            "pre_hunt_end",
                            "ng_start",
                            "gun_start",
                            "gun_end",
                            "ng_end",
                            "post_hunt_start",
                            "yr_end")

#saving for aah_disease test
save(d_fit_season, file = paste0(filepath, "d_fit_season.Rdata"))

d_fit_season$year <- 1994:2021

###
### Preliminaries gun season (not used)
###
# Z_cause_ng <- rep(0,nT_period_collar)
# Z_cause_gun <- rep(0,nT_period_collar)
# for(i in 1:5){
#     Z_cause_ng[startng[i]:endng[i]] <- 1
#     Z_cause_gun[startgun[i]:endgun[i]] <- 1
# }


# Z_overall_ng <- rep(0,nT_period_overall)
# Z_overall_gun <- rep(0,nT_period_overall)
# for(i in 1:n_year){
#     Z_overall_ng[d_fit_season$ng_start[i]:d_fit_season$ng_end[i]] <- 1
#     Z_overall_gun[d_fit_season$gun_start[i]:d_fit_season$gun_end[i]] <- 1
# }


################################################################
###
### preliminaries for hunter harvest Z at a monthly time step
###
################################################################


# startdate <- min(df_cap$date_cap)

# # interval(df_cap$date_cap[!is.na(df_cap$recap_cwd)],df_cap$recapdate_cap[!is.na(df_cap$recap_cwd)]) %/% months(1)
# year  <- 2017:2021
# startng <- c(interval(startdate,"2017-09-16") %/% months(1)+1,
#              interval(startdate,"2018-09-15") %/% months(1)+1,
#              interval(startdate,"2019-09-14") %/% months(1)+1,
#              interval(startdate,"2020-09-12") %/% months(1)+1,
#              interval(startdate,"2021-09-18") %/% months(1)+1
#             )

# endng <- c(interval(startdate,"2018-01-07") %/% months(1)+1,
#            interval(startdate,"2019-01-06") %/% months(1)+1,
#            interval(startdate,"2020-01-05") %/% months(1)+1,
#            interval(startdate,"2021-01-03") %/% months(1)+1,
#            interval(startdate,"2022-01-09") %/% months(1)+1
#             )

# startgun <- c(interval(startdate,"2017-11-18") %/% months(1)+1,
#               interval(startdate,"2018-11-17") %/% months(1)+1,
#               interval(startdate,"2019-11-23") %/% months(1)+1,
#               interval(startdate,"2020-11-21") %/% months(1)+1,
#               interval(startdate,"2021-11-20") %/% months(1)+1
#             )

# endgun <- c(interval(startdate,"2017-11-26") %/% months(1)+1,
#             interval(startdate,"2018-11-25") %/% months(1)+1,
#             interval(startdate,"2019-12-01") %/% months(1)+1,
#             interval(startdate,"2020-11-29") %/% months(1)+1,
#             interval(startdate,"2021-11-28") %/% months(1)+1
#             )


# Z_ng <- rep(0,nT_period_month)
# Z_gun <- rep(0,nT_period_month)
# for(i in 1:5){
#     Z_ng[startng[i]:endng[i]] <- 1
#     Z_gun[startgun[i]:endgun[i]] <- 1
# }

# Z_ng
# Z_gun