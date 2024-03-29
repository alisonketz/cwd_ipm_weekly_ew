#######################################
### Data for Model Fitting
#######################################

nimData <- list(
    Z_period = Z_period,
    Z_age = Z_age,
    num_period = num_period,
    adj_period = adj_period,
    weights_period = weights_period,
    age_lookup_f = age_lookup_f,
    age_lookup_m = age_lookup_m,
    # period_effect_survival = rep(NA,nT_period_overall),
    # y_hunt = 1,
    # hunt_ageweeks = d_fit_hunt$ageweeks,
    # hunt_sex = d_fit_hunt$sex,
    # hunt_age2date = d_fit_hunt$birthweek - 1,
    # hunt_n_cases = d_fit_hunt$n_cases,
    # hunt_test_status = d_fit_hunt$teststatus,
    y_hunt_pos = 1,
    hunt_pos_ageweeks = d_fit_hunt_pos$ageweeks,
    hunt_pos_sex = d_fit_hunt_pos$sex,
    hunt_pos_age2date = d_fit_hunt_pos$birthweek - 1,
    hunt_pos_n_cases = d_fit_hunt_pos$n_cases,
    y_hunt_neg = 1,
    hunt_neg_ageweeks = d_fit_hunt_neg$ageweeks,
    hunt_neg_sex = d_fit_hunt_neg$sex,
    hunt_neg_age2date = d_fit_hunt_neg$birthweek - 1,
    hunt_neg_n_cases = d_fit_hunt_neg$n_cases,
    y_sus_cens_posttest = 1,
    sus_cens_posttest_left_age_e = d_fit_sus_cens_posttest$left_age_e,
    sus_cens_posttest_right_age_r = d_fit_sus_cens_posttest$right_age_r,
    sus_cens_posttest_sex = d_fit_sus_cens_posttest$sex,
    sus_cens_posttest_age2date = sus_cens_posttest_age2date,
    y_sus_cens_postno = 1,
    sus_cens_postno_left_age_e = c(d_fit_sus_cens_postno$left_age_e, d_fit_endlive$left_age_e),
    sus_cens_postno_right_age_r = c(d_fit_sus_cens_postno$right_age_r, d_fit_endlive$right_age_r),
    sus_cens_postno_sex = c(d_fit_sus_cens_postno$sex, d_fit_endlive$sex),
    sus_cens_postno_age2date = c(sus_cens_postno_age2date, endlive_age2date),
    y_sus_mort_posttest = 1,
    sus_mort_posttest_left_age_e = d_fit_sus_mort_posttest$left_age_e,
    sus_mort_posttest_right_age_r = d_fit_sus_mort_posttest$right_age_r,
    sus_mort_posttest_right_age_s = d_fit_sus_mort_posttest$right_age_s,
    sus_mort_posttest_sex = d_fit_sus_mort_posttest$sex,
    sus_mort_posttest_age2date = sus_mort_posttest_age2date,
    y_sus_mort_postno = 1,
    sus_mort_postno_left_age_e = d_fit_sus_mort_postno$left_age_e,
    sus_mort_postno_right_age_r = d_fit_sus_mort_postno$right_age_r,
    sus_mort_postno_right_age_s = d_fit_sus_mort_postno$right_age_s,
    sus_mort_postno_dn1 = d_fit_sus_mort_postno$left_age_e,
    sus_mort_postno_sex = d_fit_sus_mort_postno$sex,
    sus_mort_postno_age2date = sus_mort_postno_age2date,
    y_icap_cens = 1,
    icap_cens_left_age_e = d_fit_icap_cens$left_age_e,
    icap_cens_right_age_r = d_fit_icap_cens$right_age_r,
    icap_cens_sex = d_fit_icap_cens$sex,
    icap_cens_age2date = icap_cens_age2date,
    y_icap_mort = 1,
    icap_mort_left_age_e = d_fit_icap_mort$left_age_e,
    icap_mort_right_age_r = d_fit_icap_mort$right_age_r,
    icap_mort_right_age_s = d_fit_icap_mort$right_age_s,
    icap_mort_sex = d_fit_icap_mort$sex,
    icap_mort_age2date = icap_mort_age2date,
    y_rec_neg_cens_posttest = 1,
    rec_neg_cens_posttest_left_age_e = d_fit_rec_neg_cens_posttest$left_age_e,
    rec_neg_cens_posttest_right_age_r = d_fit_rec_neg_cens_posttest$right_age_r,
    rec_neg_cens_posttest_sex = d_fit_rec_neg_cens_posttest$sex,
    rec_neg_cens_posttest_age2date = rec_neg_cens_posttest_age2date,
    y_rec_neg_mort = 1,
    rec_neg_mort_left_age_e = d_fit_rec_neg_mort$left_age_e,
    rec_neg_mort_right_age_r = d_fit_rec_neg_mort$right_age_r,
    rec_neg_mort_right_age_s = d_fit_rec_neg_mort$right_age_s,
    rec_neg_mort_sex = d_fit_rec_neg_mort$sex,
    rec_neg_mort_age2date = rec_neg_mort_age2date,
    y_rec_pos_mort = 1,
    rec_pos_mort_left_age_e = d_fit_rec_pos_mort$left_age_e,
    rec_pos_mort_right_age_r = d_fit_rec_pos_mort$right_age_r,
    rec_pos_mort_right_age_s = d_fit_rec_pos_mort$right_age_s,
    rec_pos_mort_dn1 = d_fit_rec_pos_mort$left_age_e,
    rec_pos_mort_dn = d_fit_rec_pos_mort$ageweek_recap,
    rec_pos_mort_sex = d_fit_rec_pos_mort$sex,
    rec_pos_mort_age2date = rec_pos_mort_age2date,
    y_rec_pos_cens = 1,
    rec_pos_cens_left_age_e = d_fit_rec_pos_cens$left_age_e,
    rec_pos_cens_right_age_r = d_fit_rec_pos_cens$right_age_r,
    rec_pos_cens_dn1 = d_fit_rec_pos_cens$ageweek_recap,
    rec_pos_cens_dn = d_fit_rec_pos_cens$right_age_r,
    rec_pos_cens_sex = d_fit_rec_pos_cens$sex,
    rec_pos_cens_age2date = rec_pos_cens_age2date,
    y_idead = 1,
    idead_left_age_e = d_fit_idead$left_age_e,
    idead_right_age_r = d_fit_idead$right_age_r,
    idead_right_age_s = d_fit_idead$right_age_s,
    idead_dn1 = d_fit_idead$left_age_e,
    idead_dn = d_fit_idead$right_age_s,
    idead_sex = d_fit_idead$sex,
    idead_age2date = idead_age2date,
    y_aah = 1,
    aah_ageweeks = d_fit_notest$ageweeks,
    aah_sex = d_fit_notest$sexnum,
    aah_age2date = d_fit_notest$age2date_weeks,
    aah_n = d_fit_notest$n,
    mort_hh = d_fit_hh$mort_h,
    sex_cause = d_fit_hh$sex,
    Z_cause_gun = Z_cause_gun,
    Z_cause_ng = Z_cause_ng,
    # Cage_less = Cage_less,
    # Cage_ant = Cage_ant,
    # O = Ototal,
    lobs = log(Ototal),
    f_logpop_sus = f_logpop_sus,
    f_logpop_inf = f_logpop_inf,
    m_logpop_sus = m_logpop_sus,
    m_logpop_inf = m_logpop_inf,
    obs_ct_fd_alpha = obs_ct_fd_alpha,
    obs_ct_fd_beta = obs_ct_fd_beta,
    Nfawn = fawndoe_df$overall_fawn,
    Ndoe = fawndoe_df$overall_doe
    )

#######################################
### Constants for MCMC
#######################################

nimConsts <- list(n_year = n_year,
    n_year_precollar = n_year_precollar,
    n_year_precollar_ext = n_year_precollar_ext,
    n_year_prestudy_ext = n_year_prestudy_ext,
    n_study_area = n_study_area,
    n_sex = n_sex,
    n_agef = n_agef,
    n_agem = n_agem,
    n_ageclassf = n_ageclassf,
    n_ageclassm = n_ageclassm,
    n_sex = n_sex,
    # cal = d_fit_season$ng_end - d_fit_season$yr_start,
    # sizeCage_f = sizeCage_f,
    # sizeCage_m = sizeCage_m,
    report_hyp_all = report_hyp_all,
    report_hyp_y = report_hyp_y,
    nT_period_overall = nT_period_overall,
    nT_period_overall_ext = nT_period_overall_ext,
    nT_period_precollar = nT_period_precollar,
    nT_period_precollar_ext = nT_period_precollar_ext,
    nT_period_collar = nT_period_collar,
    nT_period_prestudy_ext = nT_period_prestudy_ext,
    nT_age_surv = nT_age_surv,
    nT_age_surv_aah_f = nT_age_surv_aah_f,
    nT_age_surv_aah_m = nT_age_surv_aah_m,
    nT_age_short_f = nT_age_short_f,
    nT_age_short_m = nT_age_short_m,
    n_year_fec_early = n_year_fec_early,
    nknots_age = nknots_age,
    nknots_period = nknots_period,
    n_adj_period = n_adj_period,
    period_lookup_foi = period_lookup_foi,
    period_lookup_foi_study = period_lookup_foi_study,
    ng_start = d_fit_season$ng_start,
    gun_start = d_fit_season$gun_start,
    gun_end = d_fit_season$gun_end,
    ng_end = d_fit_season$ng_end,
    yr_start = d_fit_season$yr_start,
    yr_end = d_fit_season$yr_end,
    # nHarvest = nrow(d_fit_hunt),
    nInfHarvest = nrow(d_fit_hunt_pos),
    nSusHarvest = nrow(d_fit_hunt_neg),
    nSusCensTest = nrow(d_fit_sus_cens_posttest),
    nSusCensNo = nrow(d_fit_sus_cens_postno) +
                  nrow(d_fit_endlive),
    nSusMortTest = nrow(d_fit_sus_mort_posttest),
    nSusMortNoTest = nrow(d_fit_sus_mort_postno),
    nIcapCens = nrow(d_fit_icap_cens),
    nIcapMort = nrow(d_fit_icap_mort),
    nRecNegCensTest = nrow(d_fit_rec_neg_cens_posttest),
    nRecNegMort = nrow(d_fit_rec_neg_mort),
    nRecPosMort = nrow(d_fit_rec_pos_mort),
    nNegCapPosMort = nrow(d_fit_idead),
    nAAH = nrow(d_fit_notest),
    # sect_hunt = d_fit_hunt$ew,
    sect_hunt_pos = d_fit_hunt_pos$ew,
    sect_hunt_neg = d_fit_hunt_neg$ew,
    sect_sus_cens_posttest = d_fit_sus_cens_posttest$study_area,
    sect_sus_cens_postno = c(d_fit_sus_cens_postno$study_area,
                            d_fit_endlive$study_area),
    sect_sus_mort_posttest = d_fit_sus_mort_posttest$study_area,
    sect_sus_mort_postno = d_fit_sus_mort_postno$study_area,
    sect_icap_cens = d_fit_icap_cens$study_area,
    sect_icap_mort = d_fit_icap_mort$study_area,
    sect_rec_neg_cens_posttest = d_fit_rec_neg_cens_posttest$study_area,
    sect_rec_neg_mort = d_fit_rec_neg_mort$study_area,
    sect_rec_pos_mort = d_fit_rec_pos_mort$study_area,
    sect_rec_pos_cens = d_fit_rec_pos_cens$study_area,
    sect_idead = d_fit_idead$study_area,
    sect_aah = d_fit_notest$study_area,
    records_cause = records_cause,
    interval_cause = d_fit_hh$right_period_s - 1,
    # indx_mat_pe_surv = indx_mat_pe_surv,
    # intvl_step_yr = intvl_step_yr_weekly,
    m_period_foi = rep(0,n_year),
    f_period_foi = rep(0,n_year),
    # age_effect_survival = rep(0, nT_age_surv),
    # period_annual_survival = rep(0,n_year_precollar + 1),
    period_effect_survival = rep(0,nT_period_overall_ext),
    Z_foi_cgam = Z_foi_cgam,
    nknots_foi_cgam = nknots_foi_cgam,
    Z_foi_spline = Z_foi_spline,
    nknots_foi_spline = nknots_foi_spline,
    Z_collar_gun = Z_collar_gun,
    Z_collar_ng = Z_collar_ng
)


#######################################
### Initial Values for MCMC
#######################################

initsFun <- function()list(beta_male = rnorm(1, -.2, .01),
    beta0_sus_temp = rnorm(1, -5.75, 0.0001),
    # beta0_survival_sus = rnorm(1, -8.5, 0.0001),
    sus_mix = 1,
    beta0_inf_temp = rnorm(1, -3, 0.0001),
    # beta0_survival_inf = rnorm(1, -8, 0.0001),
    inf_mix = 1,
    ln_b_age_survival = rnorm(nknots_age) * 10^-4,
    # tau_age_survival = runif(1, .1, 1),
    # b_period_survival = rnorm(nknots_period) * 10^-4,
    # tau_period_survival = runif(1, .1, 1),
    # beta_harvest_gun = rnorm(1, 0, sd = 1),
    # beta_harvest_ng = rnorm(1, 0, sd = 1),
    # tau_period_precollar = rgamma(1, 1, 1),
    # period_annual_survival = rnorm(n_year_precollar + 1, 0, .1),
    # period_int_survival = rnorm(1, 0, .1),
    # m_period_foi = seq(-.25, .25, length = n_year),
    # f_period_foi = seq(-.25, .25, length = n_year),
    # m_period_foi_temp = seq(-.25, .25, length = n_year),
    # f_period_foi_temp = seq(-.25, .25, length = n_year),
    # ln_b_foi_cgam_m = rnorm(nknots_foi_cgam) * 10^-4,
    # ln_b_foi_cgam_f = rnorm(nknots_foi_cgam) * 10^-4,
    # f_b_foi_spline = rnorm(nknots_foi_spline) * 10^-4,
    # m_b_foi_spline = rnorm(nknots_foi_spline) * 10^-4,
    # tau_foi_spline_male = runif(1, .1, 1),
    # tau_foi_spline_female = runif(1, .1, 1),
    # tau_period_foi_male = runif(1, 4.2, 6.8),
    # tau_period_foi_female = runif(1, 2.27, 3.44),
    tau_age_foi_male = runif(1, 1.5, 1.7),
    tau_age_foi_female = runif(1, 2.7, 4.2),
    m_age_foi = c(rnorm(1, -10.7, sd = .1),
                  rnorm(1, -7.55, sd = .1),
                  rnorm(1, -6.15, sd = .1),
                  rnorm(1, -5.05, sd = .1),
                  rnorm(1, -4.6, sd = .1),
                  rnorm(1, -4.4, sd = .1)),
    f_age_foi = c(rnorm(1, -10.4, sd = .1),
                  rnorm(1, -8.5, sd = .1),
                  rnorm(1, -7.2, sd = .1),
                  rnorm(1, -6.7, sd = .1),
                  rnorm(1, -6.7, sd = .1),
                  rnorm(1, -6.5, sd = .1),
                  rnorm(1, -6.4, sd = .1)),
    beta0_cause = rnorm(1, -2.8, .1),
    beta_cause_male = rnorm(1, 0, .1),
    beta_cause_gun = rnorm(1, 1.5, .1),
    beta_cause_ng = rnorm(1, 3, .1),
    # tau_obs = matrix(runif(4, 1, 3), 2, 2),
    # tau_obs = runif(1, .01, .5),
    # tau_pop = runif(2, .05, .25),
    # report_overall = report_overall_init,
    # report = report_init,
    # fec_epsilon = fec_eps_init,
    # mu_fec = rnorm(1, mu_fec_init, .01),
    # fec_prec_eps = runif(1, 5, 10),
    space_temp = rnorm(1, -.55, .01),
    space_mix = 1
    )
nimInits <- initsFun()

########################################################
### Build and compile model in R
########################################################

# start_Rmodel <- Sys.time()
Rmodel <- nimbleModel(code = modelcode,
                      constants = nimConsts,
                      data = nimData,
                      inits = initsFun(),
                      calculate = FALSE,
                      check = FALSE
                      )
# end_Rmodel <- Sys.time() - start_Rmodel
Rmodel$initializeInfo()

Cnim <- compileNimble(Rmodel)

for(i in 1:10){beepr::beep(1)}

#######################################
### Parameters to trace in MCMC
#######################################

parameters <- c(
              "beta_male",
              "tau_age_foi_male",
              "tau_age_foi_female",
            #   "tau1_age_foi_male",
            #   "tau1_age_foi_female",
              "m_age_foi",
              "f_age_foi",
              "m_age_foi_mu",
              "f_age_foi_mu",
              # "tau_period_foi_male",
              # "tau_period_foi_female",
              # "f_period_foi",
              # "m_period_foi",
              "space",
              "space_temp",
              "space_mix",
              "beta0_survival_inf",
              "beta0_survival_sus",
              # "beta_harvest_gun",
              # "beta_harvest_ng",
            #   "beta0_sus_temp",
            #   "sus_mix",
            #   "beta0_inf_temp",
            #   "inf_mix",
              # "tau_age_survival",
              "age_effect_survival",
              "ln_b_age_survival",
              # "b_period_survival",
              # "tau_period_survival",
              # "tau_period_precollar",
              # "period_effect_survival",
              # "period_int_survival",
              "beta0_cause",
              "beta_cause_gun",
              "beta_cause_ng",
              "beta_cause_male",
              "p_nogun_f",
              "p_gun_f",
              "p_nogun_m",
              "p_gun_m",
            #   "report",
            #   "fec",
            #   "mu_fec",
            #   "fec_prec_eps",
            #   "fec_epsilon",
              "psi",
              "sn_inf",
              "sn_sus",
              "sh_inf",
              "sh_sus"#,
            #   "mu_obs",
            #   "tau_obs",
            #   "tau_pop"
               )

confMCMC <- configureMCMC(Rmodel,
                         monitors = parameters,
                         thin = 1,
                         # enableWAIC = TRUE,
                         useConjugacy = FALSE)
nimMCMC <- buildMCMC(confMCMC)
CnimMCMC <- compileNimble(nimMCMC,
                         project = Rmodel)
for(i in 1:10){beepr::beep(1)}

reps <- 1000
bin <- 0 # reps*.5
n_chains <- 1

set.seed(1001)
starttime <- Sys.time()
mcmcout <- runMCMC(CnimMCMC,
                  niter = reps,
                  nburnin = bin,
                  nchains = n_chains,
                  inits = initsFun,
                  samplesAsCodaMCMC = TRUE,
                  summary = TRUE
                  )
runtime <- difftime(Sys.time(),
                    starttime,
                    units = "min")
runtime
for (i in 1:10) {beepr::beep(1)}

# mcmcout$summary
# end_Rmodel
# endtime_rmodel_compile
# endtime_mcmc
# runtime

# sink("runtime_allsteps.txt")
# cat("Rmodel:\n")
# end_Rmodel
# cat("\nCompile Rmodel:\n")
# endtime_rmodel_compile
# cat("\nCompile MCMC:\n")
# endtime_mcmc
# cat("\nRun MCMC 100 iter: ",runtime)
# sink()
#############################################################
###
### Running in parallel, restartable
###
#############################################################

# reps  <- 10
# bin <- reps * .5
# n_thin <- 1
# n_chains <- 3
# starttime <- Sys.time()
# cl <- makeCluster(n_chains, timeout = 5184000)

# clusterExport(cl, c("modelcode",
#                     "initsFun",
#                     "nimData",
#                     "nimConsts",
#                     "parameters",
#                     "reps",
#                     "bin",
#                     "n_thin",
#                     "set_period_effects_constant",
#                     "dInfHarvest",
#                     "dSusHarvest",
#                     "dSusCensTest",
#                     "dSusCensNo",
#                     "dSusMortTest",
#                     "dSusMortNoTest",
#                     "dIcapCens",
#                     "dIcapMort",
#                     "dRecNegCensTest",
#                     "dRecNegMort",
#                     "dRecPosMort",
#                     "dRecPosCens",
#                     "dNegCapPosMort",
#                     "dAAH",
#                     "calc_surv_aah",
#                     "calc_surv_harvest",
#                     "calc_infect_prob"
#                     ))
# for (j in seq_along(cl)) {
#   set.seed(j + 1000)
#   init <- initsFun()
#   clusterExport(cl[j], "init")
# }
# for (i in 1:10) {beepr::beep(1)}

# # starttime <- Sys.time()
# # mcmcout1 <-  mcmc.list(clusterEvalQ(cl, {
# #   library(nimble)
# #   library(coda)

#   compile(dInfHarvest)
#   compile(dSusHarvest)
#   compile(dSusCensTest)
#   compile(dSusCensNo)
#   compile(dSusMortTest)
#   compile(dSusMortNoTest)
#   compile(dIcapCens)
#   compile(dIcapMort)
#   compile(dRecNegCensTest)
#   compile(dRecNegMort)
#   compile(dRecPosMort)
#   compile(dRecPosCens)
#   compile(dNegCapPosMort)
#   compile(dAAH)
#   compile(set_period_effects_constant)
#   compile(calc_surv_aah)
#   compile(calc_surv_harvest)
#   compile(calc_infect_prob)

# #   ##############################################################
# #   ###
# #   ### Execute MCMC
# #   ###
# #   ##############################################################

# #   Rmodel <- nimbleModel(code = modelcode,
# #                         name = "modelcode",
# #                         constants = nimConsts,
# #                         data = nimData,
# #                         inits = initsFun)
# #   Cnim <- compileNimble(Rmodel)
# #   confMCMC <- configureMCMC(Rmodel,
# #                             monitors = parameters,
# #                             thin = n_thin,
# #                             useConjugacy = FALSE)
# #   nimMCMC <- buildMCMC(confMCMC)
# #   CnimMCMC <- compileNimble(nimMCMC,
# #                             project = Rmodel)

# #   CnimMCMC$run(reps, reset = FALSE)

# #   return(as.mcmc(as.matrix(CnimMCMC$mvSamples)))
# # }))
# # runtime1 <- difftime(Sys.time(), starttime, units = "min")

# # for(chn in 1:nc) { # nc must be > 1
# #   ind.keep <- c()
# #   for(p in 1:length(parameters)) ind.keep <-
# #       c(ind.keep, which(str_detect(dimnames(out1[[chn]])[[2]], parameters[p]))) %>% unique()
# #   out1[[chn]] <- out1[[chn]][,ind.keep]
# # }

# # ## Check convergence ##
# # out2 <- out1
# # ni.saved <- nrow(out2[[1]])
# # for(chn in 1:nc) { # nc must be > 1
  
# #   if(nb < 1) {
# #     nb.real <- (round(ni.saved * nb)+1)
# #   } else {
# #     nb.real <- (round(nb/nt)+1)
# #   }
# #   out2[[chn]] <- out2[[chn]][nb.real:ni.saved,]
# # }
# # out.mcmc <- coda::as.mcmc.list(lapply(out2, coda::as.mcmc))
# # stopCluster(cl)

# # save(mcmcout1, file = "mcmcout1.Rdata")
# # save(runtime1, file = "runtime1.Rdata")
# # save(endtime_rmodel_compile, file = "endtime_rmodel_compile.Rdata")
# # save(endtime_mcmc, file = "endtime_mcmc.Rdata")

# #not calculating waic, because too many params would need to be traced
# # posteriorSamplesMatrix <- rbind(mcmcout[[1]], mcmcout[[2]], mcmcout[[3]])
# # CnimMCMC$run(5)   ## non-zero number of iterations
# # nimble:::matrix2mv(posteriorSamplesMatrix, CnimMCMC$mvSamples)
# # # CnimMCMC$enableWAIC <- TRUE
# # waic_spline <- calculateWAIC(posteriorSamplesMatrix,Rmodel)


# # waic_spline_covs <- mcmcout$WAIC
# # save(waic_spline, file = "waic_spline.Rdata")



# ###
# ### save model run
# ###

# # save(runtime,file="results/runtime.Rdata")
# # save(mcmcout,file="results/mcmcout.Rdata")