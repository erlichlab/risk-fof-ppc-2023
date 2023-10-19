## get the behavior data

# get bilateral fof opto silencing data
get_bi_opto_data <- function(git_root) {
  df <- read.csv(file.path(git_root,"csv","opto_bilateral.csv"))
  df <- df %>%
    filter(subj_choice != "violation") %>%
    filter(!is.na(lottery_poke), !is.na(surebet_poke)) %>%
    group_by(sessid) %>%
    mutate(lottery_norm = lottery_mag / max(lottery_mag),
           sb_norm = sb_mag / max(lottery_mag),
           dose=isopto,
           chose_lottery=ifelse(subj_choice=="lottery",1,0)) %>% ungroup()
  return(df)
}

# get unilateral fof opto silencing data
get_uni_opto_data <- function(git_root) {
  df <- read.csv(file.path(git_root,"csv","opto_unilateral.csv"))
  df <- df %>%
    filter(subj_choice != "violation") %>%
    filter(!is.na(lottery_poke), !is.na(surebet_poke)) %>%
    group_by(sessid) %>%
    mutate(lottery_norm = lottery_mag / max(lottery_mag),
           sb_norm = sb_mag / max(lottery_mag),
           dose=isopto,
           chose_lottery=ifelse(subj_choice=="lottery",1,0)) %>% ungroup()
  return(df)
}

# get bilateral fof muscimol silencing data
get_bifof_muscimol_data <- function(git_root) {
  control_df <- read.csv(file.path(git_root,"csv", "figure_behavior_population.csv")) %>% mutate(dose=0)
  df <- read.csv(file.path(git_root,"csv","figure_infusion_bi_fof.csv")) %>%
    filter(dosage == 0.3) %>% mutate(dose=1)
  df <- rbind(df, control_df)
  df <- df %>%
    filter(subj_choice != "violation") %>%
    filter(!is.na(lottery_poke), !is.na(surebet_poke)) %>%
    group_by(sessid) %>%
    mutate(lottery_norm = lottery_mag / max(lottery_mag),
           sb_norm = sb_mag / max(lottery_mag),
           dose=case_when(infusion == 'control' & volume == 0 ~0,
                          infusion == 'infusion' & volume == 0.3 ~ 1),
           chose_lottery=ifelse(subj_choice=="lottery",1,0)) %>% ungroup()
  return(df)
}

# get unilateral fof muscimol silencing data
get_unifof_muscimol_data <- function(git_root) {
  control_df <- read.csv(file.path(git_root,"csv", "figure_behavior_population.csv")) %>% mutate(dose=0)
  df <- read.csv(file.path(git_root,"csv","figure_infusion_uni_fof.csv")) %>%
    filter(dosage == 0.3) %>% mutate(dose=1)

  df <- rbind(df, control_df)
  df <- df %>%
    filter(subj_choice != "violation") %>%
    filter(!is.na(lottery_poke), !is.na(surebet_poke)) %>%
    group_by(sessid) %>%
    mutate(lottery_norm = lottery_mag / max(lottery_mag),
           sb_norm = sb_mag / max(lottery_mag),
           dose=case_when(infusion == 'control' & volume == 0 ~0,
                          infusion == 'infusion' & volume == 0.3 ~ 1),
           chose_lottery=ifelse(subj_choice=="lottery",1,0)) %>% ungroup()
  return(df)
}

# get bilateral ppc muscimol silencing data
get_bippc_muscimol_data <- function(git_root) {
  control_df <- read.csv(file.path(git_root,"csv", "figure_behavior_population.csv")) %>% mutate(dose=0)
  df <- read.csv(file.path(git_root,"csv","figure_infusion_bi_ppc.csv")) %>%
    filter(dosage == 0.3) %>% mutate(dose=1)

  df <- rbind(df, control_df)
  df <- df %>%
    filter(subj_choice != "violation") %>%
    filter(!is.na(lottery_poke), !is.na(surebet_poke)) %>%
    group_by(sessid) %>%
    mutate(lottery_norm = lottery_mag / max(lottery_mag),
           sb_norm = sb_mag / max(lottery_mag),
           dose=case_when(infusion == 'control' & volume == 0 ~0,
                          infusion == 'infusion' & volume == 0.3 ~ 1),
           chose_lottery=ifelse(subj_choice=="lottery",1,0)) %>% ungroup()
  return(df)
}

# get unilateral ppc muscimol silencing data
get_unippc_muscimol_data <- function(git_root) {
  control_df <- read.csv(file.path(git_root,"csv", "figure_behavior_population.csv")) %>% mutate(dose=0)
  df <- read.csv(file.path(git_root,"csv","figure_infusion_uni_ppc.csv")) %>%
    filter(dosage == 0.3) %>% mutate(dose=1)

  df <- rbind(df, control_df)
  df <- df %>%
    filter(subj_choice != "violation") %>%
    filter(!is.na(lottery_poke), !is.na(surebet_poke)) %>%
    group_by(sessid) %>%
    mutate(lottery_norm = lottery_mag / max(lottery_mag),
           sb_norm = sb_mag / max(lottery_mag),
           dose=case_when(infusion == 'control' & volume == 0 ~0,
                          infusion == 'infusion' & volume == 0.3 ~ 1),
           chose_lottery=ifelse(subj_choice=="lottery",1,0)) %>% ungroup()
  return(df)
}

# get ephys data and randomly pick 10% trials as drug
get_ephys_data <- function(git_root) {
  df <- read.csv('csv/ephys_behavior.csv', na.strings=c("","NA"))
  df$dose <- ifelse(runif(nrow(df))<0.1, 1, 0)
  df <- df %>%
    filter(subj_choice != "violation") %>%
    filter(!is.na(lottery_poke), !is.na(surebet_poke)) %>%
    group_by(sessid) %>%
    mutate(lottery_norm = lottery_mag / max(lottery_mag),
           sb_norm = sb_mag / max(lottery_mag),
           chose_lottery=ifelse(subj_choice=="lottery",1,0)) %>% ungroup()
  return(df)
}

## functions

# binomialize the data ahead
binomialize_brms <- function(df) {
  bdf <- df %>%
    group_by(subjid, lottery_norm, sb_norm, dose, lottery_prob) %>%
    summarise(total = n(), prop_lott = sum(chose_lottery == 1))
  return(bdf)
}

inv_logit <- function(x) {
  return(1 / (1 + exp(-x)))
}

normal_cdf <- pnorm
Phi_approx <- function(x) {
  return(pnorm(x, 0, 1))
}

## brms model fitting
f_drug_phi <- bf(
  prop_lott | trials(total) ~ inv_logit(omega1) *
    (Phi_approx(
      (lottery_prob * (lottery_norm) ^ exp(logrho) -
         sb_norm ^ exp(logrho)) /
        (sqrt(2)*exp(noise)))) +
    (1-inv_logit(omega1)) * inv_logit(omega2),
  logrho ~ dose + (1 | subjid),
  noise ~ dose + (1 | subjid),
  omega1 ~ dose + (1 | subjid),
  omega2 ~ dose + (1 | subjid),
  nl = TRUE)

priors <- prior(normal(0, 0.5), nlpar = 'logrho') +
  prior(normal(3, 1), nlpar = 'omega1') +
  prior(normal(0, 1), nlpar = 'omega2') +
  prior(normal(-3, 0.3), nlpar = 'noise') +
  prior(normal(0,0.35), lb=0, class = sd, group = 'subjid', nlpar = 'omega1') +
  prior(normal(0,0.35), lb=0, class = sd, group = 'subjid', nlpar = 'omega2') +
  prior(normal(0,0.35), lb=0, class = sd, group = 'subjid', nlpar = 'noise') +
  prior(normal(0,0.35), lb=0, class = sd, group = 'subjid', nlpar = 'logrho') +
  prior(normal(0,0.5), class = b, nlpar = 'logrho', coef = 'dose') +
  prior(normal(0,0.5), class = b, nlpar = 'noise', coef = 'dose') +
  prior(normal(0,1), class = b, nlpar = 'omega1', coef = 'dose') +
  prior(normal(0,1), class = b, nlpar = 'omega2', coef = 'dose')

## Synthetic data
make_agents <- function(n) {
  agents <- data.frame(
    subjid = stri_rand_strings(n, 3, pattern = "[A-Z]"),
    log_rho = runif(n, min = -0.5, max = 0.2),
    beta = runif(n) * -2 - 2,
    omega_1 = runif(n, min = 1, max = 3),
    omega_2 = runif(n, min = -2, max = 2)
  )
  return(agents)
}

make_trials <- function(n_trials= 5000, rew_multi = 8) {
  df <- data.frame(
    lottery_prob = 0.5,
    lottery_mag = rew_multi * sample(c(0,2,4,8,16,32), n_trials, replace=TRUE),
    sb_mag = rew_multi * 3
  )

  df$lottery_norm <- df$lottery_mag / max(df$lottery_mag)
  df$sb_norm <- df$sb_mag / max(df$lottery_mag)

  return(df)
}

make_agent_choice_cdf <- function(trials, agent) {
  df <- trials
  df$chose_lottery <- ifelse(
    runif(nrow(df)) < (
      inv_logit(agent$omega_1)*(1 - normal_cdf(0,
                                               df$lottery_prob * (df$lottery_norm) ^ exp(agent$log_rho) -
                                                 df$sb_norm ^ exp(agent$log_rho),
                                               sqrt(2)*exp(agent$beta))) +
        (1-inv_logit(agent$omega_1)) * inv_logit(agent$omega_2)),
    1, 0)
  df$subjid <- agent$subjid
  return(df)
}

drug_agent <- function(agent, deltaRho=-0.5, deltaOmega1=0, deltaOmega2=0, deltaBeta=0) {
  agent$log_rho <- agent$log_rho + deltaRho
  agent$omega_1 <- agent$omega_1 + deltaOmega1
  agent$omega_2 <- agent$omega_2 + deltaOmega2
  agent$beta <- agent$beta + deltaBeta
  return(agent)
}

drug_sim <- function(agents,  deltaRho=-0.5, deltaOmega1=0, deltaOmega2=0, deltaBeta=0)
{
  tdf <- make_trials(500)
  df <- NULL
  for (i in 1:nrow(agents)) {
    df0 <- make_agent_choice_cdf(tdf, agents[i,])
    df0$dose <- 0
    df1 <- make_agent_choice_cdf(tdf, drug_agent(agents[i,], deltaRho, deltaOmega1, deltaOmega2, deltaBeta))
    df1$dose <- 1
    df <- rbind(df, df0, df1)
  }
  return(df)

}

##################



make_choices <- function(trials,
    subjid="A", log_rho=-0.1, beta=15, omega_1=2, omega_2=1) {
    # make data
    print(c(subjid, log_rho, beta, omega_1, omega_2))
    df <- trials

    # make choices

    df$chose_lottery <- ifelse(
        runif(nrow(df)) < (
            inv_logit(omega_1)*inv_logit(
                beta*(df$lottery_prob * (df$lottery_norm) ^ exp(log_rho)  - (df$sb_norm) ^ exp(log_rho))) +
                (1-inv_logit(omega_1)) * inv_logit(omega_2)
            ),
        1, 0)
    df$subjid <- subjid
    return(df)

}

make_agent_choice <- function(trials, agent) {
    df <- trials
    df$chose_lottery <- ifelse(
        runif(nrow(df)) < (
            inv_logit(agent$omega_1)*inv_logit(
                agent$beta*(df$lottery_prob * (df$lottery_norm) ^ exp(agent$log_rho) -
                (df$sb_norm) ^ exp(agent$log_rho))) +
            (1-inv_logit(agent$omega_1)) * inv_logit(agent$omega_2)),
        1, 0)
    df$subjid <- agent$subjid
    return(df)
}










plot_comparison <- function(fit, agents){
    fit_param = agents %>% mutate(type='fit')
    data_param = agents %>% mutate(type='data')

    df = rbind(fit_param, data_param)
}




plot_drug <- function(drugdf){
    g = ggplot(drugdf, aes(x=lottery_norm, y=chose_lottery, group=dose, color=factor(dose))) +
    stat_summary() +
    geom_smooth(method='glm', method.args=list(family='binomial')) +
    facet_wrap(~subjid) +
    theme(text=element_text(size=16,  family="sans"))
    return(g)
}

f_id <- bf(
    prop_lott | trials(total) ~ inv_logit(omega1) *
        inv_logit(beta*(lottery_prob * (lottery_norm) ^ exp(logrho)
            - (sb_norm) ^ exp(logrho))) +
        (1-inv_logit(omega1)) * inv_logit(omega2),
    logrho ~ 1 + (1|id|subjid),
    beta ~ 1 + (1|id|subjid),
    omega1 ~ 1 + (1|id|subjid),
    omega2 ~ 1 + (1|id|subjid),
    nl = TRUE)
    # This model gives noisier estimates.
    # I think there are also more parameters to estimate.


f <- bf(
    prop_lott | trials(total) ~ inv_logit(omega1) *
        inv_logit(
            (lottery_prob * (lottery_norm) ^ exp(logrho) -
                (sb_norm) ^ exp(logrho)) / exp(noise)
            ) +
        (1-inv_logit(omega1)) * inv_logit(omega2),
    logrho ~ 1 + (1 | subjid),
    noise ~ 1 + (1 | subjid),
    omega1 ~ 1 + (1 | subjid),
    omega2 ~ 1 + (1 | subjid),
    nl = TRUE)


f_cdf <- bf(
    prop_lott | trials(total) ~ inv_logit(omega1) *
        (1 - normal_cdf(0,
                lottery_prob * (lottery_norm) ^ exp(logrho) -
                sb_norm ^ exp(logrho),
                sqrt(2)*noise)) +
        (1-inv_logit(omega1)) * inv_logit(omega2),
    logrho ~ 1 + (1 | subjid),
    noise ~ 1 + (1 | subjid),
    omega1 ~ 1 + (1 | subjid),
    omega2 ~ 1 + (1 | subjid),
    nl = TRUE)



f_drug_cdf <- bf(
    prop_lott | trials(total) ~ inv_logit(omega1) *
        (1 - normal_cdf(0,
                lottery_prob * (lottery_norm) ^ exp(logrho) -
                sb_norm ^ exp(logrho),
                sqrt(2)*exp(noise))) +
        (1-inv_logit(omega1)) * inv_logit(omega2),
    logrho ~ dose + (1 | subjid),
    noise ~ dose + (1 | subjid),
    omega1 ~ dose + (1 | subjid),
    omega2 ~ dose + (1 | subjid),
    nl = TRUE)





f_drug_phi_rho <- bf(
    prop_lott | trials(total) ~ inv_logit(omega1) *
        (Phi_approx(
                (lottery_prob * (lottery_norm) ^ exp(logrho) -
                    sb_norm ^ exp(logrho)) /
                (sqrt(2)*exp(noise)))) +
        (1-inv_logit(omega1)) * inv_logit(omega2),
    logrho ~ dose + (1 | subjid),
    noise ~ dose + (1 | subjid),
    omega1 ~ 1 + (1 | subjid),
    omega2 ~ 1 + (1 | subjid),
    nl = TRUE)


f_drug_phi_omega <- bf(
    prop_lott | trials(total) ~ inv_logit(omega1) *
        (Phi_approx(
                (lottery_prob * (lottery_norm) ^ exp(logrho) -
                    sb_norm ^ exp(logrho)) /
                (sqrt(2)*exp(noise)))) +
        (1-inv_logit(omega1)) * inv_logit(omega2),
    logrho ~ 1 + (1 | subjid),
    noise ~ dose + (1 | subjid),
    omega1 ~ dose + (1 | subjid),
    omega2 ~ dose + (1 | subjid),
    nl = TRUE)


f_drug <- bf(
    prop_lott | trials(total) ~ inv_logit(omega1) *
        inv_logit(
             (lottery_prob * (lottery_norm) ^ exp(logrho) -
                (sb_norm) ^ exp(logrho)) / exp(noise)
            ) +
        (1-inv_logit(omega1)) * inv_logit(omega2),
    logrho ~ dose + (1 | subjid),
    noise ~ dose + (1 | subjid),
    omega1 ~ dose + (1 | subjid),
    omega2 ~ dose + (1 | subjid),
    nl = TRUE)


f_drug_rho_only <- bf(
    prop_lott | trials(total) ~ inv_logit(omega1) *
        inv_logit(
             (lottery_prob * (lottery_norm) ^ exp(logrho) -
                (sb_norm) ^ exp(logrho)) / exp(noise)
            ) +
        (1-inv_logit(omega1)) * inv_logit(omega2),
    logrho ~ dose + (1 | subjid),
    noise ~ dose + (1 | subjid),
    omega1 ~ 1 + (1 | subjid),
    omega2 ~ 1 + (1 | subjid),
    nl = TRUE)


f_drug_omega_only <- bf(
    prop_lott | trials(total) ~ inv_logit(omega1) *
        inv_logit(
             (lottery_prob * (lottery_norm) ^ exp(logrho) -
                (sb_norm) ^ exp(logrho)) / exp(noise)
            ) +
        (1-inv_logit(omega1)) * inv_logit(omega2),
    logrho ~ 1 + (1 | subjid),
    noise ~ dose + (1 | subjid),
    omega1 ~ dose + (1 | subjid),
    omega2 ~ dose + (1 | subjid),
    nl = TRUE)

# These work for deltaRho but not deltaOmega
bad_priors <- prior(normal(0, 0.5), nlpar = 'logrho') +
            prior(normal(3, 2), nlpar = 'omega1') +
            prior(normal(0, 1), nlpar = 'omega2') +
            prior(normal(0, 3), nlpar = 'noise') +
            prior(student_t(3,0,1), class = sd, group = 'subjid', nlpar = 'omega1') +
            prior(student_t(3,0,1), class = sd, group = 'subjid', nlpar = 'omega2') +
            prior(student_t(3,0,1), class = sd, group = 'subjid', nlpar = 'noise') +
            prior(student_t(3,0,1), class = sd, group = 'subjid', nlpar = 'logrho') +
            prior(normal(0,0.2), class = b, nlpar = 'logrho', coef = 'dose') +
            prior(normal(0,0.2), class = b, nlpar = 'noise', coef = 'dose') +
            prior(normal(0,0.2), class = b, nlpar = 'omega1', coef = 'dose') +
            prior(normal(0,0.2), class = b, nlpar = 'omega2', coef = 'dose')

loose_priors <- prior(normal(0, 0.5), nlpar = 'logrho') +
            prior(normal(3, 2), nlpar = 'omega1') +
            prior(normal(0, 1), nlpar = 'omega2') +
            prior(normal(0, 3), nlpar = 'noise') +
            prior(student_t(3,0,3), class = sd, group = 'subjid', nlpar = 'omega1') +
            prior(student_t(3,0,3), class = sd, group = 'subjid', nlpar = 'omega2') +
            prior(student_t(3,0,3), class = sd, group = 'subjid', nlpar = 'noise') +
            prior(student_t(3,0,3), class = sd, group = 'subjid', nlpar = 'logrho') +
            prior(student_t(3,0,3), class = sd, nlpar = 'omega1') +
            prior(student_t(3,0,3), class = sd, nlpar = 'omega2') +
            prior(student_t(3,0,3), class = sd, nlpar = 'noise') +
            prior(student_t(3,0,3), class = sd, nlpar = 'logrho') +
            prior(normal(0,0.2), class = b, nlpar = 'logrho', coef = 'dose') +
            prior(normal(0,0.2), class = b, nlpar = 'noise', coef = 'dose') +
            prior(normal(0,0.2), class = b, nlpar = 'omega1', coef = 'dose') +
            prior(normal(0,0.2), class = b, nlpar = 'omega2', coef = 'dose')



priors_rho <- prior(normal(0, 0.5), nlpar = 'logrho') +
            prior(normal(3, 1), nlpar = 'omega1') +
            prior(normal(0, 1), nlpar = 'omega2') +
            prior(normal(0, 3), nlpar = 'noise') +
            prior(normal(0,0.35), lb=0, class = sd, group = 'subjid', nlpar = 'omega1') +
            prior(normal(0,0.35), lb=0, class = sd, group = 'subjid', nlpar = 'omega2') +
            prior(normal(0,0.35), lb=0, class = sd, group = 'subjid', nlpar = 'noise') +
            prior(normal(0,0.35), lb=0, class = sd, group = 'subjid', nlpar = 'logrho') +
            prior(normal(0,.5), class = b, nlpar = 'logrho', coef = 'dose') +
            prior(normal(0,1), class = b, nlpar = 'noise', coef = 'dose')

priors_omega <- prior(normal(0, 0.5), nlpar = 'logrho') +
            prior(normal(3, 2), nlpar = 'omega1') +
            prior(normal(0, 1), nlpar = 'omega2') +
            prior(normal(0, 3), nlpar = 'noise') +
            prior(normal(0,0.35), lb=0, class = sd, group = 'subjid', nlpar = 'omega1') +
            prior(normal(0,0.35), lb=0, class = sd, group = 'subjid', nlpar = 'omega2') +
            prior(normal(0,0.35), lb=0, class = sd, group = 'subjid', nlpar = 'noise') +
            prior(normal(0,0.35), lb=0, class = sd, group = 'subjid', nlpar = 'logrho') +
            prior(normal(0,1), class = b, nlpar = 'noise', coef = 'dose') +
            prior(normal(0,1), class = b, nlpar = 'omega1', coef = 'dose') +
            prior(normal(0,1), class = b, nlpar = 'omega2', coef = 'dose')
