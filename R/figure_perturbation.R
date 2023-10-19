## fit the bilateral perturbation data with brms and do the plotting
# Chaofei 27/3/2023

# library(brms)
# library(dplyr)
# library(bayesplot)
# library(stringi)
# library(ggplot2)
# library(patchwork)
# library(rprojroot)
# library(modelr)
# library(tidyverse)
# library(tidybayes)
# library(loo)
# library(GGally)
# library(stringr)

# set whether you want to refit the model
# if you want to refit the models yourself, set "do_the_fit_this_time" to TRUE,
# then it will take 12 min to fit all these models in a maching with Intel(R) Core(TM) i7-12700K,
# Or you can skip the fitting and load it from the file
do_the_fit_this_time = FALSE

# load the dataset
git_root <- find_root(is_git_root)
source(file.path(git_root, "R", "brms_functions.R"))

# bilateral opto silencing data
bidf_opto <- get_bi_opto_data(git_root)
bino_bidf_opto <- binomialize_brms(bidf_opto)

#unilateral opto silencing data
unidf_opto <- get_uni_opto_data(git_root)
bino_unidf_opto <- binomialize_brms(unidf_opto)

# bilateral fof muscimol inhibition data
bidf_fof_muscimol <- get_bifof_muscimol_data(git_root) %>% drop_na(dose)
bino_bidf_fof_m <- binomialize_brms(bidf_fof_muscimol)

# unilateral fof muscimol inhibition data
unidf_fof_muscimol <- get_unifof_muscimol_data(git_root) %>% drop_na(dose)
bino_unidf_fof_m <- binomialize_brms(unidf_fof_muscimol)

# bilateral fof muscimol inhibition data only first 40 trials
bidf_fof_muscimol_first40 <- get_bifof_muscimol_data(git_root) %>% drop_na(dose) %>% filter(trialnum<=40)
bino_bidf_fof_m_40 <- binomialize_brms(bidf_fof_muscimol_first40)

# unilateral fof muscimol inhibition data only first 40 trials
unidf_fof_muscimol_first40 <- get_unifof_muscimol_data(git_root) %>% drop_na(dose) %>% filter(trialnum<=40)
bino_unidf_fof_m_40 <- binomialize_brms(unidf_fof_muscimol_first40)

# bilateral ppc muscimol inhibition data
bidf_ppc_muscimol <- get_bippc_muscimol_data(git_root) %>% drop_na(dose)
bino_bidf_ppc_m <- binomialize_brms(bidf_ppc_muscimol)

# unilateral ppc muscimol inhibition data
unidf_ppc_muscimol <- get_unippc_muscimol_data(git_root) %>% drop_na(dose)
bino_unidf_ppc_m <- binomialize_brms(unidf_ppc_muscimol)

# bilateral ppc muscimol inhibition data only first 40 trials
bidf_ppc_muscimol_first40 <- get_bippc_muscimol_data(git_root) %>% drop_na(dose) %>% filter(trialnum<=40)
bino_bidf_ppc_m_40 <- binomialize_brms(bidf_ppc_muscimol_first40)

# unilateral ppc muscimol inhibition data only first 40 trials
unidf_ppc_muscimol_first40 <- get_unippc_muscimol_data(git_root) %>% drop_na(dose) %>% filter(trialnum<=40)
bino_unidf_ppc_m_40 <- binomialize_brms(unidf_ppc_muscimol_first40)

ephys <- get_ephys_data() %>% drop_na(dose)
bino_ephys <- binomialize_brms(ephys)

# fit the opto data with brms
if (do_the_fit_this_time) {
  iter = 13000
  warmup = 8000
  seed_opto = 3
  seed_muscimol = 33
  chains = 6

  bi_opto_cdf <- brm(f_drug_phi,
                     data = bino_bidf_opto,
                     family = binomial(link = "identity"),
                     prior = priors, iter = iter, warmup = warmup,
                     chains = chains, cores = chains, seed = seed_opto,
                     sample_prior = "yes", save_pars = save_pars(all = TRUE)
  )

  uni_opto_cdf <- brm(f_drug_phi,
                     data = bino_unidf_opto,
                     family = binomial(link = "identity"),
                     prior = priors, iter = iter, warmup = warmup,
                     chains = chains, cores = chains, seed = seed_opto,
                     sample_prior = "yes", save_pars = save_pars(all = TRUE)
  )

  # fit the muscimol data with brms
  bi_fof_muscimol_cdf <- brm(f_drug_phi,
                             data = bino_bidf_fof_m,
                             family = binomial(link = "identity"),
                             prior = priors, iter = iter, warmup = warmup,
                             chains = chains, cores = chains, seed = seed_muscimol,
                             sample_prior = "yes", save_pars = save_pars(all = TRUE)
  )

  uni_fof_muscimol_cdf <- brm(f_drug_phi,
                             data = bino_unidf_fof_m,
                             family = binomial(link = "identity"),
                             prior = priors, iter = iter, warmup = warmup,
                             chains = chains, cores = chains, seed = seed_muscimol,
                             sample_prior = "yes", save_pars = save_pars(all = TRUE)
  )

  ## fit the fof first 40 trials
  bi_fof_muscimol_cdf_40 <- brm(f_drug_phi,
                             data = bino_bidf_fof_m_40,
                             family = binomial(link = "identity"),
                             prior = priors, iter = iter, warmup = warmup,
                             chains = chains, cores = chains, seed = seed_muscimol,
                             sample_prior = "yes", save_pars = save_pars(all = TRUE)
  )

  uni_fof_muscimol_cdf_40 <- brm(f_drug_phi,
                              data = bino_unidf_fof_m_40,
                              family = binomial(link = "identity"),
                              prior = priors, iter = iter, warmup = warmup,
                              chains = chains, cores = chains, seed = seed_muscimol,
                              sample_prior = "yes", save_pars = save_pars(all = TRUE)
  )

  ## fit the ppc data
  bi_ppc_muscimol_cdf <- brm(f_drug_phi,
                             data = bino_bidf_ppc_m,
                             family = binomial(link = "identity"),
                             prior = priors, iter = iter, warmup = warmup,
                             chains = chains, cores = chains, seed = seed_muscimol,
                             sample_prior = "yes", save_pars = save_pars(all = TRUE)
  )

  uni_ppc_muscimol_cdf <- brm(f_drug_phi,
                              data = bino_unidf_ppc_m,
                              family = binomial(link = "identity"),
                              prior = priors, iter = iter, warmup = warmup,
                              chains = chains, cores = chains, seed = seed_muscimol,
                              sample_prior = "yes", save_pars = save_pars(all = TRUE)
  )


  ## fit the ppc first 40 trials
  bi_ppc_muscimol_cdf_40 <- brm(f_drug_phi,
                             data = bino_bidf_ppc_m_40,
                             family = binomial(link = "identity"),
                             prior = priors, iter = iter, warmup = warmup,
                             chains = chains, cores = chains, seed = seed_muscimol,
                             sample_prior = "yes", save_pars = save_pars(all = TRUE)
  )

  uni_ppc_muscimol_cdf_40 <- brm(f_drug_phi,
                              data = bino_unidf_ppc_m_40,
                              family = binomial(link = "identity"),
                              prior = priors, iter = iter, warmup = warmup,
                              chains = chains, cores = chains, seed = seed_muscimol,
                              sample_prior = "yes", save_pars = save_pars(all = TRUE)
  )

  ephys_cdf <- brm(f_drug_phi,
                     data = bino_ephys,
                     family = binomial(link = "identity"),
                     prior = priors, iter = iter, warmup = warmup,
                     chains = chains, cores = chains, seed = seed_opto,
                     sample_prior = "yes", save_pars = save_pars(all = TRUE)
  )

  # fit the synthetic data
  agent_seed = 1234
  set.seed(agent_seed)
  agents <- make_agents(20)
  rho_syntheticdf <- drug_sim(agents, deltaRho=-0.50, deltaOmega1=0, deltaOmega2=0, deltaBeta=0)
  omega_syntheticdf <- drug_sim(agents, deltaRho=0, deltaOmega1=-1, deltaOmega2=-3, deltaBeta=0)

  bino_rho_synthetic <- binomialize_brms(rho_syntheticdf)
  bino_omega_synthetic <- binomialize_brms(omega_syntheticdf)

  # fit the synthetic data using brms model
  rho_synthetic_cdf <- brm(f_drug_phi,
                           data = bino_rho_synthetic,
                           family = binomial(link = "identity"),
                           prior = priors, iter = 8000, warmup = 6000,
                           chains = 6, cores = 6, seed = 231,
                           sample_prior = "yes", save_pars = save_pars(all = TRUE)
  )
  #rho_synthetic_cdf <- add_criterion(rho_synthetic_cdf, c("loo"), moment_match=TRUE, pointwise=TRUE)

  omega_synthetic_cdf <- brm(f_drug_phi,
                             data = bino_omega_synthetic,
                             family = binomial(link = "identity"),
                             prior = priors, iter = 8000, warmup = 6000,
                             chains = 6, cores = 6, seed = 231,
                             sample_prior = "yes", save_pars = save_pars(all = TRUE)
  )
  #omega_synthetic_cdf <- add_criterion(omega_synthetic_cdf, c("loo"), moment_match=TRUE, pointwise=TRUE)
}else{
  load(file.path(git_root,"brms_model_fits.RData"))
}

# find the most frequent task parameter given the df
get_freq_task_param = function(df, param = 'sb_mag'){
  return(as.numeric(names(sort(table(df[param]), decreasing = TRUE)[1])))
}

bino = function(x){
  out = binom.test(sum(x), length(x))
  df = data.frame(y = mean(x), ymin=out$conf.int[1], ymax=out$conf.int[2])
  return(df)
}

inv_logit <- function(x) {
  return(1 / (1 + exp(-x)))
}

normal_cdf <- pnorm

colfunc <- colorRampPalette(c("purple4", "white"))
colfunc <- colorRampPalette(c("#D766C9", "white"))


## plot the sample subject data
plot_pred_sample_individual <- function(this_subjid, df, model, x_lab_on = TRUE, y_lab_on = TRUE, legend_on = TRUE){
  df1 <- df %>% filter(subjid == this_subjid)
  bino_df <- binomialize_brms(df1)
  bino_df_pred <- data.frame(df1) %>%
    data_grid(subjid = subjid, sb_norm = get_freq_task_param(df, 'sb_norm'),
              lottery_norm = seq_range(c(-0.5, 1.1), by=0.005), dose, lottery_prob = get_freq_task_param(df, 'lottery_prob'))

  bino_df_pred$total = 100
  bino_df_pred <- bino_df_pred %>%
    add_epred_draws(model, re_formula = NULL)

  p <- bino_df_pred %>%
    ggplot(aes(x = (lottery_norm*lottery_prob-sb_norm), y=.epred/100, color = factor(dose, levels = c(0, 1), labels = c('control', 'bi-fof')), fill = factor(dose, levels = c(0, 1), labels = c('control', 'bi-fof')))) +
    theme_classic(base_size = 10) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed',  alpha = 0.4) +
    stat_lineribbon(aes(y= .epred/100, color = factor(dose, levels = c(0, 1), labels = c('control', 'bi-fof'))), .width = c(.8, .5), alpha = 1/4) +
    stat_summary_bin(mapping = aes(x = (lottery_norm*lottery_prob-sb_norm), y = chose_lottery, color = factor(dose, levels = c(0, 1), labels = c('control', 'bi-fof'))), data = df1, fun.data = bino,
                     geom = 'pointrange', bins = 7, position = position_jitterdodge(dodge.width = 0.01)) +
    scale_x_continuous(breaks=c(0,0.4)) +
    scale_y_continuous(limits = c(-0.05, 1.05), breaks=c(0,0.5,1)) +
    scale_fill_brewer() +
    annotate("text", label = this_subjid, x = .2, y = .1, size = 4) +
    scale_color_manual(values = c('control' = 'azure4','bi-fof' = 'purple3')) +
    scale_fill_manual(values = c('control' = 'azure4','bi-fof' = 'purple3')) +
    theme(plot.margin = margin(5,5,5,5))
  if(x_lab_on){
    p <- p +labs(x = bquote(EV[lottery]-EV[surebet]))
  }else{
    p <- p + theme(axis.title.x = element_blank())
  }
  if(!legend_on){
    p <- p + theme(legend.position = "none")
  }
  if(y_lab_on){
    p <- p + labs(y = 'P(Chose lottery)', color='', fill='')
  }else{
    p <- p + theme(axis.title.y = element_blank())
  }
  return(p)
}

## plot the individual data
plot_pred_individual <- function(df, model){
  df_name <- deparse(substitute(df))
  if (str_detect(df_name, 'bi')){
    this_labels = c('Control', 'Bi-fof silencing')
    clr_values = c('Control' = 'azure4','Bi-fof silencing' = 'purple3')
  }else{
    this_labels = c('Control', 'Uni-fof silencing')
    clr_values = c('Control' = 'azure4','Uni-fof silencing' = '#D766C9')
  }
  bino_df <- binomialize_brms(df)
  subjids = unique(bino_df$subjid)
  bino_df_pred <- data.frame()
  for (i in 1:length(subjids)){
    df_grid <- data.frame(df) %>%
      filter(subjid == subjids[i]) %>%
      data_grid(subjid = subjid, sb_norm = get_freq_task_param(df %>% filter(subjid == subjids[i]), 'sb_norm'),
                lottery_norm = seq_range(c(-0.5, 1.1), by=0.005), dose, lottery_prob = get_freq_task_param(df %>% filter(subjid == subjids[i]), 'lottery_prob'))
    bino_df_pred = rbind(bino_df_pred, df_grid)
  }

  bino_df_pred$total = 100
  bino_df_pred <- bino_df_pred %>%
    add_epred_draws(model, re_formula = NULL)

  p <- bino_df_pred %>%
    ggplot(aes(x = (lottery_norm*lottery_prob-sb_norm), y=.epred/100, color = factor(dose, levels = c(0, 1), labels = this_labels), fill = factor(dose, levels = c(0, 1), labels = this_labels))) +
    theme_classic(base_size = 15) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed',  alpha = 0.4) +
    stat_lineribbon(aes(y= .epred/100, color = factor(dose, levels = c(0, 1), labels = this_labels)), .width = c(.8, .5), alpha = 1/4) +
    stat_summary_bin(mapping = aes(x = (lottery_norm*lottery_prob-sb_norm), y = chose_lottery, color = factor(dose, levels = c(0, 1), labels = this_labels)), data = df, fun.data = bino,
                     geom = 'pointrange', bins = 7, position = position_jitterdodge(dodge.width = 0.01)) +
    scale_x_continuous(breaks=c(0,0.4)) +
    scale_y_continuous(breaks=c(0,0.5,1)) +
    scale_fill_brewer() +
    labs(y = 'P(Chose lottery)', x = bquote(EV[lottery]-EV[surebet]), color='', fill='') +
    scale_color_manual(values = clr_values) +
    scale_fill_manual(values = clr_values) +
    facet_wrap(.~subjid)

  if (str_detect(df_name, 'opto')){
    p <- p + theme(legend.position="bottom")
  }else{
    p <- p + theme(legend.position="none")
  }
  return(p)
}


## plot the opto and muscimol posterior together
plot_delta_rho_psd <- function(bi_opto_model, uni_opto_model, bi_muscimol_model, uni_muscimol_model, y_lab = TRUE){
  color_scheme_set(scheme = rev(colfunc(10)[1:6]))
  ps_bi_opto <- posterior_samples(bi_opto_model, pars = '^b', as.matrix = T)
  ps_bi_opto <- as.data.frame(ps_bi_opto) %>% mutate(delta_rho_bi_opto = exp(b_logrho_Intercept + b_logrho_dose) - exp(b_logrho_Intercept)) %>% select(c('delta_rho_bi_opto'))
  ps_uni_opto <- posterior_samples(uni_opto_model, pars = '^b', as.matrix = T)
  ps_uni_opto <- as.data.frame(ps_uni_opto) %>% mutate(delta_rho_uni_opto = exp(b_logrho_Intercept + b_logrho_dose) - exp(b_logrho_Intercept)) %>% select(c('delta_rho_uni_opto'))
  ps_bi_muscimol <- posterior_samples(bi_muscimol_model, pars = '^b', as.matrix = T)
  ps_bi_muscimol <- as.data.frame(ps_bi_muscimol) %>% mutate(delta_rho_bi_muscimol = exp(b_logrho_Intercept + b_logrho_dose) - exp(b_logrho_Intercept)) %>% select(c('delta_rho_bi_muscimol'))
  ps_uni_muscimol <- posterior_samples(uni_muscimol_model, pars = '^b', as.matrix = T)
  ps_uni_muscimol <- as.data.frame(ps_uni_muscimol) %>% mutate(delta_rho_uni_muscimol = exp(b_logrho_Intercept + b_logrho_dose) - exp(b_logrho_Intercept)) %>% select(c('delta_rho_uni_muscimol'))

  ps_delta_rho <- cbind(ps_bi_opto, ps_uni_opto, ps_bi_muscimol, ps_uni_muscimol) %>% as.mcmc()
  p <- mcmc_areas(ps_delta_rho, prob=.7, prob_outer = 0.9999, point_est = 'mean', area_method = 'equal height') +
    vline_0(size = .2, color = 'gray') +
    theme_classic(base_size = 10) +
    scale_x_continuous(limits=c(-0.6, 0.2), breaks = c(-0.4, 0), labels = c('-0.4', '0')) +
    scale_y_discrete(limits = c('delta_rho_uni_muscimol','delta_rho_bi_muscimol','delta_rho_uni_opto','delta_rho_bi_opto'),
                     labels = c(
                       'delta_rho_bi_opto' = 'Bi-Opto',
                       'delta_rho_uni_opto' = 'Uni-Opto',
                       'delta_rho_bi_muscimol'= 'Bi-Muscimol',
                       'delta_rho_uni_muscimol'= 'Uni-Muscimol'),
                     expand = expansion(add = c(.3, .3))) +
    xlab(expression(paste(Delta, rho))) +
    theme(axis.text.y = element_text(angle = 60))

  if (!y_lab){
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}

plot_delta_noise_psd <- function(bi_opto_model, uni_opto_model, bi_muscimol_model, uni_muscimol_model, y_lab = TRUE){
  ps_bi_opto <- posterior_samples(bi_opto_model, pars = '^b', as.matrix = T)
  ps_bi_opto <- as.data.frame(ps_bi_opto) %>% mutate(delta_noise_bi_opto = exp(b_noise_Intercept + b_noise_dose) - exp(b_noise_Intercept)) %>% select(c('delta_noise_bi_opto'))
  ps_uni_opto <- posterior_samples(uni_opto_model, pars = '^b', as.matrix = T)
  ps_uni_opto <- as.data.frame(ps_uni_opto) %>% mutate(delta_noise_uni_opto = exp(b_noise_Intercept + b_noise_dose) - exp(b_noise_Intercept)) %>% select(c('delta_noise_uni_opto'))
  ps_bi_muscimol <- posterior_samples(bi_muscimol_model, pars = '^b', as.matrix = T)
  ps_bi_muscimol <- as.data.frame(ps_bi_muscimol) %>% mutate(delta_noise_bi_muscimol = exp(b_noise_Intercept + b_noise_dose) - exp(b_noise_Intercept)) %>% select(c('delta_noise_bi_muscimol'))
  ps_uni_muscimol <- posterior_samples(uni_muscimol_model, pars = '^b', as.matrix = T)
  ps_uni_muscimol <- as.data.frame(ps_uni_muscimol) %>% mutate(delta_noise_uni_muscimol = exp(b_noise_Intercept + b_noise_dose) - exp(b_noise_Intercept)) %>% select(c('delta_noise_uni_muscimol'))


  ps_delta_rho <- cbind(ps_bi_opto, ps_uni_opto, ps_bi_muscimol, ps_uni_muscimol) %>% as.mcmc()
  p <- mcmc_areas(ps_delta_rho, prob=.7, prob_outer = 0.9999, point_est = 'mean', area_method = 'equal height') +
    vline_0(size = .2, color = 'gray') +
    theme_classic(base_size = 10) +
    scale_x_continuous(limits=c(-0.05, 0.15), breaks = c(-0.05, 0, 0.1), labels = c('-0.05', '0', '0.1')) +
    scale_y_discrete(limits = c('delta_noise_uni_muscimol','delta_noise_bi_muscimol','delta_noise_uni_opto','delta_noise_bi_opto'),
                     labels = c(
                       'delta_noise_bi_opto' = 'Bi-Opto',
                       'delta_noise_uni_opto' = 'Uni-Opto',
                       'delta_noise_bi_muscimol'= 'Bi-Muscimol',
                       'delta_noise_uni_muscimol'= 'Uni-Muscimol'),
                     expand = expansion(add = c(.3, .3))) +
    xlab(expression(paste(Delta, sigma)))
  if (!y_lab){
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}

plot_delta_omega_rational_psd <- function(bi_opto_model, uni_opto_model, bi_muscimol_model, uni_muscimol_model, y_lab = TRUE){
  ps_bi_opto <- posterior_samples(bi_opto_model, pars = '^b', as.matrix = T)
  ps_bi_opto <- as.data.frame(ps_bi_opto) %>% mutate(delta_omega_rational_bi_opto = inv_logit(b_omega1_Intercept + b_omega1_dose) - inv_logit(b_omega1_Intercept)) %>% select(c('delta_omega_rational_bi_opto'))
  ps_uni_opto <- posterior_samples(uni_opto_model, pars = '^b', as.matrix = T)
  ps_uni_opto <- as.data.frame(ps_uni_opto) %>% mutate(delta_omega_rational_uni_opto = inv_logit(b_omega1_Intercept + b_omega1_dose) - inv_logit(b_omega1_Intercept)) %>% select(c('delta_omega_rational_uni_opto'))
  ps_bi_muscimol <- posterior_samples(bi_muscimol_model, pars = '^b', as.matrix = T)
  ps_bi_muscimol <- as.data.frame(ps_bi_muscimol) %>% mutate(delta_omega_rational_bi_muscimol = inv_logit(b_omega1_Intercept + b_omega1_dose) - inv_logit(b_omega1_Intercept)) %>% select(c('delta_omega_rational_bi_muscimol'))
  ps_uni_muscimol <- posterior_samples(uni_muscimol_model, pars = '^b', as.matrix = T)
  ps_uni_muscimol <- as.data.frame(ps_uni_muscimol) %>% mutate(delta_omega_rational_uni_muscimol = inv_logit(b_omega1_Intercept + b_omega1_dose) - inv_logit(b_omega1_Intercept)) %>% select(c('delta_omega_rational_uni_muscimol'))

  ps_delta_rho <- cbind(ps_bi_opto, ps_uni_opto, ps_bi_muscimol, ps_uni_muscimol) %>% as.mcmc()
  p <- mcmc_areas(ps_delta_rho, prob=.7, prob_outer = 0.9999, point_est = 'mean', area_method = 'equal height') +
    vline_0(size = .2, color = 'gray') +
    theme_classic(base_size = 10) +
    scale_x_continuous(limits=c(-0.55, 0.3), breaks = c(-0.4, 0, 0.2), labels = c('-0.4', '0', '0.2')) +
    scale_y_discrete(limits = c('delta_omega_rational_uni_muscimol','delta_omega_rational_bi_muscimol','delta_omega_rational_uni_opto','delta_omega_rational_bi_opto'),
                     labels = c(
                       'delta_omega_rational_bi_opto' = 'Bi-Opto',
                       'delta_omega_rational_uni_opto' = 'Uni-Opto',
                       'delta_omega_rational_bi_muscimol'= 'Bi-Muscimol',
                       'delta_omega_rational_uni_muscimol'= 'Uni-Muscimol'),
                     expand = expansion(add = c(.3, .3))) +
    xlab(expression(paste(Delta, omega['rational'])))

  if (!y_lab){
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}

plot_delta_omega_lottery_psd <- function(bi_opto_model, uni_opto_model, bi_muscimol_model, uni_muscimol_model, y_lab = TRUE){
  ps_bi_opto <- posterior_samples(bi_opto_model, pars = '^b', as.matrix = T)
  ps_bi_opto <- as.data.frame(ps_bi_opto) %>% mutate(delta_omega_lottery_bi_opto = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*inv_logit(b_omega2_Intercept + b_omega2_dose) -
                                                       (1- inv_logit(b_omega1_Intercept))*inv_logit(b_omega2_Intercept)) %>% select(c('delta_omega_lottery_bi_opto'))
  ps_uni_opto <- posterior_samples(uni_opto_model, pars = '^b', as.matrix = T)
  ps_uni_opto <- as.data.frame(ps_uni_opto) %>% mutate(delta_omega_lottery_uni_opto = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*inv_logit(b_omega2_Intercept + b_omega2_dose) -
                                                         (1- inv_logit(b_omega1_Intercept))*inv_logit(b_omega2_Intercept)) %>% select(c('delta_omega_lottery_uni_opto'))
  ps_bi_muscimol <- posterior_samples(bi_muscimol_model, pars = '^b', as.matrix = T)
  ps_bi_muscimol <- as.data.frame(ps_bi_muscimol) %>% mutate(delta_omega_lottery_bi_muscimol = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*inv_logit(b_omega2_Intercept + b_omega2_dose) -
                                                               (1- inv_logit(b_omega1_Intercept))*inv_logit(b_omega2_Intercept)) %>% select(c('delta_omega_lottery_bi_muscimol'))
  ps_uni_muscimol <- posterior_samples(uni_muscimol_model, pars = '^b', as.matrix = T)
  ps_uni_muscimol <- as.data.frame(ps_uni_muscimol) %>% mutate(delta_omega_lottery_uni_muscimol = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*inv_logit(b_omega2_Intercept + b_omega2_dose) -
                                                                 (1- inv_logit(b_omega1_Intercept))*inv_logit(b_omega2_Intercept)) %>% select(c('delta_omega_lottery_uni_muscimol'))

  ps_delta_rho <- cbind(ps_bi_opto, ps_uni_opto, ps_bi_muscimol, ps_uni_muscimol) %>% as.mcmc()
  p <- mcmc_areas(ps_delta_rho, prob=.7, prob_outer = 0.9999, point_est = 'mean', area_method = 'equal height') +
    vline_0(size = .2, color = 'gray') +
    theme_classic(base_size = 10) +
    scale_x_continuous(limits=c(-0.2, 0.1), breaks = c(-0.2, 0, 0.1), labels = c('-0.2', '0', '0.1')) +
    scale_y_discrete(limits = c('delta_omega_lottery_uni_muscimol','delta_omega_lottery_bi_muscimol','delta_omega_lottery_uni_opto','delta_omega_lottery_bi_opto'),
                     labels = c(
                       'delta_omega_lottery_bi_opto' = 'Bi-Opto',
                       'delta_omega_lottery_uni_opto' = 'Uni-Opto',
                       'delta_omega_lottery_bi_muscimol'= 'Bi-Muscimol',
                       'delta_omega_lottery_uni_muscimol'= 'Uni-Muscimol'),
                     expand = expansion(add = c(.3, .3))) +
    xlab(expression(paste(Delta, omega['lottery'])))

  if (!y_lab){
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}


plot_delta_omega_surebet_psd <- function(bi_opto_model, uni_opto_model, bi_muscimol_model, uni_muscimol_model, y_lab = TRUE){
  ps_bi_opto <- posterior_samples(bi_opto_model, pars = '^b', as.matrix = T)
  ps_bi_opto <- as.data.frame(ps_bi_opto) %>% mutate(delta_omega_surebet_bi_opto = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*(1-inv_logit(b_omega2_Intercept + b_omega2_dose)) -
                                                       (1- inv_logit(b_omega1_Intercept))*(1-inv_logit(b_omega2_Intercept))) %>% select(c('delta_omega_surebet_bi_opto'))
  ps_uni_opto <- posterior_samples(uni_opto_model, pars = '^b', as.matrix = T)
  ps_uni_opto <- as.data.frame(ps_uni_opto) %>% mutate(delta_omega_surebet_uni_opto = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*(1-inv_logit(b_omega2_Intercept + b_omega2_dose)) -
                                                         (1- inv_logit(b_omega1_Intercept))*(1-inv_logit(b_omega2_Intercept))) %>% select(c('delta_omega_surebet_uni_opto'))
  ps_bi_muscimol <- posterior_samples(bi_muscimol_model, pars = '^b', as.matrix = T)
  ps_bi_muscimol <- as.data.frame(ps_bi_muscimol) %>% mutate(delta_omega_surebet_bi_muscimol = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*(1-inv_logit(b_omega2_Intercept + b_omega2_dose)) -
                                                               (1- inv_logit(b_omega1_Intercept))*(1-inv_logit(b_omega2_Intercept))) %>% select(c('delta_omega_surebet_bi_muscimol'))
  ps_uni_muscimol <- posterior_samples(uni_muscimol_model, pars = '^b', as.matrix = T)
  ps_uni_muscimol <- as.data.frame(ps_uni_muscimol) %>% mutate(delta_omega_surebet_uni_muscimol = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*(1-inv_logit(b_omega2_Intercept + b_omega2_dose)) -
                                                                 (1- inv_logit(b_omega1_Intercept))*(1-inv_logit(b_omega2_Intercept))) %>% select(c('delta_omega_surebet_uni_muscimol'))

  ps_delta_rho <- cbind(ps_bi_opto, ps_uni_opto, ps_bi_muscimol, ps_uni_muscimol) %>% as.mcmc()
  p <- mcmc_areas(ps_delta_rho, prob=.7, prob_outer = 0.9999, point_est = 'mean', area_method = 'equal height') +
    vline_0(size = .2, color = 'gray') +
    theme_classic(base_size = 10) +
    scale_x_continuous(limits=c(-0.12, 0.6), breaks = c(0, 0.3, 0.6), labels = c('0', '0.3', '0.6')) +
    scale_y_discrete(limits = c('delta_omega_surebet_uni_muscimol','delta_omega_surebet_bi_muscimol','delta_omega_surebet_uni_opto','delta_omega_surebet_bi_opto'),
                     labels = c(
                       'delta_omega_surebet_bi_opto' = 'Bi-Opto',
                       'delta_omega_surebet_uni_opto' = 'Uni-Opto',
                       'delta_omega_surebet_bi_muscimol'= 'Bi-Muscimol',
                       'delta_omega_surebet_uni_muscimol'= 'Uni-Muscimol'),
                     expand = expansion(add = c(.3, .3))) +
    xlab(expression(paste(Delta, omega['surebet'])))

  if (!y_lab){
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}








## plot the individual postiers draw from the brms
get_population_pars_from_brms <- function(model){
  model_name <- deparse(substitute(model))

  pars <- model %>% spread_draws(b_logrho_Intercept, b_logrho_dose, r_subjid__logrho[subjid,],
                                 b_noise_Intercept, b_noise_dose, r_subjid__noise[subjid,],
                                 b_omega1_Intercept, b_omega1_dose, r_subjid__omega1[subjid,],
                                 b_omega2_Intercept, b_omega2_dose, r_subjid__omega2[subjid,]) %>%
    mutate(logrho = b_logrho_Intercept + r_subjid__logrho, logrho_treat = b_logrho_Intercept + r_subjid__logrho + b_logrho_dose,
           lognoise = b_noise_Intercept + r_subjid__noise, lognoise_treat = b_noise_Intercept + r_subjid__noise + b_noise_dose,
           omega1 = b_omega1_Intercept + r_subjid__omega1, omega1_treat = b_omega1_Intercept + r_subjid__omega1 + b_omega1_dose,
           omega2 = b_omega2_Intercept + r_subjid__omega2, omega2_treat = b_omega2_Intercept + r_subjid__omega2 + b_omega2_dose,
           rho_control = exp(logrho), rho_treat = exp(logrho_treat),
           noise_control = exp(lognoise), noise_treat = exp(lognoise_treat),
           omega_rational_control = inv_logit(omega1), omega_rational_treat = inv_logit(omega1_treat),
           omega_lottery_control = (1-inv_logit(omega1))*inv_logit(omega2), omega_lottery_treat = (1-inv_logit(omega1_treat))*inv_logit(omega2_treat),
           omega_surebet_control = (1-inv_logit(omega1))*(1-inv_logit(omega2)), omega_surebet_treat = (1-inv_logit(omega1_treat))*(1-inv_logit(omega2_treat)))
  return(list(pars = pars, model_name = model_name))
}

plot_individual_draw_dens <- function(df, dataset, xlab_on = TRUE){

  if (str_detect(dataset, 'bi')){
    this_labels = c('Control', 'Bi-fof silencing')
    clr_values = c('Control' = 'azure4','Bi-fof silencing' = 'purple3')
  }else{
    this_labels = c('Control', 'Uni-fof silencing')
    clr_values = c('Control' = 'azure4','Uni-fof silencing' = '#D766C9')
  }
  p1 <- ggplot(df %>% select(c('subjid', 'rho_control','rho_treat')) %>% gather(type, value, 2:3)) +
    theme_classic(base_size = 14) +
    geom_density(aes(x = value, y = after_stat(scaled), fill = factor(type, levels = c('rho_control', 'rho_treat'), labels = this_labels),
                     color = factor(type, levels = c('rho_control', 'rho_treat'), labels = this_labels)), alpha = 0.5, linewidth = 0.5) +
    scale_color_manual(values = clr_values) +
    scale_fill_manual(values = clr_values) +
    scale_x_continuous(limits=c(0.2, 1.8), breaks = c(0.2,1.0,1.8)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    theme(legend.position = 'none',
          axis.title.y = element_blank())

  p2 <- ggplot(df %>% select(c('subjid', 'noise_control','noise_treat')) %>% gather(type, value, 2:3)) +
    theme_classic(base_size = 14) +
    geom_density(aes(x = value, y = after_stat(scaled), fill = factor(type, levels = c('noise_control', 'noise_treat'), labels = this_labels),
                     color = factor(type, levels = c('noise_control', 'noise_treat'), labels = this_labels)), alpha = 0.5, linewidth = 0.5) +
    scale_color_manual(values = clr_values) +
    scale_fill_manual(values = clr_values) +
    scale_x_continuous(limits=c(0, 0.2), breaks = c(0,0.1,0.2)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    theme(legend.position = 'none',
          axis.title.y = element_blank())

  p3 <- ggplot(df %>% select(c('subjid', 'omega_rational_control','omega_rational_treat')) %>% gather(type, value, 2:3)) +
    theme_classic(base_size = 14) +
    geom_density(aes(x = value, y = after_stat(scaled), fill = factor(type, levels = c('omega_rational_control', 'omega_rational_treat'), labels = this_labels),
                     color = factor(type, levels = c('omega_rational_control', 'omega_rational_treat'), labels = this_labels)), alpha = 0.5, linewidth = 0.5) +
    scale_color_manual(values = clr_values) +
    scale_fill_manual(values = clr_values) +
    scale_x_continuous(limits=c(0.75, 1), breaks = c(0.8,0.9,1.0)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    theme(legend.position = 'none',
          axis.title.y = element_blank())


  p4 <- ggplot(df %>% select(c('subjid', 'omega_lottery_control','omega_lottery_treat')) %>% gather(type, value, 2:3)) +
    theme_classic(base_size = 14) +
    geom_density(aes(x = value, y = after_stat(scaled), fill = factor(type, levels = c('omega_lottery_control', 'omega_lottery_treat'), labels = this_labels),
                     color = factor(type, levels = c('omega_lottery_control', 'omega_lottery_treat'), labels = this_labels)), alpha = 0.5, linewidth = 0.5) +
    scale_color_manual(values = clr_values) +
    scale_fill_manual(values = clr_values) +
    scale_x_continuous(limits=c(0, .15), breaks = c(0,.1)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    theme(legend.position = 'none',
          axis.title.y = element_blank())

  p5 <- ggplot(df %>% select(c('subjid', 'omega_surebet_control','omega_surebet_treat')) %>% gather(type, value, 2:3)) +
    theme_classic(base_size = 14) +
    geom_density(aes(x = value, y = after_stat(scaled), fill = factor(type, levels = c('omega_surebet_control', 'omega_surebet_treat'), labels = this_labels),
                     color = factor(type, levels = c('omega_surebet_control', 'omega_surebet_treat'), labels = this_labels)), alpha = 0.5, linewidth = 0.5) +
    scale_color_manual(values = clr_values) +
    scale_fill_manual(values = clr_values) +
    scale_x_continuous(limits=c(0, .15), breaks = c(0,.1)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    theme(legend.position = 'none',
          axis.title.y = element_blank())

  if (str_detect(dataset, 'bi')){
    if (unique(df$subjid) %in% c(2152,2153,2154,2155,2156,2160,2165,2166)){
      p1 <- p1 + scale_x_continuous(limits=c(0.2,1), breaks = c(0.2,0.6,1.0))
      p3 <- p3 + scale_x_continuous(limits=c(0.2,1), breaks = c(0.2,0.6,1.0))
    }
    if (unique(df$subjid) %in% c(2152, 2153, 2156, 2160)){
      p4 <- p4 + scale_x_continuous(limits=c(0, .3), breaks = c(0,.15,.3))
    }
    if (unique(df$subjid) %in% c(2152, 2153, 2160)){
      p5 <- p5 + scale_x_continuous(limits=c(0, .1), breaks = c(0,.05,.1))
    }
    if (unique(df$subjid) %in% c(2154,2155)){
      p5 <- p5 + scale_x_continuous(limits=c(0, .5), breaks = c(0,.25,.5))
    }
    if (unique(df$subjid) %in% c(2156)){
      p5 <- p5 + scale_x_continuous(limits=c(0, .22), breaks = c(0,.1,.2))
    }
  }else{
    if (unique(df$subjid) %in% c(2172,2176,2177,2180,2182,2225,2228,2240)){
      p3 <- p3 + scale_x_continuous(limits=c(0.85,1), breaks = c(0.9,1.0))
      p4 <- p4 + scale_x_continuous(limits=c(0,0.03), breaks = c(0,0.03))
    }
    if (unique(df$subjid) %in% c(2176,2240)){
      p1 <- p1 + scale_x_continuous(limits=c(0.9,1.6), breaks = c(1,1.5))
    }else if (unique(df$subjid) %in% c(2172,2177,2180,2182)){
      p1 <- p1 + scale_x_continuous(limits=c(0.5,1.3), breaks = c(0.5,1))
    }else if (unique(df$subjid) %in% c(2225,2228)){
      p1 <- p1 + scale_x_continuous(limits=c(0.4,1), breaks = c(0.5,1))
    }
    if (unique(df$subjid) %in% c(2172,2177,2182)){
      p2 <- p2 + scale_x_continuous(limits=c(0.06,0.18), breaks = c(0.06,0.18))
    }else if (unique(df$subjid) %in% c(2176,2180,2225,2228,2240)){
      p2 <- p2 + scale_x_continuous(limits=c(0.01,0.08), breaks = c(0.03,0.06))
    }

    if (unique(df$subjid) %in% c(2152,2153,2154,2155,2156,2160,2165,2166)){
      p2 <- p2 + scale_x_continuous(limits=c(0.02,0.13), breaks = c(0.05,1.0))
      p3 <- p3 + scale_x_continuous(limits=c(0.6,1), breaks = c(0.6,1))
    }
    if (unique(df$subjid) %in% c(2152,2153)){
      p1 <- p1 + scale_x_continuous(limits=c(0.5,1), breaks = c(0.5,1.0))
    }else if(unique(df$subjid) %in% c(2160)){
      p1 <- p1 + scale_x_continuous(limits=c(0.4,0.65), breaks = c(0.4,0.6))
    }else if (unique(df$subjid) %in% c(2154,2155,2156,2165,2166)){
      p1 <- p1 + scale_x_continuous(limits=c(0.3,0.5), breaks = c(0.3,0.5))
    }
    if (unique(df$subjid) %in% c(2152,2153,2156,2160)){
      p4 <- p4 + scale_x_continuous(limits=c(0.05,0.35), breaks = c(0.1,0.3))
    }
  }


  if (xlab_on){
    p1 <- p1 + xlab(expression(rho))
    p2 <- p2 + xlab(expression(sigma))
    p3 <- p3 + xlab(expression(omega[rational]))
    p4 <- p4 + xlab(expression(omega[lottery]))
    p5 <- p5 + xlab(expression(omega[surebet]))
  }else{
    p1 <- p1 + theme(axis.title.x = element_blank())
    p2 <- p2 + theme(axis.title.x = element_blank())
    p3 <- p3 + theme(axis.title.x = element_blank())
    p4 <- p4 + theme(axis.title.x = element_blank())
    p5 <- p5 + theme(axis.title.x = element_blank())
  }

  p_individual <- plot_spacer() + grid::grid.text(sprintf('%d',unique(df$subjid)), x = 0.4, y = 0.5, gp = grid::gpar(fontsize = 7, fontface = "bold")) +
    p1 + p2 + p3 + p4 + p5 + plot_layout(widths = c(0.1,0.7,2,2,2,2,2),nrow = 1)
  return(p_individual)
}

plot_population_draw_dens <- function(model, dataset = 'bi_opto'){
  pars_list <- get_population_pars_from_brms(model)
  pars_df <- pars_list$pars %>% select(-contains(c('b_', 'r_')))
  subjids <- unique(pars_df$subjid)
  for (i in 1:length(subjids)){
    this_subjid <- subjids[i]
    this_pars_df <- pars_df %>% filter(subjid == this_subjid)
    this_p <- plot_individual_draw_dens(this_pars_df, dataset, i == length(subjids))
    if (i == 1){
      P = this_p
    }else{
      P <- P/this_p
    }
  }
  if (dataset == 'bi_opto'){
    plot_height = 2.4
  }else{
    plot_height = 3.7
  }
  #ggsave(sprintf("~/Desktop/individual_pars_dens_%s.png", dataset), plot = P, width = 4, units = "in", height = plot_height, scale = 1.7)
  #ggsave(sprintf("~/Desktop/individual_pars_dens_%s.pdf", dataset), plot = P, width = 4, units = "in", height = plot_height, scale = 1.7)
  return(P)
}


## plots for the synthetic data
get_mean_parameters_from_brms <- function(original_agents, deltaRho=-0.50, deltaOmega1=0, deltaOmega2=0, deltaBeta=0, synthetic_model){
  # get the estimated parameters
  synthetic_ps <- synthetic_model %>% spread_draws(b_logrho_Intercept, b_logrho_dose, r_subjid__logrho[subjid,],
                                                   b_noise_Intercept, b_noise_dose, r_subjid__noise[subjid,],
                                                   b_omega1_Intercept, b_omega1_dose, r_subjid__omega1[subjid,],
                                                   b_omega2_Intercept, b_omega2_dose, r_subjid__omega2[subjid,]) %>%
    mutate(subj_logrho = b_logrho_Intercept + r_subjid__logrho, subj_logrho_dose = b_logrho_Intercept + r_subjid__logrho + b_logrho_dose,
           subj_noise = b_noise_Intercept + r_subjid__noise, subj_noise_dose = b_noise_Intercept + r_subjid__noise + b_noise_dose,
           subj_omega1 = b_omega1_Intercept + r_subjid__omega1, subj_omega1_dose = b_omega1_Intercept + r_subjid__omega1 + b_omega1_dose,
           subj_omega2 = b_omega2_Intercept + r_subjid__omega2, subj_omega2_dose = b_omega2_Intercept + r_subjid__omega2 + b_omega2_dose) %>%
    mean_qi() %>% select(starts_with('subj')) %>% select(-contains('.lower')) %>% select(-contains('.upper'))

  synthetic_pe_control <- synthetic_ps %>% select(-contains('_dose')) %>% mutate(dose = 0)
  synthetic_pe_dose <- synthetic_ps %>% select(contains(c('subjid','_dose'))) %>% mutate(dose = 1) %>% rename(subj_logrho = subj_logrho_dose,
                                                                                                              subj_noise = subj_noise_dose,
                                                                                                              subj_omega1 = subj_omega1_dose,
                                                                                                              subj_omega2 = subj_omega2_dose)
  synthetic_pe <- rbind(synthetic_pe_control, synthetic_pe_dose) %>% mutate(param = 'e') %>% gather(Type, Value, 2:5)

  # get the original parameters
  synthetic_po_control <- original_agents %>% rename(subj_logrho = log_rho,
                                                     subj_noise = beta,
                                                     subj_omega1 = omega_1,
                                                     subj_omega2 = omega_2) %>% mutate(dose = 0)
  synthetic_po_dose <- synthetic_po_control %>% mutate(subj_logrho = subj_logrho + deltaRho,
                                                       subj_noise = subj_noise + deltaBeta,
                                                       subj_omega1 = subj_omega1 + deltaOmega1,
                                                       subj_omega2 = subj_omega2 + deltaOmega2,
                                                       dose = 1)
  synthetic_po <- rbind(synthetic_po_control, synthetic_po_dose) %>% mutate(param = 'o') %>% mutate(subj_noise = subj_noise) %>% gather(Type, Value, 2:5)
  parameters_df <- synthetic_po %>% inner_join(synthetic_pe, by = c('subjid', 'dose', 'Type')) %>% rename(original = Value.x, estimated = Value.y)
  return(parameters_df)
}

# plot the original vs. estimated parameters from synthetic data
plot_parameter_pairs <- function(df, legend_on = TRUE){
  p <- ggplot(df,aes(x = original, y = estimated, color = factor(dose, levels = c(0, 1), labels = c('control', 'perturbation')), fill = factor(dose, levels = c(0, 1), labels = c('control', 'perturbation')))) +
    theme_classic(base_size = 15) +
    geom_point(size = 2.3) +
    geom_smooth(method = 'lm') +
    geom_abline(linetype = 'dashed', color = 'gray37') +
    scale_color_manual(values = c('control' = 'azure4','perturbation' = 'purple3')) +
    scale_fill_manual(values = c('control' = 'gray70','perturbation' = 'purple3')) +
    theme(legend.title=element_blank())
  if (!legend_on){
    p <- p + theme(legend.position = 'none')
  }
  return(p)
}

make_parameters_comparation_brms <- function(){
  rho_change_parameters_df <- get_mean_parameters_from_brms(agents, deltaRho=-0.50, deltaOmega1=0, deltaOmega2=0, deltaBeta=0, rho_synthetic_cdf)
  omega_change_parameters_df <- get_mean_parameters_from_brms(agents, deltaRho=0, deltaOmega1=-1, deltaOmega2=-3, deltaBeta=0, omega_synthetic_cdf)
  rho_only_omega1_df <- rho_change_parameters_df %>% filter(Type == 'subj_omega1') %>% rename(omega1_o = original, omega1_e = estimated)
  rho_only_omega2_df <- rho_change_parameters_df %>% filter(Type == 'subj_omega2') %>% rename(omega2_o = original, omega2_e = estimated)
  omega_only_omega1_df <- omega_change_parameters_df %>% filter(Type == 'subj_omega1') %>% rename(omega1_o = original, omega1_e = estimated)
  omega_only_omega2_df <- omega_change_parameters_df %>% filter(Type == 'subj_omega2') %>% rename(omega2_o = original, omega2_e = estimated)
  rho_only_omega_df <- inner_join(rho_only_omega1_df, rho_only_omega2_df, by = c('subjid','dose')) %>% mutate(omega_rational_o = inv_logit(omega1_o),
                                                                                                              omega_rational_e = inv_logit(omega1_e),
                                                                                                              omega_lottery_o = (1-inv_logit(omega1_o)) * inv_logit(omega2_o),
                                                                                                              omega_lottery_e = (1-inv_logit(omega1_e)) * inv_logit(omega2_e))
  omega_only_omega_df <- inner_join(omega_only_omega1_df, omega_only_omega2_df, by = c('subjid','dose')) %>% mutate(omega_rational_o = inv_logit(omega1_o),
                                                                                                                    omega_rational_e = inv_logit(omega1_e),
                                                                                                                    omega_lottery_o = (1-inv_logit(omega1_o)) * inv_logit(omega2_o),
                                                                                                                    omega_lottery_e = (1-inv_logit(omega1_e)) * inv_logit(omega2_e))

  p1 <- plot_parameter_pairs(rho_change_parameters_df %>% filter(Type == 'subj_logrho') %>% mutate(original = exp(original), estimated = exp(estimated)), legend_on = FALSE) +
    scale_x_continuous(limits=c(0, 1.25), breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1.0')) +
    scale_y_continuous(limits=c(0, 1.25), breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1.0')) + xlab(expression(rho)) + ylab('Model estimate')
  p2 <- plot_parameter_pairs(rho_change_parameters_df %>% filter(Type == 'subj_noise') %>% mutate(original = exp(original), estimated = exp(estimated)), legend_on = FALSE) +
    scale_x_continuous(limits=c(0, 0.14), breaks = c(0, .07, 0.14), labels = c('0', '0.07', '0.14')) +
    scale_y_continuous(limits=c(0, 0.14), breaks = c(0, .07, 0.14), labels = c('0', '0.07', '0.14')) +
    xlab(expression(sigma)) + theme(axis.title.y = element_blank())
  p3 <- plot_parameter_pairs(rho_change_parameters_df %>% filter(Type == 'subj_omega1'), legend_on = FALSE) +
    scale_x_continuous(limits=c(0.8, 3.2), breaks = c(1, 2, 3), labels = c('1.0', '2.0', '3.0')) +
    scale_y_continuous(limits=c(0.8, 3.2), breaks = c(1, 2, 3), labels = c('1.0', '2.0', '3.0')) + xlab(expression(paste(omega, '1'))) + theme(axis.title.y = element_blank())
  p4 <- plot_parameter_pairs(rho_change_parameters_df %>% filter(Type == 'subj_omega2'), legend_on = FALSE) +
    scale_x_continuous(limits=c(-1.2, 2.2), breaks = c(-1, 0, 2), labels = c('-1.0', '0', '2.0')) +
    scale_y_continuous(limits=c(-1.2, 2.2), breaks = c(-1, 0, 2), labels = c('-1.0', '0', '2.0')) + xlab(expression(paste(omega, '2'))) + theme(axis.title.y = element_blank())
  p5 <- plot_parameter_pairs(omega_change_parameters_df %>% filter(Type == 'subj_logrho') %>% mutate(original = exp(original), estimated = exp(estimated)), legend_on = FALSE) +
    scale_x_continuous(limits=c(0, 1.25), breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1.0')) +
    scale_y_continuous(limits=c(0, 1.25), breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1.0')) + xlab(expression(rho)) + ylab('Model estimate')
  p6 <- plot_parameter_pairs(omega_change_parameters_df %>% filter(Type == 'subj_noise') %>% mutate(original = exp(original), estimated = exp(estimated)), legend_on = FALSE) +
    scale_x_continuous(limits=c(0, 0.14), breaks = c(0, .07, 0.14), labels = c('0', '0.07', '0.14')) +
    scale_y_continuous(limits=c(0, 0.14), breaks = c(0, .07, 0.14), labels = c('0', '0.07', '0.14')) + xlab(expression(sigma)) + theme(axis.title.y = element_blank())
  p7 <- plot_parameter_pairs(omega_change_parameters_df %>% filter(Type == 'subj_omega1'), legend_on = FALSE) +
    scale_x_continuous(limits=c(-0.5, 3.2), breaks = c(0, 1, 2, 3), labels = c('0', '1.0', '2.0', '3.0')) +
    scale_y_continuous(limits=c(-0.5, 3.2), breaks = c(0, 1, 2, 3), labels = c('0', '1.0', '2.0', '3.0')) + xlab(expression(paste(omega, '1'))) + theme(axis.title.y = element_blank())
  p8 <- plot_parameter_pairs(omega_change_parameters_df %>% filter(Type == 'subj_omega2')) +
    scale_x_continuous(limits=c(-4.5, 2.5), breaks = c(-4, -2, 0, 2), labels = c('-4.0', '-2.0', '0', '2.0')) +
    scale_y_continuous(limits=c(-4.5, 2.5), breaks = c(-4, -2, 0, 2), labels = c('-4.0', '-2.0', '0', '2.0')) + xlab(expression(paste(omega, '2'))) + theme(axis.title.y = element_blank())

  layout <- "
  ABCD
  EFGH
  "
  p <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + guide_area() + plot_layout(design = layout)
  return(p)
}

## plot the original model fitted posterior together
plot_delta_phi_psd <- function(bi_opto_model, uni_opto_model, bi_fof_muscimol_model, uni_fof_muscimol_model, y_lab = TRUE){
  ps_bi_opto <- posterior_samples(bi_opto_model, pars = '^b', as.matrix = T)
  ps_bi_opto <- as.data.frame(ps_bi_opto) %>% mutate(delta_phi_bi_opto = b_logrho_dose) %>% select(c('delta_phi_bi_opto'))
  ps_uni_opto <- posterior_samples(uni_opto_model, pars = '^b', as.matrix = T)
  ps_uni_opto <- as.data.frame(ps_uni_opto) %>% mutate(delta_phi_uni_opto = b_logrho_dose) %>% select(c('delta_phi_uni_opto'))
  ps_bi_muscimol <- posterior_samples(bi_fof_muscimol_model, pars = '^b', as.matrix = T)
  ps_bi_muscimol <- as.data.frame(ps_bi_muscimol) %>% mutate(delta_phi_bi_muscimol = b_logrho_dose) %>% select(c('delta_phi_bi_muscimol'))
  ps_uni_muscimol <- posterior_samples(uni_fof_muscimol_model, pars = '^b', as.matrix = T)
  ps_uni_muscimol <- as.data.frame(ps_uni_muscimol) %>% mutate(delta_phi_uni_muscimol = b_logrho_dose) %>% select(c('delta_phi_uni_muscimol'))

  ps_delta_rho <- cbind(ps_bi_opto, ps_uni_opto, ps_bi_muscimol, ps_uni_muscimol) %>% as.mcmc()
  p <- mcmc_areas(ps_delta_rho, prob=.7, prob_outer = 0.9999, point_est = 'mean', area_method = 'equal height') +
    vline_0(size = .2, color = 'gray') +
    theme_classic(base_size = 10) +
    scale_x_continuous(limits=c(-0.75, 0.25), breaks = c(-0.6, 0, 2), labels = c('-0.6', '0', '0.2')) +
    scale_y_discrete(limits = c('delta_phi_uni_muscimol','delta_phi_bi_muscimol','delta_phi_uni_opto','delta_phi_bi_opto'),
                     labels = c(
                       'delta_phi_bi_opto' = 'Bi-Opto',
                       'delta_phi_uni_opto' = 'Uni-Opto',
                       'delta_phi_bi_muscimol'= 'Bi-Muscimol',
                       'delta_phi_uni_muscimol'= 'Uni-Muscimol'),
                     expand = expansion(add = c(.3, .3))) +
    xlab(expression(paste(Delta, phi))) +
    theme(axis.text.y = element_text(angle = 60))

  if (!y_lab){
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}

plot_delta_psi_psd <- function(bi_opto_model, uni_opto_model, bi_fof_muscimol_model, uni_fof_muscimol_model, y_lab = TRUE){
  ps_bi_opto <- posterior_samples(bi_opto_model, pars = '^b', as.matrix = T)
  ps_bi_opto <- as.data.frame(ps_bi_opto) %>% mutate(delta_psi_bi_opto = b_noise_dose) %>% select(c('delta_psi_bi_opto'))
  ps_uni_opto <- posterior_samples(uni_opto_model, pars = '^b', as.matrix = T)
  ps_uni_opto <- as.data.frame(ps_uni_opto) %>% mutate(delta_psi_uni_opto = b_noise_dose) %>% select(c('delta_psi_uni_opto'))
  ps_bi_muscimol <- posterior_samples(bi_fof_muscimol_model, pars = '^b', as.matrix = T)
  ps_bi_muscimol <- as.data.frame(ps_bi_muscimol) %>% mutate(delta_psi_bi_muscimol = b_noise_dose) %>% select(c('delta_psi_bi_muscimol'))
  ps_uni_muscimol <- posterior_samples(uni_fof_muscimol_model, pars = '^b', as.matrix = T)
  ps_uni_muscimol <- as.data.frame(ps_uni_muscimol) %>% mutate(delta_psi_uni_muscimol = b_noise_dose) %>% select(c('delta_psi_uni_muscimol'))

  ps_delta_rho <- cbind(ps_bi_opto, ps_uni_opto, ps_bi_muscimol, ps_uni_muscimol) %>% as.mcmc()
  p <- mcmc_areas(ps_delta_rho, prob=.7, prob_outer = 0.9999, point_est = 'mean', area_method = 'equal height') +
    vline_0(size = .2, color = 'gray') +
    theme_classic(base_size = 10) +
    scale_x_continuous(limits=c(-0.9, 1.3), breaks = c(-0.5, 0, 0.5, 1.0)) +
    scale_y_discrete(limits = c('delta_psi_uni_muscimol','delta_psi_bi_muscimol','delta_psi_uni_opto','delta_psi_bi_opto'),
                     labels = c(
                       'delta_psi_bi_opto' = 'Bi-Opto',
                       'delta_psi_uni_opto' = 'Uni-Opto',
                       'delta_psi_bi_muscimol'= 'Bi-Muscimol',
                       'delta_psi_uni_muscimol'= 'Uni-Muscimol'),
                     expand = expansion(add = c(.3, .3))) +
    xlab(expression(paste(Delta, psi)))
  if (!y_lab){
    p <- p + theme(axis.text.y = element_blank())
  }
}

plot_delta_omega1_psd <- function(bi_opto_model, uni_opto_model, bi_fof_muscimol_model, uni_fof_muscimol_model, y_lab = TRUE){
  ps_bi_opto <- posterior_samples(bi_opto_model, pars = '^b', as.matrix = T)
  ps_bi_opto <- as.data.frame(ps_bi_opto) %>% mutate(delta_omega1_bi_opto = b_omega1_dose) %>% select(c('delta_omega1_bi_opto'))
  ps_uni_opto <- posterior_samples(uni_opto_model, pars = '^b', as.matrix = T)
  ps_uni_opto <- as.data.frame(ps_uni_opto) %>% mutate(delta_omega1_uni_opto = b_omega1_dose) %>% select(c('delta_omega1_uni_opto'))
  ps_bi_muscimol <- posterior_samples(bi_fof_muscimol_model, pars = '^b', as.matrix = T)
  ps_bi_muscimol <- as.data.frame(ps_bi_muscimol) %>% mutate(delta_omega1_bi_muscimol = b_omega1_dose) %>% select(c('delta_omega1_bi_muscimol'))
  ps_uni_muscimol <- posterior_samples(uni_fof_muscimol_model, pars = '^b', as.matrix = T)
  ps_uni_muscimol <- as.data.frame(ps_uni_muscimol) %>% mutate(delta_omega1_uni_muscimol = b_omega1_dose) %>% select(c('delta_omega1_uni_muscimol'))

  ps_delta_rho <- cbind(ps_bi_opto, ps_uni_opto, ps_bi_muscimol, ps_uni_muscimol) %>% as.mcmc()
  p <- mcmc_areas(ps_delta_rho, prob=.7, prob_outer = 0.9999, point_est = 'mean', area_method = 'equal height') +
    vline_0(size = .2, color = 'gray') +
    theme_classic(base_size = 10) +
    scale_x_continuous(limits=c(-3, 5), breaks = c(-2, 0, 2, 4)) +
    scale_y_discrete(limits = c('delta_omega1_uni_muscimol','delta_omega1_bi_muscimol','delta_omega1_uni_opto','delta_omega1_bi_opto'),
                     labels = c(
                       'delta_omega1_bi_opto' = 'Bi-Opto',
                       'delta_omega1_uni_opto' = 'Uni-Opto',
                       'delta_omega1_bi_muscimol'= 'Bi-Muscimol',
                       'delta_omega1_uni_muscimol'= 'Uni-Muscimol'),
                     expand = expansion(add = c(.3, .3))) +
    xlab(expression(paste(Delta, omega['1'])))

  if (!y_lab){
    p <- p + theme(axis.text.y = element_blank())
  }
}

plot_delta_omega2_psd <- function(bi_opto_model, uni_opto_model, bi_fof_muscimol_model, uni_fof_muscimol_model, y_lab = TRUE){
  ps_bi_opto <- posterior_samples(bi_opto_model, pars = '^b', as.matrix = T)
  ps_bi_opto <- as.data.frame(ps_bi_opto) %>% mutate(delta_omega2_bi_opto = b_omega2_dose) %>% select(c('delta_omega2_bi_opto'))
  ps_uni_opto <- posterior_samples(uni_opto_model, pars = '^b', as.matrix = T)
  ps_uni_opto <- as.data.frame(ps_uni_opto) %>% mutate(delta_omega2_uni_opto = b_omega2_dose) %>% select(c('delta_omega2_uni_opto'))
  ps_bi_muscimol <- posterior_samples(bi_fof_muscimol_model, pars = '^b', as.matrix = T)
  ps_bi_muscimol <- as.data.frame(ps_bi_muscimol) %>% mutate(delta_omega2_bi_muscimol = b_omega2_dose) %>% select(c('delta_omega2_bi_muscimol'))
  ps_uni_muscimol <- posterior_samples(uni_fof_muscimol_model, pars = '^b', as.matrix = T)
  ps_uni_muscimol <- as.data.frame(ps_uni_muscimol) %>% mutate(delta_omega2_uni_muscimol = b_omega2_dose) %>% select(c('delta_omega2_uni_muscimol'))

  ps_delta_rho <- cbind(ps_bi_opto, ps_uni_opto, ps_bi_muscimol, ps_uni_muscimol) %>% as.mcmc()
  p <- mcmc_areas(ps_delta_rho, prob=.7, prob_outer = 0.9999, point_est = 'mean', area_method = 'equal height') +
    vline_0(size = .2, color = 'gray') +
    theme_classic(base_size = 10) +
    scale_x_continuous(limits=c(-6.5, 3.5), breaks = c(-5, -2.5, 0, 2.5)) +
    scale_y_discrete(limits = c('delta_omega2_uni_muscimol','delta_omega2_bi_muscimol','delta_omega2_uni_opto','delta_omega2_bi_opto'),
                     labels = c(
                       'delta_omega2_bi_opto' = 'Bi-Opto',
                       'delta_omega2_uni_opto' = 'Uni-Opto',
                       'delta_omega2_bi_muscimol'= 'Bi-Muscimol',
                       'delta_omega2_uni_muscimol'= 'Uni-Muscimol'),
                     expand = expansion(add = c(.3, .3))) +
    xlab(expression(paste(Delta, omega['2'])))

  if (!y_lab){
    p <- p + theme(axis.text.y = element_blank())
  }
}

make_brms_original_psd <- function(){
  p1 <- plot_delta_phi_psd(bi_opto_cdf, uni_opto_cdf, bi_fof_muscimol_cdf, uni_fof_muscimol_cdf)
  p2 <- plot_delta_psi_psd(bi_opto_cdf, uni_opto_cdf, bi_fof_muscimol_cdf, uni_fof_muscimol_cdf, y_lab = F)
  p3 <- plot_delta_omega1_psd(bi_opto_cdf, uni_opto_cdf, bi_fof_muscimol_cdf, uni_fof_muscimol_cdf, y_lab = F)
  p4 <- plot_delta_omega2_psd(bi_opto_cdf, uni_opto_cdf, bi_fof_muscimol_cdf, uni_fof_muscimol_cdf, y_lab = F)
  P <- p1 + p2 + p3 + p4 + plot_layout(nrow = 1)
  return(P)
}

# plots the mcmc pairs for the opto and muscimol brms model
make_mcmc_pairs_plots <- function(model, dataset = 'opto'){
  bayesplot_theme_set()
  color_scheme_set(scheme = rev(colfunc(10)[1:6]))
  bayesplot_theme_update(text = element_text(size = 15))
  delta_pars <- model %>% as_draws_array(variable = '^b_.*_dose', regex = TRUE) #%>% rename(delta_rho = b_logrho_dose, delta_noise = b_noise_dose, delta_omega1 = b_omega1_dose, delta_omega2 = b_omega2_dose)
  dimnames(delta_pars)[[3]] <- c('Δϕ', 'Δψ', 'Δω₁', 'Δω₂')
  p <- mcmc_pairs(delta_pars, diag_fun = 'dens', off_diag_args = list(size = 1, alpha = 0.5))
  bayesplot_theme_update(text = element_text(size = 20, family = "sans"))
  #ggsave(sprintf("~/Desktop/parameter_mcmc_pairs_%s.png", dataset), plot = p, width = 2, units = "in", height = 2, scale = 4.5)
  #ggsave(sprintf("~/Desktop/parameter_mcmc_pairs_%s.pdf", dataset), plot = p, width = 2, units = "in", height = 2, scale = 4.5)
  return(p)
}

##  version 2 to plot the mcmc pairs for the opto and muscimol brms model
my_under_fun <- function(data, mapping, ...){
  # Using default ggplot density function

  p <- ggplot(data = data, mapping = mapping) +
    stat_density_2d(aes(fill = as.factor(..level..)), geom = "polygon", n=100, contour_var = "ndensity") +
    scale_fill_manual(values = rev(colfunc(16)[1:14]))
  return(p)
}

my_diag_fun <- function(data, mapping, ...){
  # Using default ggplot density function

  p <- ggplot(data = data, mapping = mapping) +
    geom_density(color='#D766C9', fill=colfunc(13)[6]) +
    geom_vline(xintercept = 0, color = "grey40", size=0.5, linetype="dashed")
  p
}

make_mcmc_pairs_plots_verson2 <- function(model, dataset = 'bi_opto'){
  colrs = rev(colfunc(10)[1:6])
  pars <- posterior_samples(model, pars = '^b', as.matrix = T)
  df_pars <- data.frame(pars) %>% select(c(b_logrho_dose, b_noise_dose, b_omega1_dose, b_omega2_dose)) %>%
    rename('delta_phi' = 'b_logrho_dose', 'delta_psi' = 'b_noise_dose', 'delta_omega1'='b_omega1_dose', 'delta_omega2'='b_omega2_dose')
  #delta_pars <- model %>% as_draws_array(variable = '^b_.*_dose', regex = TRUE) #%>% rename(delta_rho = b_logrho_dose, delta_noise = b_noise_dose, delta_omega1 = b_omega1_dose, delta_omega2 = b_omega2_dose)
  #dimnames(delta_pars)[[3]] <- c('Δϕ', 'Δψ', 'Δω₁', 'Δω₂')
  p <- ggpairs(df_pars, lower=list(continuous=my_under_fun), diag=list(continuous=my_diag_fun),
               columnLabels = c("Delta*phi", "Delta*psi", "Delta*omega[1]", "Delta*omega[2]"),
               labeller = "label_parsed") + theme_classic() + theme(text = element_text(size = 20))
  p2 <- p
  if (dataset == 'bi_opto'){
    pair_lims_l = c(-0.5, -0.2, -2.5, -5)
    pair_lims_r = c(0.1, 0.7, 1.5, 3)
    pair_lims_lb = c(-0.3, 0, -1.5, -3)
    pair_lims_rb = c(0, 0.5, 0, 2)
    pair_lims_ylb = c(0, 0, -1, -2)
    pair_lims_yrb = c(4, 0.4, 0, 0)
    for(i in 1:p$nrow) {
      for(j in 1:i) {
        p2[i,j] <- p[i,j] +
          scale_x_continuous(limits = c(pair_lims_l[j], pair_lims_r[j]), breaks = c(pair_lims_lb[j], 0, pair_lims_rb[j]))
        if (i>1 & j<i){
          p2[i,j] <- p2[i,j] + scale_y_continuous(breaks = c(pair_lims_ylb[i], 0, pair_lims_yrb[i]))
        }
      }
    }
  }else if (dataset == 'bi_fof_muscimol'){
    pair_lims_l = c(-0.7, -0.6, -3, -6)
    pair_lims_r = c(0.3, 1.3, 4, 4)
    pair_lims_lb = c(-0.5, 0, -2, -4)
    pair_lims_rb = c(0, 1, 2, 4)
    pair_lims_ylb = c(0, 0, -2, -4)
    pair_lims_yrb = c(4, 0.5, 2, 0)
    for(i in 1:p$nrow) {
      for(j in 1:i) {
        p2[i,j] <- p[i,j] +
          scale_x_continuous(limits = c(pair_lims_l[j], pair_lims_r[j]), breaks = c(pair_lims_lb[j], 0, pair_lims_rb[j]))
        if (i>1 & j<i){
          p2[i,j] <- p2[i,j] + scale_y_continuous(breaks = c(pair_lims_ylb[i], 0, pair_lims_yrb[i]))
        }
      }
    }
  }else if (dataset == 'uni_opto'){
    pair_lims_l = c(-0.3, -0.2, -1.5, -4)
    pair_lims_r = c(0.1, 0.4, 1, 2)
    pair_lims_lb = c(-0.2, 0, -1, -2)
    pair_lims_rb = c(0, 0.2, 0, 2)
    pair_lims_ylb = c(0, 0, -0.5, -2)
    pair_lims_yrb = c(0.2, 0.2, 0, 0)
    for(i in 1:p$nrow) {
      for(j in 1:i) {
        p2[i,j] <- p[i,j] +
          scale_x_continuous(limits = c(pair_lims_l[j], pair_lims_r[j]), breaks = c(pair_lims_lb[j], 0, pair_lims_rb[j]))
        if (i>1 & j<i){
          p2[i,j] <- p2[i,j] + scale_y_continuous(breaks = c(pair_lims_ylb[i], 0, pair_lims_yrb[i]))
        }
      }
    }
  }else{
    pair_lims_l = c(-0.2, -0.1, -1, -2)
    pair_lims_r = c(0.03, 0.6, 0.5, 0.5)
    pair_lims_lb = c(-0.2, 0, -1, -1)
    pair_lims_rb = c(0, 0.4, 0, 0)
    pair_lims_ylb = c(0, 0.1, -0.4, -1)
    pair_lims_yrb = c(0.2, 0.3, 0, -0.5)
    for(i in 1:p$nrow) {
      for(j in 1:i) {
        p2[i,j] <- p[i,j] +
          scale_x_continuous(limits = c(pair_lims_l[j], pair_lims_r[j]), breaks = c(pair_lims_lb[j], 0, pair_lims_rb[j]))
        if (i>1 & j<i){
          p2[i,j] <- p2[i,j] + scale_y_continuous(breaks = c(pair_lims_ylb[i], pair_lims_yrb[i]))
        }
      }
    }
  }

  #ggsave(sprintf("~/Desktop/parameter_mcmc_pairs_%s.png", dataset), plot = p2, width = 2, units = "in", height = 2, scale = 3)
  #ggsave(sprintf("~/Desktop/parameter_mcmc_pairs_%s.pdf", dataset), plot = p2, width = 2, units = "in", height = 2, scale = 3)
  return(p2)
}


## plot the PPC muscimol posterior together
plot_delta_rho_psd_ppc <- function(bi_opto_model, uni_opto_model, bi_muscimol_model, uni_muscimol_model, y_lab = TRUE){
  color_scheme_set(scheme = rev(colfunc(10)[1:6]))
  ps_bi_opto <- posterior_samples(bi_opto_model, pars = '^b', as.matrix = T)
  ps_bi_opto <- as.data.frame(ps_bi_opto) %>% mutate(delta_rho_bi_opto = exp(b_logrho_Intercept + b_logrho_dose) - exp(b_logrho_Intercept)) %>% select(c('delta_rho_bi_opto'))
  ps_uni_opto <- posterior_samples(uni_opto_model, pars = '^b', as.matrix = T)
  ps_uni_opto <- as.data.frame(ps_uni_opto) %>% mutate(delta_rho_uni_opto = exp(b_logrho_Intercept + b_logrho_dose) - exp(b_logrho_Intercept)) %>% select(c('delta_rho_uni_opto'))
  ps_bi_muscimol <- posterior_samples(bi_muscimol_model, pars = '^b', as.matrix = T)
  ps_bi_muscimol <- as.data.frame(ps_bi_muscimol) %>% mutate(delta_rho_bi_muscimol = exp(b_logrho_Intercept + b_logrho_dose) - exp(b_logrho_Intercept)) %>% select(c('delta_rho_bi_muscimol'))
  ps_uni_muscimol <- posterior_samples(uni_muscimol_model, pars = '^b', as.matrix = T)
  ps_uni_muscimol <- as.data.frame(ps_uni_muscimol) %>% mutate(delta_rho_uni_muscimol = exp(b_logrho_Intercept + b_logrho_dose) - exp(b_logrho_Intercept)) %>% select(c('delta_rho_uni_muscimol'))

  ps_delta_rho <- cbind(ps_bi_opto, ps_uni_opto, ps_bi_muscimol, ps_uni_muscimol) %>% as.mcmc()
  p <- mcmc_areas(ps_delta_rho, prob=.7, prob_outer = 0.9999, point_est = 'mean', area_method = 'equal height') +
    vline_0(size = .2, color = 'gray') +
    theme_classic(base_size = 10) +
    scale_x_continuous(limits=c(-0.2, 0.2), breaks = c(-0.2, 0), labels = c('-0.2', '0')) +
    scale_y_discrete(limits = c('delta_rho_uni_muscimol','delta_rho_bi_muscimol','delta_rho_uni_opto','delta_rho_bi_opto'),
                     labels = c(
                       'delta_rho_bi_opto' = 'Bi-PPC',
                       'delta_rho_uni_opto' = 'Uni-PPC',
                       'delta_rho_bi_muscimol'= 'Bi-PPC-40',
                       'delta_rho_uni_muscimol'= 'Uni-PPC-40'),
                     expand = expansion(add = c(.3, .3))) +
    xlab(expression(paste(Delta, rho))) +
    theme(axis.text.y = element_text(angle = 60))

  if (!y_lab){
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}

plot_delta_noise_psd_ppc <- function(bi_opto_model, uni_opto_model, bi_muscimol_model, uni_muscimol_model, y_lab = TRUE){
  ps_bi_opto <- posterior_samples(bi_opto_model, pars = '^b', as.matrix = T)
  ps_bi_opto <- as.data.frame(ps_bi_opto) %>% mutate(delta_noise_bi_opto = exp(b_noise_Intercept + b_noise_dose) - exp(b_noise_Intercept)) %>% select(c('delta_noise_bi_opto'))
  ps_uni_opto <- posterior_samples(uni_opto_model, pars = '^b', as.matrix = T)
  ps_uni_opto <- as.data.frame(ps_uni_opto) %>% mutate(delta_noise_uni_opto = exp(b_noise_Intercept + b_noise_dose) - exp(b_noise_Intercept)) %>% select(c('delta_noise_uni_opto'))
  ps_bi_muscimol <- posterior_samples(bi_muscimol_model, pars = '^b', as.matrix = T)
  ps_bi_muscimol <- as.data.frame(ps_bi_muscimol) %>% mutate(delta_noise_bi_muscimol = exp(b_noise_Intercept + b_noise_dose) - exp(b_noise_Intercept)) %>% select(c('delta_noise_bi_muscimol'))
  ps_uni_muscimol <- posterior_samples(uni_muscimol_model, pars = '^b', as.matrix = T)
  ps_uni_muscimol <- as.data.frame(ps_uni_muscimol) %>% mutate(delta_noise_uni_muscimol = exp(b_noise_Intercept + b_noise_dose) - exp(b_noise_Intercept)) %>% select(c('delta_noise_uni_muscimol'))


  ps_delta_rho <- cbind(ps_bi_opto, ps_uni_opto, ps_bi_muscimol, ps_uni_muscimol) %>% as.mcmc()
  p <- mcmc_areas(ps_delta_rho, prob=.7, prob_outer = 0.9999, point_est = 'mean', area_method = 'equal height') +
    vline_0(size = .2, color = 'gray') +
    theme_classic(base_size = 10) +
    scale_x_continuous(limits=c(-0.05, 0.08), breaks = c(-0.05, 0, 0.05), labels = c('-0.05', '0', '0.05')) +
    scale_y_discrete(limits = c('delta_noise_uni_muscimol','delta_noise_bi_muscimol','delta_noise_uni_opto','delta_noise_bi_opto'),
                     labels = c(
                       'delta_noise_bi_opto' = 'Bi-PPC',
                       'delta_noise_uni_opto' = 'Uni-PPC',
                       'delta_noise_bi_muscimol'= 'Bi-PPC-40',
                       'delta_noise_uni_muscimol'= 'Uni-PPC-40'),
                     expand = expansion(add = c(.3, .3))) +
    xlab(expression(paste(Delta, sigma)))
  if (!y_lab){
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}

plot_delta_omega_rational_psd_ppc <- function(bi_opto_model, uni_opto_model, bi_muscimol_model, uni_muscimol_model, y_lab = TRUE){
  ps_bi_opto <- posterior_samples(bi_opto_model, pars = '^b', as.matrix = T)
  ps_bi_opto <- as.data.frame(ps_bi_opto) %>% mutate(delta_omega_rational_bi_opto = inv_logit(b_omega1_Intercept + b_omega1_dose) - inv_logit(b_omega1_Intercept)) %>% select(c('delta_omega_rational_bi_opto'))
  ps_uni_opto <- posterior_samples(uni_opto_model, pars = '^b', as.matrix = T)
  ps_uni_opto <- as.data.frame(ps_uni_opto) %>% mutate(delta_omega_rational_uni_opto = inv_logit(b_omega1_Intercept + b_omega1_dose) - inv_logit(b_omega1_Intercept)) %>% select(c('delta_omega_rational_uni_opto'))
  ps_bi_muscimol <- posterior_samples(bi_muscimol_model, pars = '^b', as.matrix = T)
  ps_bi_muscimol <- as.data.frame(ps_bi_muscimol) %>% mutate(delta_omega_rational_bi_muscimol = inv_logit(b_omega1_Intercept + b_omega1_dose) - inv_logit(b_omega1_Intercept)) %>% select(c('delta_omega_rational_bi_muscimol'))
  ps_uni_muscimol <- posterior_samples(uni_muscimol_model, pars = '^b', as.matrix = T)
  ps_uni_muscimol <- as.data.frame(ps_uni_muscimol) %>% mutate(delta_omega_rational_uni_muscimol = inv_logit(b_omega1_Intercept + b_omega1_dose) - inv_logit(b_omega1_Intercept)) %>% select(c('delta_omega_rational_uni_muscimol'))

  ps_delta_rho <- cbind(ps_bi_opto, ps_uni_opto, ps_bi_muscimol, ps_uni_muscimol) %>% as.mcmc()
  p <- mcmc_areas(ps_delta_rho, prob=.7, prob_outer = 0.9999, point_est = 'mean', area_method = 'equal height') +
    vline_0(size = .2, color = 'gray') +
    theme_classic(base_size = 10) +
    scale_x_continuous(limits=c(-0.3, 0.3), breaks = c(-0.2, 0, 0.2), labels = c('-0.2', '0', '0.2')) +
    scale_y_discrete(limits = c('delta_omega_rational_uni_muscimol','delta_omega_rational_bi_muscimol','delta_omega_rational_uni_opto','delta_omega_rational_bi_opto'),
                     labels = c(
                       'delta_omega_rational_bi_opto' = 'Bi-PPC',
                       'delta_omega_rational_uni_opto' = 'Uni-PPC',
                       'delta_omega_rational_bi_muscimol'= 'Bi-PPC-40',
                       'delta_omega_rational_uni_muscimol'= 'Uni-PPC-40'),
                     expand = expansion(add = c(.3, .3))) +
    xlab(expression(paste(Delta, omega['rational'])))

  if (!y_lab){
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}

plot_delta_omega_lottery_psd_ppc <- function(bi_opto_model, uni_opto_model, bi_muscimol_model, uni_muscimol_model, y_lab = TRUE){
  ps_bi_opto <- posterior_samples(bi_opto_model, pars = '^b', as.matrix = T)
  ps_bi_opto <- as.data.frame(ps_bi_opto) %>% mutate(delta_omega_lottery_bi_opto = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*inv_logit(b_omega2_Intercept + b_omega2_dose) -
                                                       (1- inv_logit(b_omega1_Intercept))*inv_logit(b_omega2_Intercept)) %>% select(c('delta_omega_lottery_bi_opto'))
  ps_uni_opto <- posterior_samples(uni_opto_model, pars = '^b', as.matrix = T)
  ps_uni_opto <- as.data.frame(ps_uni_opto) %>% mutate(delta_omega_lottery_uni_opto = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*inv_logit(b_omega2_Intercept + b_omega2_dose) -
                                                         (1- inv_logit(b_omega1_Intercept))*inv_logit(b_omega2_Intercept)) %>% select(c('delta_omega_lottery_uni_opto'))
  ps_bi_muscimol <- posterior_samples(bi_muscimol_model, pars = '^b', as.matrix = T)
  ps_bi_muscimol <- as.data.frame(ps_bi_muscimol) %>% mutate(delta_omega_lottery_bi_muscimol = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*inv_logit(b_omega2_Intercept + b_omega2_dose) -
                                                               (1- inv_logit(b_omega1_Intercept))*inv_logit(b_omega2_Intercept)) %>% select(c('delta_omega_lottery_bi_muscimol'))
  ps_uni_muscimol <- posterior_samples(uni_muscimol_model, pars = '^b', as.matrix = T)
  ps_uni_muscimol <- as.data.frame(ps_uni_muscimol) %>% mutate(delta_omega_lottery_uni_muscimol = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*inv_logit(b_omega2_Intercept + b_omega2_dose) -
                                                                 (1- inv_logit(b_omega1_Intercept))*inv_logit(b_omega2_Intercept)) %>% select(c('delta_omega_lottery_uni_muscimol'))

  ps_delta_rho <- cbind(ps_bi_opto, ps_uni_opto, ps_bi_muscimol, ps_uni_muscimol) %>% as.mcmc()
  p <- mcmc_areas(ps_delta_rho, prob=.7, prob_outer = 0.9999, point_est = 'mean', area_method = 'equal height') +
    vline_0(size = .2, color = 'gray') +
    theme_classic(base_size = 10) +
    scale_x_continuous(limits=c(-0.2, 0.1), breaks = c(-0.2, 0, 0.1), labels = c('-0.2', '0', '0.1')) +
    scale_y_discrete(limits = c('delta_omega_lottery_uni_muscimol','delta_omega_lottery_bi_muscimol','delta_omega_lottery_uni_opto','delta_omega_lottery_bi_opto'),
                     labels = c(
                       'delta_omega_lottery_bi_opto' = 'Bi-PPC',
                       'delta_omega_lottery_uni_opto' = 'Uni-PPC',
                       'delta_omega_lottery_bi_muscimol'= 'Bi-PPC-40',
                       'delta_omega_lottery_uni_muscimol'= 'Uni-PPC-40'),
                     expand = expansion(add = c(.3, .3))) +
    xlab(expression(paste(Delta, omega['lottery'])))

  if (!y_lab){
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}


plot_delta_omega_surebet_psd_ppc <- function(bi_opto_model, uni_opto_model, bi_muscimol_model, uni_muscimol_model, y_lab = TRUE){
  ps_bi_opto <- posterior_samples(bi_opto_model, pars = '^b', as.matrix = T)
  ps_bi_opto <- as.data.frame(ps_bi_opto) %>% mutate(delta_omega_surebet_bi_opto = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*(1-inv_logit(b_omega2_Intercept + b_omega2_dose)) -
                                                       (1- inv_logit(b_omega1_Intercept))*(1-inv_logit(b_omega2_Intercept))) %>% select(c('delta_omega_surebet_bi_opto'))
  ps_uni_opto <- posterior_samples(uni_opto_model, pars = '^b', as.matrix = T)
  ps_uni_opto <- as.data.frame(ps_uni_opto) %>% mutate(delta_omega_surebet_uni_opto = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*(1-inv_logit(b_omega2_Intercept + b_omega2_dose)) -
                                                         (1- inv_logit(b_omega1_Intercept))*(1-inv_logit(b_omega2_Intercept))) %>% select(c('delta_omega_surebet_uni_opto'))
  ps_bi_muscimol <- posterior_samples(bi_muscimol_model, pars = '^b', as.matrix = T)
  ps_bi_muscimol <- as.data.frame(ps_bi_muscimol) %>% mutate(delta_omega_surebet_bi_muscimol = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*(1-inv_logit(b_omega2_Intercept + b_omega2_dose)) -
                                                               (1- inv_logit(b_omega1_Intercept))*(1-inv_logit(b_omega2_Intercept))) %>% select(c('delta_omega_surebet_bi_muscimol'))
  ps_uni_muscimol <- posterior_samples(uni_muscimol_model, pars = '^b', as.matrix = T)
  ps_uni_muscimol <- as.data.frame(ps_uni_muscimol) %>% mutate(delta_omega_surebet_uni_muscimol = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*(1-inv_logit(b_omega2_Intercept + b_omega2_dose)) -
                                                                 (1- inv_logit(b_omega1_Intercept))*(1-inv_logit(b_omega2_Intercept))) %>% select(c('delta_omega_surebet_uni_muscimol'))

  ps_delta_rho <- cbind(ps_bi_opto, ps_uni_opto, ps_bi_muscimol, ps_uni_muscimol) %>% as.mcmc()
  p <- mcmc_areas(ps_delta_rho, prob=.7, prob_outer = 0.9999, point_est = 'mean', area_method = 'equal height') +
    vline_0(size = .2, color = 'gray') +
    theme_classic(base_size = 10) +
    scale_x_continuous(limits=c(-0.12, 0.4), breaks = c(0, 0.2, 0.4), labels = c('0', '0.2', '0.4')) +
    scale_y_discrete(limits = c('delta_omega_surebet_uni_muscimol','delta_omega_surebet_bi_muscimol','delta_omega_surebet_uni_opto','delta_omega_surebet_bi_opto'),
                     labels = c(
                       'delta_omega_surebet_bi_opto' = 'Bi-PPC',
                       'delta_omega_surebet_uni_opto' = 'Uni-PPC',
                       'delta_omega_surebet_bi_muscimol'= 'Bi-PPC-40',
                       'delta_omega_surebet_uni_muscimol'= 'Uni-PPC-40'),
                     expand = expansion(add = c(.3, .3))) +
    xlab(expression(paste(Delta, omega['surebet'])))

  if (!y_lab){
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}

make_main_PPC_fig <- function(){
  layout2 <- "
  ACDEFG
  BCDEFG
  "
  p1 <- plot_pred_sample_individual(2156, bidf_ppc_muscimol, bi_ppc_muscimol_cdf, legend_on = FALSE, x_lab_on = FALSE)
  p2 <- plot_pred_sample_individual(2156, bidf_ppc_muscimol_first40, bi_ppc_muscimol_cdf_40, legend_on = FALSE)
  p3 <- plot_delta_rho_psd_ppc(bi_ppc_muscimol_cdf, uni_ppc_muscimol_cdf, bi_ppc_muscimol_cdf_40, uni_ppc_muscimol_cdf_40)
  p4 <- plot_delta_noise_psd_ppc(bi_ppc_muscimol_cdf, uni_ppc_muscimol_cdf, bi_ppc_muscimol_cdf_40, uni_ppc_muscimol_cdf_40, y_lab = F)
  p5 <- plot_delta_omega_rational_psd_ppc(bi_ppc_muscimol_cdf, uni_ppc_muscimol_cdf, bi_ppc_muscimol_cdf_40, uni_ppc_muscimol_cdf_40, y_lab = F)
  p6 <- plot_delta_omega_lottery_psd_ppc(bi_ppc_muscimol_cdf, uni_ppc_muscimol_cdf, bi_ppc_muscimol_cdf_40, uni_ppc_muscimol_cdf_40, y_lab = F)
  p7 <- plot_delta_omega_surebet_psd_ppc(bi_ppc_muscimol_cdf, uni_ppc_muscimol_cdf, bi_ppc_muscimol_cdf_40, uni_ppc_muscimol_cdf_40, y_lab = F)
  psum2 <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + guide_area() + plot_layout(design = layout2, widths = c(1,.8,.8,.8,.8,.8))
  ggsave("~/Desktop/brms_PPC.png", plot = psum2, width = 2.5, units = "in", height = 0.9, scale = 3.5)
  ggsave("~/Desktop/brms_PPC.pdf", plot = psum2, width = 2.5, units = "in", height = 0.9, scale = 3.5)
}

make_main_FOF_first40 <- function(){
  layout2 <- "
  ACDEFG
  BCDEFG
  "
  p1 <- plot_pred_sample_individual(2152, bidf_fof_muscimol, bi_fof_muscimol_cdf, legend_on = FALSE, x_lab_on = FALSE)
  p2 <- plot_pred_sample_individual(2152, bidf_fof_muscimol_first40, bi_fof_muscimol_cdf_40, legend_on = FALSE)
  p3 <- plot_delta_rho_psd(bi_fof_muscimol_cdf, uni_fof_muscimol_cdf, bi_fof_muscimol_cdf_40, uni_fof_muscimol_cdf_40)
  p4 <- plot_delta_noise_psd(bi_fof_muscimol_cdf, uni_fof_muscimol_cdf, bi_fof_muscimol_cdf_40, uni_fof_muscimol_cdf_40, y_lab = F)
  p5 <- plot_delta_omega_rational_psd(bi_fof_muscimol_cdf, uni_fof_muscimol_cdf, bi_fof_muscimol_cdf_40, uni_fof_muscimol_cdf_40, y_lab = F)
  p6 <- plot_delta_omega_lottery_psd(bi_fof_muscimol_cdf, uni_fof_muscimol_cdf, bi_fof_muscimol_cdf_40, uni_fof_muscimol_cdf_40, y_lab = F)
  p7 <- plot_delta_omega_surebet_psd(bi_fof_muscimol_cdf, uni_fof_muscimol_cdf, bi_fof_muscimol_cdf_40, uni_fof_muscimol_cdf_40, y_lab = F)
  psum2 <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + guide_area() + plot_layout(design = layout2, widths = c(1,.8,.8,.8,.8,.8))
  ggsave("~/Desktop/brms_ppc.png", plot = psum2, width = 2.5, units = "in", height = 0.9, scale = 3.5)
  ggsave("~/Desktop/brms_ppc.pdf", plot = psum2, width = 2.5, units = "in", height = 0.9, scale = 3.5)
}

###


# plot the posterior of brms model
plot_posterior_dose_effect <- function(model, x_lab_on = TRUE){
  p <- mcmc_areas(model, regex_pars = 'b_.*_dose', prob=.7, prob_outer = 0.9999, point_est = 'mean') +
    vline_0(size = .2, color = 'gray') +
    #vline_at(model, colMeans, size=.2, color='gray') +
    theme_classic(base_size = 15) +
    scale_x_continuous(limits=c(-4,2), breaks = c(-4, 0, 2), labels = c('-4', '0', '2')) +
    scale_y_discrete(limits = c('b_omega2_dose','b_omega1_dose','b_noise_dose','b_logrho_dose'),
                     labels = c(
                       'b_omega2_dose' = expression(paste('inv_logit(', omega, '2)')),
                       'b_omega1_dose'= expression(paste('inv_logit(', omega, '1)')),
                       'b_noise_dose' = 'noise',
                       'b_logrho_dose' = expression(paste('log(', rho, ')'))),
                     expand = expansion(add = c(.3, .3)))
  if(x_lab_on){
    p <- p + xlab('Posterior density')
  }
  return(p)
}

## plot the posteriors of brms model seperately, make the logrho and lognoise one plot and the logit_omega another plot
plot_posterior_rho_dose_effect <- function(model, x_lab_on = TRUE, y_tick_label = TRUE){
  p <- mcmc_areas(model, regex_pars = 'logrho_dose|noise_dose', prob=.7, prob_outer = 0.9999, point_est = 'mean') +
    vline_0(size = .2, color = 'gray') +
    #vline_at(model, colMeans, size=.2, color='gray') +
    theme_classic(base_size = 15) +
    scale_x_continuous(limits=c(-1.2, 0.6), breaks = c(-1.2, 0, 0.6), labels = c('-1.2', '0', '0.6')) +
    scale_y_discrete(limits = c('b_noise_dose','b_logrho_dose'),
                     labels = c(
                       'b_noise_dose' = expression(paste(Delta, sigma)),
                       'b_logrho_dose' = expression(paste(Delta ,'log(', rho, ')'))),
                     expand = expansion(add = c(.3, .3)))
  if(x_lab_on){
    p <- p + xlab('Posterior density')
  }
  if (!y_tick_label){
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}

plot_posterior_omega_dose_effect <- function(model, x_lab_on = TRUE, y_tick_label = TRUE){
  p <- mcmc_areas(model, regex_pars = 'omega1_dose|omega2_dose', prob=.7, prob_outer = 0.9999, point_est = 'mean') +
    vline_0(size = .2, color = 'gray') +
    #vline_at(model, colMeans, size=.2, color='gray') +
    theme_classic(base_size = 15) +
    scale_x_continuous(limits=c(-6, 3), breaks = c(-6, 0, 3), labels = c('-6', '0', '3')) +
    scale_y_discrete(limits = c('b_omega2_dose','b_omega1_dose'),
                     labels = c(
                       'b_omega2_dose' = expression(paste(Delta, f, '(',omega, '2)')),
                       'b_omega1_dose'= expression(paste(Delta, f, '(',omega, '1)'))),
                     expand = expansion(add = c(.3, .3)))
  if(x_lab_on){
    p <- p + xlab('Posterior density')
  }
  if (!y_tick_label){
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}

## plot the posteriors of brms model seperately, covert rho noise and omega back, make the rho and noise one plot and the omega another plot
plot_posterior_rho_dose_effect_transformed <- function(model, x_lab_on = TRUE, y_tick_label = TRUE){
  ps <- posterior_samples(model, pars = '^b', as.matrix = T)
  ps <- as.data.frame(ps) %>% mutate(delta_rho = exp(b_logrho_Intercept + b_logrho_dose) - exp(b_logrho_Intercept),
                                     delta_noise = b_noise_dose) %>%
    select(c('delta_rho', 'delta_noise')) %>% as.mcmc()

  p <- mcmc_areas(ps, prob=.7, prob_outer = 0.9999, point_est = 'mean') +
    vline_0(size = .2, color = 'gray') +
    #vline_at(model, colMeans, size=.2, color='gray') +
    theme_classic(base_size = 15) +
    scale_x_continuous(limits=c(-0.8, 0.8), breaks = c(-0.8, 0, 0.8), labels = c('-0.8', '0', '0.8')) +
    scale_y_discrete(limits = c('delta_noise','delta_rho'),
                     labels = c(
                       'delta_noise' = expression(paste(Delta, sigma)),
                       'delta_rho' = expression(paste(Delta, rho))),
                     expand = expansion(add = c(.3, .3)))
  if(x_lab_on){
    p <- p + xlab('Posterior density')
  }
  if (!y_tick_label){
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}

plot_posterior_omega_dose_effect_transformed <- function(model, x_lab_on = TRUE, y_tick_label = TRUE){
  ps <- posterior_samples(model, pars = '^b', as.matrix = T)
  ps <- as.data.frame(ps) %>% mutate(delta_omega_rational = inv_logit(b_omega1_Intercept + b_omega1_dose) - inv_logit(b_omega1_Intercept),
                                     delta_omega_lottery = (1- inv_logit(b_omega1_Intercept + b_omega1_dose))*inv_logit(b_omega2_Intercept + b_omega2_dose) -
                                       (1- inv_logit(b_omega1_Intercept))*inv_logit(b_omega2_Intercept)) %>%
    select(c('delta_omega_rational', 'delta_omega_lottery')) %>% as.mcmc()

  p <- mcmc_areas(ps, prob=.7, prob_outer = 0.9999, point_est = 'mean') +
    vline_0(size = .2, color = 'gray') +
    #vline_at(model, colMeans, size=.2, color='gray') +
    theme_classic(base_size = 15) +
    scale_x_continuous(limits=c(-0.3, 0.3), breaks = c(-0.3, 0, 0.3), labels = c('-0.3', '0', '0.3')) +
    scale_y_discrete(limits = c('delta_omega_lottery','delta_omega_rational'),
                     labels = c(
                       'delta_omega_lottery' = expression(paste(Delta, omega['lottey'])),
                       'delta_omega_rational'= expression(paste(Delta, omega['rational']))),
                     expand = expansion(add = c(.3, .3)))
  if(x_lab_on){
    p <- p + xlab('Posterior density')
  }
  if (!y_tick_label){
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}



## generate csv from brms tabel
get_individual_pars_int <- function(model){
  model_name <- deparse(substitute(model))

  pars <- model %>% spread_draws(b_logrho_Intercept, b_logrho_dose, r_subjid__logrho[subjid,],
                                 b_noise_Intercept, b_noise_dose, r_subjid__noise[subjid,],
                                 b_omega1_Intercept, b_omega1_dose, r_subjid__omega1[subjid,],
                                 b_omega2_Intercept, b_omega2_dose, r_subjid__omega2[subjid,]) %>%
    mutate(logrho = b_logrho_Intercept + r_subjid__logrho, logrho_treat = b_logrho_Intercept + r_subjid__logrho + b_logrho_dose,
           lognoise = b_noise_Intercept + r_subjid__noise, lognoise_treat = b_noise_Intercept + r_subjid__noise + b_noise_dose,
           omega1 = b_omega1_Intercept + r_subjid__omega1, omega1_treat = b_omega1_Intercept + r_subjid__omega1 + b_omega1_dose,
           omega2 = b_omega2_Intercept + r_subjid__omega2, omega2_treat = b_omega2_Intercept + r_subjid__omega2 + b_omega2_dose,
           rho_control = exp(logrho), rho_treat = exp(logrho_treat),
           noise_control = exp(lognoise), noise_treat = exp(lognoise_treat),
           omega_rational_control = inv_logit(omega1), omega_rational_treat = inv_logit(omega1_treat),
           omega_lottery_control = (1-inv_logit(omega1))*inv_logit(omega2), omega_lottery_treat = (1-inv_logit(omega1_treat))*inv_logit(omega2_treat),
           omega_surebet_control = (1-inv_logit(omega1))*(1-inv_logit(omega2)), omega_surebet_treat = (1-inv_logit(omega1_treat))*(1-inv_logit(omega2_treat)),
           delta_rho = rho_treat - rho_control,
           delta_noise = noise_treat - noise_control,
           delta_omega_rational = omega_rational_treat - omega_rational_control,
           delta_omega_lottery = omega_lottery_treat - omega_lottery_control,
           delta_omega_surebet = omega_surebet_treat - omega_surebet_control) %>%
    mean_qi() %>%
    select(1,62:106)
  return(list(pars = pars, model_name = model_name))
}

generate_individual_pars_perturb_csv <- function(){
  bi_fof_opto <- get_individual_pars_int(bi_opto_cdf)$pars %>% mutate(type='bi_fof_opto')
  uni_fof_opto <- get_individual_pars_int(uni_opto_cdf)$pars %>% mutate(type='uni_fof_opto')
  bi_fof_muscimol <- get_individual_pars_int(bi_fof_muscimol_cdf)$pars %>% mutate(type='bi_fof_muscimol')
  uni_fof_muscimol <- get_individual_pars_int(uni_fof_muscimol_cdf)$pars %>% mutate(type='uni_fof_muscimol')
  bi_ppc_muscimol <- get_individual_pars_int(bi_ppc_muscimol_cdf)$pars %>% mutate(type='bi_ppc_muscimol')
  uni_ppc_muscimol <- get_individual_pars_int(uni_ppc_muscimol_cdf)$pars %>% mutate(type='uni_ppc_muscimol')
  individual_pars <- rbind(bi_fof_opto,uni_fof_opto,bi_fof_muscimol,uni_fof_muscimol,bi_ppc_muscimol,uni_ppc_muscimol) %>%
    mutate(rho = sprintf("%.2f [%.2f, %.2f]",rho_control, rho_control.lower, rho_control.upper),
           rho_pert = sprintf("%.2f [%.2f, %.2f]",rho_treat, rho_treat.lower, rho_treat.upper),
           noise = sprintf("%.2f [%.2f, %.2f]",noise_control, noise_control.lower, noise_control.upper),
           noise_pert = sprintf("%.2f [%.2f, %.2f]",noise_treat, noise_treat.lower, noise_treat.upper),
           omega_rational = sprintf("%.2f [%.2f, %.2f]",omega_rational_control, omega_rational_control.lower, omega_rational_control.upper),
           omega_rational_pert = sprintf("%.2f [%.2f, %.2f]",omega_rational_treat, omega_rational_treat.lower, omega_rational_treat.upper),
           omega_lottery = sprintf("%.2f [%.2f, %.2f]",omega_lottery_control, omega_lottery_control.lower, omega_lottery_control.upper),
           omega_lottery_pert = sprintf("%.2f [%.2f, %.2f]",omega_lottery_treat, omega_lottery_treat.lower, omega_lottery_treat.upper),
           omega_surebet = sprintf("%.2f [%.2f, %.2f]",omega_surebet_control, omega_surebet_control.lower, omega_surebet_control.upper),
           omega_surebet_pert = sprintf("%.2f [%.2f, %.2f]",omega_surebet_treat, omega_surebet_treat.lower, omega_surebet_treat.upper),
           delta_rhos = sprintf("%.2f [%.2f, %.2f]",delta_rho, delta_rho.lower, delta_rho.upper),
           delta_noises = sprintf("%.2f [%.2f, %.2f]",delta_noise, delta_noise.lower, delta_noise.upper),
           delta_omega_rationals = sprintf("%.2f [%.2f, %.2f]",delta_omega_rational, delta_omega_rational.lower, delta_omega_rational.upper),
           delta_omega_lotterys = sprintf("%.2f [%.2f, %.2f]",delta_omega_lottery, delta_omega_lottery.lower, delta_omega_lottery.upper),
           delta_omega_surebets = sprintf("%.2f [%.2f, %.2f]",delta_omega_surebet, delta_omega_surebet.lower, delta_omega_surebet.upper)) %>%
    select(1,47:48,50,52,54,56,58:62)
  return(individual_pars)
}

generate_individual_pars_ephys_csv <- function(){
  ephys <- get_individual_pars_int(ephys_cdf)$pars %>%
    mutate(rho = sprintf("%.2f [%.2f, %.2f]",rho_control, rho_control.lower, rho_control.upper),
           rho_pert = sprintf("%.2f [%.2f, %.2f]",rho_treat, rho_treat.lower, rho_treat.upper),
           noise = sprintf("%.2f [%.2f, %.2f]",noise_control, noise_control.lower, noise_control.upper),
           noise_pert = sprintf("%.2f [%.2f, %.2f]",noise_treat, noise_treat.lower, noise_treat.upper),
           omega_rational = sprintf("%.2f [%.2f, %.2f]",omega_rational_control, omega_rational_control.lower, omega_rational_control.upper),
           omega_rational_pert = sprintf("%.2f [%.2f, %.2f]",omega_rational_treat, omega_rational_treat.lower, omega_rational_treat.upper),
           omega_lottery = sprintf("%.2f [%.2f, %.2f]",omega_lottery_control, omega_lottery_control.lower, omega_lottery_control.upper),
           omega_lottery_pert = sprintf("%.2f [%.2f, %.2f]",omega_lottery_treat, omega_lottery_treat.lower, omega_lottery_treat.upper),
           omega_surebet = sprintf("%.2f [%.2f, %.2f]",omega_surebet_control, omega_surebet_control.lower, omega_surebet_control.upper),
           omega_surebet_pert = sprintf("%.2f [%.2f, %.2f]",omega_surebet_treat, omega_surebet_treat.lower, omega_surebet_treat.upper),
           delta_rhos = sprintf("%.2f [%.2f, %.2f]",delta_rho, delta_rho.lower, delta_rho.upper),
           delta_noises = sprintf("%.2f [%.2f, %.2f]",delta_noise, delta_noise.lower, delta_noise.upper),
           delta_omega_rationals = sprintf("%.2f [%.2f, %.2f]",delta_omega_rational, delta_omega_rational.lower, delta_omega_rational.upper),
           delta_omega_lotterys = sprintf("%.2f [%.2f, %.2f]",delta_omega_lottery, delta_omega_lottery.lower, delta_omega_lottery.upper),
           delta_omega_surebets = sprintf("%.2f [%.2f, %.2f]",delta_omega_surebet, delta_omega_surebet.lower, delta_omega_surebet.upper)) %>%
    select(1,47,49,51,53,55,57:61)
  return(ephys)
}
