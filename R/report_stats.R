# list of functions used to report stats
library(car)

bi_fof_less_trials = function(region = 'bi_fof'){
  # rats did fewer trials on 300ng bilateral FOF infusion compared to control
  df1 = read.csv('csv/figure_behavior_population.csv') %>% group_by(subjid, sessid) %>% tally()
  df2 = read.csv(sprintf('csv/figure_infusion_%s.csv', region)) %>% 
    filter(dosage == 0.3) %>% filter(as.Date(trialtime) < '2020-11-19') %>% 
    group_by(subjid, sessid) %>% tally()
  return(t.test(df1$n, df2$n))
}

stat_bi_infusion = function(region = 'FOF'){
  # report stats for bilateral infusions compared to control
  if (region == 'FOF'){
    load('fits/Bilateral FOF_glmm.RData')
    df = read.csv('csv/figure_infusion_bi_fof.csv') %>% filter(dosage %in% c(0.075, 0.3))
    control_df = read.csv('csv/figure_behavior_population.csv')
  } else{
    load('fits/Bilateral PPC_glmm.RData')
    df = read.csv('csv/figure_infusion_bi_ppc.csv') %>% filter(dosage == 0.3)
    control_df = read.csv('csv/figure_behavior_population.csv') %>% filter(subjid != 2155)
  }
  df = df %>% filter(as.Date(trialtime) < '2020-11-19')
  df = rbind(control_df, df) 
  # report GLMM stats
  cf = coef(summary(m1))
  cat('[GLMM] \n')
  cat(sprintf('dosage: %.2f +- %.2f, p = %.3f \ndelta_ev:dosage: %.2f +- %.2f, p = %.3f \n',
              cf['dosage', 'Estimate'], cf['dosage', 'Std. Error'], cf['dosage', 'Pr(>|z|)'], 
              cf['delta_ev:dosage', 'Estimate'], cf['delta_ev:dosage', 'Std. Error'], cf['delta_ev:dosage', 'Pr(>|z|)']))
  
  # compare the change in indifference points
  pred_df = data_grid(df, delta_ev = seq_range(delta_ev, by = 0.05), dosage = dosage, subjid = subjid)
  pred_df$pred = predict(m1, pred_df)
  indiff_df = pred_df %>% filter(pred < 0.51 & pred > 0.49) %>% # get approximate 0.5
    group_by(subjid, dosage) %>% summarise(delta_ev = mean(delta_ev))
  control = indiff_df %>% filter(dosage == 0) %>% pull(delta_ev) # list of indifference points in control
  inactivation = indiff_df %>% filter(dosage == 0.3) %>% pull(delta_ev) # list of indifference points in 0.3ng
  t = t.test(control, inactivation)
  cat('\n')
  cat('[Indifference points] \n')
  cat(sprintf('Average control indifference point is %.2f +- %.2f, \nAverage 300ng indifference point is %.2f +- %.2f, \nT{%.0f} = %.2f, p = %.2f \n',
              mean(control), se(control), mean(inactivation), se(inactivation), length(control), t$statistic, t$p.value))
  
  # report infusion effect in individual animals
  cat('\n')
  cat('[Individual significance] \n')
  subj_list = unique(df$subjid)
  for (subj in subj_list){
    subj_df = df %>% filter(subjid == subj)
    m1 = glm('choice ~ delta_ev*dosage', subj_df, family = binomial)
    cf = coef(summary(m1))
    cat(sprintf('%d \ndosage: %.2f +- %.2f, p = %.3f \ndelta_ev:dosage: %.2f +- %.2f, p = %.3f \n', subj,
                cf['dosage', 'Estimate'], cf['dosage', 'Std. Error'], cf['dosage', 'Pr(>|z|)'], 
                cf['delta_ev:dosage', 'Estimate'], cf['delta_ev:dosage', 'Std. Error'], cf['delta_ev:dosage', 'Pr(>|z|)']))
  }
}

stat_uni_infusion = function(region = 'FOF'){
  # report stats for unilateral infusions, NOT COMPLETE
  if (region == 'FOF'){
    load('fits/Unilateral FOF_glmm.RData')
    df = read.csv('csv/figure_infusion_uni_fof.csv') %>% filter(dosage != 0)
  } else if (region == 'PPC'){
    load('fits/Unilateral PPC_glmm.RData')
    df = read.csv('csv/figure_infusion_uni_ppc.csv') %>% filter(dosage != 0)
  }
  control_df = read.csv('csv/figure_behavior_population.csv')
  df = df %>% filter(as.Date(trialtime) < '2020-11-19')
  df = rbind(control_df, df) 
  df$infusion_side = factor(df$infusion_side)
  pred_df = data_grid(df, pro_right_delta_ev = seq_range(pro_right_delta_ev, by = 1), infusion_side = infusion_side, subjid = subjid)
  # report GLMM stats
  cf = coef(summary(m1))
  cat('[GLMM] \n')
  cat(sprintf('infusion_sideL: %.2f +- %.2f, p = %.3f \npro_right_delta_ev:infusion_sideL: %.2f +- %.2f, p = %.3f \ninfusion_sideR: %.2f +- %.2f, p = %.3f \npro_right_delta_ev:infusion_sideR: %.2f +- %.2f, p = %.3f \n ',
              cf['infusion_sideL', 'Estimate'], cf['infusion_sideL', 'Std. Error'], cf['infusion_sideL', 'Pr(>|z|)'],
              cf['pro_right_delta_ev:infusion_sideL', 'Estimate'], cf['pro_right_delta_ev:infusion_sideL', 'Std. Error'], cf['pro_right_delta_ev:infusion_sideL', 'Pr(>|z|)'],
              cf['infusion_sideR', 'Estimate'], cf['infusion_sideR', 'Std. Error'], cf['infusion_sideR', 'Pr(>|z|)'],
              cf['pro_right_delta_ev:infusion_sideR', 'Estimate'], cf['pro_right_delta_ev:infusion_sideR', 'Std. Error'], cf['pro_right_delta_ev:infusion_sideR', 'Pr(>|z|)']))
  
  # report infusion effect in individual animals
  cat('\n')
  cat('[Individual significance] \n')
  subj_list = unique(df$subjid)
  for (subj in subj_list){
    subj_df = df %>% filter(subjid == subj)
    m1 = glm('choose_right ~ pro_right_delta_ev*infusion_side', subj_df, family = binomial)
    cf = coef(summary(m1))
    cat(sprintf('%d \ninfusion_sideL: %.2f +- %.2f, p = %.3f \npro_right_delta_ev:infusion_sideL: %.2f +- %.2f, p = %.3f \ninfusion_sideR: %.2f +- %.2f, p = %.3f \npro_right_delta_ev:infusion_sideR: %.2f +- %.2f, p = %.3f \n', subj,
                cf['infusion_sideL', 'Estimate'], cf['infusion_sideL', 'Std. Error'], cf['infusion_sideL', 'Pr(>|z|)'], 
                cf['pro_right_delta_ev:infusion_sideL', 'Estimate'], cf['pro_right_delta_ev:infusion_sideL', 'Std. Error'], cf['pro_right_delta_ev:infusion_sideL', 'Pr(>|z|)'],
                cf['infusion_sideR', 'Estimate'], cf['infusion_sideR', 'Std. Error'], cf['infusion_sideR', 'Pr(>|z|)'], 
                cf['pro_right_delta_ev:infusion_sideR', 'Estimate'], cf['pro_right_delta_ev:infusion_sideR', 'Std. Error'], cf['pro_right_delta_ev:infusion_sideR', 'Pr(>|z|)']))
  }
}

stat_RT_bi_infusion = function(region = 'FOF'){
  # report stats regarding reaction in bilateral infusions
  if (region == 'FOF'){
    load('fits/bi_fof_RT_glmm.RData')
    df = read.csv('csv/figure_infusion_bi_fof.csv') %>% filter(dosage %in% c(0.075, 0.3))
    control_df = read.csv('csv/figure_behavior_population.csv')
  } else{
    df = read.csv('csv/figure_infusion_bi_ppc.csv') %>% filter(dosage == 0.3)
    control_df = read.csv('csv/figure_behavior_population.csv') %>% filter(subjid != 2155)
  }
  df = df %>% filter(as.Date(trialtime) < '2020-11-19')
  df = rbind(control_df, df) 
  
  # report LMM stats
  cf = coef(summary(m1))
  p = Anova(m1)
  cat('[LMM] \n')
  cat(sprintf('dosage: %.2f +- %.2f, p = %.3f \nchoice: %.2f +- %.2f, p = %.3f \ndelta_ev:dosage: %.2f +- %.2f, p = %.3f \ndelta_ev:choice: %.2f +- %.2f, p = %.3f \n',
              cf['dosage', 'Estimate'], cf['dosage', 'Std. Error'], p$`Pr(>Chisq)`[2], 
              cf['choice', 'Estimate'], cf['choice', 'Std. Error'], p$`Pr(>Chisq)`[3], 
              cf['delta_ev:dosage', 'Estimate'], cf['delta_ev:dosage', 'Std. Error'], p$`Pr(>Chisq)`[4],
              cf['delta_ev:choice', 'Estimate'], cf['delta_ev:choice', 'Std. Error'], p$`Pr(>Chisq)`[5])) 
  
  # report RT effect in individual animals
  cat('\n')
  cat('[Individual significance] \n')
  subj_list = unique(df$subjid)
  for (subj in subj_list){
    subj_df = df %>% filter(subjid == subj)
    m1 = lm('log_RT ~ delta_ev*dosage*choice', subj_df)
    cf = coef(summary(m1))
    cat(sprintf('%d \ndosage: %.2f +- %.2f, p = %.3f \ndelta_ev:dosage: %.2f +- %.2f, p = %.3f \n', subj,
                cf['dosage', 'Estimate'], cf['dosage', 'Std. Error'], cf['dosage', 'Pr(>|t|)'], 
                cf['delta_ev:dosage', 'Estimate'], cf['delta_ev:dosage', 'Std. Error'], cf['delta_ev:dosage', 'Pr(>|t|)']))
  }
}

stat_RT_uni_infusion = function(region = 'FOF', forced = FALSE){
  # report stats regarding reaction in unilateral infusions
  if (region == 'FOF'){
    df = read.csv('csv/figure_infusion_uni_fof.csv') %>% filter(dosage == 0.3)
    control_df = read.csv('csv/figure_behavior_population.csv')
  } else{
    df = read.csv('csv/figure_infusion_uni_ppc.csv') %>% filter(dosage == 0.3)
    control_df = read.csv('csv/figure_behavior_population.csv') %>% filter(subjid != 2155)
  }
  df = df %>% filter(as.Date(trialtime) < '2020-11-19') %>% 
    filter(RT < 3) %>% filter(!is.na(RT))
  df = rbind(control_df, df) 
  m1 = lmer('log_RT ~ pro_right_delta_ev*infusion_side*choice + (pro_right_delta_ev*infusion_side*choice | subjid)', df %>% filter(infusion_side != 'R'))
  m2 = lmer('log_RT ~ pro_right_delta_ev*infusion_side*choice + (pro_right_delta_ev*infusion_side*choice | subjid)', df %>% filter(infusion_side != 'L'))
  
  # report LMM stats
  models = list(m1, m2)
  for (i in c(1,2)){
    cf = coef(summary(models[[i]]))
    p = Anova(models[[i]])
    cat(sprintf('[LMM] for %s \n', c('L', 'R')[i]))
    cat(sprintf('pro_right_delta_ev:infusion: %.3f +- %.3f, p = %.3f \ninfusion_side: %.3f +- %.3f, p = %.3f \nchoice: %.3f +- %.3f, p = %.3f \npro_right_delta_ev:infusion_side: %.3f +- %.3f, p = %.3f \ninfusion_side:choice: %.3f +- %.3f, p = %.3f \n',
                cf[2, 'Estimate'], cf[2, 'Std. Error'], p$`Pr(>Chisq)`[1], 
                cf[3, 'Estimate'], cf[3, 'Std. Error'], p$`Pr(>Chisq)`[2],
                cf[4, 'Estimate'], cf[4, 'Std. Error'], p$`Pr(>Chisq)`[3], 
                cf[5, 'Estimate'], cf[5, 'Std. Error'], p$`Pr(>Chisq)`[4],
                cf[7, 'Estimate'], cf[7, 'Std. Error'], p$`Pr(>Chisq)`[6]))
    
  }
  
  # report RT effect in individual animals
  cat('\n')
  cat('[Individual significance] \n')
  subj_list = unique(df$subjid)
  for (subj in infusion_animals){
    subj_df = df %>% filter(subjid == subj)
    m1 = lm('log_RT ~ pro_right_delta_ev*infusion_side*choice', subj_df %>% filter(infusion_side != 'R')) # L and control
    m2 = lm('log_RT ~ pro_right_delta_ev*infusion_side*choice', subj_df %>% filter(infusion_side != 'L')) # R and control
    models = list(m1, m2)
    for (i in c(1,2)){
      cf = coef(summary(models[[i]]))
      p = Anova(models[[i]])
      cat(sprintf('%d, %s \n', subj, c('L', 'R')[i]))
      cat(sprintf('pro_right_delta_ev:infusion: %.3f +- %.3f, p = %.3f \ninfusion_side: %.3f +- %.3f, p = %.3f \nchoice: %.3f +- %.3f, p = %.3f \npro_right_delta_ev:infusion_side: %.3f +- %.3f, p = %.3f \ninfusion_side:choice: %.3f +- %.3f, p = %.3f \n',
                  cf[2, 'Estimate'], cf[2, 'Std. Error'], cf[2, 'Pr(>|t|)'], 
                  cf[3, 'Estimate'], cf[3, 'Std. Error'], cf[3, 'Pr(>|t|)'],
                  cf[4, 'Estimate'], cf[4, 'Std. Error'], cf[4, 'Pr(>|t|)'], 
                  cf[5, 'Estimate'], cf[5, 'Std. Error'], cf[5, 'Pr(>|t|)'],
                  cf[7, 'Estimate'], cf[7, 'Std. Error'], cf[7, 'Pr(>|t|)']))
    }
  }
}

stat_free_choice = function(choice_type = 'side'){
  # GLMM fits using only infusions opposite to the animal's preferred side
  data_df = read.csv(sprintf('csv/figure_control_%s_df.csv', choice_type))
  bias_df = read.csv('csv/figure_control_bias_PPC.csv')
  data_df = data_df %>% mutate(chose_ipsi = case_when(choice_side == 'left' & infusion_side == 'L' ~ 1,
                                                      choice_side == 'right' & infusion_side == 'L' ~ 0,
                                                      choice_side == 'left' & infusion_side == 'R' ~ 0,
                                                      choice_side == 'right' & infusion_side == 'R' ~ 1))
  data_df$infusion_bino = factor(data_df$infusion_bino)
  df = merge(data_df, bias_df, by = c('subjid', 'infusion_side'))
  m1 = glmer('chose_ipsi ~ infusion_bino + (infusion_bino | subjid)', df, family = 'binomial')
  cf = coef(summary(m1))
  cat(sprintf('infusion_bino: %.3f +- %.3f, p = %.3f \n',
              cf[2, 'Estimate'], cf[2, 'Std. Error'], cf[2, 'Pr(>|z|)'])) 
}

permutation_test = function(model = 'rho-and-omega', type = 'individual'){
  # do permutation test for individual samples or meta-rat samples for each parameter
  control_df = read.csv(sprintf('csv/all_%s_fits.csv', model)) %>% filter(dataset == 'Control')
  fof_df = read.csv(sprintf('csv/all_%s_fits.csv', model)) %>% filter(dataset == 'Bilateral-FOF')
  ppc_df = read.csv(sprintf('csv/all_%s_fits.csv', model)) %>% filter(dataset == 'Bilateral-PPC')
  if (type == 'individual'){
    for (i in 1:8){
      for (param in params){
        fof_deviation = fof_df[[sprintf('%s%d', param, i)]] - control_df[[sprintf('%s%d', param, i)]]
        ppc_deviation = ppc_df[[sprintf('%s%d', param, i)]] - control_df[[sprintf('%s%d', param, i)]]
        a = permutation.test(fof_deviation, rep(0, 1000))
        b = permutation.test(ppc_deviation, rep(0, 1000))
        cat(sprintf('%d \n', infusion_animals[i]))
        cat(sprintf('%s Bilateral FOF - Control, p = %.6f %s \n', param, a, ifelse(a < 0.05, ifelse(a < 0.01, '**', '*'), '')))
        cat(sprintf('%s Bilateral PPC - Control, p = %.6f %s \n', param, b, ifelse(b < 0.05, ifelse(b < 0.01, '**', '*'), '')))
        cat('-------------------------------------------------------------\n')
      }
    }
  } else if (type == 'deviation'){
    for (param in params){
      control = control_df %>% select(contains(param)) %>% pivot_longer(cols = contains(param)) %>% pull(value)
      fof = fof_df %>% select(contains(param)) %>% pivot_longer(cols = contains(param)) %>% pull(value)
      ppc = ppc_df %>% select(contains(param)) %>% pivot_longer(cols = contains(param)) %>% pull(value)
      fof_deviation = fof - control
      ppc_deviation = ppc - control
      a = permutation.test(fof_deviation, n = 10000)
      b = permutation.test(ppc_deviation, n = 10000)
      cat(sprintf('%s Bilateral FOF - Control, p = %.3f \n', param, a))
      cat(sprintf('%s Bilateral PPC - Control, p = %.3f \n', param, b))
      cat('-------------------------------------------------------------\n')
    }
  }
}

rho_and_RT = function(){
  # any correlations between change in RT and change in rho?
  # get the median deviation in rho
  df_control = read.csv('csv/all_rho-and-omega_fits.csv') %>% filter(dataset == 'Control')
  df_fof = read.csv('csv/all_rho-and-omega_fits.csv') %>% filter(dataset == 'Bilateral-FOF')
  df_deviation = df_fof - df_control
  x = df_deviation %>% select(contains('rho')) %>% summarise(across(c(1:8), median))
  # get the beta main effects of RT
  df = read.csv('csv/figure_infusion_bi_fof.csv') %>%
    filter(as.Date(trialtime) < '2020-11-19') %>%
    filter(dosage != 0.15 & dosage != 0) %>% filter(!is.na(RT))
  control_df = read.csv('csv/figure_behavior_population.csv')
  df = rbind(df, control_df)
  effect_df = data.frame()
  for (subj in infusion_animals){
    subj_df = df %>% filter(subjid == subj)
    m1 = lm('log_RT ~ delta_ev*dosage*choice', subj_df) 
    cf = coef(summary(m1))
    effect_df = rbind(effect_df, data.frame(subjid = subj, main_effect = cf['dosage', 'Estimate'], 
                                            main_p = cf['dosage', 'Pr(>|t|)'], interaction = cf['delta_ev:dosage', 'Estimate'], 
                                            interaction_p = cf['delta_ev:dosage', 'Pr(>|t|)']))
  }
  effect_df$median_rho = as.numeric(x)
  return(cor.test(effect_df$interaction, effect_df$median_rho))
}

delta_p_value = function(model = 'rho-and-omega'){
  # p value is given by the fraction of samples having the opposite sign as the median sample. 
  df = read.csv(sprintf('csv/%s_delta_fits.csv', model))
  fof_params = c('rho_delta_fof', 'omega_1_delta_fof', 'omega_2_delta_fof')
  ppc_params = c('rho_delta_ppc', 'omega_1_delta_ppc', 'omega_2_delta_ppc')
  for (i in 1:8){
    for (j in 1:3){
      fof_deviation = df[[sprintf('%s%d', fof_params[j], i)]]
      ppc_deviation = df[[sprintf('%s%d', ppc_params[j], i)]]
      a = ifelse(median(fof_deviation) > 0, 
                 sum(fof_deviation < 0) / length(fof_deviation), 
                 sum(fof_deviation > 0) / length(fof_deviation))
      b = ifelse(median(ppc_deviation) > 0, 
                 sum(ppc_deviation < 0) / length(ppc_deviation), 
                 sum(ppc_deviation > 0) / length(ppc_deviation))
      cat(sprintf('%d \n', infusion_animals[i]))
      cat(sprintf('%s, p = %.6f %s \n', fof_params[j], a, ifelse(a < 0.05, ifelse(a < 0.01, '**', '*'), '')))
      cat(sprintf('%s, p = %.6f %s \n', ppc_params[j], b, ifelse(b < 0.05, ifelse(b < 0.01, '**', '*'), '')))
      cat('-------------------------------------------------------------\n')
    }
  }
  
}

get_latex_table = function(model = 'rho-and-omega'){
  # take subject fit parameters and format into the latex table
  control_df = read.csv(sprintf('csv/all_%s_fits.csv', model)) %>% filter(dataset == 'Control')
  fof_df = read.csv(sprintf('csv/all_%s_fits.csv', model)) %>% filter(dataset == 'Bilateral-FOF')
  ppc_df = read.csv(sprintf('csv/all_%s_fits.csv', model)) %>% filter(dataset == 'Bilateral-PPC')
  for (i in 1:8){
    control_rho = control_df[[sprintf('rho%d', i)]]
    control_sigma = control_df[[sprintf('sigma%d', i)]]
    control_rational = control_df[[sprintf('omega_rational%d',i)]]
    control_lott = control_df[[sprintf('omega_lottery%d', i)]]
    control_sb = control_df[[sprintf('omega_surebet%d', i)]]
    
    fof_rho = fof_df[[sprintf('rho%d', i)]]
    fof_sigma = fof_df[[sprintf('sigma%d', i)]]
    fof_rational = fof_df[[sprintf('omega_rational%d', i)]]
    fof_lott = fof_df[[sprintf('omega_lottery%d', i)]]
    fof_sb = fof_df[[sprintf('omega_surebet%d', i)]]
    
    if (i == 4){ # 2155 NA place holders
      ppc_rho = 0
      ppc_sigma = 0
      ppc_rational = 0
      ppc_lott = 0
      ppc_sb = 0
      ppc_inx = 1
    } else{
      ppc_rho = ppc_df[[sprintf('rho%d', i)]]
      ppc_sigma = ppc_df[[sprintf('sigma%d', i)]]
      ppc_rational = ppc_df[[sprintf('omega_rational%d', i)]]
      ppc_lott = ppc_df[[sprintf('omega_lottery%d', i)]]
      ppc_sb = ppc_df[[sprintf('omega_surebet%d', i)]]
      ppc_inx = which.max(ppc_df[[sprintf('lp__%d', i)]])
    }
    
    control_inx = which.max(control_df[[sprintf('lp__%d', i)]]) # for MAP estimation
    fof_inx = which.max(fof_df[[sprintf('lp__%d', i)]])
    
    cat(sprintf('\\multirow{3}{4em}{\\textbf{%d}} & Control & %.2f [%.2f, %.2f] & %.2f [%.2f, %.2f] & %.2f [%.2f, %.2f] & %.2f [%.2f, %.2f] & %.2f [%.2f, %.2f] \\\\
                & FOF & %.2f [%.2f, %.2f] & %.2f [%.2f, %.2f] & %.2f [%.2f, %.2f] & %.2f [%.2f, %.2f] & %.2f [%.2f, %.2f] \\\\
                & PPC & %.2f [%.2f, %.2f] & %.2f [%.2f, %.2f] & %.2f [%.2f, %.2f] & %.2f [%.2f, %.2f] & %.2f [%.2f, %.2f] \\\\', infusion_animals[i], 
                control_rho[control_inx],  quantile(control_rho, 0.025), quantile(control_rho, 0.975), 
                control_sigma[control_inx],  quantile(control_sigma, 0.025), quantile(control_sigma, 0.975), 
                control_rational[control_inx],  quantile(control_rational, 0.025), quantile(control_rational, 0.975), 
                control_lott[control_inx],  quantile(control_lott, 0.025), quantile(control_lott, 0.975), 
                control_sb[control_inx],  quantile(control_sb, 0.025), quantile(control_sb, 0.975), 
                fof_rho[fof_inx],  quantile(fof_rho, 0.025), quantile(fof_rho, 0.975), 
                fof_sigma[fof_inx],  quantile(fof_sigma, 0.025), quantile(fof_sigma, 0.975), 
                fof_rational[fof_inx],  quantile(fof_rational, 0.025), quantile(fof_rational, 0.975), 
                fof_lott[fof_inx],  quantile(fof_lott, 0.025), quantile(fof_lott, 0.975), 
                fof_sb[fof_inx],  quantile(fof_sb, 0.025), quantile(fof_sb, 0.975), 
                ppc_rho[ppc_inx],  quantile(ppc_rho, 0.025), quantile(ppc_rho, 0.975), 
                ppc_sigma[ppc_inx],  quantile(ppc_sigma, 0.025), quantile(ppc_sigma, 0.975), 
                ppc_rational[ppc_inx],  quantile(ppc_rational, 0.025), quantile(ppc_rational, 0.975), 
                ppc_lott[ppc_inx],  quantile(ppc_lott, 0.025), quantile(ppc_lott, 0.975), 
                ppc_sb[ppc_inx],  quantile(ppc_sb, 0.025), quantile(ppc_sb, 0.975)), 
        '\n \\midrule \n')
  }
}