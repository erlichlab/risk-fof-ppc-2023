# functions for the Extended data plots and supplementary plots of infusion data

#' infusion_timeline
#'
#'
#' @param animal: a subjid in the list of infusion_animals
#'
#' @return ggplot object
plot_infusion_timeline = function(animal = 2156){
  # timeline plot of all the experimental manipulations performed on an animal
  timeline_df = read.csv(file.path(git_root,"csv","infusion_timeline.csv"))
  df = timeline_df %>% filter(subjid == animal)
  df$area = factor(df$area, c('control' ,'FOF', 'PPC' ,'Both'))
  df$side = factor(df$side, c('control' ,'Left', 'Right', 'Both'))
  p = ggplot(df, aes(x = days_post_surgery, y = mean_chose_lott*100, color = side, shape = area)) + theme_classic(base_size = 13) +
    geom_point(aes(size = !is.na(event))) + ylab('% Chose Lottery') +
    annotate("text", label = animal, x = min(df$days_post_surgery)*1.01, y = 1) +
    scale_shape_manual(values = c(16, 18, 15, 17)) + # round, diamond, square, triangle
    scale_color_manual(values = c('black', 'mediumseagreen', 'brown2', 'dodgerblue')) +
    scale_size_manual(values = c(2, 4)) +
    scale_x_continuous(sec.axis = sec_axis(~., name = 'Days after surgery'),
                       breaks = df$days_post_surgery[!is.na(df$event)],
                       labels = df$event[!is.na(df$event)],
                       limits = c(min(df$days_post_surgery), max(df$days_post_surgery)*1.05)) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
    geom_hline(yintercept = 50, linetype = 'dotted', alpha = 0.4) +
    geom_vline(xintercept = df$days_post_surgery[!is.na(df$sb_change)], linetype = 'solid', color = 'deepskyblue') +
    theme(legend.position = 'none',
          axis.text.x.bottom = element_text(angle = -35, hjust = -0.05),
          axis.title.x.bottom = element_blank(), axis.title.x.top = element_blank())
  return(p)
}

#' Extended_Data_infusion
#'
#'
#' @param region: 'bi_fof', 'bi_ppc', 'uni_fof', 'uni_ppc'
#'
#' @return ggplot object
plot_individual_infusion = function(region = 'bi_fof'){
  # GLMM fits and plots for Bilateral FOF/PPC, Unilateral FOF/PPC
  df_m = read.csv(file.path(git_root,"csv", sprintf('figure_infusion_%s.csv', region))) %>%
    filter(dosage != 0.15 & dosage != 0)
  control_df = read.csv(file.path(git_root,"csv","figure_behavior_population.csv"))
  df = rbind(df_m, control_df)

  if (str_detect(region, 'bi')){
    if (region == 'bi_fof'){colors = c('azure4', 'rosybrown1', 'purple4')} else{ colors = c('azure4', 'gold2')}
    df$infusion_side = factor(df$infusion_side)
    df$infusion_bino = factor(df$infusion_bino)
    # load GLMM or fit GLMM
    fname = file.path(git_root,"fits",sprintf('%s_glmm.RData', region))
    if (file.exists(fname)){
      load(fname)
    } else{
      cat('Fitting GLMM...\n')
      m1 = glmer('choice ~ delta_ev*dosage + (delta_ev*dosage | subjid)', df, family = binomial)
      save(m1, file = fname)
    }
    ind_pred_df = data_grid(df, delta_ev = seq_range(delta_ev, by = 1), dosage = dosage, subjid = subjid)
    ind_pred_df$pred = predict(m1, ind_pred_df, type = 'response', allow.new.levels=TRUE)
    qt = quantile(df$delta_ev, probs = seq(0, 1, 0.25))
    p = ggplot(df) +
      #stat_summary_bin(mapping = aes(x = delta_ev, y = choice, color = dosage), breaks = qt, fun.data = bino, geom = 'pointrange', position = position_dodge(10)) +
      stat_summary_bin(mapping = aes(x = delta_ev, y = choice, color = as.factor(dosage)), fun.data = bino, bins = 7, geom = 'pointrange', position = position_dodge(10)) +
      geom_line(ind_pred_df, mapping = aes(x = delta_ev, y = pred, color = as.factor(dosage)), alpha = 0.5, size = 1) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      labs(fill = 'Dose(µg)', color = 'Dose(µg)') +
      scale_x_continuous(limits = c(-100, 300), breaks = c(-100, 0, 200)) +
      xlab(expression(EV[lottery]-EV[surebet])) + ylab("P(Chose Lottery)")
    if (region == 'bi_fof'){
      p <- p + ggtitle('Muscimol bilateral FOF')
    }else{
      p <- p + ggtitle('Muscimol bilateral PPC')
    }
  } else if (str_detect(region, 'uni')){
    df$dosage = factor(df$dosage)
    df$infusion_side = factor(df$infusion_side, )
    df$infusion_bino = factor(df$infusion_bino)
    # load GLMM or fit GLMM
    fname = file.path(git_root, "fits", sprintf('%s_glmm.RData', region))
    if (file.exists(fname)){
      load(fname)
    } else{
      cat('Fitting GLMM...\n')
      m1 = glmer('choose_right ~ pro_right_delta_ev*infusion_side + (pro_right_delta_ev*infusion_side | subjid)', df, family = binomial)
      save(m1, file = fname)
    }
    ind_pred_df = data_grid(df, pro_right_delta_ev = seq_range(pro_right_delta_ev, by = 1), infusion_side = infusion_side, subjid = subjid)
    ind_pred_df$pred = predict(m1, ind_pred_df, type = 'response', allow.new.levels=TRUE)
    qt = quantile(df$delta_ev, probs = seq(0, 1, 0.25))
    p = ggplot(df) +
      #stat_summary_bin(mapping = aes(x = pro_right_delta_ev, y = choose_right, color = infusion_side), breaks = qt, fun.data = bino, geom = 'pointrange', position = position_dodge(10)) +
      stat_summary_bin(mapping = aes(x = pro_right_delta_ev, y = choose_right, color = infusion_side), fun.data = bino, bins = 7, geom = 'pointrange', position = position_dodge(10)) +
      geom_line(ind_pred_df, mapping = aes(x = pro_right_delta_ev, y = pred, color = infusion_side), alpha = 0.5, size = 1) +
      scale_color_manual(values = side_colors) + scale_fill_manual(values = side_colors) +
      labs(fill = 'Infusion Side', color = 'Infusion Side') +
      scale_x_continuous(limits = c(-300, 300), breaks = c(-200, 0, 200)) +
      xlab(expression(EV[right]-EV[left])) + ylab("P(Chose Right)")
    if (region == 'uni_fof'){
      p <- p + ggtitle('Muscimol unilateral FOF')
    }else{
      p <- p + ggtitle('Muscimol unilateral PPC')
    }
  }
  # plot
  ntrials_df = df_m %>% group_by(subjid) %>% tally() %>%
    mutate(label = paste0('n=', n), x = max(df$delta_ev)*0.75,  y = 0.1)
  p = p + theme_classic(base_size = BASE_SIZE) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed',  alpha = 0.4) +
    #scale_x_continuous(limits = c(-100, 300), breaks = c(-100, 0, 200)) +
    scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1)) +
    geom_text(data = ntrials_df, mapping = aes(x = x, y = y, label = label)) +
    facet_wrap(~subjid) +
    theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5))
  return(p)
}

plot_individual_infusion_ev_uni = function(region = 'bi_fof'){
  # figureS4 & S5: GLMM fits and plots for Bilateral FOF/PPC, Unilateral FOF/PPC
  df_m = read.csv(file.path(git_root,"csv",sprintf('figure_infusion_%s.csv', region))) %>%
    filter(as.Date(trialtime) < '2020-11-19') %>%
    filter(dosage == 0.3)
  control_df = read.csv(file.path(git_root,"csv","figure_behavior_population.csv"))
  df = rbind(df_m, control_df)

  if (str_detect(region, 'bi')){
    if (region == 'bi_fof'){colors = c('azure4', 'rosybrown1', 'purple4')} else{ colors = c('azure4', 'gold2')}
    df$infusion_side = factor(df$infusion_side)
    df$infusion_bino = factor(df$infusion_bino)
    # load GLMM or fit GLMM
    fname = file.path(git_root,"csv",sprintf("%s_glmm.RData", region))
    if (file.exists(fname)){
      load(fname)
    } else{
      cat('Fitting GLMM...\n')
      m1 = glmer('choice ~ delta_ev*dosage + (delta_ev*dosage | subjid)', df, family = binomial)
      save(m1, file = fname)
    }
    ind_pred_df = data_grid(df, delta_ev = seq_range(delta_ev, by = 1), dosage = dosage, subjid = subjid)
    ind_pred_df$pred = predict(m1, ind_pred_df, type = 'response', allow.new.levels=TRUE)
    qt = quantile(df$delta_ev, probs = seq(0, 1, 0.25))
    p = ggplot(df) +
      #stat_summary_bin(mapping = aes(x = delta_ev, y = choice, color = dosage), breaks = qt, fun.data = bino, geom = 'pointrange', position = position_dodge(10)) +
      stat_summary_bin(mapping = aes(x = delta_ev, y = choice, color = as.factor(dosage)), fun.data = bino, bins = 7, geom = 'pointrange', position = position_dodge(10)) +
      geom_line(ind_pred_df, mapping = aes(x = delta_ev, y = pred, color = as.factor(dosage)), alpha = 0.5, size = 1) +
      scale_color_manual(values = colors) + scale_fill_manual(values = colors) +
      labs(fill = 'Dose(µg)', color = 'Dose(µg)') +
      scale_x_continuous(limits = c(-100, 300), breaks = c(-100, 0, 200)) +
      xlab(expression(EV[lottery]-EV[surebet])) + ylab("P(Chose Lottery)")
    if (region == 'bi_fof'){
      p <- p + ggtitle('Muscimol bilateral FOF')
    }else{
      p <- p + ggtitle('Muscimol bilateral PPC')
    }
  } else if (str_detect(region, 'uni')){
    df$dosage = factor(df$dosage)
    df$infusion_side = factor(df$infusion_side, )
    df$infusion_bino = factor(df$infusion_bino)

    cat('Fitting GLMM...\n')
    m1 = glmer('choice ~ delta_ev*dosage + (delta_ev*dosage | subjid)', df, family = binomial)
    ind_pred_df = data_grid(df, delta_ev = seq_range(delta_ev, by = 1), dosage = dosage, subjid = subjid)
    ind_pred_df$pred = predict(m1, ind_pred_df, type = 'response', allow.new.levels=TRUE)
    qt = quantile(df$delta_ev, probs = seq(0, 1, 0.25))
    p = ggplot(df) +
      #stat_summary_bin(mapping = aes(x = pro_right_delta_ev, y = choose_right, color = infusion_side), breaks = qt, fun.data = bino, geom = 'pointrange', position = position_dodge(10)) +
      stat_summary_bin(mapping = aes(x = delta_ev, y = choice, color = dosage), fun.data = bino, bins = 7, geom = 'pointrange', position = position_dodge(10)) +
      geom_line(ind_pred_df, mapping = aes(x = delta_ev, y = pred, color = dosage), alpha = 0.5, size = 1) +
      scale_color_manual('Dose(µg)',values = c('0' = 'azure4', '0.3' = '#D766C9'), labels = c("0", "0.3")) +
      scale_fill_manual('Dose(µg)',values = c('0' = 'azure4', '0.3' = '#D766C9'), labels = c("0", "0.3")) +
      labs(fill = 'Infusion Side', color = 'Infusion Side') +
      scale_x_continuous(limits = c(-100, 300), breaks = c(-100, 0, 200)) +
      xlab(expression(EV[lottery]-EV[surebet])) + ylab("P(Chose Lottery)")
    if (region == 'uni_fof'){
      p <- p + ggtitle('Muscimol unilateral FOF')
    }else{
      p <- p + ggtitle('Muscimol unilateral PPC')
    }
  }
  # plot
  ntrials_df = df_m %>% group_by(subjid) %>% tally() %>%
    mutate(label = paste0('n=', n), x = max(df$delta_ev)*0.75,  y = 0.1)
  p = p + theme_classic(base_size = BASE_SIZE) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed',  alpha = 0.4) +
    #scale_x_continuous(limits = c(-100, 300), breaks = c(-100, 0, 200)) +
    scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1)) +
    geom_text(data = ntrials_df, mapping = aes(x = x, y = y, label = label)) +
    facet_wrap(~subjid) +
    theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5))
  return(p)
}

#' plot_supp_RT
#'
#'
#' @param region: 'bi_fof', 'bi_ppc', 'uni_fof', 'uni_ppc'
#'
#' @return ggplot object
plot_supp_RT = function(region = 'bi_fof'){
  # LMM fits and plots for reaction time in Bilateral FOF/PPC, Unilateral FOF/PPC
  choice_title = c('Choose surebet', 'Choose lottery')
  choices = c(0, 1)
  df_m = read.csv(file.path(git_root,"csv",sprintf('figure_infusion_%s.csv', region))) %>%
    filter(dosage != 0.15 & dosage != 0) %>%
    filter(RT < 3) %>% filter(!is.na(RT))
  control_df = read.csv(file.path(git_root,"csv", "figure_behavior_population.csv"))
  df = rbind(df_m, control_df)
  fname = file.path(git_root,"fits", sprintf('%s_RT_glmm.RData', region))
  # load or fit LMM
  if (str_detect(region, 'bi')){
    # load saved model filts
    if (file.exists(fname)){
      load(fname)
    } else{
      cat('Fitting GLMM...\n')
      m1 = lmer('log_RT ~ delta_ev*dosage*choice + (delta_ev*dosage*choice|subjid)', df)
      save(m1, file = fname)
    }
    x = 'delta_ev'
    c = 'dosage'
    if (region == 'bi_fof'){colors = c('azure4', 'rosybrown1', 'purple4')} else{ colors = c('azure4', 'gold2')}
    lab = 'Dosage'
    xlab = expression(EV[lottery]-EV[surebet])
  } else if (str_detect(region, 'uni')){
    df$infusion_side = factor(df$infusion_side, levels = c('control', 'L', 'R'))
    if (file.exists(fname)){
      load(fname)
    } else{
      cat('Fitting GLMM...\n')
      m1 = lmer('log_RT ~ delta_ev*infusion_side*choice + (pro_right_delta_ev*infusion_side*choice | subjid)', df)
      save(m1, file = fname)
    }
    x = 'pro_right_delta_evx'
    c = 'infusion_side'
    colors = side_colors
    lab = 'Infusion Side'
    xlab = expression(EV[right]-EV[left])
    df$pro_right_delta_evx = df$delta_ev
  }
  df$dosage = factor(df$dosage)
  df$choice = factor(df$choice)
  df$pred = predict(m1)
  for (i in 1:2){
    p = ggplot(df %>% filter(choice == choices[i]), aes_string(x = x, y = 'log_RT', color = c, fill = c)) +
      theme_classic(base_size = BASE_SIZE) + xlab(xlab) + ylab("log(RT)") +
      stat_summary(fun.data = mean_se, position = position_dodge(10)) +
      ggtitle(choice_title[i]) +
      geom_line(mapping = aes_string(x = x, y = 'pred', color = c)) +
      scale_color_manual(values = colors) + scale_fill_manual(values = colors) +
      labs(fill = lab, color = lab) + facet_wrap(subjid~.) +
      scale_x_continuous(limits = c(-50, 200), breaks = c(-50, 0, 150)) +
      scale_y_continuous(limits = c(-1.5, 1.5), breaks = c(-1.5, 0, 1.5)) +
      theme(legend.position = 'bottom',
            axis.text.x = element_text(size = 11)) + facet_wrap(~subjid, scales = 'fixed')
    if (i == 1){P = p} else{P = P + p}
  }
  return(P)
}


# make the "Early Trial Cutoff" vs. "P-value" for muscimol bilateral ppc data
generate_trials_cutoff_bippc <- function(){
  trial_cutoff = c()
  p_value = c()
  m_bi_control_df = read.csv(file.path(git_root, 'csv','figure_behavior_population.csv'))
  for (i in 1:20){
    trial_num = 10+i*5
    m_bi_ppc_df = read.csv(file.path(git_root, 'csv','figure_infusion_bi_ppc.csv')) %>% filter(dosage != 0, trialnum<=trial_num)
    m_bi_ppc_df = rbind(m_bi_control_df, m_bi_ppc_df)
    m_bi_ppc_df$subjid <- factor(m_bi_ppc_df$subjid)
    mf_m_bi_ppc = glmer(choice ~ delta_ev*dosage + (dosage*delta_ev|subjid), m_bi_ppc_df, family = binomial, nAGQ = 0)
    trial_cutoff <- append(trial_cutoff, trial_num)
    p_value <- append(p_value, coef(summary(mf_m_bi_ppc))[3,4])
  }
  df <- data.frame(trial_cutoff, p_value)
  return(df)
}

plot_ppc_cutoff_p <- function(){
  df <- read.csv('csv/PPC_cutoff_pvalue.csv')
  p <- ggplot(df, aes(trial_cutoff, p_value)) + theme_classic(base_size = BASE_SIZE) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0.05, linetype = 'dashed', alpha = 0.4) +
    scale_y_continuous(trans = 'log10', breaks = trans_breaks('log10', function(x) 10^x),
                       labels = trans_format('log10', math_format(10^.x))) +
    xlab('Early trial cutoff') + ylab('p value') + theme(axis.text = element_text(size = 16))
  return(p)
}


# check the 'early trial cut off' for the bi_fof data
generate_trials_cutoff_bifof <- function(){
  trial_cutoff = c()
  p_value = c()
  m_bi_control_df = read.csv(file.path(git_root, 'csv','figure_behavior_population.csv'))
  for (i in 1:20){
    trial_num = 10+i*5
    m_bi_fof_df = read.csv(file.path(git_root, 'csv','figure_infusion_bi_fof.csv')) %>% filter(dosage != 0, trialnum<=trial_num)
    m_bi_fof_df = rbind(m_bi_control_df, m_bi_fof_df)
    m_bi_fof_df$subjid <- factor(m_bi_fof_df$subjid)
    mf_m_bi_fof = glmer(choice ~ delta_ev*dosage + (dosage*delta_ev|subjid), m_bi_fof_df, family = binomial, nAGQ = 0)
    trial_cutoff <- append(trial_cutoff, trial_num)
    p_value <- append(p_value, coef(summary(mf_m_bi_fof))[3,4])
  }
  df <- data.frame(trial_cutoff, p_value)
  return(df)
}

plot_fof_cutoff_p <- function(){
  df <- generate_trials_cutoff_bifof()
  p <- ggplot(df, aes(trial_cutoff, p_value)) + theme_classic(base_size = BASE_SIZE) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0.05, linetype = 'dashed', alpha = 0.4) +
    scale_y_continuous(trans = 'log10', breaks = trans_breaks('log10', function(x) 10^x),
                       labels = trans_format('log10', math_format(10^.x))) +
    xlab('Early trial cutoff') + ylab('p value') + theme(axis.text = element_text(size = 16))
return(p)
}
