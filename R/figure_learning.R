#' figure_learning_schematic
#'
#' schematic showing how shifting the surebet magnitude shifts the decision boundary
#'
#'
#' @return ggplot object
#'
plot_learning_schematic = function(){
  # schematic showing how shifting the surebet magnitude shifts the decision boundary
  syn_df = data.frame()
  sb_mags = c(3, 5, 7)
  true_params = list('rho' = 1, 'sigma' = 3, 'omega' = c(1, 0, 0))
  for (i in 1:3){
    df = data.frame(lottery_mag = seq(1, 32, by = 0.01),
                    sb_mag = sb_mags[i], total_rew_multi = 8, lottery_prob = 0.5) %>%
      mutate(delta_ev = total_rew_multi*lottery_mag*lottery_prob - total_rew_multi*sb_mag)
    df$choice_prob = rho_sigma_agent(prob = TRUE, params = true_params, df$sb_mag, df$lottery_mag,
                                     df$lottery_prob, df$total_rew_multi)
    syn_df = rbind(syn_df, df)
  }
  syn_df$sb_mag = factor(syn_df$sb_mag)
  p = ggplot(syn_df, aes(x = lottery_mag, y = choice_prob, color = as.factor(sb_mag))) +
    theme_classic(base_size = BASE_SIZE) + geom_line(size = 2) +
    geom_line(alpha = 0.4) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
    ylab('P(choose lottery)') + xlab('Lottery magnitude') +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_color_manual(values = c('deepskyblue1', 'azure4', 'deepskyblue4')) +
    labs(color = 'Surebet magnitude') + theme(legend.position = 'none')
  return(p)
}

#' figure_learning_scatter
#'
#' scatter plot of model predicted change vs. actual change in probability choosing lottery, color by infusion and control
#'
#'
#' @return ggplot object
plot_learning_scatter = function(){
  # scatter plot of model predicted change vs. actual change in probability choosing lottery, color by infusion and control
  mean_lott_df = read.csv(file.path(git_root,"csv", "figure_learning_scatter.csv"))
  corr = cor.test(mean_lott_df$pred_shift, mean_lott_df$actual_shift)
  p = ggplot(mean_lott_df,  mapping = aes(x = pred_shift, y = actual_shift)) +
    theme_bw(base_size = BASE_SIZE) +
    xlab('Predicted Shift') + ylab('Actual Shift') +
    scale_x_continuous(limits = c(-0.45, 0.4), breaks = c(-0.3, 0, 0.3)) +
    scale_y_continuous(limits = c(-0.45, 0.4), breaks = c(-0.3, 0, 0.3)) +
    annotate("text", label = sprintf('R = %.3f\n p < 0.001', corr$estimate), x = -0.3, y = 0.2) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    theme(legend.position = 'none') + labs(color = ' ') +
    geom_point(mapping = aes(fill = infusion), shape = 21, size = 2) +
    geom_smooth(method = lm, formula = y ~ x, color = 'deepskyblue') +
    geom_smooth(method = lm, formula = y ~ x, color = 'deepskyblue') +
    scale_fill_manual(values = sb_change_colors, labels = c('Control', 'Infusion'))
  return(p)
}

genSynSmooth = function(df, sb_mags = NULL, extend = 'delta_ev'){
  # generate smooth synthetic data points for every task parameter combinations
  if (is.null(sb_mags)){
    sb_mags = get_freq_task_param(df, 'sb_mag')
  }
  if (extend == 'delta_ev'){
    pred_df = df %>% group_by(subjid) %>%
      data_grid(delta_ev = seq_range(c(min(delta_ev), max(delta_ev)*1.1), by = 0.5),
                species = 2, individual = 1, prev_outcome,
                lottery_prob = get_freq_task_param(df, 'lottery_prob'),
                sb_mag = sb_mags,
                total_rew_multi = get_freq_task_param(df, 'total_rew_multi'),
                n_trials = 200) %>%
      mutate(lottery_mag = (delta_ev + sb_mag*total_rew_multi) / (lottery_prob*total_rew_multi))
  } else if (extend == 'lottery_mag'){
    pred_df = df %>% group_by(subjid) %>%
      data_grid(lottery_mag = seq_range(lottery_mag, by = 0.5),
                species = 2, individual = 1, prev_outcome = 0, # if extend based on lottery_mag, we assume history is not included
                lottery_prob = get_freq_task_param(df, 'lottery_prob'),
                sb_mag = sb_mags,
                total_rew_multi = get_freq_task_param(df, 'total_rew_multi'),
                n_trials = 200) %>%
      mutate(delta_ev = lottery_mag*lottery_prob*total_rew_multi - sb_mag*total_rew_multi)
  }
  pred_df = pred_df %>% mutate(prev_outcome_s = case_when(prev_outcome == -1 ~ 'lottery_lose',
                                                          prev_outcome == 1 ~ 'lottery_win',
                                                          prev_outcome == 0 ~ 'surebet')) %>%
    filter(lottery_mag >= 0) %>% ungroup()
  colnames(pred_df) = paste0('pred_', colnames(pred_df))
  pred_df$pred_prev_outcome_s = factor(pred_df$pred_prev_outcome_s, levels = c('lottery_lose', 'surebet', 'lottery_win'))
  return(pred_df)
}

plot_learning_individual_summary = function(){
  # psychometric curves plotting lottery magnitudes vs. p(chose lottery), colored by before / after changing sure-bet value
  # and time-series showing % choose lottery over two weeks of the surebet change experiment
  # this one cannot be run as it requires model fitting
  individual_pred_model = stan_model(file = file.path(git_root,"stan","individual_pred.stan"))
  df = read.csv(file.path(git_root,"csv", "figure_learning.csv"))
  df$sessiondate = as.Date(df$sessiondate)
  df = df %>% filter(sessiondate < '2020-12-16')
  for (animal in PPC_sb_animals){
    this_df = df %>% filter(subjid == animal)
    change_dates = unique(this_df %>% filter(infusion %in% c('sb_change','infusion')) %>% pull(sessiondate)) # find all infusion dates
    sig_days = c()
    for (i in 1:length(change_dates)){
      if (i > 1){
        before_df = this_df %>% filter(sessiondate >= change_dates[i-1] & sessiondate < change_dates[i])
      } else if (i == 1){
        before_df = this_df %>% filter(sessiondate < change_dates[i])
      }
      if (i == length(change_dates)){
        after_df = this_df %>% filter(sessiondate >= change_dates[i])
      } else{
        after_df = this_df %>% filter(sessiondate >= change_dates[i] & sessiondate < change_dates[i+1])
      }
      old_sb = unique(before_df %>% filter(sessiondate == change_dates[i] - 1) %>% pull(sb_mag))
      new_sb = unique(this_df %>% filter(sessiondate == change_dates[i]) %>% pull(sb_mag))
      # fit model to this_df, predict using task parameters from after_df
      bino_df = binomialize(this_df)
      pred_df = genSynSmooth(this_df, sb_mags = c(old_sb, new_sb))
      data = append(as.list(bino_df), as.list(pred_df))
      data$T = dim(bino_df)[1]
      data$P = dim(pred_df)[1]
      all_fit = sampling(individual_pred_model, data, refresh = 0)
      # extract simulated n_chose_lott
      draws = rstan::extract(all_fit)
      ncl_df = as.data.frame(t(draws$pred_n_chose_lott))
      pred_df = pred_df %>% mutate(y = rowMeans(ncl_df)/pred_n_trials,
                                   ymin = apply(ncl_df, 1, quantile, 0.1)/pred_n_trials,
                                   ymax = apply(ncl_df, 1, quantile, 0.9)/pred_n_trials)
      # factor old and new sb_mag so that it's not just based on value
      pred_df = pred_df %>% mutate(sb_mag_s = case_when(pred_sb_mag == old_sb ~ 'old',
                                                        pred_sb_mag == new_sb ~ 'new'))
      pred_df$sb_mag_s = factor(pred_df$sb_mag_s, c('old', 'new'))
      before_df$sb_mag_s = 'old'
      after_df$sb_mag_s = 'new'
      both_df = rbind(before_df, after_df)
      both_df$sb_mag_s = factor(both_df$sb_mag_s, c('old', 'new'))
      infusion_or_change = unique(this_df %>% filter(sessiondate == change_dates[i]) %>% pull(infusion))
      after_color = ifelse(infusion_or_change == 'infusion', 'gold2', 'deepskyblue')
      # plot
      p = ggplot(both_df) + theme_classic(base_size = BASE_SIZE) +
        geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
        scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1)) +
        annotate("text", label = sprintf('%.2f to %.2f', old_sb, new_sb), x = max(both_df$lottery_mag)*0.70, y = 0.2) +
        xlab('Lottery magnitude') + ylab("P(Chose Lottery)") +
        stat_summary_bin(mapping = aes(x = lottery_mag, y = choice, color = sb_mag_s), fun.data = bino, geom = 'pointrange') +
        geom_ribbon(pred_df, mapping = aes(x = pred_lottery_mag, ymin = ymin, ymax = ymax, fill = sb_mag_s), alpha = 0.5) +
        scale_color_manual(values = c('azure4', after_color), labels = c('Before', 'After')) +
        scale_fill_manual(values = c('azure4', after_color), labels = c('Before', 'After')) +
        labs(fill = ' ', color = ' ') + theme(legend.position = 'none')
      if (i == 1){
        P = p
      } else{
        p = p + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
        P = P + p
      }
      # when can we detect the shift? store the day number
      both_df = both_df %>% mutate(days_relative_to_shift = as.numeric(sessiondate - change_dates[i]),
                                   shift = as.numeric(sessiondate >= change_dates[i]))
      for (dx in 0:6){
        this_both_df = both_df %>% filter(days_relative_to_shift <= dx)
        m1 = glm(choice ~ lottery_mag*shift, data = this_both_df)
        p_value = summary(m1)$coefficients["shift","Pr(>|t|)"]
        if (p_value < 0.05){
          break
        }
        if (dx == 6){
          dx = 100 # impossible sessiondate to show asterik
        }
      }
      sig_days = append(sig_days, dx)
    }
    # plot % choosing lottery over days
    temp_df = this_df %>% group_by(sessiondate) %>%
      summarise(infusion = unique(infusion), y = bino(choice)$y, ymin = bino(choice)$ymin, ymax = bino(choice)$ymax) %>%
      mutate(day_of_change = infusion != 'control')
    day_change_detected = temp_df$sessiondate[temp_df$day_of_change] + sig_days
    p3 = ggplot(temp_df, aes(x = sessiondate, y = y*100, ymin = ymin*100, ymax = ymax*100)) +
      theme_classic(base_size = BASE_SIZE) + geom_pointrange() +
      geom_pointrange(data = temp_df %>% filter(infusion == 'infusion'), colour = "gold2") +
      geom_pointrange(data = temp_df %>% filter(infusion == 'sb_change'), colour = "deepskyblue") +
      annotate('text', x = temp_df$sessiondate[temp_df$sessiondate %in% day_change_detected],
               y = (temp_df$ymax[temp_df$sessiondate %in% day_change_detected]+0.02)*100, label = '*', size = 7) +
      scale_y_continuous(breaks = c(0, 50, 100), limits = c(-2, 102)) +
      geom_hline(yintercept = 50, linetype = 'dashed', alpha = 0.4) +
      xlab(' ') + ylab("% Chose Lottery")
    P = P / p3 + plot_layout(heights = c(1, 0.5)) + plot_annotation(title = animal)
    scale_save(P, sprintf('figureS_learning_%d', animal), 20, 12, 1)
  }
}
