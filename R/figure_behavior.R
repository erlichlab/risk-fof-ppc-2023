# functions for figure1 (figure_behavior)

#' figure_behavior_timeline
#'
#' figure1b: timeline of trials in one example session.
#' @param n: the number of trials to show
#'
#' @return ggplot object
#'
plot_behavior_timeline = function(n = 30){
  # figure1b: timeline of trials in one example session.
  df = read.csv(file.path(git_root,"csv", "figure_behavior_timeline.csv"))
  s = sample(unique(df$sessid), 1)
  df = df %>% filter(sessid == s) %>% # sample one random session
    filter(viol != 1 & trialnum < n) %>%
    mutate(trialtype = case_when(!forced ~ 1, forced_lottery ~ 2, forced_surebet ~ 3))
  df$freq = df %>% group_indices(lottery_mag)
  df = df %>% mutate(image = sprintf('sinewaves/%d.png', freq))
  n_trials = dim(df)[1]
  df$trialnum = seq(1, n_trials)

  reward_df = df %>% filter(reward_received > 0)
  no_reward_df = df %>% filter(reward_received == 0)
  p = ggplot(df, aes(x = trialnum)) + theme_classic(BASE_SIZE) +
    geom_point(aes(y = rep(3, n_trials), shape = as.factor(trialtype), fill = as.factor(trialtype)), size = 4) +
    geom_image(aes(y = rep(2, n_trials), image = image), size = 0.03, asp = 3.5) +
    geom_point(aes(y = rep(1, n_trials), color = as.factor(choice)), shape = 18, size = 4) +
    geom_point(reward_df, mapping = aes(x = trialnum, y = rep(0, dim(reward_df)[1]), size = as.factor(reward_received)), shape = 21, fill = 'dodgerblue', alpha = 0.5) +
    geom_point(no_reward_df, mapping = aes(x = trialnum, y = rep(0, dim(no_reward_df)[1])), shape = 4, color = 'indianred', size = 4) +
    scale_y_continuous(limits = c(-0.5, 3.5), breaks = c(0, 1, 2, 3), labels = c(expression(paste('Reward(',mu,'l)')), 'Choice', 'Tone Freq', 'Trial Type')) +
    scale_shape_manual(values = c(23, 24, 25)) + # trialtype: choice, forced_lottery, forced_surebet
    scale_fill_manual(values = c('white', 'gold3', 'dodgerblue')) + # trialtype: choice, forced_lottery, forced_surebet
    scale_color_manual(values = c('dodgerblue', 'gold3')) + # choice: surebet, lottery
    xlab('Trials in session') + ylab(' ') +
    theme(legend.position = 'none', text = element_text(color = 'black'), axis.text = element_text(size = 16),
          axis.text.y = element_text(angle = 45))
  return(p)
}

#' figure_behavior_example
#'
#' figure1c: example animal psychometric curves by GLMM, thin lines for each session, thick line for combined
#'
#' @return ggplot object
plot_behavior_example = function(){
  # figure1c: example animal psychometric curves by GLMM, thin lines for each session, thick line for combined
  risk_df = read.csv('csv/figure_behavior_population.csv') %>% filter(subjid == 2154)
  # GLMM fit
  m0 = glm('choice ~ delta_ev', risk_df, family = binomial)
  m1 = glmer('choice ~ delta_ev + (delta_ev | sessid)', risk_df, family = binomial)
  avg_pred_df = data_grid(risk_df, delta_ev = seq_range(c(min(delta_ev)*1.3, max(delta_ev)*1.1), by = 1))
  sess_pred_df = data_grid(risk_df, delta_ev = seq_range(c(min(delta_ev)*1.3, max(delta_ev)*1.1), by = 1), sessid = sessid)
  avg_pred_df$pred = predict(m0, avg_pred_df, type = 'response', allow.new.levels=TRUE)
  sess_pred_df$pred = predict(m1, sess_pred_df, type = 'response', allow.new.levels=TRUE)
  p = ggplot(risk_df) + theme_classic(base_size = BASE_SIZE) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed',  alpha = 0.4) +
    scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1)) +
    scale_x_continuous(breaks = c(round(min(risk_df$delta_ev)/10)*10, 0, round(max(risk_df$delta_ev)/10)*10)) +
    annotate("text", label = '1 animal', x = max(risk_df$delta_ev)*0.70, y = 0.2, size = ANNOTATION_SIZE) +
    annotate("text", label = sprintf('%d sessions', length(unique(risk_df$sessid))), x = max(risk_df$delta_ev)*0.70, y = 0.1, size = ANNOTATION_SIZE) +
    xlab(expression(EV[lottery]-EV[surebet])) + ylab("P(Chose Lottery)") +
    geom_line(avg_pred_df, mapping = aes(x = delta_ev, y = pred), color = 'black', alpha = 0.5, size = 2) +
    theme(axis.text = element_text(size = 16))
  for (sessidx in unique(risk_df$sessid)){
    p = p + geom_line(sess_pred_df %>% filter(sessid == sessidx),
                      mapping = aes(x = delta_ev, y = pred), color = 'darkgray', alpha = 0.3)
  }
  p = p + stat_summary_bin(mapping = aes(x = delta_ev, y = choice), fun.data = bino, geom = 'pointrange', color = 'black')
  return(p)
}

#' figure_behavior_population
#'
#' figure1d: population psychometric curves by GLMM, thin colored lines per subject, thick line for population
#'
#'
#' @return ggplot object
plot_behavior_population = function(){
  # figure1d: population psychometric curves by GLMM, thin colored lines per subject, thick line for population
  df <- read.csv(file.path(git_root,"csv", "figure_behavior_population.csv"), na.strings=c("","NA"))
  df_ephys <- read.csv(file.path(git_root,"csv", "ephys_behavior.csv"),na.strings=c("","NA"))
  df_bi_opto <- read.csv(file.path(git_root,"csv", "opto_bilateral.csv"),na.strings=c("","NA")) %>% filter(isopto == 0)
  df_uni_opto <- read.csv(file.path(git_root,"csv", "opto_unilateral.csv"),na.strings=c("","NA")) %>% filter(isopto == 0)

  df1 <- df %>% select(c('subjid','sessid','choice','delta_ev')) %>% mutate(exper = 'muscimol', exper_n = 1)
  df_ephys1 <- preprocessRisk(df_ephys) %>% select(c('subjid','sessid','choice','delta_ev')) %>% mutate(exper = 'ephys', exper_n = 2)
  df_bi_opto1 <- preprocessRisk(df_bi_opto) %>% select(c('subjid','sessid','choice','delta_ev')) %>% mutate(exper = 'opto', exper_n = 3)
  df_uni_opto1 <- preprocessRisk(df_uni_opto) %>% select(c('subjid','sessid','choice','delta_ev')) %>% mutate(exper = 'opto', exper_n = 3)

  df <- rbind(df1, df_ephys1, df_bi_opto1, df_uni_opto1)
  df$subjid <- factor(df$subjid, levels = unique(df$subjid))
  df$sessid <- factor(df$sessid, levels = unique(df$sessid))

  # GLMM fit
  m0 = glm('choice ~ delta_ev', df, family = binomial)
  m1 = glmer('choice ~ delta_ev + (delta_ev | subjid)', df, family = binomial)
  pop_pred_df = data_grid(df, delta_ev = seq_range(c(min(delta_ev)*1.3, max(delta_ev)*1.1), by = 1))
  ind_pred_df = data_grid(df, delta_ev = seq_range(c(min(delta_ev)*1.3, max(delta_ev)*1.1), by = 1), subjid = subjid)
  pop_pred_df$pred = predict(m0, pop_pred_df, type = 'response', allow.new.levels=TRUE)
  ind_pred_df$pred = predict(m1, ind_pred_df, type = 'response', allow.new.levels=TRUE)
  ind_pred_df <- ind_pred_df %>% mutate(experiment = case_when(subjid %in% c(2152, 2153, 2154, 2155, 2156, 2160, 2165, 2166) ~ 'Muscimol',
                                                subjid %in% c(2176, 2182, 2225, 2228, 2240, 2172, 2177, 2180) ~ 'Opto',
                                                subjid %in% c(2224, 2238, 2244, 2261, 2263, 2264) ~ 'Ephys'))

  p = ggplot(df) + theme_classic(base_size = BASE_SIZE) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed',  alpha = 0.4) +
    scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1)) +
    scale_x_continuous(breaks = c(round(min(df$delta_ev)/10)*10, 0, round(max(df$delta_ev)/10)*10)) +
    xlab(expression(EV[lottery]-EV[surebet])) + ylab("P(Chose Lottery)") +
    geom_line(ind_pred_df, mapping = aes(x = delta_ev, y = pred, group = subjid, color = experiment), alpha = 0.7, size = .6) +
    theme(legend.title = element_blank(), legend.position = c(0.2, 0.9),
          legend.text=element_text(size=6), axis.text = element_text(size = 16),
          legend.key = element_rect(fill = "transparent", color = NA),
          legend.spacing.y = unit(0.1, 'cm'))
  return(p)
}


#' figure_behavior_softfixation
#'
#' analyze the softfixation pattern for all the softfixation animals
#'
#'
#' @return ggplot object
plot_soft_fixation <- function(){
  fixation_df <- read.csv(file.path(git_root, "csv", "fixation_bin_tb.csv"))
  p <- ggplot(fixation_df, aes(time, inpoke)) +
    theme_classic(base_size = 15) +
    geom_area() +
    scale_x_continuous(breaks=c(0,0.5,1)) +
    scale_y_continuous(breaks=c(0,0.5,1)) +
    labs(y = 'P(in port)', x = 'Time during fixation (s)') +
    facet_wrap(.~subjid, nrow = 3)
  return(p)
}
