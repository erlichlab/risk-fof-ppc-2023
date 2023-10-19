# utility functions
options(warn = -1)
BASE_SIZE = 15
SEED = 12345 # random seed
ANNOTATION_SIZE = 5

infusion_animals = c(2152, 2153, 2154, 2155, 2156, 2160, 2165, 2166)
PPC_sb_animals = c(2153, 2154, 2156, 2160)
FOF_sb_animals = c(2153, 2154, 2155, 2156, 2160, 2165, 2166)
side_animals = c(2154, 2155, 2156, 2160, 2165, 2166)

FOF_dosage_colors = c('azure4', 'rosybrown1', 'purple4') # control, 0.075, 0.3
PPC_dosage_colors = c('azure4', 'gold2') # control, 0.3
side_colors = c('control' = 'azure4', 'L' = '#2a6a99', 'R' = '#d88546') # control, Left, Right
sb_change_colors = c('gold2', 'deepskyblue')

#' rho_sigma_agent
#'
#' returns choice_prob (probability choosing lottery) or choice (binary) from the simulated three-agent rho-sigma agent
#'
#' @param prob: Boolean, if true returns choice_prob
#' @param params: a list as in list('rho' = 1, 'sigma' = 0.4, 'omega' = c(0.8, 0.1, 0.1))
#' @param sb_mag: a vector of n_trials with surebet magnitude
#' @param lottery_mag: a vector of n_trials with lottery magnitude
#' @param lottery_prob: a vector of n_trials with lottery probability
#' @param total_rew_multi: a vector of n_trials with total reward multiplier
#'
#' @return a vector of n_trials with choice_prob or choice
rho_sigma_agent = function(prob = FALSE, params, sb_mag, lottery_mag, lottery_prob, total_rew_multi){
  # returns choice_prob from the simulated three-agent rho-sigma agent
  rho = params$rho # rho < 1 is risk-averse, rho = 1 is risk-neutral, rho > 1 is risk-seeking
  sigma = params$sigma # noise associated with the internal representation of utility
  omega = params$omega # vector of size 3, baseline fractions of being rational, lottery and surebet agent
  n_trials = length(sb_mag)
  choice_prob = vector(length = n_trials)
  for (tx in 1:n_trials){
    u_l = lottery_prob[tx]*(total_rew_multi[tx]*lottery_mag[tx])^rho
    u_sb = (total_rew_multi[tx]*sb_mag[tx])^rho
    p_rational = 1 - pnorm(0, u_l - u_sb, sqrt(2)*sigma)
    P = c(p_rational, 1, 0)
    choice_prob[tx] = omega %*% P
  }
  choice = rbinom(n_trials, 1, choice_prob)
  if (prob){
    return(choice_prob)
  } else{
    return(choice)
  }
}

#' binomialize
#'
#' format df into binomial version for model fitting
#'
#' @param df: dataframe containing risky trials
#'
#' @return dataframe
binomialize = function(df){
  # format df into binomial version for model fitting
  bino_df = df %>% group_by(lottery_mag, lottery_prob, sb_mag, total_rew_multi) %>%
    add_tally() %>% summarise(n_trials = mean(n), n_chose_lott = sum(choice)) %>%
    mutate(delta_ev = lottery_mag*lottery_prob*total_rew_multi - sb_mag*total_rew_multi) %>%
    ungroup()
  return(bino_df)
}

#' bino
#'
#' get binomial mean and confidence intervals
#'
#' @param x: a vector of 0s and 1s
#'
#' @return a dataframe with y, ymin and ymax
bino <-function(x){
  # get binomial mean and confidence intervals
  out = binom.test(sum(x), length(x))
  df = data.frame(y = mean(x), ymin=out$conf.int[1], ymax=out$conf.int[2])
  return(df)
}

#' se
#'
#' standard error
#'
#' @param x: a vector of values
#'
#' @return numeric
se = function(x){
  # standard error
  return(sd(x) / sqrt(length(x)))
}



logit2prob = function(logit){
  odds = exp(logit)
  prob = odds / (1 + odds)
  return(prob)
}

#' get_freq_task_param
#'
#' find the most frequent task parameter given the df
#'
#' @param df: dataframe with risky trials
#' @param param: task parameter name, 'sb_mag', 'lottery_mag', 'lottery_prob' and 'total_rew_multi'
#'
#' @return numeric
get_freq_task_param = function(df, param = 'sb_mag'){
  # find the most frequent task parameter given the df
  return(as.numeric(names(sort(table(df[param]), decreasing = TRUE)[1])))
}


#' scale_save
#'
#' ggplot wrapper to save the plot
#'
#' @param p: a ggplot object to be saved
#' @param name: the name of the file
#' @param width: the width of the plot
#' @param height: the height of the plot
#' @param scale: scale of the plot, 1 is the default
#'
#' @return saves the plot in the specified path
scale_save = function(p, name, width = 4, height = 4, scale = 1){
  # ggsave wrapper
  p = p + theme(text = element_text(size = BASE_SIZE),
                axis.text = element_text(size = BASE_SIZE),
                legend.text = element_text(size = BASE_SIZE))
  fname = file.path(git_root,"plots", sprintf('%s.pdf', name))
  ggsave(filename = fname, device = "pdf", width = width, height = height, scale = scale, unit = 'cm')
}


preprocessRisk = function(df, remove_forced = TRUE, remove_viol = TRUE, remove_slow = TRUE){
  # preprocess raw risk data to be fitting-friendly
  # df may contain single or multiple subjects, however the region must be the same
  df$individual = df %>% group_indices(subjid)
  df = df %>%
    arrange(subjid, sessid, trialid) %>%
    mutate(species_name = case_when(subjid > 1000 & subjid < 2000 ~ 'mouse',
                                    subjid > 2000 ~ 'rat'),
           species = ifelse(species_name == 'mouse', 1, 2)) %>%
    mutate(forced = is.na(lottery_poke) | is.na(surebet_poke), # forced trials
           forced_surebet = is.na(lottery_poke),
           forced_lottery = is.na(surebet_poke)) %>%
    mutate(is_side_choice = is.na(lottery_poke) & is.na(surebet_poke)) %>%  # side choice trials
    group_by(sessid) %>% # find session-specific total_rew_multiplier
    mutate(total_rew_multi = case_when(subj_choice == 'violation' ~ 0,
                                       subj_choice == 'lottery' ~ 0,
                                       subj_choice == 'surebet' ~ reward_received / sb_mag)) %>%
    mutate(total_rew_multi = sort(unique(total_rew_multi))[2]) %>% ungroup() %>%
    mutate(delta_ev = (total_rew_multi * lottery_mag * lottery_prob) - (total_rew_multi * sb_mag)) %>%
    mutate(prev_reward = lag(reward)) %>%  # add previous reward
    mutate(choice = case_when(subj_choice == 'lottery' ~ 1,
                              subj_choice == 'surebet' ~ 0,
                              subj_choice == 'violation' ~ 9)) %>%
    mutate(prev_choice = lag(choice)) %>% # add previous choice
    mutate(prev_outcome_s = case_when(prev_reward > 0 & prev_choice == 1 ~ 'lottery_win',
                                      prev_reward == 0 & prev_choice == 1 ~ 'lottery_lose',
                                      prev_choice == 0 ~ 'surebet',
                                      prev_choice == 9 ~ 'surebet')) %>% # trials preceded by a violation
    mutate(prev_outcome = case_when(prev_outcome_s == 'lottery_win' ~ 1,
                                    prev_outcome_s == 'lottery_lose' ~ -1,
                                    prev_outcome_s == 'surebet' ~ 0)) %>%
    filter(trialnum != 1) %>% # remove the first trials
    filter(!is_side_choice) %>% # remove side choice trials
    filter(!is.na(total_rew_multi)) %>%
    filter(!is.na(prev_outcome_s)) # fix me!
  df$prev_outcome_s = factor(df$prev_outcome_s, levels = c('lottery_lose', 'surebet', 'lottery_win'))
  #filter(!is.na(prev_outcome))
  if (remove_forced){
    df = df %>% filter(forced == 0) # remove forced trials
  }
  if (remove_viol){
    df = df %>% filter(viol == 0) # remove violation trials
  }
  if (remove_slow){
    df = df %>% filter(RT < 5) # remove slow trials, assuming no attention is paid to this trial
  }
  df = df %>% group_by(subjid) %>% mutate(log_RT = log(RT/mean(RT, na.rm = TRUE))) # normalize RT within subjects
  #df = df %>% mutate(log2_RT = log2(RT/mean(RT))) # normalize RT across subjects
  # some unilateral-specific preprocessing
  #if (any(unique(df$infusion_side) %in% c('L','R'))){
  df = df %>% group_by(subjid) %>%
    mutate(lottery_side = unique(lottery_poke)[!is.na(unique(lottery_poke))]) %>%
    ungroup() %>%
    mutate(pro_right_delta_ev = case_when(lottery_side == 'BotL' ~ -delta_ev,
                                          lottery_side == 'BotR' ~ delta_ev),
           choose_right = case_when(lottery_side == 'BotL' ~ 1 - choice,
                                    lottery_side == 'BotR' ~ choice))
  return(df)
}
