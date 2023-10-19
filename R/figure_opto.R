# plot bilateral and unilateral opto inhibition FOF data with GLM

preprocessOpto = function(df){
  # this function preprocesses opto risky choice data
  df = df %>%
    mutate(sessiondate = as.Date(trialtime)) %>%
    mutate(forced = is.na(lottery_poke) | is.na(surebet_poke)) %>%
    # find reward multiplier
    group_by(sessid) %>%
    mutate(total_rew_multi = case_when(subj_choice == "violation" ~ 0,
                                       subj_choice == "lottery" ~ 0,
                                       subj_choice == "surebet" ~ reward / sb_mag)) %>%
    mutate(total_rew_multi = sort(unique(total_rew_multi))[2]) %>% ungroup() %>%
    mutate(delta_ev = total_rew_multi * lottery_mag * lottery_prob - total_rew_multi * sb_mag) %>%
    # define outcome
    mutate(outcome_s = case_when(subj_choice == "lottery" & lottery_outcome == "win" ~ "lottery_win",
                                 subj_choice == "lottery" & lottery_outcome == "lose" ~ "lottery_lose",
                                 subj_choice == "surebet" ~ "surebet")) %>%
    mutate(outcome = case_when(outcome_s == "lottery_win" ~ 1,
                               outcome_s == "lottery_lose" ~ -1,
                               outcome_s == "surebet" ~ 0)) %>%
    mutate(choice = case_when(subj_choice == "lottery" ~ 1,
                              subj_choice == "surebet" ~ 0)) %>%
    # add data about previous trial
    mutate(pre_choice = lag(choice)) %>%
    mutate(prev_reward = lag(reward_received)) %>%
    mutate(prev_outcome_s = lag(outcome_s)) %>%
    mutate(prev_outcome = lag(outcome)) %>%
    mutate(region = toupper(region)) %>%
    mutate(opto_out = case_when(isopto == 1 & region == 'LEFT FOF' ~ "opto_left",
                                isopto == 1 & region == 'RIGHT FOF' ~ "opto_right",
                                isopto == 1 & region == 'BILATERAL FOF' ~ "opto_bilateral",
                                isopto == 0 ~ "no_opto")) %>%
    mutate(opto_contr_ipsi = case_when(opto_out == "opto_left" & lottery_poke == "BotL" ~ "ipsilateral",
                                       opto_out == "opto_left" & lottery_poke == "BotR" ~ "contralateral",
                                       opto_out == "opto_right" & lottery_poke == "BotL" ~ "contralateral",
                                       opto_out == "opto_right" & lottery_poke == "BotR" ~ "ipsilateral",
                                       opto_out == "opto_bilateral" ~ "bilateral",
                                       opto_out == "no_opto" ~ "no_opto")) %>%
    mutate(prev_opto = lag(opto_contr_ipsi)) %>%
    mutate(prev_prev_opto = lag(prev_opto)) %>%
    # remove first/slow trials
    filter(trialnum != 1 & RT_choice < 2) %>%
    arrange(subjid, sessid, trialid) %>%
    filter(!is.na(outcome))
  df$prev_outcome_f = factor(df$prev_outcome_s, levels = c("lottery_win", "surebet", "lottery_lose"))
  df$opto_out = factor(df$opto_out, levels = unique(df$opto_out))
  df$opto_contr_ipsi = factor(df$opto_contr_ipsi, levels = unique(df$opto_contr_ipsi))
  df$subjid = factor(df$subjid, levels = unique(df$subjid))
  df$sessid = factor(df$sessid, levels = unique(df$sessid))
  return(df)
}

plot_opto_bi <- function(){
  out = read.csv(file.path(git_root,"csv","opto_bilateral.csv"))
  df <- preprocessOpto(out)
  free_df <- df %>% filter(forced == 0)
  #model prediction
  m = glm(choice ~ delta_ev * opto_contr_ipsi, data = free_df, family = binomial)
  pred_df = data_grid(free_df, delta_ev = seq(min(delta_ev) - 20, max(delta_ev) + 20, by = 1), opto_contr_ipsi = opto_contr_ipsi)
  pred = predict(m, pred_df, type = 'response', allow.new.levels=TRUE, se.fit = TRUE)
  pred_df$pred = pred$fit
  pred_df$ymin = pred$fit - pred$se.fit
  pred_df$ymax = pred$fit + pred$se.fit

  p = ggplot(free_df) + theme_classic(base_size = 15) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed',  alpha = 0.4) +
    scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1)) +
    scale_x_continuous(limits = c(-70, 240), breaks = c(round(min(free_df$delta_ev)/10)*10, 0, round(max(free_df$delta_ev)/10)*10)) +
    annotate("text", label = sprintf('%d trials \n %d sessions', dim(free_df)[1], length(unique(free_df$sessid))), x = max(free_df$delta_ev)*0.65, y = 0.1, size = ANNOTATION_SIZE) +
    #annotate("text", label = paste0('Control n=', dim(free_df %>% filter(opto_contr_ipsi == 'no_opto'))[1]), x = max(free_df$delta_ev)*0.8, y = 0.2) +
    #annotate("text", label = paste0('Bilateral n=', dim(free_df %>% filter(opto_contr_ipsi == 'bilateral'))[1]), x = max(free_df$delta_ev)*0.8, y = 0.1) +
    #annotate("text", label = unique(subj_df$subjid), fontface = 2, x = max(subj_df$delta_ev)*0.8, y = 0.2) +
    xlab(expression(EV[lottery]-EV[surebet])) + ylab("P(Chose Lottery)") +
    stat_summary_bin(mapping = aes(x = delta_ev, y = choice, color = opto_contr_ipsi), bins = 6, fun.data = bino, geom = 'pointrange', position = position_dodge(10)) +
    #geom_line(pred_df, mapping = aes(x = delta_ev, y = pred, color = opto_contr_ipsi)) +
    geom_ribbon(pred_df, mapping = aes(x = delta_ev, ymin = ymin, ymax = ymax, fill = opto_contr_ipsi), alpha = 0.5) +
    scale_color_manual('Laser',values = c(no_opto = 'azure4', bilateral = 'purple3'), labels = c('Off', 'On')) +
    scale_fill_manual('Laser',values = c(no_opto = 'azure4', bilateral = 'purple3'),, labels = c('Off', 'On')) +
    ggtitle('Opto bilateral FOF') +
    theme(legend.position = 'bottom', axis.text = element_text(size = 16), plot.title = element_text(hjust = 0.5))
  return(p)
}


plot_opto_uni <- function(){
  out = read.csv('csv/opto_unilateral.csv', na.strings=c("","NA"))
  df <- preprocessOpto(out)
  free_df <- df %>% filter(forced == 0) %>% group_by(subjid) %>%
    mutate(lottery_side = unique(lottery_poke)[!is.na(unique(lottery_poke))]) %>%
    ungroup() %>%
    mutate(pro_right_delta_ev = case_when(lottery_side == 'BotL' ~ -delta_ev,
                                          lottery_side == 'BotR' ~ delta_ev),
           choose_right = case_when(lottery_side == 'BotL' ~ 1 - choice,
                                    lottery_side == 'BotR' ~ choice))
  #model prediction
  m = glm(choose_right ~ pro_right_delta_ev * opto_out, data = free_df, family = binomial)
  pred_df = data_grid(free_df, pro_right_delta_ev = seq(min(pro_right_delta_ev) - 20, max(pro_right_delta_ev) + 20, by = 1), opto_out)
  pred = predict(m, pred_df, type = 'response', allow.new.levels=TRUE, se.fit = TRUE)
  pred_df$pred = pred$fit
  pred_df$ymin = pred$fit - pred$se.fit
  pred_df$ymax = pred$fit + pred$se.fit
  p = ggplot(free_df) + theme_classic(base_size = 15) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed',  alpha = 0.4) +
    scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1)) +
    scale_x_continuous(breaks = c(round(min(free_df$pro_right_delta_ev)/10)*10, 0, round(max(free_df$pro_right_delta_ev)/10)*10)) +
    annotate("text", label = sprintf('%d trials \n %d sessions', dim(free_df)[1], length(unique(df$sessid))), x = max(free_df$pro_right_delta_ev)*0.6, y = 0.1, size = ANNOTATION_SIZE) +
    xlab(expression(EV[right]-EV[left])) + ylab("P(Chose Right)") +
    stat_summary_bin(mapping = aes(x = pro_right_delta_ev, y = choose_right, color = opto_out), bins = 4, fun.data = bino, geom = 'pointrange', position = position_dodge(10)) +
    #geom_line(pred_df, mapping = aes(x = pro_right_delta_ev, y = pred, color = opto_contr_ipsi)) +
    geom_ribbon(pred_df, mapping = aes(x = pro_right_delta_ev, ymin = ymin, ymax = ymax, fill = opto_out), alpha = 0.5) +
    scale_color_manual('Laser side', values = c(no_opto = 'azure4', opto_left = '#2a6a99', opto_right='#d88546'), labels = c('Off', 'L', 'R')) +
    scale_fill_manual('Laser side', values = c(no_opto = 'azure4', opto_left = '#2a6a99', opto_right='#d88546'), labels = c('Off', 'L', 'R')) +
    ggtitle('Opto unilateral FOF') +
    theme(legend.position = 'bottom', axis.text = element_text(size = 16), plot.title = element_text(hjust = 0.5))
  return(p)
}

plot_opto_uni_mixside <- function(){
  out = read.csv(file.path(git_root,"csv", "opto_unilateral.csv"), na.strings=c("","NA"))
  df <- preprocessOpto(out)
  free_df <- df %>% filter(forced == 0) %>% group_by(subjid) %>%
    mutate(lottery_side = unique(lottery_poke)[!is.na(unique(lottery_poke))]) %>%
    ungroup() %>%
    mutate(pro_right_delta_ev = case_when(lottery_side == 'BotL' ~ -delta_ev,
                                          lottery_side == 'BotR' ~ delta_ev),
           choose_right = case_when(lottery_side == 'BotL' ~ 1 - choice,
                                    lottery_side == 'BotR' ~ choice))
  free_df$isopto = factor(free_df$isopto)
  #model prediction
  m = glm(choice ~ delta_ev * isopto, data = free_df, family = binomial)
  pred_df = data_grid(free_df, delta_ev = seq(min(delta_ev) - 20, max(delta_ev) + 20, by = 1), isopto)
  pred = predict(m, pred_df, type = 'response', allow.new.levels=TRUE, se.fit = TRUE)
  pred_df$pred = pred$fit
  pred_df$ymin = pred$fit - pred$se.fit
  pred_df$ymax = pred$fit + pred$se.fit
  p = ggplot(free_df) + theme_classic(base_size = 15) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed',  alpha = 0.4) +
    scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1)) +
    scale_x_continuous(breaks = c(round(min(free_df$delta_ev)/10)*10, 0, round(max(free_df$delta_ev)/10)*10)) +
    annotate("text", label = sprintf('%d trials \n %d sessions', dim(free_df)[1], length(unique(df$sessid))), x = max(free_df$delta_ev)*0.6, y = 0.1, size = ANNOTATION_SIZE) +
    xlab(expression(EV[lottery]-EV[surebet])) + ylab("P(Chose Lottey)") +
    stat_summary_bin(mapping = aes(x = delta_ev, y = choice, color = isopto), bins = 4, fun.data = bino, geom = 'pointrange', position = position_dodge(10)) +
    #geom_line(pred_df, mapping = aes(x = pro_right_delta_ev, y = pred, color = opto_contr_ipsi)) +
    geom_ribbon(pred_df, mapping = aes(x = delta_ev, ymin = ymin, ymax = ymax, fill = isopto), alpha = 0.5) +
    scale_color_manual('Laser', values = c('0' = 'azure4', '1' = '#D766C9'), labels = c("Off", "On")) +
    scale_fill_manual('Laser', values = c('0' = 'azure4', '1' = '#D766C9'), labels = c("Off", "On")) +
    ggtitle('Opto unilateral FOF') +
    theme(legend.position = 'bottom', axis.text = element_text(size = 16), plot.title = element_text(hjust = 0.5))
  return(p)
}

## check RT in bilateral opto data
# Chaofei 2023.4.12
plot_optoRT_bi <- function(){
  git_root <- find_root(is_git_root)
  out <- read.csv(file.path(git_root, "csv", "opto_bilateral.csv"), na.strings=c("","NA"))
  df <- preprocessOpto(out)
  free_df <- df %>%
    mutate(log_rt = log(RT), subjid = factor(df$subjid, levels = unique(df$subjid)),
                                                   opto_contr_ipsi = factor(df$opto_contr_ipsi, levels = unique(df$opto_contr_ipsi)),
                                                   choice = factor(df$choice)) %>%
    filter(forced == 0 & RT < 3 & !is.na(RT))

  m_rt_opto <- lmer(log_rt ~ delta_ev * opto_contr_ipsi * choice + (delta_ev * opto_contr_ipsi * choice|subjid), data = free_df)
  free_df$pred = predict(m_rt_opto)

  for (i in 1:2){
  p = ggplot(free_df %>% filter(choice == c(0,1)[i]), aes_string(x = 'delta_ev', y = 'log_rt', color = 'opto_contr_ipsi', fill = 'opto_contr_ipsi')) + theme_classic(base_size = 15) +
    scale_y_continuous(limits = c(-1.5, 1.5), breaks = c(-1.5, 0, 1.5)) +
    scale_x_continuous(breaks = c(0, round(max(free_df$delta_ev)/10)*10)) +
    xlab(expression(EV[lottery]-EV[surebet])) + ylab("log(RT)") +
    ggtitle(c('Choose surebet', 'Choose lottery')[i]) +
    stat_summary_bin(fun.data = mean_se, geom = 'pointrange', bins = 5, position = position_dodge(10)) +
    geom_line(mapping = aes_string(x = 'delta_ev', y = 'pred', color = 'opto_contr_ipsi')) +
    scale_color_manual('Laser',values = c(no_opto = 'azure4', bilateral = 'purple3'), labels = c('Off', 'On')) +
    scale_fill_manual('Laser',values = c(no_opto = 'azure4', bilateral = 'purple3'), labels = c('Off', 'On')) +
    theme(legend.position = 'bottom', axis.text.x = element_text(size = 11)) +
    facet_wrap(~subjid)
  if (i == 1){P = p} else{P = P + p}
  }
  return(P)
}

## check RT in unilateral opto data
plot_optoRT_uni <- function(){
  git_root <- find_root(is_git_root)
  out <- read.csv(file.path(git_root, "csv", "opto_unilateral.csv"), na.strings=c("","NA"))
  df <- preprocessOpto(out)
  free_df <- df %>% group_by(subjid) %>%
    mutate(lottery_side = unique(lottery_poke)[!is.na(unique(lottery_poke))]) %>%
    ungroup() %>%
    mutate(pro_right_delta_ev = case_when(lottery_side == 'BotL' ~ -delta_ev,
                                          lottery_side == 'BotR' ~ delta_ev),
           choose_right = case_when(lottery_side == 'BotL' ~ 1 - choice,
                                    lottery_side == 'BotR' ~ choice)) %>%
    mutate(log_rt = log(RT), subjid = factor(df$subjid, levels = unique(df$subjid)),
           opto_contr_ipsi = factor(df$opto_contr_ipsi, levels = unique(df$opto_contr_ipsi)),
           choice = factor(df$choice)) %>%
    filter(forced == 0 & RT < 3 & !is.na(RT))

  m_rt_opto <- lmer(log_rt ~ pro_right_delta_ev * opto_contr_ipsi * choice + (pro_right_delta_ev * opto_contr_ipsi * choice|subjid), data = free_df)
  free_df$pred = predict(m_rt_opto)

  for (i in 1:2){
    p = ggplot(free_df %>% filter(choice == c(0,1)[i]), aes_string(x = 'pro_right_delta_ev', y = 'log_rt', color = 'opto_out', fill = 'opto_out')) + theme_classic(base_size = 15) +
      scale_y_continuous(limits = c(-1.5, 1.5), breaks = c(-1.5, 0, 1.5)) +
      scale_x_continuous(breaks = c(round(min(free_df$pro_right_delta_ev)/10)*10, 0, round(max(free_df$pro_right_delta_ev)/10)*10)) +
      xlab(expression(EV[right]-EV[left])) + ylab("log(RT)") +
      ggtitle(c('Choose surebet', 'Choose lottery')[i]) +
      stat_summary_bin(fun.data = mean_se, geom = 'pointrange', bins = 5, position = position_dodge(10)) +
      geom_line(mapping = aes_string(x = 'pro_right_delta_ev', y = 'pred', color = 'opto_out')) +
      scale_color_manual('Laser side', values = c(no_opto = 'azure4', opto_left = '#2a6a99', opto_right='#d88546'), labels = c('Off', 'L', 'R')) +
      scale_fill_manual('Laser side', values = c(no_opto = 'azure4', opto_left = '#2a6a99', opto_right='#d88546'), labels = c('Off', 'L', 'R')) +
      theme(legend.position = 'bottom', axis.text.x = element_text(size = 11)) +
      facet_wrap(~subjid)
    if (i == 1){P = p} else{P = P + p}
  }
  # ggsave("~/Desktop/uni_rt_opto.png", plot = P, width = 3, units = "in", height = 2, scale = 4)
  # ggsave("~/Desktop/uni_rt_opto.pdf", plot = P, width = 3, units = "in", height = 2, scale = 4)
  return(P)
}

## plot individual data by using mixed linear model
plot_individual_opto_bi <- function(){
  git_root <- find_root(is_git_root)
  out = read.csv(file.path(git_root, "csv","opto_bilateral.csv"))
  df <- preprocessOpto(out)
  free_df <- df %>% filter(forced == 0)
  free_df$subjid <- factor(free_df$subjid, levels = unique(free_df$subjid))
  ntrials_df = free_df %>% group_by(subjid) %>% tally() %>%
    mutate(label = paste0('n=', n), x = max(free_df$delta_ev)*0.75,  y = 0.1)

  #model prediction
  m = glmer(choice ~ delta_ev * opto_contr_ipsi + (delta_ev * opto_contr_ipsi|subjid), data = free_df, family = binomial)
  pred_df = data_grid(free_df, delta_ev = seq(min(delta_ev) - 20, max(delta_ev) + 20, by = 1), opto_contr_ipsi = opto_contr_ipsi, subjid)
  pred = predict(m, pred_df, type = 'response', allow.new.levels=TRUE, se.fit = TRUE)
  pred_df$pred = pred

  p = ggplot(free_df) + theme_classic(base_size = 15) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed',  alpha = 0.4) +
    scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1)) +
    scale_x_continuous(limits = c(-70, 240), breaks = c(round(min(free_df$delta_ev)/10)*10-10, 0, round(max(free_df$delta_ev)/10)*10)) +
    geom_text(data = ntrials_df, mapping = aes(x = x, y = y, label = label)) +
    xlab(expression(EV[lottery]-EV[surebet])) + ylab("P(Chose Lottery)") +
    stat_summary_bin(mapping = aes(x = delta_ev, y = choice, color = opto_contr_ipsi), bins = 6, fun.data = bino, geom = 'pointrange', position = position_dodge(10)) +
    geom_line(pred_df, mapping = aes(x = delta_ev, y = pred, color = opto_contr_ipsi), alpha = 0.5, size = 1) +
    scale_color_manual('Laser',values = c(no_opto = 'azure4', bilateral = 'purple3'), labels = c('Off', 'On')) +
    scale_fill_manual('Laser',values = c(no_opto = 'azure4', bilateral = 'purple3'),, labels = c('Off', 'On')) +
    ggtitle('Opto bilateral FOF') +
    facet_wrap(.~subjid) +
    theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5))

  # ggsave("~/Desktop/individual_opto.png", plot = p, width = 3, units = "in", height = 2, scale = 4)
  # ggsave("~/Desktop/individual_opto.pdf", plot = p, width = 3, units = "in", height = 2, scale = 4)
  return(p)
}

## plot individual data by using mixed linear model
plot_individual_opto_uni <- function(){
  git_root <- find_root(is_git_root)
  out <- read.csv(file.path(git_root, "csv", "opto_unilateral.csv"))
  free_df <- preprocessOpto(out) %>% filter(forced == 0 & RT < 3 & !is.na(RT)) %>% group_by(subjid) %>%
    mutate(lottery_side = unique(lottery_poke)[!is.na(unique(lottery_poke))]) %>%
    ungroup() %>%
    mutate(pro_right_delta_ev = case_when(lottery_side == 'BotL' ~ -delta_ev,
                                          lottery_side == 'BotR' ~ delta_ev),
           choose_right = case_when(lottery_side == 'BotL' ~ 1 - choice,
                                    lottery_side == 'BotR' ~ choice))
  free_df$subjid <- factor(free_df$subjid, levels = unique(free_df$subjid))
  ntrials_df = free_df %>% group_by(subjid) %>% tally() %>%
    mutate(label = paste0('n=', n), x = max(free_df$delta_ev)*0.75,  y = 0.1)

  #model prediction
  m = glmer(choose_right ~ pro_right_delta_ev * opto_out + (pro_right_delta_ev * opto_out|subjid), data = free_df, family = binomial)
  pred_df = data_grid(free_df, pro_right_delta_ev = seq(min(pro_right_delta_ev) - 20, max(pro_right_delta_ev) + 20, by = 1), opto_out, subjid)
  pred = predict(m, pred_df, type = 'response')
  pred_df$pred = pred
  p = ggplot(free_df) + theme_classic(base_size = 15) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed',  alpha = 0.4) +
    scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1)) +
    scale_x_continuous(breaks = c(round(min(free_df$pro_right_delta_ev)/10)*10, 0, round(max(free_df$pro_right_delta_ev)/10)*10)) +
    geom_text(data = ntrials_df, mapping = aes(x = x, y = y, label = label)) +
    xlab(expression(EV[right]-EV[left])) + ylab("P(Chose Right)") +
    stat_summary_bin(mapping = aes(x = pro_right_delta_ev, y = choose_right, color = opto_out), bins = 6, fun.data = bino, geom = 'pointrange', position = position_dodge(10)) +
    geom_line(pred_df, mapping = aes(x = pro_right_delta_ev, y = pred, color = opto_out), alpha = 0.5, size = 1) +
    scale_color_manual('Laser side', values = c(no_opto = 'azure4', opto_left = '#2a6a99', opto_right='#d88546'), labels = c('Off', 'L', 'R')) +
    scale_fill_manual('Laser side', values = c(no_opto = 'azure4', opto_left = '#2a6a99', opto_right='#d88546'), labels = c('Off', 'L', 'R')) +
    ggtitle('Opto unilateral FOF') +
    facet_wrap(.~subjid) + theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5))

  # ggsave("~/Desktop/individual_opto.png", plot = p, width = 3, units = "in", height = 2, scale = 4)
  # ggsave("~/Desktop/individual_opto.pdf", plot = p, width = 3, units = "in", height = 2, scale = 4)
  return(p)
}

plot_individual_opto_ev_uni <- function(){
  git_root <- find_root(is_git_root)
  out <- read.csv(file.path(git_root, "csv", "opto_unilateral.csv"))
  free_df <- preprocessOpto(out) %>% filter(forced == 0 & RT < 3 & !is.na(RT)) %>% group_by(subjid) %>%
    mutate(lottery_side = unique(lottery_poke)[!is.na(unique(lottery_poke))]) %>%
    ungroup() %>%
    mutate(pro_right_delta_ev = case_when(lottery_side == 'BotL' ~ -delta_ev,
                                          lottery_side == 'BotR' ~ delta_ev),
           choose_right = case_when(lottery_side == 'BotL' ~ 1 - choice,
                                    lottery_side == 'BotR' ~ choice))
  free_df$subjid <- factor(free_df$subjid, levels = unique(free_df$subjid))
  free_df$isopto <- factor(free_df$isopto)
  ntrials_df = free_df %>% group_by(subjid) %>% tally() %>%
    mutate(label = paste0('n=', n), x = max(free_df$delta_ev)*0.75,  y = 0.1)

  #model prediction
  m = glmer(choice ~ delta_ev * isopto + (delta_ev * isopto|subjid), data = free_df, family = binomial)
  pred_df = data_grid(free_df, delta_ev = seq(min(delta_ev) - 20, max(delta_ev) + 20, by = 1), isopto, subjid)
  pred = predict(m, pred_df, type = 'response')
  pred_df$pred = pred
  p = ggplot(free_df) + theme_classic(base_size = 15) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed',  alpha = 0.4) +
    scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1)) +
    scale_x_continuous(breaks = c(round(min(free_df$delta_ev)/10)*10, 0, round(max(free_df$delta_ev)/10)*10)) +
    geom_text(data = ntrials_df, mapping = aes(x = x, y = y, label = label)) +
    xlab(expression(EV[lottery]-EV[surebet])) + ylab("P(Chose Lottery)") +
    stat_summary_bin(mapping = aes(x = delta_ev, y = choice, color = isopto), bins = 6, fun.data = bino, geom = 'pointrange', position = position_dodge(10)) +
    geom_line(pred_df, mapping = aes(x = delta_ev, y = pred, color = isopto), alpha = 0.5, size = 1) +
    scale_color_manual('Laser', values = c('0' = 'azure4', '1' = '#D766C9'), labels = c("Off", "On")) +
    scale_fill_manual('Laser', values = c('0' = 'azure4', '1' = '#D766C9'), labels = c("Off", "On")) +
    ggtitle('Opto unilateral FOF') +
    facet_wrap(.~subjid) + theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5))

  # ggsave("~/Desktop/individual_opto.png", plot = p, width = 3, units = "in", height = 2, scale = 4)
  # ggsave("~/Desktop/individual_opto.pdf", plot = p, width = 3, units = "in", height = 2, scale = 4)
  return(p)
}
