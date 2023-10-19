# functions for infusion part of figure2 (figure_infusion)

#' figure_infusion_bi
#'
#' figure2b & 2e: plots psychometric curves (model prediction + binned data), color by bilateral infusion and control
#'
#' @param region: 'FOF' or 'PPC'
#'
#' @return ggplot object
plot_infusion_bi = function(region = 'FOF', trial_num = 'all'){
  # figure2b & 2e: plots psychometric curves (model prediction + binned data), color by bilateral infusion and control
  if (region == 'FOF'){
    df_m = read.csv(file.path(git_root,"csv", "figure_infusion_bi_fof.csv"))
    dosage_colors = FOF_dosage_colors
  } else{
    if (trial_num == 'all'){
      df_m = read.csv(file.path(git_root,"csv", "figure_infusion_bi_ppc.csv"))
    }else{
      df_m = read.csv(file.path(git_root,"csv", "figure_infusion_bi_ppc.csv")) %>% filter(trialnum<=trial_num)
    }
    dosage_colors = PPC_dosage_colors
  }
  control_df = read.csv(file.path(git_root,"csv", "figure_behavior_population.csv"))
  df = rbind(control_df, df_m)

  # Use GLM to plot ribbons
  m0 = glm('choice ~ delta_ev*dosage', df, family = binomial)
  pop_pred_df = data_grid(df, delta_ev = seq(min(delta_ev) - 20, max(delta_ev) + 20, by = 1), dosage = dosage)
  pop_pred = predict(m0, pop_pred_df, type = 'response', allow.new.levels=TRUE, se.fit = TRUE)
  pop_pred_df$pred = pop_pred$fit
  pop_pred_df$ymin = pop_pred$fit - pop_pred$se.fit
  pop_pred_df$ymax = pop_pred$fit + pop_pred$se.fit

  # plot
  df$dosage = factor(df$dosage)
  pop_pred_df$dosage = factor(pop_pred_df$dosage)
  p = ggplot(df) + theme_classic(base_size = BASE_SIZE) +
    stat_summary_bin(mapping = aes(x = delta_ev, y = choice, color = dosage), bins = 7, fun.data = bino, position = position_dodge(10)) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed',  alpha = 0.4) +
    scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1)) +
    scale_x_continuous(limits = c(-80, 282), breaks = c(-60, 0, 260)) +
    annotate("text", label = sprintf('%d trials \n %d sessions', dim(df_m)[1], length(unique(df_m$sessid))), x = max(df$delta_ev)*0.65, y = 0.1, size = ANNOTATION_SIZE) +
    xlab(expression(EV[lottery]-EV[surebet])) + ylab("P(Chose Lottery)") +
    geom_ribbon(pop_pred_df, mapping = aes(x = delta_ev, ymin = ymin, ymax = ymax, fill = dosage), alpha = 0.5) +
    scale_color_manual(values = dosage_colors) + scale_fill_manual(values = dosage_colors) +
    labs(fill = 'Dose(µg)', color = 'Dose(µg)') + theme(legend.position = 'bottom', axis.text = element_text(size = 16))

  if (region == 'FOF'){
    p <- p + ggtitle('Muscimol bilateral FOF') + theme(plot.title = element_text(hjust = 0.5))
  }else{
    p <- p + ggtitle('Muscimol bilateral PPC') + theme(plot.title = element_text(hjust = 0.5))
  }
  return(p)
}

#' figure_infusion_uni
#'
#' figure2c & 2f: plots psychometric curves (model prediction + binned data), color by unilateral infusion and control
#'
#' @param region = 'FOF' or 'PPC'
#'
#' @return ggplot object
plot_infusion_uni = function(region = 'FOF', trial_num = 'all', n_bins = 7){
  # figure2c & 2f: plots psychometric curves (model prediction + binned data), color by unilateral infusion and control
  if (region == 'FOF'){
    df_m = read.csv(file.path(git_root,"csv", "figure_infusion_uni_fof.csv"))
  } else{
    if (trial_num == 'all'){
      df_m = read.csv(file.path(git_root,"csv", "figure_infusion_uni_ppc.csv"))
    }else{
      df_m = read.csv(file.path(git_root,"csv", "figure_infusion_uni_ppc.csv")) %>% filter(trialnum<=trial_num)
    }
  }
  control_df = read.csv(file.path(git_root,"csv", "figure_behavior_population.csv"))
  df = rbind(control_df, df_m)
  df$infusion_side = factor(df$infusion_side)

  # Use GLM fit for plotting ribbons
  m0 = glm('choose_right ~ pro_right_delta_ev*infusion_side', df, family = binomial)
  pop_pred_df = data_grid(df, pro_right_delta_ev = seq(min(pro_right_delta_ev) - 20, max(pro_right_delta_ev) + 20, by = 1), infusion_side = infusion_side)
  pop_pred = predict(m0, pop_pred_df, type = 'response', allow.new.levels=TRUE, se.fit = TRUE)
  pop_pred_df$pred = pop_pred$fit
  pop_pred_df$ymin = pop_pred$fit - pop_pred$se.fit
  pop_pred_df$ymax = pop_pred$fit + pop_pred$se.fit

  # plot
  p = ggplot(df) + theme_classic(base_size = BASE_SIZE) +
    stat_summary_bin(mapping = aes(x = pro_right_delta_ev, y = choose_right, color = infusion_side), bins = n_bins, fun.data = bino, position = position_dodge(10)) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed',  alpha = 0.4) +
    scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1)) +
    scale_x_continuous(limits = c(-130, 259), breaks = c(-60, 0, 240)) +
    annotate("text", label = sprintf('%d trials \n %d sessions', dim(df_m)[1], length(unique(df_m$sessid))), x = max(df$delta_ev)*0.65, y = 0.1, size = ANNOTATION_SIZE) +
    xlab(expression(EV[right]-EV[left])) + ylab("P(Chose Right)") +
    geom_ribbon(pop_pred_df, mapping = aes(x = pro_right_delta_ev, ymin = ymin, ymax = ymax, fill = infusion_side), alpha = 0.5) +
    scale_color_manual(values = side_colors, labels = c("Control", "L", "R")) + scale_fill_manual(values = side_colors, labels = c("Control", "L", "R")) +
    labs(fill = 'Side', color = 'Side') + theme(legend.position = 'bottom', axis.text = element_text(size = 16))

  if (region == 'FOF'){
    p <- p + ggtitle('Muscimol unilateral FOF') + theme(plot.title = element_text(hjust = 0.5))
  }else{
    p <- p + ggtitle('Muscimol unilateral PPC') + theme(plot.title = element_text(hjust = 0.5))
  }
  return(p)
}


#' figure_infusion_uni_ev_mixside
#'
#' figure2c & 2f: plots psychometric curves (model prediction + binned data), color by unilateral infusion and control
#'
#' @param region = 'FOF' or 'PPC'
#'
#' @return ggplot object
plot_infusion_uni_ev_mixside = function(region = 'FOF', trial_num = 'all'){
  # figure2c & 2f: plots psychometric curves (model prediction + binned data), color by unilateral infusion and control
  if (region == 'FOF'){
    df_m = read.csv(file.path(git_root,"csv", "figure_infusion_uni_fof.csv"))
  } else{
    if (trial_num == 'all'){
      df_m = read.csv(file.path(git_root,"csv", "figure_infusion_uni_ppc.csv"))
    }else{
      df_m = read.csv(file.path(git_root,"csv", "figure_infusion_uni_ppc.csv")) %>% filter(trialnum<=trial_num)
    }
  }
  control_df = read.csv(file.path(git_root,"csv", "figure_behavior_population.csv"))
  df = rbind(control_df, df_m)
  df$infusion_side = factor(df$infusion_side)
  df$drug <- tolower(df$drug)
  df$drug = factor(df$drug)

  # Use GLM fit for plotting ribbons
  m0 = glm('choice ~ delta_ev*drug', df, family = binomial)
  pop_pred_df = data_grid(df, delta_ev = seq(min(delta_ev) - 20, max(delta_ev) + 20, by = 1), drug = drug)
  pop_pred = predict(m0, pop_pred_df, type = 'response', allow.new.levels=TRUE, se.fit = TRUE)
  pop_pred_df$pred = pop_pred$fit
  pop_pred_df$ymin = pop_pred$fit - pop_pred$se.fit
  pop_pred_df$ymax = pop_pred$fit + pop_pred$se.fit

  # plot
  p = ggplot(df) + theme_classic(base_size = BASE_SIZE) +
    stat_summary_bin(mapping = aes(x = delta_ev, y = choice, color = drug), bins = 7, fun.data = bino, position = position_dodge(10)) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed',  alpha = 0.4) +
    scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1)) +
    scale_x_continuous(limits = c(-130, 259), breaks = c(-60, 0, 240)) +
    annotate("text", label = sprintf('%d trials \n %d sessions', dim(df_m)[1], length(unique(df_m$sessid))), x = max(df$delta_ev)*0.65, y = 0.1, size = ANNOTATION_SIZE) +
    xlab(expression(EV[lottery]-EV[surebet])) + ylab("P(Chose Lottery)") +
    geom_ribbon(pop_pred_df, mapping = aes(x = delta_ev, ymin = ymin, ymax = ymax, fill = drug), alpha = 0.5) +
    scale_color_manual('Dose(µg)',values = c('control' = 'azure4', 'muscimol' = '#D766C9'), labels = c("0", "0.3")) +
    scale_fill_manual('Dose(µg)',values = c('control' = 'azure4', 'muscimol' = '#D766C9'), labels = c("0", "0.3")) +
    labs(fill = 'Side', color = 'Side') + theme(legend.position = 'bottom', axis.text = element_text(size = 16))

  if (region == 'FOF'){
    p <- p + ggtitle('Muscimol unilateral FOF') + theme(plot.title = element_text(hjust = 0.5))
  }else{
    p <- p + ggtitle('Muscimol unilateral PPC') + theme(plot.title = element_text(hjust = 0.5))
  }
  return(p)
}
