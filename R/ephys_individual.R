
figure_ephys_supp = function(){
  pb_2224 <- figure_behavior_example(2224)
  pg_2224 <- figure_ephysGlm_scatter(2224)
  pb_2238 <- figure_behavior_example(2238)
  pg_2238 <- figure_ephysGlm_scatter(2238)
  pb_2244 <- figure_behavior_example(2244)
  pg_2244 <- figure_ephysGlm_scatter(2244)
  pb_2263 <- figure_behavior_example(2263)
  pg_2263 <- figure_ephysGlm_scatter(2263)
  pb_2261 <- figure_behavior_example(2261,3,4)
  pg_2261 <- figure_ephysGlm_scatter(2261)
  pb_2264 <- figure_behavior_example(2264)
  pg_2264 <- figure_ephysGlm_scatter(2264,1)

  design <- "
  1122#3344
  5566#7788
  99aa#bbcc
  "
  p <- pb_2224 + pg_2224 + pb_2238 + pg_2238 +
    pb_2244 + pg_2244 + pb_2263 + pg_2263 +
    pb_2261 + pg_2261 + pb_2264 + pg_2264 + plot_layout(design = design)
  return(p)
}

figure_ephysGlm_scatter = function(subjid1, p_legend = 0){
  cellinfo <- read.csv(file.path(git_root,"csv", "cellinfo.csv")) %>% filter(subjid == subjid1) %>% select(cellid,rate)
  choice_df <- read.csv(file.path(git_root,"csv", "choice_df_zscore_tstat.csv")) %>% filter(start == 0.5) %>% right_join(cellinfo, by = 'cellid') %>% select(-aligned_to,-start,-end,-rate)
  lottery_df <- read.csv(file.path(git_root,"csv", "lottery_df_zscore_tstat.csv")) %>% filter(start == 0.5) %>% right_join(cellinfo, by = 'cellid') %>% select(-aligned_to,-start,-end,-rate)
  glm_df <- left_join(choice_df, lottery_df, by = 'cellid')
  glm_df <- glm_df[complete.cases(glm_df),] %>%
    mutate(grp_indx = case_when(p.x>=0.05 & p.y>=0.05 ~ 0,
                                p.x<0.05 & p.y>=0.05 ~ 1,
                                p.x>=0.05 & p.y<0.05 ~ 2,
                                p.x <0.05 & p.y<0.05 ~ 3))


  p <- ggplot(glm_df) + theme_classic(base_size = BASE_SIZE) +
    geom_point(aes(x = beta.x, y = beta.y, fill = factor(grp_indx)),
               size = 3, alpha = 0.9,shape = 21, stroke = 0.3)+
    scale_fill_manual(values = c('#ECECEC','#40ACCB','#ED866B','#78C678'),
                       labels = c("non-selective","choice-only","lottery-only","both-selective"))+
    scale_color_manual(values = c('#777777','#205766','#774436','#3C623C'))+
    annotate("text", label = sprintf('%d of recorded neurons:',nrow(glm_df)), x = (max(glm_df$zscore_beta_x)-min(glm_df$zscore_beta_x))*0.05+min(glm_df$zscore_beta_x),
             y = max(glm_df$zscore_beta_y), size = ANNOTATION_SIZE, hjust=0,fontface='bold')+
    annotate("text", label = sprintf('%d choice only',nrow(glm_df[which(glm_df$grp_indx == 1),])), x = (max(glm_df$zscore_beta_x)-min(glm_df$zscore_beta_x))*0.05+min(glm_df$zscore_beta_x),
             y = max(glm_df$zscore_beta_y)*0.85, size = ANNOTATION_SIZE, hjust=0)+
    annotate("text", label = sprintf('%d lottery only',nrow(glm_df[which(glm_df$grp_indx == 2),])), x = (max(glm_df$zscore_beta_x)-min(glm_df$zscore_beta_x))*0.05+min(glm_df$zscore_beta_x),
             y = max(glm_df$zscore_beta_y)*0.7, size = ANNOTATION_SIZE, hjust=0)+
    annotate("text", label = sprintf('%d both',nrow(glm_df[which(glm_df$grp_indx == 3),])), x = (max(glm_df$zscore_beta_x)-min(glm_df$zscore_beta_x))*0.05+min(glm_df$zscore_beta_x),
             y = max(glm_df$zscore_beta_y)*0.55, size = ANNOTATION_SIZE, hjust=0)+
    geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed',  alpha = 0.4) +
    xlab(bquote(beta[choice] (t-stat))) +
    ylab(bquote(beta[lottery] (t-stat)))

  if (p_legend){
    p <- p + theme(legend.title = element_blank(),
                   legend.position = c(0.86,0.15),
                   legend.background = element_rect(fill = 'white',color = NA),
                   legend.text = element_text(size = 11),
                   legend.key.size = unit(0.2, 'cm'))+
      guides(color=guide_legend(ncol = 1))
  }else{
    p <- p + theme(legend.position = 'none')
  }
  return(p)
}


figure_behavior_example = function(subjid1,sb_low = 2.5,sb_high = 4.5){
  # figure1c: example animal psychometric curves by GLMM, thin lines for each session, thick line for combined
  risk_df = read.csv(file.path(git_root,"csv", "ephys_beh_df.csv")) %>% filter(subjid == subjid1 & sb_mag>sb_low & sb_mag<sb_high) %>% preprocessRisk()
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
    annotate("text", label = sprintf('SubjectID: %d',subjid1), x = max(risk_df$delta_ev)*0.50, y = 0.2, size = ANNOTATION_SIZE) +
    annotate("text", label = sprintf('%d sessions', length(unique(risk_df$sessid))), x = max(risk_df$delta_ev)*0.50, y = 0.1, size = ANNOTATION_SIZE) +
    xlab(expression(EV[lottery]-EV[surebet])) + ylab("P(Chose Lottery)") +
    geom_line(avg_pred_df, mapping = aes(x = delta_ev, y = pred), color = 'black', alpha = 0.5, size = 2)
  for (sessidx in unique(risk_df$sessid)){
    p = p + geom_line(sess_pred_df %>% filter(sessid == sessidx),
                      mapping = aes(x = delta_ev, y = pred), color = 'darkgray', alpha = 0.3)
  }
  p = p + stat_summary_bin(mapping = aes(x = delta_ev, y = choice), fun.data = bino, geom = 'pointrange', color = 'black')
  return(p)
}
