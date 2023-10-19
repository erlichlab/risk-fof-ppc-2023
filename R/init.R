# load the libraries and source all functions

library(rprojroot)
library(dplyr)
library(ggplot2)
#library(ggimage)
library(lme4)
library(modelr)
library(brms)
library(loo)
library(tidyr)
library(stringi)
library(stringr)
library(tidybayes)
library(bayesplot)
library(patchwork)
library(GGally)
library(scales)
library(rstan)


rstan_options(auto_write = TRUE) # cache compiled models
options(mc.cores = parallel::detectCores())  # detect cpu cores
git_root <- find_root(is_git_root)

## source functions
source(file.path(git_root, "R","utils.R"))
source(file.path(git_root, "R","figure_behavior.R"))
source(file.path(git_root, "R","figure_infusion.R"))
source(file.path(git_root, "R","figure_opto.R"))
source(file.path(git_root, "R", "brms_functions.R"))
source(file.path(git_root, "R", "figure_perturbation.R"))
source(file.path(git_root, "R", "figureS.R"))
source(file.path(git_root, "R", "figure_control.R"))
source(file.path(git_root, "R", "figure_learning.R"))
source(file.path(git_root, "R", "ephys_individual.R"))


## functions for generate figures

# layout for fig1
make_fig1 <- function(){
  #svg_fig1a <- image_read(file.path(git_root, "R", 'behavior.svg'), density = 2400)
  #fig1a <- image_ggplot(svg_fig1a)
  fig1b <- plot_behavior_timeline()
  fig1c <- plot_behavior_example()
  fig1d <- plot_behavior_population()

  fig1 <- plot_spacer() / fig1b / (fig1c + fig1d) + plot_layout(heights = c(1.1, 1, 1), guides = 'keep') #+ plot_annotation(tag_levels = 'a')
  ggsave(file.path(git_root,"plots", "fig1.png"), plot = fig1, width = 2.5, units = "in", height = 3.75, scale = 3.2)
  return(fig1)
}

# layout for fig2
make_fig2 <- function(){
  #svg_fig2a <- image_read_svg(file.path(git_root, "R",'brain.svg'), density = 2400)
  fig2b <- plot_infusion_bi('PPC')
  fig2c <- plot_infusion_bi('FOF')
  fig2d <- plot_opto_bi()
  fig2e <- plot_infusion_uni_ev_mixside('PPC')
  fig2f <- plot_infusion_uni_ev_mixside('FOF')
  fig2g <- plot_opto_uni_mixside()
  fig2h <- plot_infusion_uni('PPC')
  fig2i <- plot_infusion_uni('FOF')
  fig2j <- plot_opto_uni()
  fig2 <- plot_spacer() + fig2b + fig2c + fig2d +
    plot_spacer() + fig2e + fig2f + fig2g +
    plot_spacer() + fig2h + fig2i + fig2j + plot_layout(nrow = 3)
  ggsave(file.path(git_root,"plots", "fig2.png"), plot = fig2, width = 4, units = "in", height = 4, scale = 3.7)
  return(fig2)
}

# plots for fig3
make_fig3 <- function(){
  layout <- "
  ACDEFG
  BCDEFG
  "
  p1 <- plot_pred_sample_individual(2228, bidf_opto, bi_opto_cdf, legend_on = FALSE, x_lab_on = FALSE)
  p2 <- plot_pred_sample_individual(2160, bidf_fof_muscimol, bi_fof_muscimol_cdf, legend_on = FALSE)
  p3 <- plot_delta_rho_psd(bi_opto_cdf, uni_opto_cdf, bi_fof_muscimol_cdf, uni_fof_muscimol_cdf)
  p4 <- plot_delta_noise_psd(bi_opto_cdf, uni_opto_cdf, bi_fof_muscimol_cdf, uni_fof_muscimol_cdf, y_lab = F)
  p5 <- plot_delta_omega_rational_psd(bi_opto_cdf, uni_opto_cdf, bi_fof_muscimol_cdf, uni_fof_muscimol_cdf, y_lab = F)
  p6 <- plot_delta_omega_lottery_psd(bi_opto_cdf, uni_opto_cdf, bi_fof_muscimol_cdf, uni_fof_muscimol_cdf, y_lab = F)
  p7 <- plot_delta_omega_surebet_psd(bi_opto_cdf, uni_opto_cdf, bi_fof_muscimol_cdf, uni_fof_muscimol_cdf, y_lab = F)
  fig3 <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + guide_area() + plot_layout(design = layout, widths = c(1,.8,.8,.8,.8,.8))
  ggsave(file.path(git_root,"plots", "fig3.png"), plot = fig3, width = 2.5, units = "in", height = 0.9, scale = 3.5)
  return(fig3)
}

# layout for fig4 is in julia folder

# layout for fig5 is in matlab folder

# plots for fig6
make_fig6g <- function(){
  p1 <- plot_control_dot_line()
  p2 <- plot_control_risky_ppc()
  return(p1+p2)
}


## Plots for Extended Data Figures

# plots for extended data fig 1 is in matlab and julia folder

# plots for extended data fig2
make_ED_fig2 <- function(){
  pa <- plot_individual_infusion('bi_ppc')
  pb <- plot_individual_infusion('bi_fof')
  pc <- plot_individual_opto_bi()
  pd <- plot_individual_infusion_ev_uni('uni_ppc')
  pe <- plot_individual_infusion_ev_uni('uni_fof')
  pf <- plot_individual_opto_ev_uni()
  pg <- plot_individual_infusion('uni_ppc')
  ph <- plot_individual_infusion('uni_fof')
  pi <- plot_individual_opto_uni()
  ED_fig2 <- pa + pb + pc + pd + pe + pf + pg + ph + pi + plot_layout(ncol = 3, guides = 'keep')
  return(ED_fig2)
}

# plots for extended data fig3
make_ED_fig3a_left <- function(){
  plot_bi_opto_individual <- plot_pred_individual(bidf_opto, bi_opto_cdf)
  ggsave(file.path(git_root,"plots", "ED_fig3a_left.png"), plot = plot_bi_opto_individual, width = 2.2, units = "in", height = 2, scale = 2)
  return(plot_bi_opto_individual)
}

make_ED_fig3a_right <- function(){
  bi_opto_fof_posterior <- plot_population_draw_dens(bi_opto_cdf, 'bi_opto')
  return(bi_opto_fof_posterior)
}

make_ED_fig3b_left <- function(){
  plot_uni_opto_individual <- plot_pred_individual(unidf_opto, uni_opto_cdf)
  ggsave(file.path(git_root,"plots", "ED_fig3b_left.png"), plot = plot_uni_opto_individual, width = 2.2, units = "in", height = 2.7, scale = 2)
  return(plot_uni_opto_individual)
}

make_ED_fig3b_right <- function(){
  uni_opto_fof_posterior <- plot_population_draw_dens(uni_opto_cdf, 'uni_opto')
  return(uni_opto_fof_posterior)
}

# plots for extended data fig4
make_ED_fig4a_left <- function(){
  plot_bi_fof_muscimol_individual <- plot_pred_individual(bidf_fof_muscimol, bi_fof_muscimol_cdf)
  ggsave(file.path(git_root,"plots", "ED_fig4a_left.png"), plot = plot_bi_fof_muscimol_individual, width = 2.2, units = "in", height = 2.7, scale = 2)
  return(plot_bi_fof_muscimol_individual)
}

make_ED_fig4a_right <- function(){
  bi_muscimol_fof_posterior <- plot_population_draw_dens(bi_fof_muscimol_cdf, 'bi_muscimol')
  return(bi_muscimol_fof_posterior)
}

make_ED_fig4b_left <- function(){
  plot_uni_fof_muscimol_individual <- plot_pred_individual(unidf_fof_muscimol, uni_fof_muscimol_cdf)
  ggsave(file.path(git_root,"plots", "ED_fig4b_left.png"), plot = plot_uni_fof_muscimol_individual, width = 2.2, units = "in", height = 2.7, scale = 2)
  return(plot_uni_fof_muscimol_individual)
}

make_ED_fig4b_right <- function(){
  uni_muscimol_fof_posterior <- plot_population_draw_dens(uni_fof_muscimol_cdf, 'uni_muscimol')
  return(uni_muscimol_fof_posterior)
}

# plots for extended data fig8
make_ED_fig8b <- function(){
  p1 <- plot_infusion_bi('PPC', 40)
  p2 <- plot_infusion_uni_ev_mixside('PPC', 40)
  p3 <- plot_infusion_uni('PPC', 40, 4)
  p4 <- plot_ppc_cutoff_p()
  ED_fig8b <- p1 + p2 + p3 + p4 + plot_layout(nrow = 1)
  return(ED_fig8b)
}

## functions for supplementary figs

# plots for supplementary fig2
make_supp_fig2 <- function(){
  p1 <- plot_infusion_timeline(2152)
  p2 <- plot_infusion_timeline(2153)
  p3 <- plot_infusion_timeline(2154)
  p4 <- plot_infusion_timeline(2155)
  p5 <- plot_infusion_timeline(2156)
  p6 <- plot_infusion_timeline(2160)
  p7 <- plot_infusion_timeline(2165)
  p8 <- plot_infusion_timeline(2166)
  supp_fig2 <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + plot_layout(ncol = 2, widths = c(1, 1), heights = c(0.5,0.5,0.5,0.5))
  return(supp_fig2)
}

# plots for supplementary fig4
make_supp_fig4 <- function(){
  pa <- plot_supp_RT('bi_ppc') + plot_annotation(title = "Muscimol bilateral PPC", theme = theme(plot.title = element_text(hjust = 0.5, size = 16))) +
    plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
  pb <- plot_supp_RT('uni_ppc') + plot_annotation(title = "Muscimol unilateral PPC", theme = theme(plot.title = element_text(hjust = 0.5, size = 16))) +
    plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
  pc <- plot_supp_RT('bi_fof') + plot_annotation(title = "Muscimol bilateral FOF", theme = theme(plot.title = element_text(hjust = 0.5, size = 16))) +
    plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
  pd <- plot_supp_RT('uni_fof') + plot_annotation(title = "Muscimol unilateral FOF", theme = theme(plot.title = element_text(hjust = 0.5, size = 16))) +
    plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
  pe <- plot_optoRT_bi() + plot_annotation(title = "Opto bilateral FOF", theme = theme(plot.title = element_text(hjust = 0.5, size = 16))) +
    plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
  pf <- plot_optoRT_uni() + plot_annotation(title = "Opto unilateral FOF", theme = theme(plot.title = element_text(hjust = 0.5, size = 16))) +
    plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

  supp_fig4 <- (wrap_elements(pa) + wrap_elements(pb)) / (wrap_elements(pc) + wrap_elements(pd)) / (wrap_elements(pe) + wrap_elements(pf)) + plot_layout(heights = c(1, 1, 1))
}

