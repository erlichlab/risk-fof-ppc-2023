# generate the layout of figs using R for the paper
# Chaofei 2023.4.26

## Plots for Main Figures

# layout for fig1
Fig1 <- make_fig1()

# layout for fig2
Fig2 <- make_fig2()

# plots for fig3
Fig3 <- make_fig3()

# layout for fig4 is in julia folder

# layout for fig5 is in matlab folder
# run the matlab code "main_ephys_plots.m" in "matlab" folder

# plots for fig6
Fig6a <- plot_learning_schematic()
Fig6d <- plot_learning_scatter()
Fig6f <- plot_control_bias()
Fig6g <- make_fig6g()
plot_learning_individual_summary()
# for Fig6b and Fig6c, run "plot_learning_individual_summary()" and then check "plots" folder

## Plots for Extended Data Figures

# plots for extended data fig 1 is in matlab and julia folder

# plots for extended data fig2
ED_fig2 <- make_ED_fig2()

# plots for extended data fig3
ED_fig3a_left <- make_ED_fig3a_left()
ED_fig3a_right <- make_ED_fig3a_right()
ED_fig3b_left <- make_ED_fig3b_left()
ED_fig3b_right <- make_ED_fig3b_right()

# plots for extended data fig4
ED_fig4a_left <- make_ED_fig4a_left()
ED_fig4a_right <- make_ED_fig4a_right()
ED_fig4b_left <- make_ED_fig4b_left()
ED_fig4b_right <- make_ED_fig4b_right()

# plots for extended data fig5
ED_fig5a <- make_parameters_comparation_brms()
ED_fig5b <- make_brms_original_psd()
ED_fig5c <- make_mcmc_pairs_plots_verson2(bi_opto_cdf, 'bi_opto')
ED_fig5d <- make_mcmc_pairs_plots_verson2(uni_opto_cdf, 'uni_opto')
ED_fig5e <- make_mcmc_pairs_plots_verson2(bi_fof_muscimol_cdf, 'bi_fof_muscimol')''
ED_fig5f <- make_mcmc_pairs_plots_verson2(uni_fof_muscimol_cdf, 'uni_fof_muscimol')

# plots for extended data fig6 is in julia folder

# plots for extended data fig7
ED_fig7a <- figure_ephys_supp()
# To generate ephys control plots, run "lottery_sound_control.m" in "matlab" folder

# plots for extended data fig8
# for ED_fig8a, check "plots" folder
ED_fig8b <- make_ED_fig8b()
ED_fig8c <- make_main_PPC_fig()

## plots for Supplementary Figures

# plots for supplementary fig1
supp_fig1 <- plot_soft_fixation()

# plots for supplementary fig2
supp_fig2 <- make_supp_fig2()

# plots for supplementary fig4
supp_fig4 <- make_supp_fig4()

# generate the data for supplementary tables
infusion_fit_tb <- generate_individual_pars_perturb_csv()
ephys_fit_tb <- generate_individual_pars_ephys_csv()

# Statistical Appendix
# Knit the file "likelihood_ratio_tests_fig2.Rmd" to get the results

