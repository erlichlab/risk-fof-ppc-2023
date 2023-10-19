% code to generate the main figure for the ephys part
% Chaofei Bao, 2023.3.29

%% load or get the data
if ~exist(fullfile(utils.git_root_path,'matlab/main_ephys_plot_data.mat'),'file')
    dbc = db.labdb.getConnection('old_client');
    % get the spk data
    cell_ids = [21789,23347,23159];
    spkinfo = db.getCellData(cell_ids, 'full_info', 0);
    % get the behavior data
    SD = db.getSessData(spkinfo.sessid);

    % get the neuro selectivity data
    choice_glm_tb = readtable(fullfile(utils.git_root_path, "csv/choice_df_zscore_tstat.csv"));
    lottery_glm_tb = readtable(fullfile(utils.git_root_path, "csv/lottery_df_zscore_tstat.csv"));

    sessid_decode_deomo = 246697;
    model_file = fullfile(utils.git_root_path, 'matlab/lm1/', sprintf('%d_lm.mat',sessid_decode_deomo));
    f_file = fullfile(utils.git_root_path, 'matlab/features/', sprintf('%d.mat',sessid_decode_deomo));
    load(model_file);
    load(f_file);

    % load the single session decoding table
    decode_tb0 = readtable(fullfile(utils.git_root_path, 'csv/ephys_session_decode_force_b.csv'));

    % load the pseudo population decoding data
    pseudo_lines = readlines(fullfile(utils.git_root_path, 'matlab/fof_pseudo.json'));
    fof_pseudo = jsondecode(pseudo_lines);
    pseudo_correlation = readtable(fullfile(utils.git_root_path, 'csv/fof_pseudo_correlation.csv'));

    save(fullfile(utils.git_root_path,"matlab/main_ephys_plot_data.mat"));
else
    load(fullfile(utils.git_root_path,'matlab/main_ephys_plot_data.mat'))
end



%% plot the example rasters and PSTHS
FigA = figure('visible','on','Position', [0 0 1260 1782]);
%fig_title = title({sprintf('PSTHs, neurometric and firing rate tuning to deltaEV for cell %s',num2str(cell_id)); ...
%    sprintf('[%s %s] after cue onset (beta = %f, p = %f from glme for lottery)',num2str(cell_info.start),num2str(cell_info.xEnd), cell_info.beta, cell_info.p)});
%ax = gca;
%ax.Visible = "off";
%fig_title.Visible = "on";

ttest = @(x,y)ttest2(x,y,'alpha',0.01);
SD_sessids = [SD.sessid];

for i = 1:length(cell_ids)
    col_posi = [0.1,0.206]+(i-1)*0.12;
    row_posi = [0.8,0.74,0.703];

    fig_width = 0.08;
    fig_height = 0.12;
    raster_width = 0.08;
    surebet_clr = utils.linspecer(8,'blue');
    surebet_clr = surebet_clr(3:8,:);

    lottery_clr = [52,41,145;6,119,223;5,164,204;90,193,138;201,189,87;248,218,34]/255;
    choice_clr = utils.linspecer(8,'gray');
    choice_clr = {choice_clr(5,:),choice_clr(6,:)};
    lottery_clr = mat2cell(lottery_clr,[1 1 1 1 1 1]);
    lottery_clr1 = utils.linspecer(8,'red');
    neurometric_clr = [241,105,19]/255;
    scater_clr = [65,173,203;238,135,108;120,198,121]/255;

    cellid = cell_ids(i);
    sessid = spkinfo.sessid(find(spkinfo.cellid == cellid));
    this_spkinfo = spkinfo(find(spkinfo.cellid == cellid),:);

    % process the behavior data
    this_SD = SD(find(SD_sessids == sessid));
    [choice, lott_mag, events, subjid, PD] = getBehData(this_SD, 'trial_type', 'free', 'state', 'cue_target');

    % surebet = 0, lottery = 1
    choice_bino = zeros(length(choice), 1);
    choice_bino(strcmp(choice, 'surebet')) = 0;
    choice_bino(strcmp(choice, 'lottery')) = 1;

    % we are focusing on four time windows,fixation early, mid and
    % late period, and target arrival early period, each time window is 0.5s.
    new_events = zeros(length(events), 2);
    new_events(:,1) = events(:,1);
    new_events(:,2) = events(:,1)+0.5;

    beh_data = [choice_bino, lott_mag, new_events];
    beh_data_surebet = beh_data(beh_data(:,1)==0,:);
    beh_data_lottery = beh_data(beh_data(:,1)==1,:);


    % base reward for rat is 8
    base_reward = 8;

    surebet_choice = PD(beh_data(:,1)==0,:);
    surebet_choice1 = surebet_choice(1,:);
    reward_multiplier = surebet_choice1.reward_received/base_reward/surebet_choice1.sb_mag;
    delta_ev = (lott_mag*surebet_choice1.lottery_prob-surebet_choice1.sb_mag)*reward_multiplier*base_reward;

    % get the spikes counts for events*cellid*events_type
    spk_ma = get_spk_win(this_spkinfo, new_events, [0 0.5]);

    beh_tb = table(beh_data(:,1), delta_ev, 'VariableNames', {'choice', 'delta_ev'});


    m_beh = fitglm(beh_tb, 'choice ~ delta_ev', 'Link', 'logit', 'Distribution', 'binomial');
    beh_tb_predict = table();
    beh_tb_predict.delta_ev = linspace(min(delta_ev), max(delta_ev), 100)';
    [ypred_beh, yCI_beh] = predict(m_beh, beh_tb_predict);

    spk_tb1 = table(beh_data(:,1), delta_ev, spk_ma(:,1,1), 'VariableNames', {'choice', 'delta_ev', 'spk_counts'});
    spk_tb1_surebet = spk_tb1(spk_tb1.choice == 0, :);
    spk_tb1_lottery = spk_tb1(spk_tb1.choice == 1, :);

    spk_tb2 = table(beh_data(:,1), delta_ev, spk_ma(:,1,2), 'VariableNames', {'choice', 'delta_ev', 'spk_counts'});
    spk_tb2_surebet = spk_tb2(spk_tb2.choice == 0, :);
    spk_tb2_lottery = spk_tb2(spk_tb2.choice == 1, :);

    %m_spk1 = fitglme(spk_tb1,'spk_counts ~ delta_ev + (delta_ev|choice)', 'Distribution', 'Poisson','Link','identity');
    spk_predict_sb = zeros(height(beh_tb_predict),2);
    spk_predict_sb(:,2) = beh_tb_predict.delta_ev;
    spk_predict_lo = ones(height(beh_tb_predict),2);
    spk_predict_lo(:,2) = beh_tb_predict.delta_ev;
    spk_predict = [spk_predict_sb;spk_predict_lo];
    spk_predict_tb = table(spk_predict(:,1), spk_predict(:,2), 'VariableNames', {'choice', 'delta_ev'});

    %[spk_tb1_y, spk_tb1_CI] = predict(m_spk1, spk_predict_tb);
    try
        m_spk2 = fitglme(spk_tb2,'spk_counts ~ delta_ev + (delta_ev|choice)', 'Distribution', 'Poisson','Link','identity');
        [spk_tb2_y, spk_tb2_CI] = predict(m_spk2, spk_predict_tb);

        %     m_spk3 = fitglme(spk_tb3,'spk_counts ~ delta_ev + (delta_ev|choice)', 'Distribution', 'Poisson','Link','identity');
        %     [spk_tb3_y, spk_tb3_CI] = predict(m_spk3, spk_predict_tb);
    catch
        % nothing to do in this loop
        fprintf('Cellid = %s can not be fitted !!! \n',num2str(cell_info.cellid));
        continue
    end

    % Choose a window for plot
    choose_win = 2;

    switch choose_win
        case 1
            spk_tb_y = spk_tb1_y;
            spk_tb_CI = spk_tb1_CI;
            spk_tb = spk_tb1;
            choose_win_str = '0 to 0.5s after cue';
        case 2
            spk_tb_y = spk_tb2_y;
            spk_tb_CI = spk_tb2_CI;
            spk_tb = spk_tb2;
            choose_win_str = '0.5 to 1s after cue';
    end

    beh_tb_sum = grpstats(beh_tb, 'delta_ev',{'mean','sem'});
    spk_tb_sb = spk_tb(spk_tb.choice == 0,:);
    spk_tb_lott = spk_tb(spk_tb.choice == 1,:);
    spk_tb_sb_sum = grpstats(spk_tb_sb, 'delta_ev',{'mean','sem'});
    spk_tb_lott_sum = grpstats(spk_tb_lott, 'delta_ev',{'mean','sem'});

    % plot the raster and psth
    beh_data_tb = array2table(beh_data(:,1:3), 'VariableNames',{'choice','lott_mag','cue'});
    beh_data_tb.delta_ev = (beh_data_tb.lott_mag*surebet_choice1.lottery_prob-surebet_choice1.sb_mag)*reward_multiplier*base_reward;
    beh_data_tb.spk_count = spk_ma(:,1,choose_win);

    beh_sb_tb = beh_data_tb(beh_data_tb.choice == 0,:);
    grp_sb = groupsummary(beh_sb_tb,'lott_mag');
    % drop any trial type with trials num fewer than 5
    grp_sb = grp_sb(grp_sb.GroupCount>4,:);
    beh_sb_tb_new = beh_sb_tb(ismember(beh_sb_tb.lott_mag, grp_sb.lott_mag),:);
    grp_sb_sum_new = grpstats(beh_sb_tb_new, 'delta_ev',{'mean','sem'});

    beh_lo_tb = beh_data_tb(beh_data_tb.choice == 1,:);
    grp_lo = groupsummary(beh_lo_tb,'lott_mag');
    % drop any trial type with trials num fewer than 5
    grp_lo = grp_lo(grp_lo.GroupCount>4,:);
    beh_lo_tb_new = beh_lo_tb(ismember(beh_lo_tb.lott_mag, grp_lo.lott_mag),:);
    grp_lo_sum_new = grpstats(beh_lo_tb_new, 'delta_ev',{'mean','sem'});

    beh_tb_final = [beh_sb_tb_new;beh_lo_tb_new];

    if i == 1
        legend_str1 = {'surebet', 'lottery'};
        legend_pos1 = [0.43,0.8,0.1,0.03];
    else
        legend_str1 = '';
        legend_pos1 = [];
        legend_str2 = '';
        legend_pos2 = [];
    end


    [ras,R]=exampleraster_risky_2psths(beh_tb_final.cue,this_spkinfo.ts{1},'cnd',beh_tb_final.choice, ...
        'cnd2',beh_tb_final.lott_mag,'pre',.1,'post',1.1,...
        'renderer','painters','errorbars',0,'corner',[col_posi(1),row_posi(1)],...
        'x_label','Time from cue onset','clrs',choice_clr,'clrs2',lottery_clr,...
        'psth_height',fig_height*0.5,'total_height',fig_height,'ax_width',raster_width,...
        'testfunc',ttest,'cout',beh_tb_final.cue, 'font_size', 10, ...
        'krn',0.2,'TickDir','out', ...
        'font_size',8,'legend_fontsz',4);
    psth_ylim = ras(end).YLim;
    set(ras(end),'FontSize',8);
    set(ras(end-1),'FontSize',8);
    if i ~=1
        set(ras(end),'XLabel',[]);
        set(ras(end),'YLabel',[]);
        set(ras(end-1),'YLabel',[]);
    end

    % plot the fitted spkcounts into a seperat figure
    plot(draw.jaxes([col_posi(1),row_posi(2),fig_width,fig_height*0.27]),beh_tb_predict.delta_ev, spk_tb_y(1:100,:)/0.5, 'Color',choice_clr{1},'LineWidth',1,'LineStyle','-.');
    hold on;
    hes = errorbar(grp_sb_sum_new.delta_ev-2, grp_sb_sum_new.mean_spk_count/0.5, grp_sb_sum_new.sem_spk_count/0.5, ...
        'Color',choice_clr{1},'Marker','.','MarkerSize',10);
    hes.CapSize = 0;
    hes.LineWidth = 1.3;
    hes.LineStyle = 'none';
    plot(beh_tb_predict.delta_ev, spk_tb_y(101:200,:)/0.5,'Color',choice_clr{2},'LineWidth',1,'LineStyle','-');
    hel = errorbar(grp_lo_sum_new.delta_ev+2, grp_lo_sum_new.mean_spk_count/0.5, grp_lo_sum_new.sem_spk_count/0.5, ...
        'Color',choice_clr{2},'Marker','.','MarkerSize',10);
    hel.CapSize = 0;
    hel.LineWidth = 1.3;
    hel.LineStyle = 'none';

    xlim([min(beh_tb_predict.delta_ev)-10, max(beh_tb_predict.delta_ev)+10]);

    annotation('textbox',[col_posi(1),row_posi(2)+0.2, 0.13,0.02],'String', ...
        {sprintf('subjid: %d', spkinfo.subjid(spkinfo.cellid == cellid)),sprintf('cellid: %d', cellid)}, ...
        'LineStyle','none', 'FontSize',6)

    if i ==1
        xlabel('\Delta EV (\mul of water)','FontSize',8);
        ylabel('Hz','FontSize',8);
    end
    if i ~=1
        set(gca, 'XLabel',[]);
    end

    if i==3
        lla1=legend({'fr_{surebet}','','fr_{lottery}'},FontSize=6);
        set(lla1, 'Position', [0.434,0.74,0.07,0.02]);
        lla1.ItemTokenSize = [15,18];
        lla1.Box = 'off';
    end
    %set(gca, 'xticklabel',[]);

    %title({'Behavior performance & spke tuning', sprintf('for deltaEV %s\n',choose_win_str)});
    %legend({'fr for surebet','','fr for lottery'},Location="best")



    if i == 1
        %draw a fake legend for the psth plot
        fake_l = zeros(8,1);

        fake_l(1) = plot(NaN,NaN,'LineStyle','-.','Color',choice_clr{1},'LineWidth',1.2);
        fake_l(2) = plot(NaN,NaN,'LineStyle','-','Color',choice_clr{2},'LineWidth',1.2);
        for fi = 3:8
            fake_l(fi) = plot(NaN,NaN,'LineStyle','-','Color',lottery_clr{fi-2},'LineWidth',1.2);
        end
        fake_legend = legend(fake_l,{'surebet','lottery','0.5','2','4','8','16','32'},FontSize=6);
        set(fake_legend,'Position',[0.445,0.815,0.05,0.04]);
        fake_legend.ItemTokenSize = [15,18];
        fake_legend.Box = 'off';
    end
    hold off
end

%% do the scatter plot here
choice_glm_tb = choice_glm_tb(:,{'cellid','beta','p','start','R_squared'});
lottery_glm_tb = lottery_glm_tb(:,{'cellid','beta','p','start','R_squared'});

choice_glm = choice_glm_tb(choice_glm_tb.start == 0.5, {'cellid','beta', 'p','R_squared'});
lottery_glm = lottery_glm_tb(lottery_glm_tb.start == 0.5, {'cellid','beta', 'p','R_squared'});

sum_tb = outerjoin(choice_glm,lottery_glm,"Keys",{'cellid','cellid'});
sum_tb = sum_tb((~isnan(sum_tb.beta_choice_glm)) & (~isnan(sum_tb.beta_lottery_glm)),:);
grp_indx = zeros(height(sum_tb),1);
grp_indx(sum_tb.p_choice_glm<0.05 & (~sum_tb.p_lottery_glm<0.05)) = 1;
grp_indx((~sum_tb.p_choice_glm<0.05) & sum_tb.p_lottery_glm<0.05) = 2;
grp_indx(sum_tb.p_choice_glm<0.05 & sum_tb.p_lottery_glm<0.05) = 3;

sum_tb.grp_indx = grp_indx;
choice_clr = utils.linspecer(8,'blue');
choice_clr = choice_clr(4,:);
lottery_clr = utils.linspecer(8,'red');
lottery_clr = lottery_clr(4,:);
non_clr = utils.linspecer(8,'gray');
non_clr1 = non_clr(2,:);
both_clr = [0.494 0.184 0.556];

sum_tb.beta_choice_glm = sum_tb.beta_choice_glm;
sum_tb.beta_lottery_glm = sum_tb.beta_lottery_glm;

s2 = gscatter(draw.jaxes([0.58, 0.76, 0.21, 0.15]), ...
    sum_tb.beta_choice_glm, ...
    sum_tb.beta_lottery_glm, ...
    sum_tb.grp_indx,[non_clr1;scater_clr], ...
    '.',9);

xlabel('\beta _{choice} (t-stat)','FontSize',8);
ylabel('\beta _{\Delta EV} (t-stat)','FontSize',8);
set(gca,'FontSize');
set(gca,'XTick',[-15,-10,-5,0,5,10,15]);
set(gca,'ytick',[-15,-10,-5,0,5,10,15]);
xlim([-15,15]);
ylim([-15,15]);

hold on
sum_tb_hl = sum_tb(ismember(sum_tb.cellid_choice_glm,cell_ids),:);
plot(sum_tb_hl.beta_choice_glm,sum_tb_hl.beta_lottery_glm,'o','MarkerSize',5,'LineWidth',1.5,'Color',[44,127,184]/255/1.3);
plot([0,0],ylim,'LineWidth',1,'LineStyle',':','Color',non_clr(6,:));
plot(xlim,[0,0],'LineWidth',1,'LineStyle',':','Color',non_clr(6,:));
ls = legend({'non selective','choice only','lottery only','both','examples'},'FontSize',6,'Location',[0.68,0.76,0.1,0.05]);
ls.Box = 'off';
hold off

% tweak the scatter plot
current_ax = gca;
ch = current_ax.Children;
for cx = 4:7
    clr = ch(cx).Color;
    set(ch(cx), 'Marker','o', 'MarkerSize',2.5)
    ch(cx).MarkerFaceColor = clr;
    ch(cx).MarkerEdgeColor = clr/3;
end

% egh = ch(3);
% egh.Color = [0.3 0.3 0.3];

%% make the decode plots

col_posi_decode = [0.1,0.215,0.35]+0.4;
row_posi_decode = 0.58;
fig_height_decode = 0.08;
fig_width_decode = 0.09;

fit_tb = table();
good = isfinite(models(3).lottery_prediction);
fit_tb.mag = meta.lottery_mag(good);
fit_tb.pred = models(3).lottery_prediction(good);

pred_stats = grpstats(fit_tb,'mag','mean','DataVars','pred');
fit_tb.norm_mag = fit_tb.mag/max(fit_tb.mag);
fit_tb.norm_pred = models(3).lottery_prediction(good)/(pred_stats.mean_pred(pred_stats.mag==32));

f = @(p) mseutility(p, fit_tb.norm_mag,fit_tb.norm_pred);
[x,fval, ex, ou] = fmincon(f, [1,1], [],[]);

axes('Position',[col_posi_decode(1),row_posi_decode,fig_width_decode,fig_height_decode*0.9])
violins = violinplot(fit_tb.norm_pred,fit_tb.norm_mag,'width',0.02, 'ShowBox', false, ...
    'ShowWhiskers',false, 'LineWidth',1.2, 'MarkerSize',6, 'MedianMarkerSize', 12);
%plot(draw.jaxes([col_posi_decode(1),row_posi_decode,fig_width_decode,fig_height_decode*0.9]),fit_tb.norm_mag, fit_tb.norm_pred,'.','color',[1,1,1]*0.7,'MarkerSize',10);
hold on
%plot(pred_stats.mag/32,pred_stats.mean_pred/pred_stats.mean_pred(pred_stats.mag==32),'b.','MarkerSize',20)
xax  = 0:0.001:1;
plot(xax, x(2).*xax.^x(1),'k-','LineWidth',1.3);
xlim([-0.06,1.06])
%ylim([0,1]);
xlabel('Lottery magnitude');
ylabel('Estimated lottery magnitude')
subtitle(sprintf('subjid: %d, sessid: %d', decode_tb0.subjid(decode_tb0.sessid == sessid_decode_deomo), decode_tb0.sessid(decode_tb0.sessid == sessid_decode_deomo)))
annotation('textbox',[col_posi_decode(1)+0.0061,row_posi_decode-0.002,0.13,0.02],'String', ...
    {sprintf('r = %.2f, p < 0.01', decode_tb0.nl_r(decode_tb0.sessid == sessid_decode_deomo)), ...
    sprintf('%d trials, %d cells', size(data(1).features,1), size(data(1).features,2))}, ...
    'LineStyle','none', 'FontSize',6)
box off
ax = gca;
ax.TitleHorizontalAlignment = 'left';
%plot(xax, m.a*exp(m.b*xax)+1,'r-.');
%pause
hold off

%% plot the fitted curve for lottery value vs. predected value

decode_tb0 = decode_tb0(decode_tb0.nl_rho>0.01,:);
decode_tb = sortrows(decode_tb0,"nl_r");
xax  = 0:0.001:1;
color_scale = decode_tb.nl_r.^2;
color_scale = color_scale/max(color_scale);
plot(draw.jaxes([col_posi_decode(2),row_posi_decode,fig_width_decode,fig_height_decode*0.9]),xax, decode_tb.nl_scale(1).*xax.^decode_tb.nl_rho(1),'LineWidth',.2,'Color',[1,1,1]-color_scale(1))
hold on
for i = 2:height(decode_tb)
    plot(xax, decode_tb.nl_scale(i).*xax.^decode_tb.nl_rho(i),'LineWidth',.2,'Color',[1,1,1]-color_scale(i));
end
xlim([-0.06,1.06])
xlabel('Lottery magnitude')
%ylabel('Estimated lottery magnitude')
set(gca,'FontSize',8);
%subtitle('non-linear fit for 56 sessions')
annotation('textbox',[col_posi_decode(2)+0.027,row_posi_decode,0.1,0.02],'String', ...
    sprintf('%d sessions', height(decode_tb0)), ...
    'LineStyle','none', 'FontSize',6)
hold off
% hist graph
histsig(decode_tb.lott_cor,decode_tb.lott_p<0.05,'origin_col',col_posi_decode(3),'origin_row',row_posi_decode,'width',fig_width_decode,'height',fig_height_decode*0.8);
xlabel({"Single-session lottery decoding","(Pearson's r)"});
ylabel('# of sessions');
set(gca,'FontSize',8);
hold off

%% plot the pseudo population decoding
col_posi_decode = [0.1,0.215,0.35];
row_posi_pseudo_decode = 0.58;

pseudo_norm_lott = repmat([0.5, 2, 4, 8, 16, 32]/32, 20, 1);
pseudo_norm_lott = pseudo_norm_lott(:);
pseudo_hat_lott = reshape(fof_pseudo.x128, 120, []);
pseudo_nl_r = nan(size(pseudo_hat_lott,2),1);
pseudo_nl_p = pseudo_nl_r;
pseudo_nl_rho = pseudo_nl_r;
pseudo_nl_scale = pseudo_nl_r;

for jj = 1:length(pseudo_nl_p)
    this_pseudo_hat_lott = pseudo_hat_lott(:,jj)/mean(pseudo_hat_lott(101:120,jj));
    f = @(p) mseutility(p, pseudo_norm_lott,this_pseudo_hat_lott);
    [x,fval, ex, ou] = fmincon(f, [1,1], [],[]);
    [r,p]=corr(this_pseudo_hat_lott, x(2)*pseudo_norm_lott.^x(1));
    pseudo_nl_p(jj,1) = p;
    pseudo_nl_r(jj,1) = r;
    pseudo_nl_rho(jj,1) = x(1);
    pseudo_nl_scale(jj,1) = x(2);
end
pseudo_nltb = array2table([1:size(pseudo_hat_lott,2)]',"VariableNames",{'cellid'});
pseudo_nltb.nl_p = pseudo_nl_p;
pseudo_nltb.nl_r = pseudo_nl_r;
pseudo_nltb.nl_rho = pseudo_nl_rho;
pseudo_nltb.nl_scale = pseudo_nl_scale;
pseudo_nltb = [pseudo_nltb, pseudo_correlation(pseudo_correlation.n == 128,2)];

axes('Position',[col_posi_decode(1),row_posi_pseudo_decode,fig_width_decode,fig_height_decode*0.9])
violins = violinplot(pseudo_hat_lott(:,1)/mean(pseudo_hat_lott(101:120,1)),pseudo_norm_lott, ...
    'width',0.02, 'ShowBox', false, 'ShowWhiskers',false, 'LineWidth',1.2, 'MarkerSize',6, 'MedianMarkerSize', 12);
hold on
xax  = 0:0.001:1;
plot(xax, pseudo_nltb.nl_scale(1).*xax.^pseudo_nltb.nl_rho(1),'k-','LineWidth',1.3);
xlim([-0.06,1.06])
ylim([-0.12,1.55]);
xlabel('Lottery magnitude');
ylabel('Estimated lottery magnitude')
subtitle('Example pseudosession')
annotation('textbox',[col_posi_decode(1)-0.0016,row_posi_pseudo_decode+0.0527,0.13,0.02],'String', ...
    {sprintf('r = %.2f, p < 0.01', pseudo_nltb.nl_r(1,1)),'120 trials, 128 cells'}, ...
    'LineStyle','none', 'FontSize',6)
box off
ax = gca;
ax.TitleHorizontalAlignment = 'left';
%plot(xax, m.a*exp(m.b*xax)+1,'r-.');
%pause
hold off

% plot the 50 pseudo sessions
pseudo_nltb_sort = sortrows(pseudo_nltb, "nl_r");
color_scale = pseudo_nltb_sort.nl_r.^2;

color_scale = color_scale/max(color_scale)*0.8;
plot(draw.jaxes([col_posi_decode(2),row_posi_pseudo_decode,fig_width_decode,fig_height_decode*0.9]), ...
    xax, pseudo_nltb_sort.nl_scale(1).*xax.^pseudo_nltb_sort.nl_rho(1),'LineWidth',0.2,'Color',[1,1,1]-color_scale(1))
hold on
for i = 2:height(pseudo_nltb_sort)
    plot(xax, pseudo_nltb_sort.nl_scale(i).*xax.^pseudo_nltb_sort.nl_rho(i),'LineWidth',0.2,'Color',[1,1,1]-color_scale(i));
end
xlim([-0.06,1.06])
ylim([-0.12,1.55]);
xlabel('Lottery magnitude')
%ylabel('Estimated lottery magnitude')
subtitle('50 pseudosessions')
annotation('textbox',[col_posi_decode(2),row_posi_pseudo_decode+0.0527,0.13,0.02],'String', ...
    ('120 trials, 128 cells'), ...
    'LineStyle','none', 'FontSize',6)
hold off

% plot the boxchart
axes('Position',[col_posi_decode(3),row_posi_pseudo_decode,fig_width_decode,fig_height_decode*0.9])
boxchart(log2(pseudo_correlation.n), pseudo_correlation.Pearson, 'Notch','on','MarkerSize',3.5)
%set(gca, 'XScale', 'log')
xlim([4.2,11.8])
xticks([5,7,9,11])
xticklabels({'2^{5}','2^{7}','2^{9}','2^{11}'})
ylim([0.09, 1.01])
yticks([0.1, 0.5, 1.0])
xlabel('Population size')
ylabel({'Decoding accuracy',"(Pearson's r)"})
set(gca, 'TickDir', 'out')



set(findall(gcf,'-property','FontSize'),'FontSize',6);
set(findall(gcf,'-property','FontName'),'FontName','Arial');

%F = getframe(FigA);

%% Supplement for ephys
pseudo_correlation_shuffle = readtable(fullfile(utils.git_root_path, 'csv/fof_shuffle_pseudo_correlation.csv'));
FigB = figure('visible','on','Position', [0 0 1260 1782]);
% plot the boxchart
axes('Position',[0.2,0.85,0.09,0.08*0.9])
boxchart(log2(pseudo_correlation_shuffle.n), pseudo_correlation_shuffle.Pearson, 'Notch','on','MarkerSize',3.5)
%set(gca, 'XScale', 'log')
xlim([4.2,11.8])
xticks([5,7,9,11])
xticklabels({'2^{5}','2^{7}','2^{9}','2^{11}'})
ylim([-0.98, 1.01])
yticks([-0.5,0, 0.5, 1.0])
xlabel('Population size')
ylabel({'Decoding accuracy',"(Pearson's r)"})
set(gca, 'TickDir', 'out')
set(findall(gcf,'-property','FontSize'),'FontSize',6);
set(findall(gcf,'-property','FontName'),'FontName','Arial');

%% functions
% spikes counter in time window for spkinfor table
function spk_ma = get_spk_win(spkinfo, events, timewindow)
    [~, n_type_events] = size(events);
    spk_ma = zeros(length(events), height(spkinfo), n_type_events);
    for i = 1:height(spkinfo)
        spk_ts = spkinfo(i,:).ts{1};
        for j = 1:n_type_events
            spk_ma(:, i, j) = stats.spike_count(events(:,j), spk_ts, 'pre',timewindow(1), 'post', timewindow(2));
        end
    end
end

function [choice, lott_mag, events, subjid, PD] = getBehData(SD, varargin)
iod = @utils.inputordefault;
trial_type = iod('trial_type', 'free', varargin);
state = iod('state', 'cue_target', varargin);
long_rein_cutoff = iod('long_rein_cutoff', 0.3, varargin);
long_softfix_cutoff = iod('long_softfix_cutoff', 1.5, varargin);

PD = SD.protocol_data; % behavioral data
subjid = SD.subjid;
good = PD.viol ==0;

if trial_type == "free"
    good = good & (cell2mat(cellfun(@(x) ~isempty(x), PD.lottery_poke, 'UniformOutput',false)) &...
        cell2mat(cellfun(@(x) ~isempty(x), PD.surebet_poke, 'UniformOutput',false)));
elseif trial_type == "forced"
    good = good & ~((cell2mat(cellfun(@(x) ~isempty(x), PD.lottery_poke, 'UniformOutput',false)) &...
        cell2mat(cellfun(@(x) ~isempty(x), PD.surebet_poke, 'UniformOutput',false))));
end

% there are some trials need to be exclude:
% in soft-fixation condition, we may want to exclude trials with long soft fixation duration (>1.2s)
% in hard-fixation condition, we may want to exclude trials with a long-rein state before complete the hard fixation
hard_fix_trials = arrayfun(@(x) any(contains(get_entered_states(x), 'play_lottery_sound2_hard_fixation')), SD.peh, 'UniformOutput',1);
entered_choice_poke = arrayfun(@(x) any(contains(get_entered_states(x), 'deliver_surebet_stopsound')) | any(contains(get_entered_states(x), 'deliver_lottery_stopsound')), SD.peh, 'UniformOutput',1);

if any(hard_fix_trials) % hard fixation condition
    long_rein_status = arrayfun(@(x) extract_long_rein(x, long_rein_cutoff), SD.peh, 'UniformOutput',1);
    good = (good & entered_choice_poke & ~long_rein_status);

    event_cue = getStateChange(SD.peh, 'choice_wait_for_','prefix',1)-1;
    event_cue = event_cue(good);
    event_choice_poke = getStateChange(SD.peh, 'deliver_surebet_stopsound', 'prefix', 1);
    event_lottery_poke = getStateChange(SD.peh, 'deliver_lottery_stopsound', 'prefix', 1);
    event_lottery_poke_indx = arrayfun(@(x) any(contains(get_entered_states(x), 'deliver_lottery_stopsound')), SD.peh, 'UniformOutput',1);
    event_choice_poke(event_lottery_poke_indx) = event_lottery_poke(~isnan(event_lottery_poke));

    event_target = event_choice_poke(good);

else %soft fixation condition
    long_soft_status = arrayfun(@(x) extract_long_softfix(x,long_softfix_cutoff), SD.peh, 'UniformOutput',1);
    good = (good & entered_choice_poke & ~long_soft_status);

    event_cue = getStateChange(SD.peh, 'play_lottery_sound2','prefix',0);
    event_cue = event_cue(good);
    event_choice_poke = getStateChange(SD.peh, 'deliver_surebet_stopsound', 'prefix', 1);
    event_lottery_poke = getStateChange(SD.peh, 'deliver_lottery_stopsound', 'prefix', 1);
    event_lottery_poke_indx = arrayfun(@(x) any(contains(get_entered_states(x), 'deliver_lottery_stopsound')), SD.peh, 'UniformOutput',1);
    event_choice_poke(event_lottery_poke_indx) = event_lottery_poke(~isnan(event_lottery_poke));

    event_target = event_choice_poke(good);
end

choice = PD.subj_choice(good);
lott_mag = PD.lottery_mag(good);

switch state
    case 'cue'
        events = event_cue;
    case 'target'
        events = event_target;
    case 'cue_target'
        events = [event_cue, event_target];
end
PD = PD(good,:);
end


% return all entered states in one trial
function entered_states = get_entered_states(peh)
all_states = fields(peh.States);
entered_states_index = structfun(@(x) ~isnan(x(1)), peh.States);
entered_states = all_states(entered_states_index);
end

% extract trials with a long rein state before complete the hard fixation
function long_rein = extract_long_rein(peh, threshold)

if nargin < 2
    threshold = 0.3;
end

entered_states = get_entered_states(peh);

if ~any(contains(entered_states, 'play_lottery_sound2_hard_fixation'))
    long_rein = 0;
elseif ~any(contains(entered_states, 'fixation_rein'))
    long_rein = 0;
else
    rein_dur = peh.States.fixation_rein(:, 2) - peh.States.fixation_rein(:, 1);
    %rein_dur(rein_dur>0.99)=[];
    rein_dur = [rein_dur; peh.States.play_lottery_sound2_hard_fixation(2) - peh.States.play_lottery_sound2_hard_fixation(1)];
    rein_dur_sort = sort(rein_dur);
    long_rein = 1;
    if rein_dur_sort(end - 1) < threshold
        long_rein = 0;
    end
end
end

% extract trials with a long soft fixation period before arriving to target
function long_softfix = extract_long_softfix(peh, threshold)
if nargin < 2
    threshold = 1.3;
end

entered_states = get_entered_states(peh);
if ~any(contains(entered_states, 'play_lottery_sound2_hard_fixation'))
    long_softfix = 0;
    if any(contains(entered_states, 'choice_wait_for_choice_BotR_BotL'))
        soft_fix_dur = peh.States.choice_wait_for_choice_BotR_BotL(1, 1) - peh.States.play_lottery_sound2(1, 1);
        if soft_fix_dur > threshold
            long_softfix = 1;
        else
            long_softfix = 0;
        end
    end
end
end