% analyse the control data for lottery sound in naive animals
% Chaofei Bao, 5/2/2023
if ~exist(fullfile(utils.git_root_path,'matlab/lottery_sound_control.mat'),'file')
    subjid_list = [2234, 2253, 2256, 2257];
    sessid_list = getSubjSess(subjid_list);
    time_win = [0.5 1];
    fitted_param_sum = [];
    fitted_param_sum_shuffle = [];
    n_suffle = 1000;


    for i = 1:numel(sessid_list)
        sessid = sessid_list.sessid(i);
        SD = db.getSessData(sessid);
        lottery_cue = get_lottery_cue(SD.peh);
        %lottery_cue = lottery_cue(randsample(1:length(lottery_cue),length(lottery_cue)),:);
        lottery_cue_t = getStateChange(SD.peh, 'PlayLotterySound_mag_','prefix',1);

        CD = getCellData(sessid);
        spk_ma = get_spk_win(CD, lottery_cue_t, time_win);
        spk_ma_zscore = zscore(spk_ma);
        %spk_ma_zscore = spk_ma;

        fitted_para_ma = zeros(CD.cell_num, 5);
        fitted_para_ma_shuffle = zeros(CD.cell_num, 5, n_suffle);

        for j = 1:CD.cell_num
            spk_tb = table(lottery_cue, spk_ma_zscore(:,j), 'VariableNames',{'lott_mag','spkc'});
            try
                m = fitlm(spk_tb,'spkc ~ lott_mag');
                %m = fitglm(spk_tb,'spkc ~ lott_mag','Distribution','poisson','Link','identity');
                fitted_para_ma(j,1:4) = get_param_lm(m, 'alpha', 0.05);
            catch
                fitted_para_ma(j,1:4) = nan(1,4);
                fprintf('cellid %d in window %d can not be fitted!!! \n', CD.info.cellid(j), k);
            end
            fitted_para_ma(j,5) = CD.info.cellid(j);

            for ii = 1:n_suffle
                lottery_cue_shuffle = lottery_cue(randsample(1:length(lottery_cue),length(lottery_cue)),:);

                spk_tb_shuffle = table(lottery_cue_shuffle, spk_ma_zscore(:,j), 'VariableNames',{'lott_mag','spkc'});
                try
                    m = fitlm(spk_tb_shuffle,'spkc ~ lott_mag');
                    %m = fitglm(spk_tb,'spkc ~ lott_mag','Distribution','poisson','Link','identity');
                    fitted_para_ma_shuffle(j,1:4,ii) = get_param_lm(m, 'alpha', 0.05);
                catch
                    fitted_para_ma_shuffle(j,1:4,ii) = nan(1,4);
                    fprintf('cellid %d in window %d can not be fitted!!! \n', CD.info.cellid(j), k);
                end
                fitted_para_ma_shuffle(j,5,ii) = CD.info.cellid(j);

            end

        end
        fitted_param_sum = [fitted_param_sum;fitted_para_ma];
        fitted_param_sum_shuffle = [fitted_param_sum_shuffle;fitted_para_ma_shuffle];
    end
else
    load('lottery_sound_control.mat')
end
%% plots
position = [0.1, 0.8, 0.3, 0.16];
FigA = figure('visible','on','Position', [0 0 1260 1782]);
histsig(fitted_param_sum(:,1),fitted_param_sum(:,3),'origin_col',position(1),'origin_row',position(2),'width',position(3), ...
    'height',position(4));
xlabel('\beta _{Lottey value} (tStat)','FontSize',8);
ylabel('# of cells');
set(gca,'FontSize',8);
hold off

% compute the shuffled significant num cells
sig_num_shuffle = sum(squeeze(fitted_param_sum_shuffle(:,3,:)),1);

% plot the boxchart
axes('Position',[position(1)+0.4, position(2), position(3)/2, position(4)])
boxplot([sig_num_shuffle',repelem(6,1,1000)'], 'Notch','on','Labels',{'Shuffled labels', 'Original labels'})
box off
yticks([0, 4, 8, 12])
ylabel('# of significant neurons')
set(gca, 'TickDir', 'out', 'FontSize',8)

set(findall(gcf,'-property','FontSize'),'FontSize',9);
set(findall(gcf,'-property','FontName'),'FontName','Arial');
%% functions
function para_ma = get_param_lm(m, varargin)
    iod = @utils.inputordefault;
    alpha = iod('alpha', 0.05, varargin);
    para_ma = zeros(1,4);
    para_ma(1,1) = m.Coefficients.tStat(2);
    para_ma(1,2) = m.Coefficients.pValue(2);
    para_ma(1,4) = m.Rsquared.Ordinary;
    if m.Coefficients.pValue(2) < alpha
        para_ma(1,3) = 1;
    end

end


function sessids = getSubjSess(subjid)
    dbc = db.labdb.getConnection('client');
    subjid_str = '';

    for i = 1:numel(subjid)
        subjid_str = sprintf('%s, %d', subjid_str, subjid(i));
    end

    subjid_str = subjid_str(2:end);
    sessids = dbc.query('select unique(sessid) from phy.cellview where subjid in (%s) and cellid > 26154', {subjid_str});
end

function entered_states = get_entered_states(peh)
    all_states = fields(peh.States);
    entered_states_index = structfun(@(x) ~isnan(x(1)), peh.States);
    entered_states = all_states(entered_states_index);
end

function lottery_cue = get_lottery_cue(peh)
    lottery_cue = nan(length(peh),1);
    lottery_cue(arrayfun(@(x) any(contains(get_entered_states(x), 'PlayLotterySound_mag_5_prob_7')), peh, 'UniformOutput',1),1)=0.5;
    lottery_cue(arrayfun(@(x) any(contains(get_entered_states(x), 'PlayLotterySound_mag_20_prob_7')), peh, 'UniformOutput',1),1)=2;
    lottery_cue(arrayfun(@(x) any(contains(get_entered_states(x), 'PlayLotterySound_mag_40_prob_7')), peh, 'UniformOutput',1),1)=4;
    lottery_cue(arrayfun(@(x) any(contains(get_entered_states(x), 'PlayLotterySound_mag_80_prob_7')), peh, 'UniformOutput',1),1)=8;
    lottery_cue(arrayfun(@(x) any(contains(get_entered_states(x), 'PlayLotterySound_mag_160_prob_7')), peh, 'UniformOutput',1),1)=16;
    lottery_cue(arrayfun(@(x) any(contains(get_entered_states(x), 'PlayLotterySound_mag_320_prob_7')), peh, 'UniformOutput',1),1)=32;
end

function CD = getCellData(sessid)
    dbc = db.labdb.getConnection('client');
    out = dbc.query('select * from phy.cellview where presence_ratio > 0.95 and firing_rate > 1 and snr_best_chan > 1.5 and sessid = %d',{sessid});
    
    CD.cell_num = numel(out.cellid);
    CD.info = out;
    for i = 1:CD.cell_num
        spkinfo = db.getCellData(out.cellid(i), 'full_info', 0);
        %CD.info(i, 3) = spkinfo.nSpikes;
        CD.ts{i} = spkinfo.ts{1};
    end
end

% spikes counter for selectivity index
function spk_ma = get_spk_win(CD, events, timewindow)
    [~, n_type_events] = size(events);
    spk_ma = zeros(length(events), CD.cell_num, n_type_events);
    for i = 1:CD.cell_num
        spk_ts = CD.ts{i};
        for j = 1:n_type_events
            spk_ma(:, i, j) = stats.spike_count(events(:,j), spk_ts, 'pre',timewindow(1), 'post', timewindow(2));
        end
    end
end