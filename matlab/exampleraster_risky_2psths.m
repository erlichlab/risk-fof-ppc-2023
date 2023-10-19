function [ras,R]=exampleraster_risky_2psths(ev, ts,varargin)
% [ax_handle,data]=rasterC(ev, ts, varargin)
% pairs={'pre'        3;...
%        'post'       3;...
%        'binsz'      0.050;...
%        'cnd'        1;...
%        'meanflg'    0;...
%        'krn'        0.25;...
%        'ax_handle'  [];...
%        'legend_str' '';...
%        'renderer', 'opengl';...
%        'ref_label', 'REF';...
%        'psth_height', 0.248;...
%        'total_height' 0.8;...
%        'corner'       [0.1 0.1];...
%        'ax_width'      0.55;...
% 	   'font_name'	   'Helvetica';...
% 	   'font_size'		9;...
% 	   'legend_pos'     [0.73 0.1 0.2 0.15];...
% 	   'clrs'	{'c','b','r','m','r','g','m'};...
% 	   'x_label','';...
%     'pre_mask', -inf;...
%     'post_mask',+inf;...
%     'cout',[];...
%     'stim_back',[];...
%     'errorbars', 0;...
%     'testfunc', [];...
%     'sortby', [];...
%

corner=[];  % this is necessary because corner is a function
iod = @utils.inputordefault;
pre = iod('pre',3,varargin);
post = iod('post',3,varargin);
binsz = iod('binsz',0.01,varargin);
cnd = iod('cnd',1,varargin);
cnd2 = iod('cnd2',[],varargin); 
meanflg = iod('meanflg',0,varargin);
krn = iod('krn',0.1,varargin);
ax_handle = iod('ax_handle',[],varargin);
legend_str = iod('legend_str','',varargin);
legend_str2 = iod('legend_str2','',varargin);
renderer =iod('renderer','painters',varargin);
ref_label=iod('ref_label','REF',varargin);
psth_height=iod('psth_height',0.248,varargin);
total_height=iod('total_height',0.8,varargin);
corner=iod('corner',[0.1 0.1],varargin);
ax_width=iod('ax_width',0.55,varargin);
font_name=iod('font_name','Helvetica',varargin);
font_size=iod('font_size',14,varargin);
legend_pos=iod('legend_pos',[0.73 0.1 0.2 0.15],varargin);
legend_pos2=iod('legend_pos2',[0.73 0.3 0.2 0.15],varargin);
clrs=iod('clrs',{'b','m','r','c','k','g','y',[1 0.5 0],[0.5 0.5 0.5]},varargin);
clrs2=iod('clrs2',{'b','m','r','c','k','g','y',[1 0.5 0],[0.5 0.5 0.5]},varargin);
x_label=iod('x_label','',varargin);
pre_mask=iod('pre_mask', -inf,varargin);
post_mask=iod('post_mask',+inf,varargin);
cout=iod('cout',[],varargin);
stim_back=iod('stim_back',[],varargin);
sb_clr=iod('sb_clr',[0.8 0.8 0.4],varargin);
errorbars=iod('errorbars',1,varargin);
testfunc=iod('testfunc',[],varargin);
show_yinfo=iod('show_yinfo',1,varargin);
sortby=iod('sortby',[],varargin);
psth_ylim=iod('psth_ylim',[],varargin);
TickDir=iod('TickDir',[],varargin);
line_alpha=iod('line_alpha',0.8,varargin);
legend_fontsz=iod('legend_fontsz',6,varargin);

line_style = {'-.','-',':','--'};


set(gcf, 'Renderer',renderer);
[ntrials,nrefs]=size(ev);


mutau=zeros(1,nrefs);
for rx=2:nrefs
    mutau(rx)=nanmedian(ev(:,rx)-ev(:,1));
end

if isscalar(pre_mask)
    pre_mask=zeros(1,ntrials)+pre_mask;
elseif numel(pre_mask)~=ntrials
    fprintf(1,'numel(pre_mask) must equal num ref events or be scalar');
    return;
end


if isscalar(post_mask)
    post_mask=zeros(1,ntrials)+post_mask;
elseif numel(post_mask)~=ntrials
    fprintf(1,'numel(post_mask) must equal num ref events or be scalar');
    return;
end

if isscalar(krn)
    dx=ceil(5*krn);
    kx=-dx:binsz:dx;
    krn=normpdf(kx,0, krn);
	if isempty(find(kx==0, 1))
		error('Your binsz needs to divide 1 second into interger # of bins');
	end
    krn(kx<0)=0;
    krn=krn/sum(krn);
end

if numel(cnd)==1
    cnd=ones(1,ntrials);
end

if isempty(cnd2)
    if iscell(cnd)
        cnd_nan = cellfun(@(x)any(isnan(x)), cnd); % use any to deal with character arrays.
        cnd(cnd_nan) = {'NaN'};
        cnd = categorical(cnd);
        n_cnd = categories(cnd);
    else
        cnd = categorical(cnd);
        n_cnd = categories(cnd);
    end
else
    cnd1_el = sort(unique(cnd));
    cnd2_el = sort(unique(cnd2));
    cnd_ma = [cnd,cnd2];
    rank_ma = zeros(length(cnd2_el),1);
    for ii = 1:length(cnd2_el)
        rank_ma(cnd2 == cnd2_el(ii)) = ii;
    end
    cnd_ma = [cnd,cnd2,rank_ma];
    for i = 2:numel(cnd1_el)
        for j = 1:numel(cnd2_el)
            cnd_ma(cnd_ma(:,1) == cnd1_el(i) & cnd_ma(:,2) == cnd2_el(j),2) = cnd2_el(j) + max(cnd2_el);
        end
    end
    cnd = categorical(cnd_ma(:,2));
    n_cnd = categories(cnd);
    cnd_indx_ma = unique(cnd_ma,'rows');
end

raster_height=total_height-psth_height;
y_ind=psth_height+corner(2)+0.005;

height_per_trial=raster_height/ntrials;
if isempty(cnd2)
    psthax=axes('Position',[corner(1) corner(2) ax_width psth_height]);
    hold(psthax,'on');
    set(psthax,'FontName',font_name);
    set(psthax,'FontSize',font_size);
else
    psth_height1 = psth_height*0.49;
    psthax1=axes('Position',[corner(1) corner(2) ax_width psth_height1]);
    hold(psthax1,'on');
    set(psthax1,'FontName',font_name);
    set(psthax1,'FontSize',font_size);

    psthax2=axes('Position',[corner(1) corner(2)+psth_height1*1.1 ax_width psth_height1]);
    hold(psthax2,'on');
    set(psthax2,'XTickLabel',[]);
    set(psthax2,'FontName',font_name);
    set(psthax2,'FontSize',font_size);
end



%[Y,x,W]=warpfilter(ev,ts,krn,'pre',pre,'post',post,'kernel_bin_size',binsz);
[Y,x]=stats.spike_filter(ev,ts,krn,'pre',pre,'post',post,'kernel_bin_size',binsz);
W = ts;

for ci=1:numel(n_cnd)
    sampz=sum(cnd==n_cnd(ci));
    
    ref=cnd==n_cnd(ci);
    if ~isempty(sortby)
        idx=1:ntrials; idx=idx(:); ref=ref(:);
        ref=sortrows([sortby(ref) idx(ref)]);
        ref=ref(:,2);
    end
    
    y=Y(ref,:);
    
    [y2,x2]=draw.rasterplot(ev(ref,1),W,pre,post,'pre_mask',pre_mask(ref),'post_mask',post_mask(ref),'plotthis',0);
    ras(ci)=axes('Position',[corner(1) y_ind ax_width height_per_trial*sampz]);
    y_ind=y_ind+height_per_trial*sampz+0.001;
    
    if ~isempty(stim_back)
        patchplot(stim_back(ref,:),'clr',sb_clr)
    end
    
    
    %% Plot the rasters
    ll=line(x2,y2);
    set(ll,'color','k');
    set(gca,'XTickLabel',[]);
    set(gca,'YTick',[]);
    set(gca,'Box','off')
    set(gca,'YLim',[0 max(y2)])
    set(gca,'XLim',[-pre post]);
    
    
    for rx=1:nrefs
        if isempty(cnd2)
            ll=line([mutau(rx) mutau(rx)],[0 max(y2)]);
            set(ll,'LineStyle','-','color',clrs{ci},'LineWidth',1);
        else
            ll1=line([max(x) max(x)],[0 max(y2)]);
            set(ll1,'LineStyle',line_style{find(cnd1_el == cnd_indx_ma(ci,1))},'color',clrs{find(cnd1_el == cnd_indx_ma(ci,1))},'LineWidth',1.3);
            ll2=line([mutau(rx) mutau(rx)],[0 max(y2)]);
            set(ll2,'LineStyle',line_style{find(cnd1_el == cnd_indx_ma(ci,1))},'color',clrs2{cnd_indx_ma(ci,3)},'LineWidth',1.3);
        end
    end
    
    if ~isempty(cout)
        hold on;
        if isempty(cnd2)
            h=plot(ras(ci),cout(ref),1:sampz,'o','Color','k','MarkerFaceColor',clrs{ci},'MarkerSize',2);
        else
            h=plot(ras(ci),cout(ref),1:sampz,'o','Color','k','MarkerFaceColor',clrs2{cnd_indx_ma(ci,3)},'MarkerSize',2);
        end
    end
    %% Calculate the mean and ci of the
    
    [y x]=draw.maskraster(x,y,pre_mask(ref),post_mask(ref));
    
    ymn(ci,:) = nanmean(y,1);
    yst(ci,:)= stats.nanstderr(y,1);
    R{ci}={y x};
    %axes(psthax);
    hold on
    %     hh=line(x/1000,ymn(ci,:));
    % 	set(hh,'LineWidth',1,'LineStyle','-','Color',clrs{ci});
    if strcmpi(renderer,'opengl')
        sh(ci)=draw.shadeplot(x,ymn(ci,:)-yst(ci,:),ymn(ci,:)+yst(ci,:),{clrs{ci},psthax,0.3});
        % lh=line(x,ymn(ci,:),'Color',clrs{ci},'LineWidth',2);
    else
        if errorbars
            hh(1)=line(x,ymn(ci,:)-yst(ci,:),'Parent',psthax);
            hh(2)=line(x,ymn(ci,:)+yst(ci,:),'Parent',psthax);
            set(hh,'LineWidth',1,'LineStyle','-','Color',clrs{ci});
            %lh=line(x,ymn(ci,:),'Color',clrs{ci},'LineWidth',1,'Parent',psthax);
            lh=hh(1);
            sh(ci)=lh;
        else
            if isempty(cnd2)
                hh(1)=line(x,ymn(ci,:),'Parent',psthax);
                set(hh,'LineWidth',2,'LineStyle','-','Color',clrs{ci});
                %lh=line(x,ymn(ci,:),'Color',clrs{ci},'LineWidth',1,'Parent',psthax);
                lh=hh(1);
                sh(ci)=lh;
            else
                hh(1)=line(x,ymn(ci,:),'Parent',psthax2);
                set(hh,'LineWidth',1.2,'LineStyle',line_style{find(cnd1_el == cnd_indx_ma(ci,1))}, ...
                    'Color',[clrs2{cnd_indx_ma(ci,3)},line_alpha]); 
                %lh=line(x,ymn(ci,:),'Color',clrs{ci},'LineWidth',1,'Parent',psthax);
                lh=hh(1);
                sh2(ci)=lh;
            end
        end
            
    end
    peaky(ci)=max(ymn(ci,:));
    if isempty(cnd2)
        set(psthax,'XLim',[-pre,post]);
    else
        set(psthax2,'XLim',[-pre,post]);
    end

    legstr{ci}=[n_cnd{ci} ', n=' num2str(sampz)];
    
end

for ci2 = 1:numel(cnd1_el)
    sampz=sum(cnd_ma(:,1)==cnd1_el(ci2));
    ref=cnd_ma(:,1)==cnd1_el(ci2);
    y=Y(ref,:);
    [y x]=draw.maskraster(x,y,pre_mask(ref),post_mask(ref));

    ymn1(ci2,:) = nanmean(y,1);
    yst1(ci2,:)= stats.nanstderr(y,1);
    R1{ci2}={y x};

    hh(1)=line(x,ymn1(ci2,:),'Parent',psthax1);
    set(hh,'LineWidth',1.2,'LineStyle',line_style{ci2}, ...
        'Color',clrs{ci2});
    %lh=line(x,ymn(ci,:),'Color',clrs{ci},'LineWidth',1,'Parent',psthax);
    lh=hh(1);
    sh1(ci2)=lh;
    set(psthax1,'XLim',[-pre,post]);
end

if isempty(cnd2)
    cur_ylim=get(psthax,'YLim');
else
    cur_ylim = get(psthax2,'yLim');
end
cur_ylim(cur_ylim<0) = 0;% y-axis not going to negative, this is a problem at somewhere else and we will fix that later
if isempty(cnd2)
    if isempty(psth_ylim)
        ylim(psthax,[cur_ylim(1)*0.8 max(peaky)*1.15]);
        height_zeroline = max(peaky)*1.15;
    else
        ylim(psthax,psth_ylim);
        height_zeroline = psth_ylim(2);
    end
else
    if isempty(psth_ylim)
        ylim(psthax2,[cur_ylim(1)*0.8 max(peaky)*1.15]);
        ylim(psthax1,[cur_ylim(1)*0.8 max(peaky)*1.15]);
        height_zeroline = max(peaky)*1.15;
    else
        ylim(psthax2,psth_ylim);
        ylim(psthax1,psth_ylim);
        height_zeroline = psth_ylim(2);
    end
end

if isempty(cnd2)
    for rx=1:nrefs
        ll=plot(psthax,[mutau(rx) mutau(rx)],[0 height_zeroline]);
        set(ll,'LineStyle','-','color',[0.5 0.5 0.5],'LineWidth',1);
    end

    ch=get(psthax,'Children');
    set(psthax,'Children',[ch(nrefs+1:end); ch(1:nrefs)]);
    if TickDir == 'out'
        set(psthax,'TickDir','out'); %set the tick direction outside
    end

    xticks=get(psthax,'XTick');
    set(psthax,'XTick',xticks);
    set(ras,'XTick',xticks);
else
    for rx=1:nrefs
        ll=plot(psthax2,[mutau(rx) mutau(rx)],[0 height_zeroline]);
        ll1=plot(psthax1,[mutau(rx) mutau(rx)],[0 height_zeroline]);
        set(ll,'LineStyle','-','color',[0.5 0.5 0.5],'LineWidth',1);
        set(ll1,'LineStyle','-','color',[0.5 0.5 0.5],'LineWidth',1);
    end

    ch=get(psthax2,'Children');
    set(psthax2,'Children',[ch(nrefs+1:end); ch(1:nrefs)]);
    if TickDir == 'out'
        set(psthax2,'TickDir','out'); %set the tick direction outside
        set(psthax1,'TickDir','out'); %set the tick direction outside
    end

    xticks=get(psthax2,'XTick');
    set(psthax2,'XTick',xticks);
    set(ras,'XTick',xticks);
end


if isempty(cnd2)
    if ~isempty(legend_pos) && ~isempty(legend_str)
        [lh,oh]=legend(sh,legend_str);
        %legend boxoff %this code will generate an unexpected legend
        % keyboard
        set(lh,'Position',legend_pos);
    end
    ras(end+1) = lh;

    hold off
    %set(gca,'FontSize',36);

    if isempty(x_label)
        xh=xlabel(psthax,[]); set(xh,'interpreter','none');
    else
        xh=xlabel(psthax,x_label);set(xh,'interpreter','none');
    end
    if show_yinfo
        ylabel(psthax,'Hz \pm SE')
    else
        set(psthax,'YTickLabel',[]);
    end
else
    if ~isempty(legend_pos) && ~isempty(legend_str)
        [lh,oh]=legend(sh1,legend_str);
        %legend boxoff %this code will generate an unexpected legend
        % keyboard
        set(lh,'Position',legend_pos);
        set(lh,'FontSize',legend_fontsz);
        if ~isempty(legend_pos2) && ~isempty(legend_str2)
            [lh1,~] = legend(sh2,legend_str2);
            set(lh1,'Position',legend_pos2);
            set(lh1,'FontSize',legend_fontsz);
        end

    end
    ras(end+1) = lh;

    hold off

    if isempty(x_label)
        xh=xlabel(psthax1,[]); set(xh,'interpreter','none');
    else
        xh=xlabel(psthax1,x_label);set(xh,'interpreter','none');
    end
    if show_yinfo
        ylabel(psthax1,'Hz')
        ylabel(psthax2,'Hz')
    else
        set(psthax1,'YTickLabel',[]);
        set(psthax2,'YTickLabel',[]);
    end
end
if isempty(cnd2)
    ras(end+1)=psthax;
else
    ras(end+1) = psthax2;
    ras(end+1) = psthax1;
end
