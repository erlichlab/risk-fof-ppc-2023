function [hx]=histsig(x, x_sig, varargin)
% [ hx]=histsig(x, x_sig,  x_lim)


inpd = @utils.inputordefault;
[origin_col, varargin]=inpd('origin_col', 0.1, varargin);
[origin_row, varargin]=inpd('origin_row', 0.1, varargin);
[wdth, varargin]=inpd('width',0.2,varargin);
[hist_h, varargin]=inpd('height',0.180,varargin);
[num_bins, varargin]=inpd('n_bins',17,varargin);
[bins, varargin]=inpd('bins',[],varargin);
[ax, varargin]=inpd('ax',[],varargin);
[x_lim, varargin]=inpd('x_lim',[],varargin);
[y_lim, varargin]=inpd('y_lim',[],varargin);
[normed, varargin]=inpd('normed',false,varargin);
[zero, varargin]=inpd('zero',0,varargin);
inpd(varargin)

gd=~isnan(x);
x=x(gd);
x_sig=x_sig(gd);
x_lim_b=[min(x)-0.1 max(x)+0.1];

if isempty(x_lim)
    x_lim=x_lim_b;
end
    


if isempty(ax)
    ax=draw.jaxes([origin_col origin_row wdth hist_h]);
end

marker_size=zeros(size(x))+12;
marker_size(x_sig==1)=24;

% make the x-axis histogram
if isempty(bins)
    bins=linspace(x_lim_b(1),x_lim_b(2),num_bins);
end
nsig=histcounts(x(x_sig==0), bins);
sig=histcounts(x(x_sig==1), bins);
if normed
    total = sum(nsig) + sum(sig);
    nsig = nsig / total;
    sig = sig / total;
end
maxy=max(nsig(:)+sig(:))*1.3;
cbins=edge2cen(bins);
[hh]=bar(ax,cbins, [sig(:) nsig(:)],'stacked');
xlim(ax,x_lim);
set(ax, 'YAxisLocation','left');
set(hh(1),'FaceColor','k')
set(hh(2),'FaceColor',[1 1 1])
if isempty(y_lim)
    y_lim = [0 maxy];
end

set(ax,'box','off','YLim',y_lim)
set(ax,'Color','none')
text(ax, getx(ax),gety(ax),sprintf('%d%% p<.05, n=%d',round(100*mean(x_sig)),sum(~isnan(x))))
%text(ax, getx(ax),gety(ax),[num2str(round(100*mean(x_sig))) '% p<0.05'])

x_mean=nanmean(x);
[xt_sig,~,B]=stats.bootmean(x-zero);

[CI]=prctile(B+zero,[2.5 97.5]);

if xt_sig<0.05
    y_lim=ylim(ax);
    y_pos=0.85 * y_lim(2);
    plot(ax,x_mean, y_pos,'.k','MarkerSize',6);
    plot(ax,[CI(1) CI(2)], [y_pos y_pos], '-k');
    
end
 

function y=edge2cen(x)
b2b_dist=x(2)-x(1);
y=x+0.5*b2b_dist;
y=y(1:end-1);

function y=getx(ax)
x=xlim(ax);
y=0.03*(x(2)-x(1))+x(1);

function y=gety(ax)
x=ylim(ax);
y=1.1*(x(2)-x(1))+x(1);



