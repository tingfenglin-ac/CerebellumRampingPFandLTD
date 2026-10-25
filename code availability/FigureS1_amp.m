% 1. Add the path of the current script to the search path
% 2. Change the current folder to "Figure S1",
% 3. Using ctrl+enter within each function to run each section

%% Section 1: Open files
clear

%% merge optic data from 2 mice
% open "1x50_optogenetics_2.mat" and "electro_CF_amp.mat" in data
% availability/Figure S1
signals=mean_m2_1_50(2:618, :); % ensure stim onset is aligned



%% plot 
co=lines;
co(6,:)=[.5 .5 .5];
cnum=[5 2];
amp{1}=ymat(1,~isnan(ymat(1,:)));


x   = (1:size(signals,1))-309;
x   = x./30.98;
window_idx=x>0 & x<0.5;
amp{2}=max(signals(window_idx,:));

y=[];
ysem=[];
matconncat=[];
concatgroup=[];
for con=1:2
    y(con)=mean(mean(amp{con},1));
    ysem(con)=nanstd(mean(amp{con},1),0,2)./sqrt(size(mean(amp{con},1),2));
    matconncat=[matconncat mean(amp{con},1)];
    concatgroup=[concatgroup con.*ones(size(mean(amp{con},1)))];
end


figure;
hold on;

MarkerSize=70;
meanSize=10;
space=0.19;
jrange=0.25;
arrayfun(@(a) errorbar(a+space,y(a),ysem(a),'o-',...
    'color',co(cnum(a),:),...
    'CapSize',30,...
    'MarkerSize',meanSize,...
    'LineWidth',3),1:2,'uni',0);
arrayfun(@(a) text(a+space,y(a)+ysem(a)+0.1,num2str(y(a),'%4.2f'),...
    'color',co(cnum(a),:),...
    'HorizontalAlignment','center',...
    'FontSize',15),1:2,'uni',0);
compgroup=[1 2 3 4];
arrayfun(@(a) scatter(a-space-0.5*jrange+jrange*rand(1,length(matconncat(concatgroup==a))),matconncat(concatgroup==a),MarkerSize,...
    'MarkerFaceColor','flat',...
    'MarkerEdgeColor','flat',...
    'MarkerFaceAlpha',0.2,...
    'MarkerEdgeAlpha',0.2,...
    'CData',co(cnum(a),:)),1:2,'uni',0);


% Two-sample t-test
data1 = amp{1}(~isnan(amp{1}));
data2 = amp{2}(~isnan(amp{2}));
[h,p.ttest2,ci,stats.ttest2] = ttest2(data1,data2);

Ysig = 1.3;
sigincre = 0.01;
% significance symbol
if p.ttest2 < 0.001
    sigtext = '***';
    plot([1 2],[Ysig Ysig],'color',[0 0 0],'linewidth',3);
    text(1.5,Ysig+sigincre,sigtext,...
        'FontSize',30,...
        'HorizontalAlignment','center',...
        'color',[0 0 0]);
elseif p.ttest2 < 0.01
    sigtext = '**';
    plot([1 2],[Ysig Ysig],'color',[0 0 0],'linewidth',3);
    text(1.5,Ysig+sigincre,sigtext,...
        'FontSize',30,...
        'HorizontalAlignment','center',...
        'color',[0 0 0]);
elseif p.ttest2 < 0.05
    sigtext = '*';
    plot([1 2],[Ysig Ysig],'color',[0 0 0],'linewidth',3);
    text(1.5,Ysig+sigincre,sigtext,...
        'FontSize',30,...
        'HorizontalAlignment','center',...
        'color',[0 0 0]);
end


xticks(1:4)
yticks(0:0.5:2)
xticklabels([]);
xticklabels({'Electro','Opto'});
ylim([0 1.4])
xlim([0.5 2.5])
ylabel('Amplitude (\DeltaF/F)')
set(gca,'FontSize',23)
set(gcf,'color',[1 1 1])