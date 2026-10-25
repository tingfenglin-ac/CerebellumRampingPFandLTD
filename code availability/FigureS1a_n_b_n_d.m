%% For Supplementary Figure S1a, S1b, and S1d
% ctrl+enter to run each session
%% Section 1 
% ctrl+enter and open all the "YYYYMMDD_pre_XXX.mat" data in "data
% availability/FigureS1"

addpath(genpath('G:\My Drive\UChicago\Hansel lab\databackup\scripts'));



clear
FN=[];
    FP=[];
FNlist={};

FPlist={};
condlist=[];


f=1;
while f
    [FileName,FolderPath] = uigetfile({'*SpikeAna.mat;*SpikeAnalist.mat'},'Select SpikeData files', 'Multiselect', 'on');
    
    if FolderPath==0;f=0;end
    
    if iscell(FileName)
        NewAddFile=size(FileName,2);
    elseif FileName~=0
        NewAddFile=1;
        if strfind(FileName,'SpikeAnalist')
            load([FolderPath,FileName])
            NewAddFile=size(FN,1);
        end
    else
        NewAddFile=0;
    end
    %%
    if NewAddFile~=0;
%         if ~size(strfind(FileName,'Spikelist'),1)|iscell(FileName)
            for fnumber=1:NewAddFile
                if iscell(FileName)
                    FNlist=cat(1,FNlist,FileName{fnumber});
                    FPlist=cat(1,FPlist,FolderPath);
                elseif strfind(FileName,'SpikeAnalist')
                    FNlist=cat(1,FNlist,FN{fnumber});
                    FPlist=cat(1,FPlist,FP{fnumber});
                else
                    FNlist=cat(1,FNlist,FileName);
                    FPlist=cat(1,FPlist,FolderPath);
                end
                condlist=cat(1,condlist,f);
            end
                            f=f+1;

%         end
    end
end


% import data
for i=1:size(FNlist,1)
data(i)=load([FPlist{i},FNlist{i}]);
end

%
% pull the data out
riseLOCS=arrayfun(@(x) data(x).riseLOCS,1:length(data),'uni',0);
upLOCS=arrayfun(@(x) data(x).upLOCS,1:length(data),'uni',0);
pksLOCS=arrayfun(@(x) data(x).pksLOCS,1:length(data),'uni',0);
PKS=arrayfun(@(x) data(x).PKS,1:length(data),'uni',0);
time=arrayfun(@(x) data(x).time,1:length(data),'uni',0);
smoothBC_signal=arrayfun(@(x) data(x).smoothBC_signal,1:length(data),'uni',0);

% categorize the accumulated signal (during bust) from individual spike
% signal
ind_rise_LOCS=cell(size(riseLOCS));
accm_rise_LOCS=cell(size(riseLOCS));
ind_pks_LOCS=cell(size(riseLOCS));
accm_pks_LOCS=cell(size(riseLOCS));
ind_PKS=cell(size(riseLOCS));
accm_PKS=cell(size(riseLOCS));
for i=1:length(data)
ind_rise_LOCS{i}=cellfun(@(x,y) x(~ismember(x,y)),riseLOCS{i},upLOCS{i},'uni',0);
accm_rise_LOCS{i}=cellfun(@(x,y) x(ismember(x,y)),riseLOCS{i},upLOCS{i},'uni',0);
ind_pks_LOCS{i}=cellfun(@(x,y,z) z(~ismember(x,y)),riseLOCS{i},upLOCS{i},pksLOCS{i},'uni',0);
accm_pks_LOCS{i}=cellfun(@(x,y,z) z(ismember(x,y)),riseLOCS{i},upLOCS{i},pksLOCS{i},'uni',0);
ind_PKS{i}=cellfun(@(x,y,z) z(~ismember(x,y)),riseLOCS{i},upLOCS{i},PKS{i},'uni',0);
accm_PKS{i}=cellfun(@(x,y,z) z(ismember(x,y)),riseLOCS{i},upLOCS{i},PKS{i},'uni',0);
end
ROIlabel=arrayfun(@(a) ['Cell ' num2str(a)],1:size(smoothBC_signal{1},2),'uni',0);

%% Section2
% open "1x50_optogenetics_2.mat" data in "data availability/FigureS1"
% ctrl+enter

signals=mean_m2_1_50(2:618, :); % ensure stim onset is aligned


%% Section 3 (Fig. S1a) normalized ave cells of responsive trace 
function Normalized_ave_responsive_overlap_trace(Marker)
%% 
figure;hold on

LOCS=pksLOCS;
ana_window=[0 0.2];
increvalue=0.1;

[co] = return_colorbrewer('Dark2', 8);
co=distinguishable_colors(50);
co=lines;
cnum=[1 2 4 5];
amps=cell(size(condlist));

%
y=[];
counts=cell(size(condlist));
for cond=1:max(condlist)
    condidx=find(condlist==cond)';
    ind_activ_cell=arrayfun(@(a) nanmean(smoothBC_signal{a},3),condidx,'UniformOutput',0);
    ind_std_cell=arrayfun(@(a) nanstd(smoothBC_signal{a},0,3),condidx,'UniformOutput',0);

    ind_activ{cond}=horzcat(ind_activ_cell{:});
    ind_std{cond}=horzcat(ind_std_cell{:});
end
T=time{1}-time{1}(153);
catmean=[];
for cond=1:max(condlist)
    condidx=find(condlist==cond)';
    mean_activ=[];
    amps=[];
    amps_idx=[];
    cellmean=[];
    for i=1:length(condidx)
        idx=condidx(i);
        incre=increvalue.*(1:size(LOCS{idx},2));
        for roi=1:size(LOCS{idx},2)
            % aveF
            % because the histcounts include left left edge, so
            % convert the sign, and then left-right flip.
            % because the histcounts include the right edge of the
            % last bin, so add one more bin at right
            counts{idx}{roi,1} = cellfun(@(a) logical(histcounts(fliplr(-T(a)),fliplr(-[ana_window(1) ana_window]))),LOCS{idx}(:,roi),'uni',0);
            counts{idx}{roi,1} = cellfun(@(a) a(1),counts{idx}{roi,1});

            activ=smoothBC_signal{idx}(:,roi,:);
            mean_activ{idx}(:,roi,:)=nanmean(activ(:,:,counts{idx}{roi,1}),3);


            ana_idx=LOCS{idx}(:,roi); %index of onset or peak depend on analysis

            amps{idx}{roi,1} = arrayfun(@(a) activ(ana_idx{a}(T(ana_idx{a})>ana_window(1) & T(ana_idx{a})<=ana_window(2)),:,a),...
                1:length(ana_idx),'uni',0);%filter out the data out of the ana_window

            amps_idx{idx}{roi,1} = arrayfun(@(a) find(amps{idx}{roi,1}{a}),...
                1:length(ana_idx),'uni',0);%find all the index of amplitude
            amps_idx{idx}{roi,1} = arrayfun(@(a) min(amps_idx{idx}{roi,1}{a}),...
                1:length(ana_idx),'uni',0);%find the first index
            amps{idx}{roi,1} = arrayfun(@(a) amps{idx}{roi,1}{a}(amps_idx{idx}{roi,1}{a}),...
                1:length(ana_idx),'uni',0);%first amplitude within the ana_window

            cellmean{idx}{roi,1} = nanmean(vertcat(amps{idx}{roi}{:}),1);
        end
    end
    cellmean=vertcat(cellmean{condidx});
    catmean{cond}=horzcat(cellmean{:});
    ind_activ{cond}=horzcat(mean_activ{:});
end

for cond=1:max(condlist)
    ave_activ{cond}=nanmean(ind_activ{cond}./catmean{1},2);
    std_activ{cond}=nanstd(ind_activ{cond}./catmean{1},0,2)./sqrt(size(ind_activ{cond},2));
end

% figure;hold on
co=lines;
co(6,:)=[.5 .5 .5];
cnum=[1 2 4 5];
cnum=[5 5 5];% when ploting WT tetanus group
delay=1.6; %delay of each spike
display_window=[-0.2 1.3]; % display of the time window for each spike
window_idx=T>display_window(1) & T<display_window(2);
zero_idx=153;
linestyle={'-' '-' '-'};


for cond=1:max(condlist)

    plot_trace=(ave_activ{cond}(window_idx)-ave_activ{cond}(zero_idx))./max(ave_activ{cond}(window_idx)-ave_activ{cond}(zero_idx));
    %plotting error
    incre=increvalue.*(1:size(ind_activ{cond},2));
    ErrArea_Smooth(T(window_idx)+(cond-1).*delay,plot_trace,...
        std_activ{cond}(window_idx)./max(ave_activ{cond}(window_idx)-ave_activ{cond}(zero_idx)),...
        [co(cnum(cond),:) .5]);

    %     subplot(1,5,cond)
    incre=increvalue.*(1:size(ind_activ{cond},2));

    hold on
    plot(T(window_idx)+(cond-1).*delay,plot_trace,linestyle{cond},...
        'color',[co(cnum(cond),:)],...
        'LineWidth',1.5);

    %%
    box off
    % b.BaseLine.Visible = 'off'

    xlabel('Time (sec)')
    % ylabel('Amplitude (dF/F0)')
    xticks(-10:10)
    set(gca,'FontSize',18,'TickDir','out')
%     set(gcf,'color','none')
%     exportgraphics(gcf, 'figure.png', ...
%     'BackgroundColor', 'none', ...
%     'Resolution', 300);
end
plot(repmat(0:0.003:0.009,2,1),repmat([-0.1 ;-0.2],1,4),'color',[co(cnum(cond),:)])



% Abby's plot 
% 1x50
m   = nanmean(signals, 2);
s   = nanstd(signals,  0, 2) ./ sqrt(size(signals, 2));
m   = m./max(m);
s   = s./max(m);

x   = (1:length(m))'-309;
x   = x./30.98;
fill([x; flipud(x)], [m - s; flipud(m + s)], co(2,:), ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot(x, m, 'Color', co(2,:), 'LineWidth', 1.5);

% plot(repmat(0:0.125:0.25,2,1),repmat([-0.1 ;-0.2]-0.15,1,3),'color',co(1,:))

stim_x=[0;0.05];
fill([stim_x; flipud(stim_x)], [-0.1;-0.1; -0.2;-0.2]-0.15, co(2,:), ...
    'EdgeColor', 'none');



 xlim([-0.25 1.35])
    ylim([-0.55 1.05])
end

function FWHM
%% Section 2 (Fig. S1d) FWHM with peak constrained to first 500 ms

% electro
window_idx = T>=0 & T<display_window(2);
data = ind_activ{1}(window_idx,2:end);
t = T(window_idx);

amp{1} = nan(1,size(data,2));

for n = 1:size(data,2)

    trace = data(:,n);

    % Peak must occur within first 500 ms
    peak_window = t >= 0 & t <= 0.5;

    [peakAmp,peakIdx_temp] = max(trace(peak_window));

    peak_indices = find(peak_window);
    peakIdx = peak_indices(peakIdx_temp);

    % Half maximum
    halfAmp = peakAmp/2;

    % Rising half-max crossing before peak
    idxRise = find(trace(1:peakIdx) >= halfAmp,1,'first');

    % Falling half-max crossing after peak
    idxFallRel = find(trace(peakIdx:end) <= halfAmp,1,'first');

    if ~isempty(idxRise) && idxRise > 1 && ...
       ~isempty(idxFallRel) && idxFallRel > 1

        idxFall = peakIdx + idxFallRel - 1;

        % Rising crossing
        tRise = interp1( ...
            trace(idxRise-1:idxRise), ...
            t(idxRise-1:idxRise), ...
            halfAmp);

        % Falling crossing
        tFall = interp1( ...
            trace(idxFall-1:idxFall), ...
            t(idxFall-1:idxFall), ...
            halfAmp);

        % FWHM
        amp{1}(n) = tFall - tRise;
    end

end


% opto 
x = (1:size(signals,1))-309;
x = x./30.98;

window_idx = x>=0 & x<display_window(2);

data = signals(window_idx,:);
t = x(window_idx);

amp{2} = nan(1,size(data,2));

for n = 1:size(data,2)

    trace = data(:,n);

    % Peak within first 500 ms
    peak_window = t >= 0 & t <= 0.5;

    [peakAmp,peakIdx_temp] = max(trace(peak_window));

    peak_indices = find(peak_window);
    peakIdx = peak_indices(peakIdx_temp);

    halfAmp = peakAmp/2;

    idxRise = find(trace(1:peakIdx) >= halfAmp,1,'first');
    idxFallRel = find(trace(peakIdx:end) <= halfAmp,1,'first');

    if ~isempty(idxRise) && idxRise > 1 && ...
       ~isempty(idxFallRel) && idxFallRel > 1

        idxFall = peakIdx + idxFallRel - 1;

        tRise = interp1( ...
            trace(idxRise-1:idxRise), ...
            t(idxRise-1:idxRise), ...
            halfAmp);

        tFall = interp1( ...
            trace(idxFall-1:idxFall), ...
            t(idxFall-1:idxFall), ...
            halfAmp);

        amp{2}(n) = tFall - tRise;
    end

end
%
co=lines;
co(6,:)=[.5 .5 .5];
cnum=[5 2];


y=[];
ysem=[];
matconncat=[];
concatgroup=[];
for con=1:2
    y(con)=nanmean(nanmean(amp{con},1));
    ysem(con)=nanstd(nanmean(amp{con},1),0,2)./sqrt(size(~isnan(nanmean(amp{con},1)),2));
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
arrayfun(@(a) text(a+space,y(a)+ysem(a)+0.07,num2str(y(a),'%4.2f'),...
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

Ysig = 0.95;
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
xticklabels({'Electro','Opto'});
yticks(0:0.4:1)
ylim([0 1.1])
xlim([0.5 2.5])
ylabel('FWHM (sec)')
set(gca,'FontSize',23)
set(gcf,'color',[1 1 1])
end

function risetime_20_80
%% Section 3 (Fig. S1b) 20-80% rise time

% electro
window_idx = T >= 0 & T < 0.5;
data = ind_activ{1}(window_idx,2:end);
t = T(window_idx);

amp{1} = nan(1,size(data,2));

for n = 1:size(data,2)

    trace = data(:,n);
    trace=trace-trace(1);

    % Find peak
    [peakAmp,peakIdx] = max(trace);

    % 10% and 90% peak levels
    amp10 = 0.2 * peakAmp;
    amp90 = 0.8 * peakAmp;

    % Search only during rising phase
    riseTrace = trace(1:peakIdx);
    riseT = t(1:peakIdx);

    % Find first crossings
    idx10 = find(riseTrace >= amp10,1,'first');
    idx90 = find(riseTrace >= amp90,1,'first');

    if ~isempty(idx10) && ~isempty(idx90) && idx10 > 1 && idx90 > 1

        % Interpolate crossing times
        t10 = interp1(riseTrace(idx10-1:idx10), ...
                      riseT(idx10-1:idx10), ...
                      amp10);

        t90 = interp1(riseTrace(idx90-1:idx90), ...
                      riseT(idx90-1:idx90), ...
                      amp90);

        % 10-90% rise time
        amp{1}(n) = t90 - t10;
    end
end





% opto 
x = (1:size(signals,1))-309;
x = x./30.98;

window_idx = x >= 0 & x < 0.5;
data = signals(window_idx,:);
t = x(window_idx);

amp{2} = nan(1,size(data,2));

for n = 1:size(data,2)

    trace = data(:,n);
    trace=trace-trace(1);

    [peakAmp,peakIdx] = max(trace);

    amp10 = 0.2 * peakAmp;
    amp90 = 0.8 * peakAmp;

    riseTrace = trace(1:peakIdx);
    riseT = t(1:peakIdx);

    idx10 = find(riseTrace >= amp10,1,'first');
    idx90 = find(riseTrace >= amp90,1,'first');

    if ~isempty(idx10) && ~isempty(idx90) && idx10 > 1 && idx90 > 1

        t10 = interp1(riseTrace(idx10-1:idx10), ...
                      riseT(idx10-1:idx10), ...
                      amp10);

        t90 = interp1(riseTrace(idx90-1:idx90), ...
                      riseT(idx90-1:idx90), ...
                      amp90);

        amp{2}(n) = t90 - t10;
    end
end


co=lines;
co(6,:)=[.5 .5 .5];
cnum=[5 2];


y=[];
ysem=[];
matconncat=[];
concatgroup=[];
for con=1:2
    y(con)=nanmean(nanmean(amp{con},1));
    ysem(con)=nanstd(nanmean(amp{con},1),0,2)./sqrt(size(nanmean(amp{con},1),2));
    matconncat=[matconncat nanmean(amp{con},1)];
    concatgroup=[concatgroup con.*ones(size(nanmean(amp{con},1)))];
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
arrayfun(@(a) text(a+space,y(a)+ysem(a)+0.015,num2str(y(a),'%4.2f'),...
    'color',co(cnum(a),:),...
    'HorizontalAlignment','center',...
    'FontSize',15),1:2,'uni',0);

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

Ysig = 0.35;
sigincre = 0.005;
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
xticklabels({'Electro','Opto'});
yticks(0:0.1:0.4)
ylim([0 0.35])
xlim([0.5 2.5])
ylabel('20%-80% rise time (sec)')
set(gca,'FontSize',23)
set(gcf,'color',[1 1 1])
end