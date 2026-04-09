% 1. Add the path of the current script to the search path
% 2. Change the current folder to "Figure 4 5 S1 S2 calcium raw traces"
% 3. analyze the data within each subfolder, which is "single pulse", "short burst", "tonic" or "ramping"
% 4. Using ctrl+enter within each function to run each section

%% Section 1: Open files
% Select all the pre data (YYYYMMDD_pre_SpikeAna.mat) and click "open",
% And then elect all the early post PF+100Hz CF data (YYYYMMDD_postCF100hz1min_early_SpikeAna.mat) and click "open",
% And then elect all the late post PF+100Hz CF data (YYYYMMDD_postCF100hz1min_late_SpikeAna.mat) and click "open",
% And then elect all the early post PF+100Hz CF data (YYYYMMDD_postCF333hz1min_early_SpikeAna.mat) and click "open",
% And then elect all the late post PF+100Hz CF data (YYYYMMDD_postCF333hz1min_late_SpikeAna.mat) and click "open",
% ctrl+enter to run this section
% if single pulse proto=1; if short burst proto=2;if tonic pulse proto=3;if single pulse proto=4;
clear
proto=3;

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
            NewAddFile=size(FileName,2);
        end
    else
        NewAddFile=0;
    end  
    if NewAddFile~=0;
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
    end
end

% import data
for i=1:size(FNlist,1)
    data(i)=load([FPlist{i},FNlist{i}]);
end

% pull the data out
pksLOCS=arrayfun(@(x) data(x).pksLOCS,1:length(data),'uni',0);
time=arrayfun(@(x) data(x).time,1:length(data),'uni',0);
smoothBC_signal=arrayfun(@(x) data(x).smoothBC_signal,1:length(data),'uni',0);

% categorize the accumulated signal (during bust) from individual spike
% signal
ROIlabel=arrayfun(@(a) ['Cell ' num2str(a)],1:size(smoothBC_signal{1},2),'uni',0);

%% Section 3: plot normalized average calcium traces of each stimulus conditions
% From left to right, plot from pre- to late-tetanus conditions
% ctrl+enter to run this section

increvalue=0.1;
amps=cell(size(condlist));
ana_window=[0.05 0.3];
T=time{1}-time{1}(data(1).stim(1));
ana_idx=T>ana_window(1) & T<=ana_window(2);

y=[];
counts=cell(size(condlist));
for cond=1:max(condlist)
    condidx=find(condlist==cond)';
    ind_activ_cell=arrayfun(@(a) nanmean(smoothBC_signal{a},3),condidx,'UniformOutput',0);
    ind_std_cell=arrayfun(@(a) nanstd(smoothBC_signal{a},0,3),condidx,'UniformOutput',0);
    ind_activ{cond}=horzcat(ind_activ_cell{:});
    ind_std{cond}=horzcat(ind_std_cell{:});
end
Nind_activ=cellfun(@(a) a./max(ind_activ{1}(ana_idx,:)),ind_activ,'uni',0);

figure; hold on
for cond=1:max(condlist)
    ave_activ{cond}=nanmean(Nind_activ{cond},2);
    std_activ{cond}=nanstd(Nind_activ{cond},0,2)./sqrt(size(Nind_activ{cond},2));
end

co=lines;
co(6,:)=[.5 .5 .5];
colortransfer=[6 5 2 4];
cnum=repmat(colortransfer(proto),1,5);
delay=1.6; 
display_window=[-0.2 1.3]; 
window_idx=T>display_window(1) & T<display_window(2);
linestyle={'-' '-' '-' '-' '-'};

stimT=[T(data(1).stim)';T(data(1).stim)'];
stimT=stimT(:);
Stim_signal=[0 2 2 0]';

anaT=[0.05 0.3;0.05 0.3];
anaT=anaT(:);
ana_signal=[-1 2 2 -1]';

for cond=1:max(condlist)
    %plotting stimulus
    b=area(stimT+(cond-1).*delay,Stim_signal,...
        'FaceColor',[255 96 0]/255,...
        'EdgeColor','none',...
        'FaceAlpha',0.4);

    %plotting analysis window
    b=area(anaT+(cond-1).*delay,Stim_signal,...
        'FaceColor',[249, 191, 69]./255,...
        'EdgeColor','none',...
        'FaceAlpha',0.4);

    %plotting error
    incre=increvalue.*(1:size(ind_activ{cond},2));
    ErrArea_Smooth(T(window_idx)+(cond-1).*delay,ave_activ{cond}(window_idx),...
        std_activ{cond}(window_idx),...
        [co(cnum(cond),:) .5]);
end

for cond=1:max(condlist)
    incre=increvalue.*(1:size(ind_activ{cond},2));

    hold on
    plot(T(window_idx)+(cond-1).*delay,ave_activ{cond}(window_idx),linestyle{cond},...
        'color',[co(cnum(cond),:)],...
        'LineWidth',1.5);
end

    box off
    b.BaseLine.Visible = 'off'
    xlim([-0.3 7.5])
    ylim([-0.02 1.5])
    xlabel('Time (sec)')
    xticks(-10:10)
    set(gca,'FontSize',18,'TickDir','out')
    set(gcf,'color',[1 1 1])
