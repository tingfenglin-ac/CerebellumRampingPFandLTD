% Figure 3a calcium signal trace
% 1. Add the path of the current script to the search path
% 2. Change the current folder to "Figure 3 PF and CF synaptic signals",
% 3. Using ctrl+enter within each function to run each section

%% Section 1: Open files 
% Select all the PF data or CF data

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
CEllnumber=arrayfun(@(a) size(data(a).signal,2),find(condlist==1));%cell number of each experiment
Cellid=arrayfun(@(a) arrayfun(@(b) [FNlist{a} '_cell' num2str(b)],1:CEllnumber(a),'uni',0),find(condlist==1),'UniformOutput',0)
Cellid=horzcat(Cellid{:});

%% Section 2: plot the traces

increvalue=1;
base=21;
LOCS=pksLOCS;

n=15;
co = return_colorbrewer('Oranges', n);
co = cat(3,co,return_colorbrewer('Blues', n));
co = cat(3,co,return_colorbrewer('BuGn', n));
co = cat(3,co,return_colorbrewer('Greys', n));
co = cat(3,co,return_colorbrewer('RdPu', n));
co = co(7:13,:,:);
co=repmat(co,1,1,7);

cnum=1;
amps=cell(size(condlist));

y=[];
counts=cell(size(condlist));
for cond=1:max(condlist)
    condidx=find(condlist==cond)';
    figure
    incre=0;
    [~,NsizeRank]=sort(cellfun(@(a) size(a,2),smoothBC_signal));
    for i=NsizeRank
        idx=condidx(i);
        T=time{idx}-time{idx}(data(idx).stim(1));
        incre=increvalue.*(1:size(LOCS{idx},2))+incre(end);
        for roi=1:size(LOCS{idx},2)
            activ=smoothBC_signal{idx}(:,roi,:);
            mean_activ=nanmean(activ,3);

            hold on
            plot(T,mean_activ+incre(roi),...
                'color',co(roi,:,cnum),...
                'LineWidth',1.5);

            signalStd=nanstd(activ,0,3)./sqrt(size(activ,3));
            ErrArea_Smooth(T,mean_activ+incre(roi),...
                signalStd,...
                co(roi,:,cnum));
        end
        cnum=cnum+1;
    end

    xlim([-5.5 5.5])
    xticks([-5:5])
    xlabel('Time (sec)')
    set(gca,'TickDir','out');
    ylim([0 40])
    set(gca,'FontSize',18)
    set(gcf,'color',[1 1 1])
    figposition=get(gcf,'position');
    figposition(4)=1200;
    figposition(2)=-500;
    set(gcf,'Position',figposition);
end