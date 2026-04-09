function Fig3statistics
% 1. Add the path of the current script to the search path
% 2. Change the current folder to "Figure 3 PF and CF synaptic signals",
% 3. Using ctrl+enter within each function to run each section

%% Section 1: Open files
% select all PF data all at once and enter, and then all the CF data

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
signal_raw=arrayfun(@(x) data(x).signal,1:length(data),'uni',0);

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
ana_window=[0.05 0.3];
%%
    function prob_state(Marker)
        %% Section 2: plotting the probability within "ana window" sec after stim
        LOCS=pksLOCS;
        counts=cell(size(condlist));
        y=[];
        for cond=1:max(condlist)
            condidx=find(condlist==cond)';
            for i=1:length(condidx)
                % average clls
                idx=condidx(i);
                T=time{idx}-time{idx}(data(idx).stim(1)+1);
                for roi=1:size(LOCS{idx},2)
                    counts{idx}{roi,1} = cellfun(@(a) logical(histcounts(fliplr(-T(a)),fliplr(-[ana_window(1) ana_window]))),LOCS{idx}(:,roi),'uni',0);
                    counts{idx}{roi,1} = cellfun(@(a) a(1),counts{idx}{roi,1},'uni',0);
                    prob{idx}{roi,1} = sum(vertcat(counts{idx}{roi}{:}),1)./length(LOCS{idx}(:,roi));
                end
            end
            catprob=vertcat(prob{condidx});
            y=cat(1,y,horzcat(catprob{:}));
        end

        figure
        x=1:max(condlist);
        hold on
        for n=1:size(y,2)
            jitterx=x'+0.*rand(size(x'));
            scatter(jitterx,y(:,n),50,...
                'MarkerFaceColor','flat',...
                'MarkerEdgeColor','flat',...
                'MarkerFaceAlpha',0.2,...
                'MarkerEdgeAlpha',0.2,...
                'MarkerFaceColor',[.5 .5 .5],...
                'MarkerEdgeColor',[.5 .5 .5]);
            h=plot(jitterx,y(:,n),...
                'color',[.5 .5 .5 .2],...
                'LineWidth',1);
            h.Color(4) = 1;
        end

        y_prob=y;
        hold on
        ymean=nanmean(y,2);
        ystd=nanstd(y,0,2)./sqrt(size(y,2));
        errorbar([0.5 2.5],ymean,ystd,'o',...
            'CapSize',15,...
            'MarkerSize',7,...
            'MarkerFaceColor',[.3 .3 .3],...
            'Color',[.3 .3 .3],...
            'LineWidth',3);
        text([0.5 2.5],ymean+0.075,compose('%.2f', ymean),'FontSize',16,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);

        % stat
        [h,p] = ttest(y(1,:)',y(2,:)');
        Ysig=0.21.*1000;Yincre=0.014.*1000;sigincre=0.0014.*1000;
        sigidx=find(p<0.05 & p>=0.01)';
        if sigidx
            for sn=sigidx
                plot(1:2,[Ysig Ysig],'color',[0 0 0],'linewidth',3)
                text(mean(1:2),Ysig+sigincre,'*','FontSize',30,...
                    'HorizontalAlignment','center',...
                    'color',[0 0 0]);
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(p<0.01 & p>=0.001)';
        if sigidx
            for sn=sigidx
                plot(1:2,[Ysig Ysig],'color',[0 0 0],'linewidth',3)
                text(mean(1:2),Ysig+sigincre,'**','FontSize',30,...
                    'HorizontalAlignment','center',...
                    'color',[0 0 0]);
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(p<0.001)';
        if sigidx
            for sn=sigidx
                plot(1:2,[Ysig Ysig],'color',[0 0 0],'linewidth',3)
                text(mean(1:2),Ysig+sigincre,'***','FontSize',30,...
                    'HorizontalAlignment','center',...
                    'color',[0 0 0]);
                Ysig=Ysig+Yincre;
            end
        end

        xticks([0.75 2.25])
        xticklabels({'PF' 'CF'});
        xlim([0 3])
        ylim([0 1.1])
        ylabel('Probability')
        set(gca,'FontSize',18)
        set(gcf,'color',[1 1 1])
        figposition=get(gcf,'position');
        figposition(3)=300;
        set(gcf,'Position',figposition);
    end

    function max_WithinTimewindow_state(Marker)
        %% Section 3: plotting the amplitude within "ana window" sec after stim
        LOCS=pksLOCS;
        amps=cell(size(condlist));
        amps_idx=cell(size(condlist));
        %
        y=[];
        for cond=1:max(condlist)
            condidx=find(condlist==cond)';
            for i=1:length(condidx)
                idx=condidx(i);
                T=time{idx}-time{idx}(data(idx).stim(1)+1);
                for roi=1:size(LOCS{idx},2)
                    activ=smoothBC_signal{idx}(:,roi,:);
                    amps{idx}{roi,1} = activ(T>ana_window(1) & T<=ana_window(2),:,:);%filter out the data out of the ana_window
                    amps{idx}{roi,1} = max(amps{idx}{roi,1},[],1);%find the max within ana_window
                    amps{idx}{roi,1} = amps{idx}{roi,1}(:);
                    cellmean{idx}{roi,1} = nanmean(amps{idx}{roi},1);
                end
            end
            catmean=vertcat(cellmean{condidx});
            y=cat(1,y,horzcat(catmean{:}));
        end

        figure
        x=1:max(condlist);
        yIDX=y(1,:)>0.3;

        subplot(1,10,3:10);hold on
        for n=1:size(y,2)
            jitterx=x'+0.*rand(size(x'))
            scatter(jitterx,y(:,n),50,...
                'MarkerFaceColor','flat',...
                'MarkerEdgeColor','flat',...
                'MarkerFaceAlpha',0.2,...
                'MarkerEdgeAlpha',0.2,...
                'MarkerFaceColor',[.5 .5 .5],...
                'MarkerEdgeColor',[.5 .5 .5]);
            h=plot(jitterx,y(:,n),...
                'color',[.5 .5 .5 .2],...
                'LineWidth',1);
            h.Color(4) = 1;
        end

        y_prob=y;
        %
        hold on
        ymean=nanmean(y,2);
        ystd=nanstd(y,0,2)./sqrt(size(y,2));
        errorbar([0.5 2.5],ymean,ystd,'o',...
            'CapSize',15,...
            'MarkerSize',7,...
            'MarkerFaceColor',[.3 .3 .3],...
            'Color',[.3 .3 .3],...
            'LineWidth',3);
        text([0.5 2.5],ymean+0.4,compose('%.2f', ymean),'FontSize',16,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);

        % stat
        [h,p] = ttest(y(1,:)',y(2,:)');
        Ysig=4.6;Yincre=0.1;sigincre=0.1;
        sigidx=find(p<0.05 & p>=0.01)';
        if sigidx
            for sn=sigidx
                plot(1:2,[Ysig Ysig],'color',[0 0 0],'linewidth',3)
                text(mean(1:2),Ysig+sigincre,'*','FontSize',30,...
                    'HorizontalAlignment','center',...
                    'color',[0 0 0]);
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(p<0.01 & p>=0.001)';
        if sigidx
            for sn=sigidx
                plot(1:2,[Ysig Ysig],'color',[0 0 0],'linewidth',3)
                text(mean(1:2),Ysig+sigincre,'**','FontSize',30,...
                    'HorizontalAlignment','center',...
                    'color',[0 0 0]);
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(p<0.001)';
        if sigidx
            for sn=sigidx
                plot(1:2,[Ysig Ysig],'color',[0 0 0],'linewidth',3)
                text(mean(1:2),Ysig+sigincre,'***','FontSize',30,...
                    'HorizontalAlignment','center',...
                    'color',[0 0 0]);
                Ysig=Ysig+Yincre;
            end
        end

        xticks([0.75 2.25])
        xticklabels({'PF' 'CF'});
        xlim([0 3])
        ylabel('Amplitude (\DeltaF/F)')
        set(gca,'FontSize',18)
        set(gcf,'color',[1 1 1])
        figposition=get(gcf,'position');
        figposition(3)=300;
        set(gcf,'Position',figposition);
    end

    function CVofAmp(Marker)
        %% Section 4: plotting the CV within "ana window" sec after stim
        LOCS=pksLOCS;
        amps=cell(size(condlist));
        y=[];
        for cond=1:max(condlist)
            condidx=find(condlist==cond)';
            for i=1:length(condidx)
                idx=condidx(i);
                T=time{idx}-time{idx}(data(idx).stim(1)+1);
                for roi=1:size(LOCS{idx},2)
                    activ=smoothBC_signal{idx}(:,roi,:);
                    pks_idx=pksLOCS{idx}(:,roi); %index of peak
                    ana_idx=LOCS{idx}(:,roi);
                    amps{idx}{roi,1} = activ(T>ana_window(1) & T<=ana_window(2),:,:);%filter out the data out of the ana_window
                    amps{idx}{roi,1} = max(amps{idx}{roi,1},[],1);%find the mean within ana_window
                    amps{idx}{roi,1} = amps{idx}{roi,1}(:);
                    cellmean{idx}{roi,1} = std(amps{idx}{roi},[],1)/nanmean(amps{idx}{roi},1);
                end
            end
            catmean=vertcat(cellmean{condidx});
            y=cat(1,y,horzcat(catmean{:}));
        end

        figure
        x=1:max(condlist);
        yIDX=y(1,:)>0.3;
        hold on
        for n=1:size(y,2)
            jitterx=x'+0.*rand(size(x'))
            scatter(jitterx,y(:,n),50,...
                'MarkerFaceColor','flat',...
                'MarkerEdgeColor','flat',...
                'MarkerFaceAlpha',0.2,...
                'MarkerEdgeAlpha',0.2,...
                'MarkerFaceColor',[.5 .5 .5],...
                'MarkerEdgeColor',[.5 .5 .5]);
            h=plot(jitterx,y(:,n),...
                'color',[.5 .5 .5 .2],...
                'LineWidth',1);
            h.Color(4) = 1;
        end

        y_prob=y;
        hold on
        ymean=nanmean(y,2);
        ystd=nanstd(y,0,2)./sqrt(size(y,2));
        errorbar([0.5 2.5],ymean,ystd,'o',...
            'CapSize',15,...
            'MarkerSize',7,...
            'MarkerFaceColor',[.3 .3 .3],...
            'Color',[.3 .3 .3],...
            'LineWidth',3);
        text([0.5 2.5],ymean+0.075,compose('%.2f', ymean),'FontSize',16,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);

        % stat
        [h,p] = ttest(y(1,:)',y(2,:)');
        Ysig=1;Yincre=0.014.*1000;sigincre=0.02;
        sigidx=find(p<0.05 & p>=0.01)';
        if sigidx
            for sn=sigidx
                plot(1:2,[Ysig Ysig],'color',[0 0 0],'linewidth',3)
                text(mean(1:2),Ysig+sigincre,'*','FontSize',30,...
                    'HorizontalAlignment','center',...
                    'color',[0 0 0]);
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(p<0.01 & p>=0.001)';
        if sigidx
            for sn=sigidx
                plot(1:2,[Ysig Ysig],'color',[0 0 0],'linewidth',3)
                text(mean(1:2),Ysig+sigincre,'**','FontSize',30,...
                    'HorizontalAlignment','center',...
                    'color',[0 0 0]);
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(p<0.001)';
        if sigidx
            for sn=sigidx
                plot(1:2,[Ysig Ysig],'color',[0 0 0],'linewidth',3)
                text(mean(1:2),Ysig+sigincre,'***','FontSize',30,...
                    'HorizontalAlignment','center',...
                    'color',[0 0 0]);
                Ysig=Ysig+Yincre;
            end
        end

        xticks([0.75 2.25])
        xticklabels({'PF' 'CF'});
        xlim([0 3])
        ylim([0 1.1])
        ylabel('CV of amp. (\DeltaF/F)')
        set(gca,'FontSize',18)
        set(gcf,'color',[1 1 1])
        figposition=get(gcf,'position');
        figposition(3)=300;
        set(gcf,'Position',figposition);
    end

    function time_to_max(Marker)
        %% Section 5: plotting the time to peak after stim
        LOCS=pksLOCS;
        amps=cell(size(condlist));
        y=[];
        for cond=1:max(condlist)
            condidx=find(condlist==cond)';
            for i=1:length(condidx)
                idx=condidx(i);
                T=time{idx}-time{idx}(data(idx).stim(1)+1);
                for roi=1:size(LOCS{idx},2)
                    activ=smoothBC_signal{idx}(:,roi,:);
                    TpostStime=T(T>ana_window(1) & T<=ana_window(2));
                    amps{idx}{roi,1} = activ(T>ana_window(1) & T<=ana_window(2),:,:);%filter out the data out of the ana_window
                    [~,amps{idx}{roi,1}] = max(amps{idx}{roi,1},[],1);
                    amps{idx}{roi,1} = TpostStime(amps{idx}{roi,1}(:));
                    cellmean{idx}{roi,1} = nanmean(amps{idx}{roi},1);
                end
            end
            catmean=vertcat(cellmean{condidx});
            y=cat(1,y,horzcat(catmean{:}));
        end


        figure
        x=1:max(condlist);
        yIDX=y(1,:)>0.3;
        hold on
        for n=1:size(y,2)
            jitterx=x'+0.*rand(size(x'))
            scatter(jitterx,y(:,n),50,...
                'MarkerFaceColor','flat',...
                'MarkerEdgeColor','flat',...
                'MarkerFaceAlpha',0.2,...
                'MarkerEdgeAlpha',0.2,...
                'MarkerFaceColor',0.5*[1 1 1],...
                'MarkerEdgeColor',0.5*[1 1 1]);
            h=plot(jitterx,y(:,n),...
                'color',[.5 .5 .5 .2],...
                'LineWidth',1);
            h.Color(4) = 1;
        end

        y_prob=y;
        hold on
        ymean=nanmean(y,2);
        ystd=nanstd(y,0,2)./sqrt(size(y,2));
        errorbar([0.5 2.5],ymean,ystd,'o',...
            'CapSize',15,...
            'MarkerSize',7,...
            'MarkerFaceColor',[.3 .3 .3],...
            'Color',[.3 .3 .3],...
            'LineWidth',3);
        text([0.5 2.5],ymean+0.012,compose('%.2f', ymean),'FontSize',16,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);

        % stat
        [h,p] = ttest(y(1,:)',y(2,:)');
        Ysig=0.265;Yincre=0.01;sigincre=0.003;
        sigidx=find(p<0.05 & p>=0.01)';
        if sigidx
            for sn=sigidx
                plot(1:2,[Ysig Ysig],'color',[0 0 0],'linewidth',3)
                text(mean(1:2),Ysig+sigincre,'*','FontSize',30,...
                    'HorizontalAlignment','center',...
                    'color',[0 0 0]);
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(p<0.01 & p>=0.001)';
        if sigidx
            for sn=sigidx
                plot(1:2,[Ysig Ysig],'color',[0 0 0],'linewidth',3)
                text(mean(1:2),Ysig+sigincre,'**','FontSize',30,...
                    'HorizontalAlignment','center',...
                    'color',[0 0 0]);
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(p<0.001)';
        if sigidx
            for sn=sigidx
                plot(1:2,[Ysig Ysig],'color',[0 0 0],'linewidth',3)
                text(mean(1:2),Ysig+sigincre,'***','FontSize',30,...
                    'HorizontalAlignment','center',...
                    'color',[0 0 0]);
                Ysig=Ysig+Yincre;
            end
        end

        xticks([0.75 2.25])
        xticklabels({'PF' 'CF'});
        xlim([0 3])
        ylim([0.125 0.275])
        yticks(0.15:0.05:0.25)
        ylabel('Time to peak (sec)')
        set(gca,'FontSize',18)
        set(gcf,'color',[1 1 1])
        figposition=get(gcf,'position');
        figposition(3)=300;
        set(gcf,'Position',figposition);
    end
end