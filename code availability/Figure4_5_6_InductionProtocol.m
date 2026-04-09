% 1. Add the path of the current script to the search path
% 2. Change the current folder to "Figure 4 5 6 induction protocol",
% 3. Using ctrl+enter within each function to run each section

%% Section 1: Open files
% Select all the files starting with "pulse_" and click "open",
% and then select all the files starting with "short_burst_" and click "open",
% and then select all the files starting with "tonic_" and click "open",
% and finally select all the files starting with "ramp_" and click "open",

clear
co=lines;
co(6,:)=[.5 .5 .5];
cnum=[6 5 2 4];

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

        %         end
    end
end

% import data
for i=1:size(FNlist,1)
    data(i)=load([FPlist{i},FNlist{i}]);
    baseline(i)=load([FPlist{i},FNlist{i}(1:end-12),'baseSignal.mat']);
end

% pull the data out
riseLOCS=arrayfun(@(x) data(x).riseLOCS,1:length(data),'uni',0);
upLOCS=arrayfun(@(x) data(x).upLOCS,1:length(data),'uni',0);
pksLOCS=arrayfun(@(x) data(x).pksLOCS,1:length(data),'uni',0);
PKS=arrayfun(@(x) data(x).PKS,1:length(data),'uni',0);
time=arrayfun(@(x) data(x).time,1:length(data),'uni',0);
smoothBC_signal=arrayfun(@(x) data(x).smoothBC_signal,1:length(data),'uni',0);
signal=arrayfun(@(x) data(x).signal,1:length(data),'uni',0);
baseSignal=arrayfun(@(x) baseline(x).signal_raw,1:length(baseline),'uni',0);

TimePhase=[1 3 5];%select the last recording during the pre phase as baseline
for i=1:length(signal)
    pcr=arrayfun(@(a) prctile(baseSignal{i}{TimePhase(a)}(:,:,end),20,1),1:3,'uni',0);
    pcr=cat(3,pcr{:});
    BC_signal{i}=(signal{i}-repmat(pcr,size(signal{i},1),1))./pcr;
end

%% Section 2: set up the induction protocol
incre=2;
PFstim=[];
CFstim=[];
% PF stim timing of each burst
% % pulse
D(1)=0.01783;
PFstim{1}=[0];% PF pulse
CFstim{1,1}=[0 0.01 0.02 0.03];protocol=2;% 100Hz CF
CFstim{1,2}=[0 0.003 0.006 0.009];protocol=3;% 333Hz CF
% short burst
D(2)=0.02549;
PFstim{2}=[0 0.02 0.04];% PF ramp
CFstim{2,1}=[0 0.01 0.02 0.03]+0.04;protocol=2;% 100Hz CF
CFstim{2,2}=[0 0.003 0.006 0.009]+0.04;protocol=3;% 333Hz CF
% % tonic
D(3)=0.0298;
PFstim{3}=[0 0.1 0.2 0.3 0.4];% PF tonic
CFstim{3,1}=[0 0.01 0.02 0.03]+0.4;protocol=2;% 100Hz CF
CFstim{3,2}=[0 0.003 0.006 0.009]+0.4;protocol=3;% 333Hz CF
% % ramp
D(4)=0.03213;
PFstim{4}=[0 0.1 0.2 0.3 0.35 0.37];% PF ramp
CFstim{4,1}=[0 0.01 0.02 0.03]+0.37;protocol=2;% 100Hz CF
CFstim{4,2}=[0 0.003 0.006 0.009]+0.37;protocol=3;% 333Hz CF

adjtime=(0:mean(diff(time{1})):time{1}(end))';

%% Section 3: plot mean +/- sem of calcium traces, alinged by the first CF stimulus (Fig. 6d)
protocol=2; % If protocol=2 plot protocol (1) in which PF stimulation was paired with 4CF at 100Hz
% If protocol=3 plot protocol (2) in which PF stimulation was paired with 4CF at 333Hz

cellactiv=cell(1,4);
plotidx=-15:30;
figure;hold on
plot([-1 1],[0 0],'k--')

for con=[1 3 4 2]
    t=adjtime-PFstim{con}(end)+D(con);
    activ=cat(2,BC_signal{condlist==con});

    for celln=1:size(activ,2)
        [~, idx] = min(abs(t - (1:8)), [], 1);
        activarray=arrayfun(@(c) activ(c+plotidx,celln,protocol),idx,'uni',0);
        cellactiv{con}(:,celln)=mean(horzcat(activarray{:}),2);
    end
    ErrArea_Smooth(t(idx(2)+plotidx)-t(idx(2)),mean(cellactiv{con},2),...
        nanstd(cellactiv{con},0,2)./sqrt(size(cellactiv{con},2)),...
        [co(cnum(con),:) .5]);hold on;
    plot(t(idx(2)+plotidx)-t(idx(2)),mean(cellactiv{con},2),...
        'color',co(cnum(con),:),...
        'LineWidth',1.5);

    %plotting the stimulus
    % y value
    PFy=repmat([0;incre*3],1,length(PFstim{con}(:)));%times busrst number
    PFy=repmat(1.8+[0;0.1]+0.2*(con-1),1,length(PFstim{con}(:)));%times busrst number

    % x value
    PFx=repmat(PFstim{con},2,1);

    plot(PFx-PFstim{con}(end),PFy,'color',[co(cnum(con),:) 1],'LineWidth',1.5);hold on
end
% plot CF stimulus
CFy=repmat([0;incre*3],1,length(CFstim{con,protocol-1}(:)));%times busrst number
CFy=repmat([0;1.7],1,length(CFstim{con,protocol-1}(:)));%times busrst number

CFx=repmat(CFstim{con,protocol-1},2,1);

plot(CFx-PFstim{con}(end),CFy,'color',[.5 .5 .5 .5],'LineWidth',1);

xlim([-0.5 1])
ylim([-0.1 2.6])
set(gca,'TickDir','out');

%% Section 4: plot mean +/- sem of calcium traces, full 10 sec recording (Fig. 4b, 4c, S1b, S1c, 5b, 5c, S2b, S2c)
protocol=2; % If protocol=2 plot protocol (1) in which PF stimulation was paired with 4CF at 100Hz
% If protocol=3 plot protocol (2) in which PF stimulation was paired with 4CF at 333Hz

cellactiv=cell(1,4);
plotidx=-15:30;

for con=1:condlist(end)
    %plotting the stimulus
    % y value
    PFy=repmat([0;incre*3],1,length(PFstim{con}(:)));%times busrst number
    PFy=repmat(1.8+[0;0.1],1,length(PFstim{con}(:)));%times busrst number
    PFy=repmat(PFy,1,10);%10 times per 10sec

    CFy=repmat([0;incre*3],1,length(CFstim{con,protocol-1}(:)));%times busrst number
    CFy=repmat([0;1.7],1,length(CFstim{con,protocol-1}(:)));%times busrst number
    CFy=repmat(CFy,1,10);%10 times per 10sec

    % x value
    TRep=repmat(0:9,length(PFstim{con}),1);
    PFx=repmat(PFstim{con},2,10)+repmat(TRep(:)',2,1);
    TRep=repmat(0:9,length(CFstim{con,protocol-1}),1);
    CFx=repmat(CFstim{con,protocol-1},2,10)+repmat(TRep(:)',2,1);

    figure
    plot(PFx-PFstim{con}(end),PFy,'color',[[238 133 74]./255 1],'linewidth',1.2);hold on
    plot(CFx-PFstim{con}(end),CFy,'color',[.3 .3 .3 .3],'linewidth',1);

    % plot the data
    t=adjtime-PFstim{con}(end)+D(con);
    activ=cat(2,BC_signal{condlist==con});

    ErrArea_Smooth(t,mean(activ(:,:,protocol),2),...
        nanstd(activ(:,:,protocol),0,2)./sqrt(size(activ(:,:,protocol),2)),...
        [co(cnum(con),:) .5]);hold on;
    plot(t,mean(activ(:,:,protocol),2),...
        'color',co(cnum(con),:),...
        'LineWidth',1.5);

    xlim([-0.5 10.5])
    %     ylim([0 1]*incre+(i-1)*incre)
    ylim([0 2])
    set(gca,'TickDir','out');
    %
    p=get(gcf,'position');p(3)=p(3)*2.6;p(4)=p(4)/2;
    set(gcf,'position',p)
end

%% Section 5: plot Delta baseline (Fig. 6f)

protocol=2;
amp=cell(1,4);
plotidx=-12:0;
for con=1:condlist(end)
    t=adjtime-PFstim{con}(end)+D(con);
    activ=cat(2,BC_signal{condlist==con});

    for c=1:8
        [~, idx] = min(abs(t - c), [], 1);
        tw=idx+plotidx;% PF
        amp{con}=[amp{con};arrayfun(@(a) activ(tw(end),a,protocol)-activ(tw(1),a,protocol),1:size(activ,2))];
    end
end

y=[];
ysem=[];
matconncat=[];
concatgroup=[];
for con=1:condlist(end)
    y(con)=mean(mean(amp{con},1));
    ysem(con)=nanstd(mean(amp{con},1),0,2)./sqrt(size(mean(amp{con},1),2));
    matconncat=[matconncat mean(amp{con},1)];
    concatgroup=[concatgroup con.*ones(size(mean(amp{con},1)))];
end

% One-way ANOVA stat
p=[];
stats=[];
data=matconncat;
g2=concatgroup;
[p.anova1,~,stats.anova1] = anova1(matconncat,g2);
c=multcompare(stats.anova1);

figure;
hold on;
plot([-1 5],[0 0],'k--')

MarkerSize=30;
meanSize=10;
space=0.01;
jrange=0.02;
interval=0.22;
arrayfun(@(a) errorbar(a*interval+space,y(a),ysem(a),'o-',...
    'color',co(cnum(a),:),...
    'CapSize',30,...
    'MarkerSize',meanSize,...
    'LineWidth',3),1:4,'uni',0);
arrayfun(@(a) text(a*interval+space,y(a)+ysem(a)+0.08,num2str(y(a),'%4.2f'),...
    'color',co(cnum(a),:),...
    'HorizontalAlignment','center',...
    'FontSize',15),1:4,'uni',0);

compgroup=[1 2 3 4];
arrayfun(@(a) scatter(a*interval-space-0.08-0.5*jrange+jrange*rand(1,length(data(g2==compgroup(a)))),data(g2==compgroup(a)),MarkerSize,...
    'MarkerFaceColor','flat',...
    'MarkerEdgeColor','flat',...
    'MarkerFaceAlpha',0.2,...
    'MarkerEdgeAlpha',0.2,...
    'CData',co(cnum(a),:)),1:4,'uni',0);

%significance
sigX=1;
sigY=0.1;
dsigY=0.04;
formatSpec = '%.3f';
FS=16;
Ysig=0.85;Yincre=0.09;sigincre=0.009;

sigidx=find(c(:,end)<0.05 & c(:,end)>=0.01)';
if sigidx
    for sn=sigidx
        plot(c(sn,1:2)*interval+space,[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(c(sn,1:2))*interval+space,Ysig+sigincre,'*','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end
sigidx=find(c(:,end)<0.01 & c(:,end)>=0.001)';
if sigidx
    for sn=sigidx
        plot(c(sn,1:2)*interval+space,[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(c(sn,1:2))*interval+space,Ysig+sigincre,'**','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end
sigidx=find(c(:,end)<0.001)';
if sigidx
    for sn=sigidx
        plot(c(sn,1:2)*interval+space,[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(c(sn,1:2))*interval+space,Ysig+sigincre,'***','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end

% mean traces alinged by CF stimulus and zero the AUC area
protocol=2;
cellactiv=cell(1,4);
for con=1:condlist(end)
    t=adjtime-PFstim{con}(end)+D(con);
    activ=cat(2,BC_signal{condlist==con});

    for celln=1:size(activ,2)
        %         tw=t>c+0.05 & t<c+0.3;
        [~, idx] = min(abs(t - (1:8)), [], 1);
        activarray=arrayfun(@(c) activ(c+plotidx,celln,protocol),idx,'uni',0);
        cellactiv{con}(:,celln)=mean(horzcat(activarray{:}),2);
    end

    %plotting the stimulus
    % y value
    PFy=repmat([0;incre*3],1,length(PFstim{con}(:)));%times busrst number
    PFy=repmat(0.8+[0;0.1]+0.2*(con-1),1,length(PFstim{con}(:)));%times busrst number

    % x value
    PFx=repmat(PFstim{con},2,1);

end

for con=1:condlist(end)
    ErrArea_Smooth(t(idx(2)+plotidx)-t(idx(2)),mean(cellactiv{con},2)-mean(cellactiv{con}(1,:)),...
        nanstd(cellactiv{con},0,2)./sqrt(size(cellactiv{con},2)),...
        [co(cnum(con),:) .5]);hold on;
    plot(t(idx(2)+plotidx)-t(idx(2)),mean(cellactiv{con},2)-mean(cellactiv{con}(1,:)),...
        'color',co(cnum(con),:),...
        'LineWidth',1.5);
end

xlim([-0.5 1])
ylim([-0.5 1.2])
set(gca,'TickDir','out');

%% Section 6: plot area under curve of 400ms window before CF input (Fig. 6g)

protocol=2;
amp=cell(1,4);
ts=mean(diff(t));
for con=1:condlist(end)
    t=adjtime-PFstim{con}(end)+D(con);
    activ=cat(2,BC_signal{condlist==con});
    for c=1:9
        tw=t>c-0.4 & t<c;
        amp{con}=[amp{con};arrayfun(@(a) sum(activ(tw,a,protocol)).*ts,1:size(activ,2))];
    end
end

y=[];
ysem=[];
matconncat=[];
concatgroup=[];
for con=1:condlist(end)
    y(con)=mean(mean(amp{con},1));
    ysem(con)=nanstd(mean(amp{con},1),0,2)./sqrt(size(mean(amp{con},1),2));
    matconncat=[matconncat mean(amp{con},1)];
    concatgroup=[concatgroup con.*ones(size(mean(amp{con},1)))];
end

% One-way ANOVA stat
p=[];
stats=[];
data=matconncat;
g2=concatgroup;
[p.anova1,~,stats.anova1] = anova1(matconncat,g2);
c=multcompare(stats.anova1);

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
    'LineWidth',3),1:4,'uni',0);
arrayfun(@(a) text(a+space,y(a)+ysem(a)+0.05,num2str(y(a),'%4.2f'),...
    'color',co(cnum(a),:),...
    'HorizontalAlignment','center',...
    'FontSize',15),1:4,'uni',0);

compgroup=[1 2 3 4];
arrayfun(@(a) scatter(a-space-0.5*jrange+jrange*rand(1,length(data(g2==compgroup(a)))),data(g2==compgroup(a)),MarkerSize,...
    'MarkerFaceColor','flat',...
    'MarkerEdgeColor','flat',...
    'MarkerFaceAlpha',0.2,...
    'MarkerEdgeAlpha',0.2,...
    'CData',co(cnum(a),:)),1:4,'uni',0);

%significance
sigX=1;
sigY=0.04;
dsigY=0.04;
formatSpec = '%.3f';
FS=16;
Ysig=0.65;Yincre=0.045;sigincre=0.005;

sigidx=find(c(:,end)<0.05 & c(:,end)>=0.01)';
if sigidx
    for sn=sigidx
        plot(c(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(c(sn,1:2)),Ysig+sigincre,'*','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end
sigidx=find(c(:,end)<0.01 & c(:,end)>=0.001)';
if sigidx
    for sn=sigidx
        plot(c(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(c(sn,1:2)),Ysig+sigincre,'**','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end
sigidx=find(c(:,end)<0.001)';
if sigidx
    for sn=sigidx
        plot(c(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(c(sn,1:2)),Ysig+sigincre,'***','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end


xticks(1:4)
xticklabels([]);
ylim([0 0.85])
xlim([0.5 4.5])
ylabel('AUC (\DeltaF/F; -400-0 ms)')
set(gca,'FontSize',23)
set(gcf,'color',[1 1 1])
%% Section 7: plot Relative Amplitude (Fig. 6h)

protocol=2;
amp=cell(1,4);
for con=1:condlist(end)
    t=adjtime-PFstim{con}(end)+D(con);
    activ=cat(2,BC_signal{condlist==con});
    for c=1:9
        tw1=find(t>=c-0.04 & t<=c);%IDX for pre CF level
        tw2=find(t>=c+0.05 & t<c+0.3);%IDX for CF amplitude

        amp{con}=[amp{con};arrayfun(@(a) max(activ(tw2,a,protocol)-activ(tw1(end),a,protocol)),1:size(activ,2))];
    end
end

y=[];
ysem=[];
matconncat=[];
concatgroup=[];
for con=1:condlist(end)
    y(con)=mean(mean(amp{con},1));
    ysem(con)=nanstd(mean(amp{con},1),0,2)./sqrt(size(mean(amp{con},1),2));
    matconncat=[matconncat mean(amp{con},1)];
    concatgroup=[concatgroup con.*ones(size(mean(amp{con},1)))];
end

% One-way ANOVA stat
p=[];
stats=[];
data=matconncat;
g2=concatgroup;
[p.anova1,~,stats.anova1] = anova1(matconncat,g2);
c=multcompare(stats.anova1);

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
    'LineWidth',3),1:4,'uni',0);
arrayfun(@(a) text(a+space,y(a)+ysem(a)+0.1,num2str(y(a),'%4.2f'),...
    'color',co(cnum(a),:),...
    'HorizontalAlignment','center',...
    'FontSize',15),1:4,'uni',0);

compgroup=[1 2 3 4];
arrayfun(@(a) scatter(a-space-0.5*jrange+jrange*rand(1,length(data(g2==compgroup(a)))),data(g2==compgroup(a)),MarkerSize,...
    'MarkerFaceColor','flat',...
    'MarkerEdgeColor','flat',...
    'MarkerFaceAlpha',0.2,...
    'MarkerEdgeAlpha',0.2,...
    'CData',co(cnum(a),:)),1:4,'uni',0);

%significance
sigX=1;
sigY=0.04;
dsigY=0.04;
formatSpec = '%.3f';
FS=16;
Ysig=1.9;Yincre=0.12;sigincre=0.0101;

sigidx=find(c(:,end)<0.05 & c(:,end)>=0.01)';
if sigidx
    for sn=sigidx
        plot(c(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(c(sn,1:2)),Ysig+sigincre,'*','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end
sigidx=find(c(:,end)<0.01 & c(:,end)>=0.001)';
if sigidx
    for sn=sigidx
        plot(c(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(c(sn,1:2)),Ysig+sigincre,'**','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end
sigidx=find(c(:,end)<0.001)';
if sigidx
    for sn=sigidx
        plot(c(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(c(sn,1:2)),Ysig+sigincre,'***','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end

xticks(1:4)
xticklabels([]);
ylim([0 2.3])
xlim([0.5 4.5])
ylabel('Rel. peak amp. (\DeltaF/F)')
set(gca,'FontSize',23)
set(gcf,'color',[1 1 1])

%% Section 8: plot absolute amplitude (Fig. 6i)

protocol=2;
amp=cell(1,4);
for con=1:condlist(end)
    t=adjtime-PFstim{con}(end)+D(con);
    activ=cat(2,BC_signal{condlist==con});
    for c=1:9
        tw=t>=c+0.05 & t<c+0.3;
        amp{con}=[amp{con};arrayfun(@(a) max(activ(tw,a,protocol)),1:size(activ,2))];
    end
end

y=[];
ysem=[];
matconncat=[];
concatgroup=[];
for con=1:condlist(end)
    y(con)=mean(mean(amp{con},1));
    ysem(con)=nanstd(mean(amp{con},1),0,2)./sqrt(size(mean(amp{con},1),2));
    matconncat=[matconncat mean(amp{con},1)];
    concatgroup=[concatgroup con.*ones(size(mean(amp{con},1)))];
end

% One-way ANOVA stat
p=[];
stats=[];
data=matconncat;
g2=concatgroup;%stimulus condition
[p.anova1,~,stats.anova1] = anova1(matconncat,g2);
c=multcompare(stats.anova1);

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
    'LineWidth',3),1:4,'uni',0);
arrayfun(@(a) text(a+space,y(a)+ysem(a)+0.3,num2str(y(a),'%4.2f'),...
    'color',co(cnum(a),:),...
    'HorizontalAlignment','center',...
    'FontSize',15),1:4,'uni',0);

compgroup=[1 2 3 4];
arrayfun(@(a) scatter(a-space-0.5*jrange+jrange*rand(1,length(data(g2==compgroup(a)))),data(g2==compgroup(a)),MarkerSize,...
    'MarkerFaceColor','flat',...
    'MarkerEdgeColor','flat',...
    'MarkerFaceAlpha',0.2,...
    'MarkerEdgeAlpha',0.2,...
    'CData',co(cnum(a),:)),1:4,'uni',0);

%significance
sigX=1;
sigY=0.04;
dsigY=0.04;
formatSpec = '%.3f';
FS=16;
Ysig=3.45;Yincre=0.25;sigincre=0.02;

sigidx=find(c(:,end)<0.05 & c(:,end)>=0.01)';
if sigidx
    for sn=sigidx
        plot(c(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(c(sn,1:2)),Ysig+sigincre,'*','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end
sigidx=find(c(:,end)<0.01 & c(:,end)>=0.001)';
if sigidx
    for sn=sigidx
        plot(c(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(c(sn,1:2)),Ysig+sigincre,'**','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end
sigidx=find(c(:,end)<0.001)';
if sigidx
    for sn=sigidx
        plot(c(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(c(sn,1:2)),Ysig+sigincre,'***','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end

xticks(1:4)
xticklabels([]);
ylim([0 4.5])
xlim([0.5 4.5])
ylabel('Abs. peak amp. (\DeltaF/F)')
set(gca,'FontSize',23)
set(gcf,'color',[1 1 1])

%% Section 9: generate LTD data
% all below is LTD vs induction protocol parameters
% 1. open the excel file "peak amplitude"
% 2. copy the data from excel file "peak amplitude" to the cell array "LTD"
%    copy the entire singl pulse data to LTD{1}
%    copy the entire short bust data to LTD{2}
%    copy the entire tonic bust data to LTD{3}
%    copy the entire ramping burst data to LTD{4}

LTD=cell(1,4);

%% Section 9: LTD vs Delta baseline (Fig. 6j)
protocol=2;
amp=cell(1,4);

for con=1:condlist(end)
    t=adjtime-PFstim{con}(end)+D(con);
    activ=cat(2,BC_signal{condlist==con});
    for c=1:8
        [~, idx] = min(abs(t - c), [], 1);
        tw=idx+plotidx;
        amp{con}=[amp{con};arrayfun(@(a) activ(tw(end),a,protocol)-activ(tw(1),a,protocol),1:size(activ,2))];
    end
end

x=[];
normy=cell(1,4);
xmean=[];
xStd=[];
ymean=[];
yStd=[];
for con=1:4
    x{con}=mean(amp{con},1);
    xmean(con)=mean(x{con});
    xStd(con)=nanstd(x{con},0,2)./sqrt(size(x{con},2));
    normy{con}=LTD{con}((protocol-1)*2+[-1:1],:)./LTD{con}((protocol-1)*2-1,:);
    ymean(:,con)=mean(normy{con},2);
    yStd(:,con)=nanstd(normy{con}-normy{con}(1,:),0,2)./sqrt(size(normy{con},2));
end

for EL=3
    figure;hold on
    MarkerSize=50;
    for con=4:-1:1
        scatter(xmean(con),ymean(EL,con)-ymean(1,con),MarkerSize,...
            'MarkerFaceColor','flat',...
            'MarkerEdgeColor','flat',...
            'MarkerFaceAlpha',0.8,...
            'MarkerEdgeAlpha',0.8,...
            'CData',co(cnum(con),:))
        errorbar(xmean(con),ymean(EL,con)-ymean(1,con),xStd(con),'horizontal',...
            "MarkerSize",MarkerSize,...
            'MarkerFaceColor',co(cnum(con),:),...
            'MarkerEdgeColor',co(cnum(con),:),...
            'LineWidth', 1.5, ...
            'Color',          co(cnum(con),:))
        errorbar(xmean(con),ymean(EL,con)-ymean(1,con),yStd(EL,con),...
            "MarkerSize",MarkerSize,...
            'MarkerFaceColor',co(cnum(con),:),...
            'MarkerEdgeColor',co(cnum(con),:),...
            'LineWidth', 1.5,...
            'Color',          co(cnum(con),:))
    end

    %linear fit pre
    mdl = fitlm(xmean,ymean(EL,:)-ymean(1,:),'Intercept',true);
    beta = mdl.Coefficients.Estimate;
    Xnew = linspace(min(xmean)-0.2, max(xmean)+0.2, 1000)';
    [ypred,yci] = predict(mdl, Xnew,'Alpha',0.05);

    hold on
    plot(Xnew, ypred, '-','LineWidth',1,'Color',[0 0 0 .5]);
    xp=-0.08;
    yp=0.1;
    marksize=40;
    [R,P] = corrcoef(xmean,ymean(EL,:)-ymean(1,:));

    xlabel('\Deltabaseline (\DeltaF/F)');
    ylabel('\Deltanorm. amplitude');

    if P(2)<0.001
        text(xp,yp,['R = ' num2str(R(2),'%.3f') '; {\it p} < 0.001'],...
            'FontSize',20);
    else
        text(xp,yp,['R = ' num2str(R(2),'%.3f') '; {\it p} = ' num2str(P(2),'%.3f')],...
            'FontSize',20);
    end
    ylim([-0.18 0.1])
    xlim([-0.1 0.21])
    axis square
    set(gca,'FontSize',23)
    set(gcf,'color',[1 1 1])
end

%% Section 9: LTD vs PF AUC (-400 - 0 ms) (Fig. 6k)

protocol=2;
amp=cell(1,4);
ts=mean(diff(t));
for con=1:condlist(end)
    t=adjtime-PFstim{con}(end)+D(con);
    activ=cat(2,BC_signal{condlist==con});
    for c=1:9
        tw=t>c-0.40 & t<c;% PF
        amp{con}=[amp{con};arrayfun(@(a) sum(activ(tw,a,protocol)).*ts,1:size(activ,2))];
    end
end

x=[];
normy=cell(1,4);
xmean=[];
xStd=[];
ymean=[];
yStd=[];
for con=1:4
    x{con}=mean(amp{con},1);
    xmean(con)=mean(x{con});
    xStd(con)=nanstd(x{con},0,2)./sqrt(size(x{con},2));
    normy{con}=LTD{con}((protocol-1)*2+[-1:1],:)./LTD{con}((protocol-1)*2-1,:);
    ymean(:,con)=mean(normy{con},2);
    yStd(:,con)=nanstd(normy{con}-normy{con}(1,:),0,2)./sqrt(size(normy{con},2));
end
xhorzcat=horzcat(x{:});
yhorzcat=horzcat(normy{:});

for EL=3
    figure;hold on
    MarkerSize=50;
    for con=4:-1:1
        scatter(xmean(con),ymean(EL,con)-ymean(1,con),MarkerSize,...
            'MarkerFaceColor','flat',...
            'MarkerEdgeColor','flat',...
            'MarkerFaceAlpha',0.8,...
            'MarkerEdgeAlpha',0.8,...
            'CData',co(cnum(con),:))
        errorbar(xmean(con),ymean(EL,con)-ymean(1,con),xStd(con),'horizontal',...
            "MarkerSize",MarkerSize,...
            'MarkerFaceColor',co(cnum(con),:),...
            'MarkerEdgeColor',co(cnum(con),:),...
            'LineWidth', 1.5, ...
            'Color',          co(cnum(con),:))
        errorbar(xmean(con),ymean(EL,con)-ymean(1,con),yStd(EL,con),...
            "MarkerSize",MarkerSize,...
            'MarkerFaceColor',co(cnum(con),:),...
            'MarkerEdgeColor',co(cnum(con),:),...
            'LineWidth', 1.5,...
            'Color',          co(cnum(con),:))
    end

    %linear fit pre
    mdl = fitlm(xmean,ymean(EL,:)-ymean(1,:),'Intercept',true);
    beta = mdl.Coefficients.Estimate;
    Xnew = linspace(min(xmean)-0.1, max(xmean)+0.1, 1000)';
    [ypred,yci] = predict(mdl, Xnew,'Alpha',0.05);

    hold on
    plot(Xnew, ypred, '-','LineWidth',1,'Color',[0 0 0 .5]);
    xp=0.02;
    yp=0.1;
    marksize=40;
    [R,P] = corrcoef(xmean,ymean(EL,:)-ymean(1,:));

    xlabel('AUC (\DeltaF/F;-400-0 ms)');
    ylabel('\Deltanorm. amplitude');

    if P(2)<0.001
        text(xp,yp,['R = ' num2str(R(2),'%.3f') '; {\it p} < 0.001'],...
            'FontSize',20);
    else
        text(xp,yp,['R = ' num2str(R(2),'%.3f') '; {\it p} = ' num2str(P(2),'%.3f')],...
            'FontSize',20);
    end
    ylim([-0.18 0.1])
    xlim([0 0.25])
    axis square
    set(gca,'FontSize',23)
    set(gcf,'color',[1 1 1])
end

%% MEAN - LTP vs rel amplitude (Fig. 6l)
protocol=2;
amp=cell(1,4);

for con=1:condlist(end)
    t=adjtime-PFstim{con}(end)+D(con);
    activ=cat(2,BC_signal{condlist==con});
    for c=1:9
        tw=find(t>=c+0.0 & t<c+0.3);
        amp{con}=[amp{con};arrayfun(@(a) max(activ(tw(2:end),a,protocol)-activ(tw(1),a,protocol)),1:size(activ,2))];
    end
end

% PF integral: plot mean+sem
x=[];
normy=cell(1,4);
xmean=[];
xStd=[];
ymean=[];
yStd=[];
for con=1:4
    x{con}=mean(amp{con},1);
    xmean(con)=mean(x{con});
    xStd(con)=nanstd(x{con},0,2)./sqrt(size(x{con},2));
    normy{con}=LTD{con}((protocol-1)*2+[-1:1],:)./LTD{con}((protocol-1)*2-1,:);
    ymean(:,con)=mean(normy{con},2);
    yStd(:,con)=nanstd(normy{con}-normy{con}(1,:),0,2)./sqrt(size(normy{con},2));
end

for EL=3
    figure;hold on
    MarkerSize=50;

    for con=4:-1:1
        scatter(xmean(con),ymean(EL,con)-ymean(1,con),MarkerSize,...
            'MarkerFaceColor','flat',...
            'MarkerEdgeColor','flat',...
            'MarkerFaceAlpha',0.8,...
            'MarkerEdgeAlpha',0.8,...
            'CData',co(cnum(con),:))
        errorbar(xmean(con),ymean(EL,con)-ymean(1,con),xStd(con),'horizontal',...
            "MarkerSize",MarkerSize,...
            'MarkerFaceColor',co(cnum(con),:),...
            'MarkerEdgeColor',co(cnum(con),:),...
            'LineWidth', 1.5, ...
            'Color',          co(cnum(con),:))
        errorbar(xmean(con),ymean(EL,con)-ymean(1,con),yStd(EL,con),...
            "MarkerSize",MarkerSize,...
            'MarkerFaceColor',co(cnum(con),:),...
            'MarkerEdgeColor',co(cnum(con),:),...
            'LineWidth', 1.5,...
            'Color',          co(cnum(con),:))
    end

    %linear fit pre
    mdl = fitlm(xmean,ymean(EL,:)-ymean(1,:),'Intercept',true);
    beta = mdl.Coefficients.Estimate;
    Xnew = linspace(min(xmean)-0.2, max(xmean)+0.2, 1000)';
    [ypred,yci] = predict(mdl, Xnew,'Alpha',0.05);

    hold on
    plot(Xnew, ypred, '-','LineWidth',1,'Color',[0 0 0 .5]);
    xp=0.43;
    yp=0.1;
    marksize=40;
    [R,P] = corrcoef(xmean,ymean(EL,:)-ymean(1,:));

    xlabel('Rel. peak amp (\DeltaF/F)');
    ylabel('\Deltanorm. amplitude');



    if P(2)<0.001
        text(xp,yp,['R = ' num2str(R(2),'%.3f') '; {\it p} < 0.001'],...
            'FontSize',20);
    else
        text(xp,yp,['R = ' num2str(R(2),'%.3f') '; {\it p} = ' num2str(P(2),'%.3f')],...
            'FontSize',20);
    end
    % xticks(0.8:0.4:1.6)
    % yticks(-0.15:0.1:0.05)
    ylim([-0.18 0.1])

    xlim([0.4 0.85])
    axis square
    set(gca,'FontSize',23)
    set(gcf,'color',[1 1 1])
end

%% Section 9: LTP vs abs amplitude (Fig. 6m)
protocol=2;
amp=cell(1,4);
for con=1:condlist(end)
    t=adjtime-PFstim{con}(end)+D(con);
    activ=cat(2,BC_signal{condlist==con});
    for c=1:9
        tw=t>=c+0.05 & t<c+0.3;
        amp{con}=[amp{con};arrayfun(@(a) max(activ(tw,a,protocol)),1:size(activ,2))];
    end
end

% PF integral: plot mean+sem
x=[];
normy=cell(1,4);
xmean=[];
xStd=[];
ymean=[];
yStd=[];
for con=1:4
    x{con}=mean(amp{con},1);
    xmean(con)=mean(x{con});
    xStd(con)=nanstd(x{con},0,2)./sqrt(size(x{con},2));
    normy{con}=LTD{con}((protocol-1)*2+[-1:1],:)./LTD{con}((protocol-1)*2-1,:);
    ymean(:,con)=mean(normy{con},2);
    yStd(:,con)=nanstd(normy{con}-normy{con}(1,:),0,2)./sqrt(size(normy{con},2));
end

for EL=3
    figure;hold on
    MarkerSize=50;

    for con=4:-1:1
        scatter(xmean(con),ymean(EL,con)-ymean(1,con),MarkerSize,...
            'MarkerFaceColor','flat',...
            'MarkerEdgeColor','flat',...
            'MarkerFaceAlpha',0.8,...
            'MarkerEdgeAlpha',0.8,...
            'CData',co(cnum(con),:))
        errorbar(xmean(con),ymean(EL,con)-ymean(1,con),xStd(con),'horizontal',...
            "MarkerSize",MarkerSize,...
            'MarkerFaceColor',co(cnum(con),:),...
            'MarkerEdgeColor',co(cnum(con),:),...
            'LineWidth', 1.5, ...
            'Color',          co(cnum(con),:))
        errorbar(xmean(con),ymean(EL,con)-ymean(1,con),yStd(EL,con),...
            "MarkerSize",MarkerSize,...
            'MarkerFaceColor',co(cnum(con),:),...
            'MarkerEdgeColor',co(cnum(con),:),...
            'LineWidth', 1.5,...
            'Color',          co(cnum(con),:))
    end

    %linear fit pre
    mdl = fitlm(xmean,ymean(EL,:)-ymean(1,:),'Intercept',true);
    beta = mdl.Coefficients.Estimate;
    Xnew = linspace(min(xmean)-0.2, max(xmean)+0.2, 1000)';
    [ypred,yci] = predict(mdl, Xnew,'Alpha',0.05);

    hold on
    plot(Xnew, ypred, '-','LineWidth',1,'Color',[0 0 0 .5]);
    xp=0.7;
    yp=0.1;
    marksize=40;
    [R,P] = corrcoef(xmean,ymean(EL,:)-ymean(1,:));

    xlabel('Abs. peak amp (\DeltaF/F)');
    ylabel('\Deltanorm. amplitude');

    if P(2)<0.001
        text(xp,yp,['R = ' num2str(R(2),'%.3f') '; {\it p} < 0.001'],...
            'FontSize',20);
    else
        text(xp,yp,['R = ' num2str(R(2),'%.3f') '; {\it p} = ' num2str(P(2),'%.3f')],...
            'FontSize',20);
    end
    xticks(0.8:0.4:1.6)
    ylim([-0.18 0.1])
    xlim([0.6 1.7])
    axis square
    set(gca,'FontSize',23)
    set(gcf,'color',[1 1 1])
end

%% Section 9: LTD vs PF AUC (-200 - 0 ms) (Fig. 6n-q)
AnaWind=[-0.4 -0.35;-0.25 -0.2;-0.15 -0.1;-0.05 0]; %time window to plot AUC-LTD relationship

for ana=1:4
    protocol=2;
    amp=cell(1,4);
    for con=1:condlist(end)
        t=adjtime-PFstim{con}(end)+D(con);
        ts=mean(diff(t));%sampling rate

        activ=cat(2,BC_signal{condlist==con});
        for c=1:9
            tw=t>c+AnaWind(ana,1) & t<c+AnaWind(ana,2);% PF
            amp{con}=[amp{con};arrayfun(@(a) sum(activ(tw,a,protocol)).*ts,1:size(activ,2))];
        end
    end

    x=[];
    normy=cell(1,4);
    xmean=[];
    xStd=[];
    ymean=[];
    yStd=[];
    for con=1:4
        x{con}=mean(amp{con},1);
        xmean(con)=mean(x{con});
        xStd(con)=nanstd(x{con},0,2)./sqrt(size(x{con},2));
        normy{con}=LTD{con}((protocol-1)*2+[-1:1],:)./LTD{con}((protocol-1)*2-1,:);
        ymean(:,con)=mean(normy{con},2);
        yStd(:,con)=nanstd(normy{con}-normy{con}(1,:),0,2)./sqrt(size(normy{con},2));
    end
    xhorzcat=horzcat(x{:});
    yhorzcat=horzcat(normy{:});

    for EL=3
        figure;hold on
        MarkerSize=50;

        for con=4:-1:1
            scatter(xmean(con),ymean(EL,con)-ymean(1,con),MarkerSize,...
                'MarkerFaceColor','flat',...
                'MarkerEdgeColor','flat',...
                'MarkerFaceAlpha',0.8,...
                'MarkerEdgeAlpha',0.8,...
                'CData',co(cnum(con),:))
            errorbar(xmean(con),ymean(EL,con)-ymean(1,con),xStd(con),'horizontal',...
                "MarkerSize",MarkerSize,...
                'MarkerFaceColor',co(cnum(con),:),...
                'MarkerEdgeColor',co(cnum(con),:),...
                'LineWidth', 1.5, ...
                'Color',          co(cnum(con),:))
            errorbar(xmean(con),ymean(EL,con)-ymean(1,con),yStd(EL,con),...
                "MarkerSize",MarkerSize,...
                'MarkerFaceColor',co(cnum(con),:),...
                'MarkerEdgeColor',co(cnum(con),:),...
                'LineWidth', 1.5,...
                'Color',          co(cnum(con),:))
        end

        %linear fit pre
        mdl = fitlm(xmean,ymean(EL,:)-ymean(1,:),'Intercept',true);
        beta = mdl.Coefficients.Estimate;
        Xnew = linspace(min(xmean)-0.1, max(xmean)+0.1, 1000)';
        [ypred,yci] = predict(mdl, Xnew,'Alpha',0.05);

        hold on
        plot(Xnew, ypred, '-','LineWidth',1,'Color',[0 0 0 .5]);
        xp=0.007;
        yp=0.1;
        marksize=40;
        [R,P] = corrcoef(xmean,ymean(EL,:)-ymean(1,:));

        xlabels = {
            'AUC (\DeltaF/F;-400 - -350 ms)'
            'AUC (\DeltaF/F;-250 - -200 ms)'
            'AUC (\DeltaF/F;-150 - -100 ms)'
            'AUC (\DeltaF/F;-50 - 0 ms)'};

        xlabel(xlabels{ana})
        ylabel('\Deltanorm. amplitude');

        if P(2)<0.001
            text(xp,yp,['R = ' num2str(R(2),'%.3f') '; {\it p} < 0.001'],...
                'FontSize',20);
        else
            text(xp,yp,['R = ' num2str(R(2),'%.3f') '; {\it p} = ' num2str(P(2),'%.3f')],...
                'FontSize',20);
        end
%         yticks([-0.1 0])
        ylim([-0.18 0.04])
        ylim([-0.18 0.1])

        xlim([0.005 0.035])
        axis square
        set(gca,'FontSize',23)
        set(gcf,'color',[1 1 1])
    end
end

%% Section 10: LTD vs PF AUC (every 50 ms) (Fig. 6r)
step=0.05
AnaWind=(-0.4:step:-step)';
AnaWind=[AnaWind AnaWind+step];
R=[];
P=[];
for ana=1:length(AnaWind)
    protocol=2;
    amp=cell(1,4);
    for con=1:condlist(end)
        t=adjtime-PFstim{con}(end)+D(con);
        ts=mean(diff(t));
        activ=cat(2,BC_signal{condlist==con});
        for c=1:9
            tw=t>c+AnaWind(ana,1) & t<c+AnaWind(ana,2);% PF
            amp{con}=[amp{con};arrayfun(@(a) sum(activ(tw,a,protocol)).*ts,1:size(activ,2))];
        end
    end

    x=[];
    normy=cell(1,4);
    xmean=[];
    ymean=[];
    for con=1:4
        x{con}=mean(amp{con},1);
        xmean(con)=mean(x{con});
        normy{con}=LTD{con}((protocol-1)*2+[-1:1],:)./LTD{con}((protocol-1)*2-1,:);
        ymean(:,con)=mean(normy{con},2);
    end
    xhorzcat=horzcat(x{:});
    yhorzcat=horzcat(normy{:});

    EL=3
    MarkerSize=50
    [r,p] = corrcoef(xmean,ymean(EL,:)-ymean(1,:));
    R(ana,:)=r(2);
    P(ana,:)=p(2);
end

figure;hold on
yyaxis left
plot(mean(AnaWind,2)*1000,R,'o-',...
    'color',[.5 .5 .5],...
    'MarkerFaceColor',[1 1 1],...
    'MarkerEdgeColor',[.5 .5 .5],...
    'linewidth',2);
ylim([-1 -0.5])
yticks([ -0.95 -0.9 -0.8 -0.6])
ylabel('R');
ax = gca;
ax.YColor = [0 0 0];

yyaxis right
plot(mean(AnaWind([1 4 6 8],:),2)*1000,P([1 4 6 8]),'o',...
    'color',[.5 .5 .5],...
    'MarkerFaceColor',[.5 .5 .5],...
    'MarkerEdgeColor',[.5 .5 .5],...
    'MarkerSize',10);
ylim([0 0.5])
yticks([ 0.05 0.1 0.2 0.4])
ylabel('{\it p}');
ax.YColor = [0 0 0];

xticks(-400:100:0)
set(gca,'FontSize',23)
set(gcf,'color',[1 1 1])
xlabel('Time before CF input (ms)')
grid on
p=get(gcf,'position');
p(4)=p(4).*1.2;
set(gcf,'position',p)



%% PF induced change plus raw traces_individual condition
protocol=2;
amp=cell(1,4);
ts=mean(diff(t));
plotidx=-12:0;
for con=1:condlist(end)
    t=adjtime-PFstim{con}(end)+D(con);

    activ=cat(2,BC_signal{condlist==con});
    %     figure;hold on


    for c=1:8
        [~, idx] = min(abs(t - c), [], 1);

        tw=idx+plotidx;% PF
        %         tw=t>c+0.00 & t<c+0.4;% PF
        amp{con}=[amp{con};arrayfun(@(a) activ(tw(end),a,protocol)-activ(tw(1),a,protocol),1:size(activ,2))];
        %         arrayfun(@(a) plot(t(tw),activ(tw,a,protocol)),1:size(activ,2));

    end
    %     figure;plot(amp{con});
end

y=[];
ysem=[];
matconncat=[];
concatgroup=[];
for con=1:condlist(end)
    y(con)=mean(mean(amp{con},1));
    ysem(con)=nanstd(mean(amp{con},1),0,2)./sqrt(size(mean(amp{con},1),2));
    matconncat=[matconncat mean(amp{con},1)];
    concatgroup=[concatgroup con.*ones(size(mean(amp{con},1)))];
end

protocol=2;
cellactiv=cell(1,4);
for con=1:condlist(end)

    figure;
    hold on;
    plot([-1 5],[0 0],'k--')

    MarkerSize=30;
    meanSize=10;
    space=0.01;
    jrange=0.02;
    interval=0.22;
    arrayfun(@(a) errorbar(a*interval+space,y(a),ysem(a),'o-',...
        'color',co(cnum(a),:),...
        'CapSize',30,...
        'MarkerSize',meanSize,...
        'LineWidth',3),1:4,'uni',0);
    arrayfun(@(a) text(a*interval+space,y(a)+ysem(a)+0.08,num2str(y(a),'%4.2f'),...
        'color',co(cnum(a),:),...
        'HorizontalAlignment','center',...
        'FontSize',15),1:4,'uni',0);

    compgroup=[1 2 3 4];
    arrayfun(@(a) scatter(a*interval-space-0.08-0.5*jrange+jrange*rand(1,length(data(g2==compgroup(a)))),data(g2==compgroup(a)),MarkerSize,...
        'MarkerFaceColor','flat',...
        'MarkerEdgeColor','flat',...
        'MarkerFaceAlpha',0.2,...
        'MarkerEdgeAlpha',0.2,...
        'CData',co(cnum(a),:)),1:4,'uni',0);

    t=adjtime-PFstim{con}(end)+D(con);
    activ=cat(2,BC_signal{condlist==con});

    for celln=1:size(activ,2)
        [~, idx] = min(abs(t - (1:8)), [], 1);
        activarray=arrayfun(@(c) activ(c+plotidx,celln,protocol),idx,'uni',0);
        cellactiv{con}(:,celln)=mean(horzcat(activarray{:}),2);
    end

    plot(t(idx(2)+plotidx)-t(idx(2)),cellactiv{con}-cellactiv{con}(1,:),...
        'color',[co(cnum(con),:) 0.5],...
        'LineWidth',0.2);

    % y value
    PFy=repmat([0;incre*3],1,length(PFstim{con}(:)));%times busrst number
    PFy=repmat(0.8+[0;0.1]+0.2*(con-1),1,length(PFstim{con}(:)));%times busrst number

    % x value
    PFx=repmat(PFstim{con},2,1);
    plot(t(idx(2)+plotidx)-t(idx(2)),mean(cellactiv{con},2)-mean(cellactiv{con}(1,:)),...
        'color',co(cnum(con),:),...
        'LineWidth',1.5);

    xlim([-0.5 1])
    ylim([-0.5 1.2])
    set(gca,'TickDir','out');
end

