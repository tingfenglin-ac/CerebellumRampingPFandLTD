% 1. Add the path of the current script to the search path
% 2. open the excel file "peak amplitude"
% 3. Using ctrl+enter within each function to run each section

%% Section 1: open files
% ctrl+enter to run this section
% copy the data from excel file "peak amplitude" to the following matices:
% con_data, ton_data, ram_data, shr_data

clear
co=lines;
co(6,:)=[.5 .5 .5];
cnum=[6 2 4 5];
colorlist=[];
markersize=30;
criteria=[0 100];
groupnum=4;

con_data=[]; % copy the first three rows (pre, early, late) of singl-pulse data to here
shr_data=[]; % copy the first three rows (pre, early, late) of short burst data to here
ton_data=[]; % copy the first three rows (pre, early, late) of tonic burst data to here
ram_data=[]; % copy the first three rows (pre, early, late) of ramping burst data to here

%% Section 2: normalization
ncon_data=con_data./con_data(1,:);
nshr_data=shr_data./shr_data(1,:);
nton_data=ton_data./ton_data(1,:);
nram_data=ram_data./ram_data(1,:);

%% Section 3: two-way ANOVA test
con_mean=nanmean(ncon_data,2);
con_sem=nanstd(ncon_data,0,2)./sqrt(size(ncon_data,2));
matconncat=[ncon_data];
concatgroup=ones(size(ncon_data));

shr_mean=nanmean(nshr_data,2);
shr_sem=nanstd(nshr_data,0,2)./sqrt(size(nshr_data,2));
matconncat=[matconncat nshr_data];
concatgroup=[concatgroup 2.*ones(size(nshr_data))];

ton_mean=nanmean(nton_data,2);
ton_sem=nanstd(nton_data,0,2)./sqrt(size(nton_data,2));
matconncat=[matconncat nton_data];
concatgroup=[concatgroup 3.*ones(size(nton_data))];

ram_mean=nanmean(nram_data,2);
ram_sem=nanstd(nram_data,0,2)./sqrt(size(nram_data,2));
matconncat=[matconncat nram_data];
concatgroup=[concatgroup 4.*ones(size(nram_data))];

conn=3;
data=matconncat(1:end);
g1=repmat(1:conn,1,length(matconncat));%time
g2=concatgroup(1:end);%genotype
[p.anovan,~,stats.anovan] = anovan(data,{g1,g2},'model','interaction','varnames',{'time','genotype'});
c.anovan=multcompare(stats.anovan,'Dimension',[1,2],'Display','off');
g2mat=repmat(1:groupnum,conn,1);
groupindex=[repmat(1:conn,1,groupnum);g2mat(1:end)];%index of combinations and the corresponding groups

%% Section 4: plot mean +/- sem with significance signs
cnum=[6 5 2 4];
comparegroup=3;
x=repmat((1:comparegroup)',1,groupnum);
y=[con_mean shr_mean ton_mean ram_mean];
ysem=[con_sem shr_sem ton_sem ram_sem];
figure;hold on;

MarkerSize=30;
meanSize=10;
space=0.06;
jrange=0.04;
interval=0.2;

plot([0 4],[0 0],'--','Color',[0 0 0]);

arrayfun(@(a) errorbar((1:comparegroup)'+(a-1)*interval,...
    y(1:comparegroup,a)-1,...
    ysem(1:comparegroup,a),'o',...
    'color',co(cnum(a),:),...
    'CapSize',10,...
    'MarkerSize',5,...
    'LineWidth',3),1:groupnum,'uni',0);

arrayfun(@(a) scatter((a-1)*interval-space-0.5*jrange+jrange*rand(size(matconncat(2:3,concatgroup(1,:)==a)))+ones(size(matconncat(2:3,concatgroup(1,:)==a))).*[2:3]',...
    matconncat(2:3,concatgroup(1,:)==a)-1,MarkerSize,...
    'MarkerFaceColor','flat',...
    'MarkerEdgeColor','flat',...
    'MarkerFaceAlpha',0.2,...
    'MarkerEdgeAlpha',0.2,...
    'CData',co(cnum(a),:)),1:4,'uni',0);

% half-violin plots with white boxplots, jittered scatter and linkline 
violin_x = [2:3;
            (2:3)+interval;
            (2:3)+2*interval;
            (2:3)+3*interval];
h=daviolinplot(matconncat(2:3,:)'-1,'groups',concatgroup(1,:)',...
    'xvalues',violin_x,...
    'violin','half',...
    'colors',co(cnum,:),...
    'violinalpha',0.5,...
    'violinwidth',0.5,...
    'outliers',0,...
    'box' ,0);
set(h.ds,'EdgeColor','none');



xticks([1:3]+2*interval-0.05)
xticklabels({'Pre','Early-Pre','Late-Pre'});
ylim([-0.95 0.65])
xlim([1.8 3.8])
ylabel('\Delta Norm. amplitude')
set(gca,'FontSize',23)
set(gcf,'color',[1 1 1])
figposition=get(gcf,'position');
figposition(3)=700;
% figposition(2)=-500;
set(gcf,'Position',figposition);


clearvars sig

gcomb=[1 2 3 4];% combinations of groups to compare
Yincre=0.07;textincre=0.001;text_xspace=-0.03;fontsize=28;
    % plot significance of anovan over time
for i=1:size(gcomb,1)
        Ysig=0.3;
    for g=1:groupnum%plot significance of each group
        idx=find(groupindex(2,:)==gcomb(i,g));%find index of combinations between time
        sig(1,:)=[1 2 c.anovan(c.anovan(:,1)==idx(1) & c.anovan(:,2)==idx(2),end)];% x1 x2 and p-value
        sig(2,:)=[1 3 c.anovan(c.anovan(:,1)==idx(1) & c.anovan(:,2)==idx(3),end)];% x1 x2 and p-value

        sigidx=find(sig(:,end)<0.05 & sig(:,end)>=0.01)';
        if sigidx
            for sn=sigidx
                plot(sig(sn,1:2)+(g-1)*interval,[Ysig Ysig],'color',co(cnum(gcomb(i,g)),:),'linewidth',1.5)
                text(mean(sig(sn,1:2))+(g-1)*interval,Ysig+textincre,'*','FontSize',fontsize,...
                    'HorizontalAlignment','center',...
                    'color',co(cnum(gcomb(i,g)),:));
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(sig(:,end)<0.01 & sig(:,end)>=0.001)';
        if sigidx
            for sn=sigidx
                plot(sig(sn,1:2)+(g-1)*interval,[Ysig Ysig],'color',co(cnum(gcomb(i,g)),:),'linewidth',1.5)
                text(mean(sig(sn,1:2))+(g-1)*interval,Ysig+textincre,'**','FontSize',fontsize,...
                    'HorizontalAlignment','center',...
                    'color',co(cnum(gcomb(i,g)),:));
                Ysig=Ysig+Yincre;
            end
        end
        sigidx=find(sig(:,end)<0.001)';
        if sigidx
            for sn=sigidx
%                 plot(sig(sn,1:2)+(g-1)*interval,[Ysig Ysig],'color',co(cnum(gcomb(i,g)),:),'linewidth',1.5)
                text(mean(sig(sn,2))+(g-1)*interval,Ysig,'***','FontSize',fontsize,...
                    'HorizontalAlignment','center',...
                    'color',co(cnum(gcomb(i,g)),:));
%                 Ysig=Ysig+Yincre;
            end
        end
    end
end

for i=1:size(gcomb,1)
    % plot significance of anovan between groups
    for t=2:comparegroup%plot significance of each time
            Ysig=1.38-1;

        idx=find(groupindex(2,:)==gcomb(i,1) & groupindex(1,:)==t);%find index the first group
        idx(2)=find(groupindex(2,:)==gcomb(i,2) & groupindex(1,:)==t);%find index the second group
        idx(3)=find(groupindex(2,:)==gcomb(i,3) & groupindex(1,:)==t);%find index the second group
        idx(4)=find(groupindex(2,:)==gcomb(i,4) & groupindex(1,:)==t);%find index the second group

        genocomp=c.anovan(c.anovan(:,1)==idx(1) | c.anovan(:,1)==idx(2) | c.anovan(:,1)==idx(3)| c.anovan(:,1)==idx(4),:);%find the first column belong to either groups
        genocomp=genocomp(genocomp(:,2)==idx(1) | genocomp(:,2)==idx(2) | genocomp(:,2)==idx(3)| genocomp(:,2)==idx(4),:);%find the second column belong to either groups

%         Xsig=0.6;Xincre=0.1;
        sigidx=find(genocomp(:,end)<0.05 & genocomp(:,end)>=0.01)';
        if sigidx
            for sn=sigidx
                plot(i+interval*groupindex(2,genocomp(sn,[1 2])),[Ysig Ysig],'color',[0 0 0],'linewidth',1.5,...
                    'color',[0 0 0 0.5])
                tx=text(mean(groupindex(2,genocomp(sn,[1 2])))+text_xspace,Ysig+textincre,'*','FontSize',fontsize,...
                    'HorizontalAlignment','center',...
                    'color',[0 0 0 0.2])
%                 tx.Rotation=90;
                Ysig=Ysig+Yincre;
%                 Xsig=Xsig+Xincre;
            end
        end
        sigidx=find(genocomp(:,end)<0.01 & genocomp(:,end)>=0.001)';
        if sigidx
            for sn=sigidx
                plot(i+interval*groupindex(2,genocomp(sn,[1 2])),[Ysig Ysig],'color',[0 0 0],'linewidth',1.5,...
                    'color',[0 0 0 0.5])
                tx=text(mean(groupindex(2,genocomp(sn,[1 2])))+text_xspace,Ysig+textincre,'**','FontSize',fontsize,...
                    'HorizontalAlignment','center',...
                    'color',[0 0 0 0.2])
%                 tx.Rotation=90;
                Ysig=Ysig+Yincre;
%                 Xsig=Xsig+Xincre;
            end
        end
        sigidx=find(genocomp(:,end)<0.001)';
        if sigidx
            for sn=sigidx
                plot(t+interval*(groupindex(2,genocomp(sn,[1 2]))-1),[Ysig Ysig],'color',[0 0 0],'linewidth',1.5,...
                    'color',[0 0 0 0.5])
                tx=text(t+interval*(mean(groupindex(2,genocomp(sn,[1 2])))-0.85)+text_xspace,Ysig+textincre,'***','FontSize',fontsize,...
                    'HorizontalAlignment','center',...
                    'color',[0 0 0 0.2])
%                 tx.Rotation=90;
                Ysig=Ysig+Yincre;
%                 Xsig=Xsig+Xincre;
            end
        end
    end
end

