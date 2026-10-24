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
legendlist=[];
markersize=30;
criteria=[0 100];
groupnum=4;

con_data=[]; % copy the entire singl-pulse data to here
shr_data=[]; % copy the entire short burst data to here
ton_data=[]; % copy the entire tonic burst data to here
ram_data=[]; % copy the entire ramping burst data to here


%%
clearvars -except co cnum legendlist colorlist markersize criteria legendlist groupnum con_data ton_data ram_data shr_data

%% Section 2: normalization
ncon_data=con_data./con_data(1,:);
nton_data=ton_data./ton_data(1,:);
nram_data=ram_data./ram_data(1,:);
nshr_data=shr_data./shr_data(1,:);

%% Section 3: rearrange data for statistics and ploting
% ctrl+enter to run this section
con_mean=nanmean(ncon_data,2);
con_sem=nanstd(ncon_data,0,2)./sqrt(size(ncon_data,2));
matconncat=[ncon_data];
concatgroup=ones(size(ncon_data));

ton_mean=nanmean(nton_data,2);
ton_sem=nanstd(nton_data,0,2)./sqrt(size(nton_data,2));
matconncat=[matconncat nton_data];
concatgroup=[concatgroup 2.*ones(size(nton_data))];

ram_mean=nanmean(nram_data,2);
ram_sem=nanstd(nram_data,0,2)./sqrt(size(nram_data,2));
matconncat=[matconncat nram_data];
concatgroup=[concatgroup 3.*ones(size(nram_data))];

shr_mean=nanmean(nshr_data,2);
shr_sem=nanstd(nshr_data,0,2)./sqrt(size(nshr_data,2));
matconncat=[matconncat nshr_data];
concatgroup=[concatgroup 4.*ones(size(nshr_data))];

%% Section 4: plot comparison of peak amplitude over time (Fig. 4 and Fig. 5)
comparegroup=3; % 3 for Fig. 4e, g and Fig. 5e, g; 5 for Fig. S1e, g and Fig. S2e, g

cnum=[6 2 4 5];
for g=1:groupnum

    % 1 way anova
    conn=size(con_data,1);
    data=matconncat(1:5,concatgroup(1,:)==g);
    [p,~,stats] = anova1(data');
    c=multcompare(stats,'Display','off');

    %plotting
    figure;hold on;

    % proportion
        SHT=0.15;
        PlastPopu(g,1)=sum(data(3,:)>SHT+1);
        PlastPopu(g,2)=sum(data(3,:)<=SHT+1 & data(3,:)>-SHT+1);
        PlastPopu(g,3)=sum(data(3,:)<-SHT+1);
    if comparegroup==3 %plot
        numcell=num2cell(PlastPopu(g,:)./size(data,2).*100)
        text(ones(1,3).*3.2, [SHT 0 -SHT].*2+1, cellfun(@(a) [num2str(a,'%.f') '%'], numcell, 'uniform', 0),'FontSize',23);
        plot([0 0;5 5],[SHT -SHT;SHT -SHT]+1,'--','Color',[.3 .3 .3],'linewidth',1.5)
    end

    % mean +/- sem
    x=(1:conn)';
    y=nanmean(data,2);
    ysem=nanstd(data,0,2)./sqrt(size(data,2));
    plot((1:comparegroup)',data((1:comparegroup),:),...
        'color',[co(cnum(g),:) 0.2],...
        'LineWidth',1.5);
    errorbar((1:comparegroup)',y(1:comparegroup,:),ysem(1:comparegroup,:),'o-',...
        'color',co(cnum(g),:),...
        'CapSize',20,...
        'MarkerSize',5,...
        'LineWidth',3);

    xticks(1:7)
    xticklabels({'Pre','Early','Late','Early','Late','Early','Late'});
    if comparegroup==3; xticklabels({'Pre','Early','Late'}); end
    xlim([0.5 comparegroup+0.5])
    if comparegroup==3; xlim([0.5 4.5]); end
    ylim([0.1 1.5]) 
    if comparegroup==3;ylim([0.1 1.35]);end
    yticks(0.4:0.3:1.3)
    ylabel('Normalized amplitude')
    set(gca,'FontSize',23)
    set(gcf,'color',[1 1 1])

    clearvars sig
    gcomb=[1 2 3 4];% combinations of groups to compare
    i=1;
    % plot significance of anovan over time
    Ysig=1.45;Yincre=-0.06;textincre=0.003;% comparing first 3 group
    if comparegroup==3;Ysig=1.3;Yincre=-0.06;textincre=0.003;end% comparing first 3 group
    ccom=c(c(:,1)==1|c(:,1)==3,:);
    if comparegroup==3; ccom=c((c(:,1)==1 & c(:,2)==2) | (c(:,1)==1 & c(:,2)==3),:); end
    sigidx=find(ccom(:,end)<0.05 & ccom(:,end)>=0.01)';
    if sigidx
        for sn=sigidx
            plot(ccom(sn,1:2),[Ysig Ysig],'color',co(cnum(gcomb(i,g)),:),'linewidth',1.5)
            text(mean(ccom(sn,1:2)),Ysig+textincre,'*','FontSize',22,...
                'HorizontalAlignment','center',...
                'color',co(cnum(gcomb(i,g)),:));
            Ysig=Ysig+Yincre;
        end
    end
    sigidx=find(ccom(:,end)<0.01 & ccom(:,end)>=0.001)';
    if sigidx
        for sn=sigidx
            plot(ccom(sn,1:2),[Ysig Ysig],'color',co(cnum(gcomb(i,g)),:),'linewidth',1.5)
            text(mean(ccom(sn,1:2)),Ysig+textincre,'**','FontSize',22,...
                'HorizontalAlignment','center',...
                'color',co(cnum(gcomb(i,g)),:));
            Ysig=Ysig+Yincre;
        end
    end
    sigidx=find(ccom(:,end)<0.001)';
    if sigidx
        for sn=sigidx
            plot(ccom(sn,1:2),[Ysig Ysig],'color',co(cnum(gcomb(i,g)),:),'linewidth',1.5)
            text(mean(ccom(sn,1:2)),Ysig+textincre,'***','FontSize',22,...
                'HorizontalAlignment','center',...
                'color',co(cnum(gcomb(i,g)),:));
            Ysig=Ysig+Yincre;
        end
    end
end

%% Section 5: plot comparison of peak amplitude over time (Fig. 6c)


[stat,Tpost]=chi_square_for_two_categories_test(PlastPopu([1,4,2,3],:))
figure
b=bar(fliplr(PlastPopu([1,4,2,3],:))./sum(PlastPopu([1,4,2,3],:),2).*100,'stacked','LineStyle','none');hold on;
b(1).FaceColor = [0.20 0.60 0.80];   % RGB 0–1
b(2).FaceColor = [249 191 69]./255;
b(3).FaceColor = [208 16 76]./255;


%significance
Ysig=102;Yincre=7.5;sigincre=1;

sigidx=find(stat(:,end)<0.05 & stat(:,end)>=0.01)';
if sigidx
    for sn=sigidx
        plot(stat(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(stat(sn,1:2)),Ysig+sigincre,'*','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end
sigidx=find(stat(:,end)<0.01 & stat(:,end)>=0.001)';
if sigidx
    for sn=sigidx
        plot(stat(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(stat(sn,1:2)),Ysig+sigincre,'**','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end
sigidx=find(stat(:,end)<0.001)';
if sigidx
    for sn=sigidx
        plot(stat(sn,1:2),[Ysig Ysig],'color',[0 0 0],'linewidth',3)
        text(mean(stat(sn,1:2)),Ysig+sigincre,'***','FontSize',30,...
            'HorizontalAlignment','center',...
            'color',[0 0 0]);
        Ysig=Ysig+Yincre;
    end
end

xticks(1:4)
xticklabels([]);
xlim([0.5 5.5])
yticks(0:50:100)
ylabel('Percentage (%)')
set(gca,'FontSize',23,'box','off')
set(gcf,'color',[1 1 1])