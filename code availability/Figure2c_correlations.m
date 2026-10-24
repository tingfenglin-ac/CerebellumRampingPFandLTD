% 1. Add the path of the current script to the search path
% 2. Change the current folder to "Figure 2 spontaneous synaptic signals",
% 3. Using ctrl+enter within each function to run each section

%% Section 1: Open files 
% Select all Dendrites for statistic analysis

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

riseLOCS=arrayfun(@(x) data(x).riseLOCS,1:length(data),'uni',0);
upLOCS=arrayfun(@(x) data(x).upLOCS,1:length(data),'uni',0);
pksLOCS=arrayfun(@(x) data(x).pksLOCS,1:length(data),'uni',0);
PKS=arrayfun(@(x) data(x).PKS,1:length(data),'uni',0);
time=arrayfun(@(x) data(x).time,1:length(data),'uni',0);
smoothBC_signal=arrayfun(@(x) data(x).smoothBC_signal,1:length(data),'uni',0);

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


%% Section 2: pull out within and between groups
[within, between] = pearson_groups(data);

% --- Collect ALL within-structure r (upper triangles only, no self-corr) ---
S = [within{:}];              % make a 1×S struct array from the cell array
T_within = vertcat(S.pairs);  % vertical concat of all 'pairs' tables
r_within_all = abs(T_within.r);    % column vector of r

% --- Collect ALL between-structure r ---
r_between_all = abs(between.r);    % column vector of r

% 2 sample ttest
[h,p] = ttest2(r_within_all,r_between_all)
sem_within=nanstd(r_within_all,0,1)./sqrt(length(r_within_all));
sem_between=nanstd(r_between_all,0,1)./sqrt(length(r_between_all));

figure;
hold on;

MarkerSize=10;
meanSize=10;
space=0.19;
jrange=0.25;

%within group
errorbar(1+space,mean(r_within_all),sem_within,'o-',...
    'color',[.5 .5 .5],...
    'CapSize',30,...
    'MarkerSize',meanSize,...
    'LineWidth',3)
text(1+space,mean(r_within_all)+sem_within+0.1,num2str(mean(r_within_all),'%4.2f'),...
    'color',[.5 .5 .5],...
    'HorizontalAlignment','center',...
    'FontSize',15);
scatter(1-space-0.5*jrange+jrange*rand(1,length(r_within_all)),r_within_all,MarkerSize,...
    'MarkerFaceColor','flat',...
    'MarkerEdgeColor','flat',...
    'MarkerFaceAlpha',0.2,...
    'MarkerEdgeAlpha',0,...
    'CData',[0 0 0]);

%between group
errorbar(2+space,mean(r_between_all),sem_between,'o-',...
    'color',[.5 .5 .5],...
    'CapSize',30,...
    'MarkerSize',meanSize,...
    'LineWidth',3)
text(2+space,mean(r_between_all)+sem_between+0.1,num2str(mean(r_between_all),'%4.2f'),...
    'color',[.5 .5 .5],...
    'HorizontalAlignment','center',...
    'FontSize',15);
scatter(2-space-0.5*jrange+jrange*rand(1,length(r_between_all)),r_between_all,MarkerSize,...
    'MarkerFaceColor','flat',...
    'MarkerEdgeColor','flat',...
    'MarkerFaceAlpha',0.2,...
    'MarkerEdgeAlpha',0,...
    'CData',[0 0 0]);

%significance
sigY=0.04;
formatSpec = '%.3f';
FS=16;
Ysig=1;Yincre=0.25;sigincre=0.02;

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

xticks(1:4)
xticklabels({'Within','Between'});
xlim([0.5 2.5])
yticks(0:0.5:1)
ylim([-0 1.1])
ylabel('|r|')
set(gca,'FontSize',23)
set(gcf,'color',[1 1 1])


%%
function [within, between] = pearson_groups(data)
% PEARSON_GROUPS
%   Computes Pearson correlations among time-series in data(k).smoothBC_signal.
%   Uses columns 2:end in each structure.
%   Returns:
%     within  - cell array (1 x numStructs). within{k} has:
%                  .rMat  : Ck x Ck correlation matrix (Ck = size(smoothBC_signal,2)-1)
%                  .pairs : table with columns (i,j,r) for upper-triangular pairs
%     between - table of cross-structure pairs:
%                  columns: s1,c1,s2,c2,r
%
% Notes:
%   - If time grids differ, signals are linearly interpolated onto the union grid.
%   - NaNs are handled pairwise (overlapping valid samples only).
%
% Example:
%   [within, between] = pearson_groups(data);

nS = numel(data);

% Extract time and the matrix of interest (cols 2:end)
T  = cell(1,nS);
X  = cell(1,nS);
C  = zeros(1,nS);
for k = 1:nS
    if ~isfield(data(k),'smoothBC_signal')
        error('data(%d) lacks field smoothBC_signal.',k);
    end
    if ~isfield(data(k),'time')
        error('data(%d) lacks field time.',k);
    end
    T{k} = data(k).time(:);
    Xi   = data(k).smoothBC_signal;
    if size(Xi,2) < 2
        error('data(%d).smoothBC_signal has < 2 columns; need at least 2.',k);
    end
    X{k} = Xi(:,2:end);             % <-- use columns 2:end
    C(k) = size(X{k},2);            % number of series used in structure k
end

%% -------- Within-structure correlations --------
within = cell(1,nS);
for k = 1:nS
    Xi = X{k};
    % z-score columns (ignores NaNs; constant columns -> NaNs)
    Xi = zscore_cols(X{k});
    % Corr within a structure: use corr with pairwise rows (handles NaNs)
    rMat = corr(Xi, Xi, 'Rows','pairwise');  % Ck x Ck
    % Upper-triangular pairs
    [I,J] = find(triu(true(size(rMat)),1));
    r     = arrayfun(@(a,b) rMat(a,b), I, J);
    within{k} = struct( ...
        'rMat', rMat, ...
        'pairs', table(I,J,r,'VariableNames',{'i','j','r'}) );
end

%% -------- Between-structure correlations --------
% Pre-allocate list (size unknown; collect and convert to table)
btw_s1 = [];
btw_c1 = [];
btw_s2 = [];
btw_c2 = [];
btw_r  = [];

for s1 = 1:nS
    for s2 = (s1+1):nS
        % Try fast path if time grids match exactly
        sameGrid = numel(T{s1})==numel(T{s2}) && isequal(T{s1}, T{s2});

        if sameGrid
            % Normalize columns (z-score with NaN handling)
            X1 = zscore_cols(X{s1});
            X2 = zscore_cols(X{s2});
            % Pairwise corr using only overlapping finite rows
            % If both have no NaNs, we can do full matrix multiply:
            if all(all(isfinite(X1))) && all(all(isfinite(X2)))
                n = size(X1,1);
                R = (X1' * X2) / (n-1);  % C1 x C2
            else
                % Fall back to pairwise nan-safe loop for columns
                R = nan(size(X1,2), size(X2,2));
                for c1 = 1:size(X1,2)
                    for c2 = 1:size(X2,2)
                        R(c1,c2) = paircorr_nan(X1(:,c1), X2(:,c2));
                    end
                end
            end
        else
            % Build a common grid = union of unique time points
            tUnion = unique([T{s1}; T{s2}]);
            % Interpolate both onto the union grid
            X1i = interp1(T{s1}, X{s1}, tUnion, 'linear', NaN);
            X2i = interp1(T{s2}, X{s2}, tUnion, 'linear', NaN);
            % z-score columns (after interpolation)
            X1 = zscore_cols(X1i);
            X2 = zscore_cols(X2i);
            % Compute pairwise corr with NaN-overlap
            R = nan(size(X1,2), size(X2,2));
            for c1 = 1:size(X1,2)
                for c2 = 1:size(X2,2)
                    R(c1,c2) = paircorr_nan(X1(:,c1), X2(:,c2));
                end
            end
        end

        % Append results to long lists
        [c1_idx, c2_idx] = ndgrid(1:size(R,1), 1:size(R,2));
        btw_s1 = [btw_s1; repmat(s1, numel(R), 1)];
        btw_c1 = [btw_c1; c1_idx(:)];
        btw_s2 = [btw_s2; repmat(s2, numel(R), 1)];
        btw_c2 = [btw_c2; c2_idx(:)];
        btw_r  = [btw_r;  R(:)];
    end
end

between = table(btw_s1,btw_c1,btw_s2,btw_c2,btw_r, ...
    'VariableNames', {'s1','c1','s2','c2','r'});

end % function pearson_groups


%% ---------- Helpers ----------
function Xz = zscore_cols(X)
% Z-score each column of X independently, ignoring NaNs.
mu = mean(X, 1, 'omitnan');
sd = std (X, 0, 'omitnan');
% Avoid divide-by-zero for constant columns: mark their sd as NaN so output becomes NaN
% (corr with constant series is undefined and will be skipped pairwise).
sd(sd==0) = NaN;
Xz = (X - mu) ./ sd;
end

function r = paircorr_nan(x,y)
% Pearson correlation using only overlapping finite samples
mask = isfinite(x) & isfinite(y);
n = sum(mask);
if n < 2
    r = NaN;
    return
end
xm = x(mask); ym = y(mask);
r = corr(xm,ym, 'Rows','pairwise');
r=r(2);
end
