function [stat,Tpost]=chi_square_for_two_categories_test(counts)
% Supports 4x3 counts: rows = conditions 1..4, columns = [Cat1 Cat2 Cat3]
% Omnibus test uses log-linear Poisson GLMs (no weights, no data expansion).
% Requires Statistics and Machine Learning Toolbox (glmfit).

%% ----- Basic checks -----
if ~ismatrix(counts) || size(counts,1)~=4 || size(counts,2)~=3
    error('counts must be a 4x3 matrix: 4 conditions (rows) x 3 categories (cols).');
end
if any(counts(:)<0) || any(mod(counts(:),1)~=0)
    error('counts must contain nonnegative integers.');
end
rowTot = sum(counts,2);
colTot = sum(counts,1);
if any(rowTot==0) || any(colTot==0)
    error('All row and column totals must be > 0 for chi-square computations.');
end

%% ----- Chi-square test of independence (4x3) -----
N      = sum(rowTot);
E      = (rowTot*colTot)/N;                    % expected counts
chi2   = sum((counts - E).^2 ./ E, 'all');
df_chi = (size(counts,1)-1)*(size(counts,2)-1); % (4-1)*(3-1)=6
p_chi  = 1 - chi2cdf(chi2, df_chi);
fprintf('Chi-square test (4x3): chi2(%d)=%.3f, p=%.4g\n', df_chi, chi2, p_chi);

%% ----- Log-linear Poisson models for omnibus LRT (Condition effect) -----
% Vectorize table to 12 rows (cells), build design matrices.
y = counts(:);                                % 12x1 cell counts
[iIdx, jIdx] = ndgrid(1:4, 1:3);              % cell coordinates
iIdx = iIdx(:);  jIdx = jIdx(:);

% Condition (4 levels) -> 3 dummy cols; Category (3 levels) -> 2 dummy cols
C = dummyvar(categorical(iIdx,1:4));  C = C(:,2:4);   % 12x3
K = dummyvar(categorical(jIdx,1:3));  K = K(:,2:3);   % 12x2

% Independence model: log(mu_ij) = intercept + alpha_i + beta_j
X_ind = [C K];                                % glmfit adds an intercept

% Saturated model adds interaction terms ( (4-1)*(3-1)=6 columns )
inter = zeros(size(y,1), size(C,2)*size(K,2));
t = 1;
for a = 1:size(C,2)
    for b = 1:size(K,2)
        inter(:,t) = C(:,a).*K(:,b);
        t = t + 1;
    end
end
X_full = [X_ind inter];

% Fit Poisson GLMs with log link
[~, dev_ind, stats_ind]  = glmfit(X_ind,  y, 'poisson', 'link','log'); %#ok<ASGLU>
[~, dev_full, stats_full]= glmfit(X_full, y, 'poisson', 'link','log'); %#ok<ASGLU>

% Likelihood-ratio test: difference in deviance ~ chi2 with df = added params
LRT  = dev_ind - dev_full;
dfL  = stats_ind.dfe - stats_full.dfe;        % should be 6
p_LRT = 1 - chi2cdf(LRT, dfL);
fprintf('Log-linear omnibus (Condition x Category): chi2(%d)=%.3f, p=%.4g\n', dfL, LRT, p_LRT);

%% ----- Helpful: per-condition percentages for each category -----
pctCat = 100*counts ./ rowTot;                % 4x3
Tsum = table((1:4)', rowTot, pctCat(:,1), pctCat(:,2), pctCat(:,3), ...
    'VariableNames', {'Condition','N','Pct_Category1','Pct_Category2','Pct_Category3'});
disp(Tsum);

%% ----- Post-hoc pairwise tests across conditions (2x3 chi-square) + Bonferroni -----
pairs = nchoosek(1:4,2);
rawP  = nan(size(pairs,1),1);
chi2_pair = nan(size(pairs,1),1);
for k = 1:size(pairs,1)
    sub = counts(pairs(k,:), :);              % 2x3 contingency for the pair
    rT = sum(sub,2); cT = sum(sub,1); Np = sum(rT);
    Epair = (rT*cT)/Np;
    chi2_pair(k) = sum((sub - Epair).^2 ./ Epair, 'all');
    df_pair = (2-1)*(3-1);                    % = 2
    rawP(k) = 1 - chi2cdf(chi2_pair(k), df_pair);
end
pAdj = min(rawP * size(pairs,1), 1);          % Bonferroni

Tpost = table(pairs(:,1), pairs(:,2), chi2_pair, rawP, pAdj, ...
    'VariableNames', {'Cond_i','Cond_j','chi2','p_raw','p_Bonf'});
disp(Tpost);

%% ----- Return stats: [i j p_raw p_Bonf] (same shape as before) -----
stat = [pairs(:,1), pairs(:,2), rawP, pAdj];
end