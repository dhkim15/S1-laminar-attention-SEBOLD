% -------------------------------------------------------------------------
% Laminar qT1 Profile and Derivative Analysis (CSF → WM)
%
% This script:
%   1) Loads subject-wise laminar profiles (qT1, dqT1, ddqT1) across 8 bins.
%   2) Flips LAYNII-order vectors (often WM→CSF) to CSF→WM for presentation.
%   3) Plots group mean ± SEM profiles (qT1, dqT1, ddqT1).
%   4) Estimates subject-wise peak location of the second derivative (ddqT1)
%      and performs statistical testing of peak depth.
%   5) Visualizes one-sample t-test p-values comparing peak locations to each
%      bin-center depth (log scale), and marks Bin 4 (depth 0.5625).
%
% Input:
%   ddqT1.mat containing variable "data" with size (nSubs x 3) cell array:
%       data{s,1} = qT1 profile (1x8)
%       data{s,2} = dqT1 profile (1x8)  (normalized)
%       data{s,3} = ddqT1 profile (1x8) (normalized)
%
% Output (saved to ./outputs):
%   - laminar_profiles.png
%   - peak_pvalues.png
%   - stats_log.txt
%
% Dependencies (optional):
%   - shadedErrorBar.m (if not found, script falls back to errorbar)
% -------------------------------------------------------------------------

clearvars;
close all;
clc;

%% ===============================
% User options
%% ===============================
matFile   = 'ddqT1.mat';
outdir    = fullfile(pwd, 'outputs');
ref_depth = 0.6;          % reference normalized depth for one-sample t-test
peak_mode = 'max';        % 'max' | 'min' | 'absmax' (depends on ddqT1 definition)

% Output options
save_png_resolution = 300;
show_bin_tests      = false;   % set true to print p-values for all bin-center depths

%% ===============================
% Setup outputs
%% ===============================
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

logFile = fullfile(outdir, 'stats_log.txt');
if exist(logFile, 'file')
    delete(logFile);
end
diary(logFile);

fprintf('Laminar qT1 analysis\n');
fprintf('MAT file      : %s\n', matFile);
fprintf('Output dir    : %s\n', outdir);
fprintf('ref_depth     : %.3f\n', ref_depth);
fprintf('peak_mode     : %s\n', peak_mode);
fprintf('show_bin_tests: %d\n\n', show_bin_tests);

%% ===============================
% Load and validate input
%% ===============================
S = load(matFile, 'data');
assert(isfield(S, 'data'), 'Missing variable "data" in %s', matFile);
data = S.data;

assert(iscell(data) && size(data,2) == 3, 'data must be an (nSubs x 3) cell array');

nSubs = size(data,1);

% LAYNII bin centers (often expressed WM→CSF). You used these previously.
% WM→CSF centers:
d_wm2csf = [0.0625, 0.1875, 0.3125, 0.4375, 0.5625, 0.6875, 0.8125, 0.9375];
nDepths  = numel(d_wm2csf);

% Presentation axis: CSF→WM
x = 1:nDepths;
layer_labels = compose('%d', x);
d_csf2wm = fliplr(d_wm2csf);   % 0.9375 near CSF → 0.0625 near WM

% Optional dependency
use_shadedErrorBar = exist('shadedErrorBar','file') == 2;
if use_shadedErrorBar
    fprintf('Dependency: shadedErrorBar found (using shaded plots).\n\n');
else
    fprintf('Dependency: shadedErrorBar NOT found (fallback to errorbar).\n\n');
end

%% ===============================
% Extract matrices from cell array
%% ===============================
qT1_all   = nan(nSubs, nDepths);
dqT1_all  = nan(nSubs, nDepths);
ddqT1_all = nan(nSubs, nDepths);

for s = 1:nSubs
    q   = data{s,1};
    dq  = data{s,2};
    ddq = data{s,3};

    assert(isvector(q)   && numel(q)   == nDepths, 'data{%d,1} must be 1x%d', s, nDepths);
    assert(isvector(dq)  && numel(dq)  == nDepths, 'data{%d,2} must be 1x%d', s, nDepths);
    assert(isvector(ddq) && numel(ddq) == nDepths, 'data{%d,3} must be 1x%d', s, nDepths);

    qT1_all(s,:)   = q(:)';
    dqT1_all(s,:)  = dq(:)';
    ddqT1_all(s,:) = ddq(:)';
end

%% ===============================
% Flip ONLY for presentation/stats: WM→CSF --> CSF→WM
%% ===============================
qT1_csf2wm   = fliplr(qT1_all);
dqT1_csf2wm  = fliplr(dqT1_all);
ddqT1_csf2wm = fliplr(ddqT1_all);

%% ===============================
% Compute mean and SEM (bin-wise effective N)
%% ===============================
% qT1
mean_qT1 = mean(qT1_csf2wm, 1, 'omitnan');
sd_qT1   = std(qT1_csf2wm, 0, 1, 'omitnan');
nEff_qT1 = sum(~isnan(qT1_csf2wm), 1);
sem_qT1  = sd_qT1 ./ sqrt(max(nEff_qT1, 1));

% dqT1
mean_dqT1 = mean(dqT1_csf2wm, 1, 'omitnan');
sd_dqT1   = std(dqT1_csf2wm, 0, 1, 'omitnan');
nEff_dqT1 = sum(~isnan(dqT1_csf2wm), 1);
sem_dqT1  = sd_dqT1 ./ sqrt(max(nEff_dqT1, 1));

% ddqT1
mean_ddqT1 = mean(ddqT1_csf2wm, 1, 'omitnan');
sd_ddqT1   = std(ddqT1_csf2wm, 0, 1, 'omitnan');
nEff_ddqT1 = sum(~isnan(ddqT1_csf2wm), 1);
sem_ddqT1  = sd_ddqT1 ./ sqrt(max(nEff_ddqT1, 1));

%% ===============================
% Plot: mean profiles (CSF → WM)
%% ===============================
fig1 = figure('Color', 'w', 'Name', 'Laminar Profiles');
set(fig1, 'Position', [100 100 1400 450]);

subplot(1,3,1);
if use_shadedErrorBar
    s1 = shadedErrorBar(x, mean_qT1, sem_qT1);
    s1.mainLine.Color = 'k';
    s1.mainLine.LineWidth = 2;
else
    errorbar(x, mean_qT1, sem_qT1, 'k-o', 'LineWidth', 1.5);
end
xlabel({'Equivolume Bin (CSF \rightarrow WM)', '(1 = near CSF, 8 = near WM)'});
ylabel('qT1 (ms)');
xticks(x);
xticklabels(layer_labels);
grid on;
axis square;

subplot(1,3,2);
if use_shadedErrorBar
    s2 = shadedErrorBar(x, mean_dqT1, sem_dqT1);
    s2.mainLine.Color = 'b';
    s2.mainLine.LineWidth = 2;
else
    errorbar(x, mean_dqT1, sem_dqT1, 'b-o', 'LineWidth', 1.5);
end
xlabel({'Equivolume Bin (CSF \rightarrow WM)', '(1 = near CSF, 8 = near WM)'});
ylabel('Normalized dqT1/dx');
xticks(x);
xticklabels(layer_labels);
grid on;
axis square;

subplot(1,3,3);
if use_shadedErrorBar
    s3 = shadedErrorBar(x, mean_ddqT1, sem_ddqT1);
    s3.mainLine.Color = 'r';
    s3.mainLine.LineWidth = 2;
else
    errorbar(x, mean_ddqT1, sem_ddqT1, 'r-o', 'LineWidth', 1.5);
end
xlabel({'Equivolume Bin (CSF \rightarrow WM)', '(1 = near CSF, 8 = near WM)'});
ylabel('Normalized d^2qT1/dx^2');
xticks(x);
xticklabels(layer_labels);
grid on;
axis square;

set(findall(fig1,'-property','FontSize'),'FontSize',14);
set(findall(fig1,'type','axes'),'XLim',[0.5 8.5]);

exportgraphics(fig1, fullfile(outdir, 'laminar_profiles.png'), ...
    'Resolution', save_png_resolution);

%% ===============================
% Stats: Peak of SECOND derivative (bin + depth)
%% ===============================
peak_bins = nan(nSubs, 1);

for s = 1:nSubs
    v = ddqT1_csf2wm(s,:);

    if all(isnan(v))
        peak_bins(s) = nan;
        continue;
    end

    switch lower(peak_mode)
        case 'max'
            [~, idx] = max(v);
        case 'min'
            [~, idx] = min(v);
        case 'absmax'
            [~, idx] = max(abs(v));
        otherwise
            error('Unknown peak_mode: %s', peak_mode);
    end

    peak_bins(s) = idx;   % 1..8 (1=CSF, 8=WM)
end

peak_locs = d_csf2wm(peak_bins);   % normalized depth coordinate in CSF→WM sense

valid = ~isnan(peak_locs);
nValid = sum(valid);

fprintf('[SECOND DERIVATIVE PEAK LOCATION STATS]\n');
fprintf('Valid N          : %d / %d\n', nValid, nSubs);

if nValid >= 2
    [~, p, ci, stats] = ttest(peak_locs(valid), ref_depth);

    m   = mean(peak_locs(valid), 'omitnan');
    sd  = std(peak_locs(valid),  'omitnan');
    sem = sd / sqrt(nValid);

    fprintf('Mean ± SD        : %.4f ± %.4f\n', m, sd);
    fprintf('Mean ± SEM       : %.4f ± %.4f\n', m, sem);
    fprintf('Reference depth  : %.4f\n', ref_depth);
    fprintf('t(%d)=%.3f, p=%.4g\n', stats.df, stats.tstat, p);
    fprintf('95%% CI (mean)    : [%.4f, %.4f]\n\n', ci(1), ci(2));
else
    fprintf('Not enough valid subjects for t-test.\n\n');
end

%% ===============================
% Visualization: p-values vs each bin depth (CSF→WM)
%% ===============================
depths = d_csf2wm;
p_values = nan(size(depths));

if nValid >= 2
    for i = 1:numel(depths)
        [~, p_values(i)] = ttest(peak_locs(valid), depths(i));
    end
end

fig2 = figure('Color', 'w', 'Name', 'Peak p-values vs Bin Depth');
set(fig2, 'Position', [200 200 700 500]);

plot(x, p_values, 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'r');
hold on;
yline(0.05, '--b', 'p = 0.05', 'LineWidth', 1.5);

% Mark Bin 4 explicitly (as stated in caption): Bin 4 center depth = 0.5625
xline(4, '--g', 'LineWidth', 1.5);

xlabel('Equivolume Bin (1 = near CSF, 8 = near WM)');
ylabel('p-value (one-sample t-test vs bin-center depth)');
set(gca, 'YScale', 'log');
xticks(x);
xticklabels(layer_labels);
grid on;
xlim([0.5 8.5]);

legend({'p-values', 'p = 0.05 threshold', 'Bin 4 (0.5625)'}, ...
    'Location', 'southwest');

set(findall(fig2,'-property','FontSize'),'FontSize',14);

exportgraphics(fig2, fullfile(outdir, 'peak_pvalues.png'), ...
    'Resolution', save_png_resolution);

%% ===============================
% Optional: print all bin-center tests
%% ===============================
if show_bin_tests
    fprintf('[P-VALUES VS EACH BIN-CENTER DEPTH]\n');
    for i = 1:numel(depths)
        fprintf('Bin %d (depth = %.4f): p = %.4g\n', i, depths(i), p_values(i));
    end
    fprintf('\n');
end

%% ===============================
% Heuristic diagnostic
% (Not used for the caption marker; retained for reproducibility)
%% ===============================
if all(~isnan(p_values))
    % Identify the bin whose depth is closest to the distribution of
    % individual peak locations (used only for diagnostics).
    [~, idx_min_p] = min(abs(p_values - 0.5));
    fprintf('[Heuristic diagnostic]\n');
    fprintf('Bin closest to p≈0.5: Bin %d (depth=%.4f), p=%.4g\n\n', ...
        idx_min_p, depths(idx_min_p), p_values(idx_min_p));
end

%% ===============================
% Done
%% ===============================
fprintf('Saved figures to: %s\n', outdir);
diary off;