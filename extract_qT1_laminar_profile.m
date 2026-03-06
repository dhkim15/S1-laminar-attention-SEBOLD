% extract_qT1_laminar_profile.m
%
% Extract a laminar qT1 profile from a discrete LAYNII layer map (1..8)
% within an ROI, while preserving the original analysis logic:
%   - optional ROI dilation
%   - layer-wise extraction using (layer == i)
%   - first derivative normalization
%   - second derivative computation
%
% This public-friendly version:
%   - removes hard-coded local paths
%   - saves outputs to ./out/
%   - keeps the original ROI dilation + layer==i logic unchanged
%
% Requirements:
%   - MATLAB
%   - Image Processing Toolbox (for imdilate)
%
% Important:
%   qT1, layer map, and ROI must be in the same voxel space and dimensions.
%
% Inputs to edit below:
%   - qT1_path
%   - layer_path
%   - roi_path
%   - subject_id
%
% Outputs:
%   ./out/<subject_id>_qT1profile.mat
%   ./out/<subject_id>_qT1profile.png

clc;
clear;
close all;

%% =========================
% USER SETTINGS
%% =========================
subject_id = 'CYH_sb01'; % e.g., 'sub01'

% Input files (edit these paths)
qT1_path  = '/path/to/qT1_map_new.nii';
layer_path = '/path/to/rim_mask_layers_equidist.nii.gz';   % discrete LAYNII labels: 1..8
roi_path   = '/path/to/rh.S1_3b_in_qT1.nii.gz';

% Number of cortical bins/layers
n_layers = 8;

% ROI dilation
% Keep this identical to the original logic if you want to reproduce results.
do_dilate = true;
dilate_kernel = ones(3,3,3);

% Flip control
% In the original code, a double flip effectively resulted in no net flip.
% To explicitly reproduce the original behavior, keep this false.
% Set true only if you intentionally want reversed profile ordering.
flip_profile = false;

% Output directory
out_dir = fullfile(pwd, 'out');

%% =========================
% CHECK INPUTS
%% =========================
assert(isfile(qT1_path),  'qT1 file not found: %s', qT1_path);
assert(isfile(layer_path), 'Layer file not found: %s', layer_path);
assert(isfile(roi_path),   'ROI file not found: %s', roi_path);

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% =========================
% LOAD DATA
%% =========================
qT1 = double(niftiread(qT1_path));
layer_map = double(niftiread(layer_path));
roi = niftiread(roi_path) ~= 0;

% Check dimensions
assert(isequal(size(qT1), size(layer_map), size(roi)), ...
    'qT1, layer map, and ROI must have identical dimensions.');

% For discrete LAYNII layer maps, robustly enforce integer labels
layer_lbl = round(layer_map);

%% =========================
% ROI DILATION
%% =========================
if do_dilate
    roi_dilated = imdilate(roi, dilate_kernel);
else
    roi_dilated = roi;
end

% Restrict analysis to valid layers inside the dilated ROI
layered_roi = (layer_lbl > 0) & roi_dilated;

%% =========================
% EXTRACT qT1 PROFILE
%% =========================
qT1_profile = nan(n_layers, 1);
nvox_layer  = zeros(n_layers, 1);

for i = 1:n_layers
    mask = (layer_lbl == i) & layered_roi;
    nvox_layer(i) = nnz(mask);

    if nvox_layer(i) > 0
        qT1_profile(i) = mean(qT1(mask), 'omitnan');
    end
end

% Optional explicit flip
if flip_profile
    qT1_profile = flipud(qT1_profile);
    nvox_layer  = flipud(nvox_layer);
end

%% =========================
% DERIVATIVES
%% =========================
% Match original style:
%   1) first derivative from qT1_profile
%   2) normalize by max absolute value
%   3) second derivative from normalized first derivative
dqT1 = gradient(qT1_profile);

m = max(abs(dqT1));
if isfinite(m) && m > 0
    dqT1_normalized = dqT1 / m;
else
    dqT1_normalized = dqT1;
end

ddqT1 = gradient(dqT1_normalized);

% Bin-center depth coordinates
depth_centers = linspace(0.5/n_layers, 1 - 0.5/n_layers, n_layers);

%% =========================
% PLOT
%% =========================
fig = figure('Color', 'w', 'Position', [200 200 650 750]);

subplot(3,1,1);
plot(depth_centers, qT1_profile, 'k-o', 'LineWidth', 2);
title('qT1 Profile', 'Interpreter', 'none');
xlabel('Normalized depth (bin centers)');
ylabel('qT1');
grid on;

subplot(3,1,2);
plot(depth_centers, dqT1_normalized, 'b-o', 'LineWidth', 2);
title('First Derivative (normalized)', 'Interpreter', 'none');
xlabel('Normalized depth (bin centers)');
ylabel('dqT1 (normalized)');
grid on;

subplot(3,1,3);
plot(depth_centers, ddqT1, 'r-o', 'LineWidth', 2);
title('Second Derivative', 'Interpreter', 'none');
xlabel('Normalized depth (bin centers)');
ylabel('ddqT1');
grid on;

sgtitle(sprintf('%s | voxels/layer: %s', ...
    subject_id, strjoin(string(nvox_layer'), ' ')), ...
    'Interpreter', 'none');

png_path = fullfile(out_dir, sprintf('%s_qT1profile.png', subject_id));
saveas(fig, png_path);

%% =========================
% SAVE OUTPUT
%% =========================
out = struct();
out.subject_id = subject_id;
out.qT1_path = qT1_path;
out.layer_path = layer_path;
out.roi_path = roi_path;
out.n_layers = n_layers;
out.do_dilate = do_dilate;
out.dilate_kernel = dilate_kernel;
out.flip_profile = flip_profile;

out.nvox_layer = nvox_layer(:);
out.depth_centers = depth_centers(:);

out.qT1_profile = qT1_profile(:);
out.dqT1 = dqT1(:);
out.dqT1_normalized = dqT1_normalized(:);
out.ddqT1 = ddqT1(:);

out.figure_png = png_path;

mat_path = fullfile(out_dir, sprintf('%s_qT1profile.mat', subject_id));
save(mat_path, 'out', '-v7.3');

fprintf('\nSaved outputs:\n');
fprintf('  %s\n', mat_path);
fprintf('  %s\n\n', png_path);