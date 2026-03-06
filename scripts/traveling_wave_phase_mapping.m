% traveling_wave_phase_mapping.m
%
% Public-friendly version of the original hybrid traveling-wave analysis.
%
% This script:
%   1. loads concatenated 4D fMRI data
%   2. trims rest periods at the beginning and end of each run
%   3. identifies tuned voxels using sinusoidal regression (R^2 threshold)
%   4. estimates phase and amplitude at the stimulus frequency using FFT
%   5. flips phase for reverse runs
%   6. averages phase across valid runs using a circular mean
%   7. saves phase, R^2, amplitude, and thresholded phase maps
%
% The analysis logic is kept as close as possible to the original script.
%
% Requirements:
%   - MATLAB with NIfTI support (niftiread, niftiinfo, niftiwrite)
%
% Inputs to edit:
%   - nii_file
%   - optional mask_file
%
% Outputs:
%   - phase_map.nii.gz
%   - r2_map.nii.gz
%   - amplitude_map.nii.gz
%   - phase_masked.nii.gz
%   - r2_mask.nii.gz

clc;
clear;
close all;

%% =========================
% USER SETTINGS
%% =========================

nii_file = 'all_runs_mc_concat.nii.gz';   % concatenated 4D BOLD data
use_mask = false;                         % set true to use an external mask
mask_file = 'gray_mask.nii.gz';           % optional binary mask

TR = 2;                    % repetition time (s)
stimCycle = 40;            % one stimulus cycle duration (s)
stimFreq = 1 / stimCycle;  % stimulus frequency (Hz)

r2_thresh = 0.04;          % threshold for tuned voxels during run-wise selection

nRuns = 6;                 % total concatenated runs
TRs_per_run = 210;         % includes rest at start and end
trim_TRs = 5;              % trim 10 s from each end when TR = 2 s

% Optional debug voxel
enable_debug = false;
debug_x = 73;
debug_y = 84;
debug_z = 15;

% Threshold used for the saved masked phase map
r2_mask_thresh = 0.06;

%% =========================
% LOAD DATA
%% =========================

bold_nii = double(niftiread(nii_file));
[nx, ny, nz, nt] = size(bold_nii);

assert(nt == nRuns * TRs_per_run, ...
    'Mismatch between data time dimension and run settings.');

if use_mask
    assert(isfile(mask_file), 'Mask file not found: %s', mask_file);
    mask_nii = logical(niftiread(mask_file));
    assert(isequal(size(mask_nii), [nx ny nz]), ...
        'Mask dimensions must match the fMRI spatial dimensions.');
else
    mask_nii = true(nx, ny, nz);
end

%% =========================
% MODEL SETUP
%% =========================

valid_TRs = TRs_per_run - 2 * trim_TRs;
t_run = (0:valid_TRs - 1) * TR;

model_sin = sin(2 * pi * stimFreq * t_run);
model_cos = cos(2 * pi * stimFreq * t_run);
X = zscore([model_sin' model_cos']);

freqs = (0:valid_TRs - 1) / (valid_TRs * TR);
[~, stimIdx] = min(abs(freqs - stimFreq));

%% =========================
% PREALLOCATE OUTPUTS
%% =========================

amplitude = nan(nx, ny, nz);
phase_deg = nan(nx, ny, nz);
r2 = nan(nx, ny, nz);

fprintf('Running hybrid regression + FFT analysis on %d runs...\n', nRuns);

%% =========================
% VOXEL LOOP
%% =========================

for x = 1:nx
    for y = 1:ny
        for z = 1:nz

            if ~mask_nii(x, y, z)
                continue;
            end

            r2_vals = [];
            amps = [];
            phases = [];

            for r = 1:nRuns

                run_start = (r - 1) * TRs_per_run + 1 + trim_TRs;
                run_end   = r * TRs_per_run - trim_TRs;

                ts = squeeze(bold_nii(x, y, z, run_start:run_end));

                if any(isnan(ts)) || std(ts) == 0
                    continue;
                end

                ts = detrend(ts);

                %% 1) sinusoidal regression
                ts_z = zscore(ts);
                b = X \ ts_z;
                y_fit = X * b;
                r2_val = corr(ts_z, y_fit)^2;

                if r2_val > r2_thresh

                    %% 2) FFT at stimulus frequency
                    Y = fft(ts_z);
                    phase_rad = angle(Y(stimIdx));
                    amp_val = abs(Y(stimIdx)) / valid_TRs;

                    %% flip phase for reverse runs (2, 4, 6)
                    is_reverse = mod(r, 2) == 0;
                    if is_reverse
                        phase_rad = mod(-phase_rad, 2*pi);
                    end

                    %% optional debug print
                    if enable_debug && x == debug_x && y == debug_y && z == debug_z
                        if is_reverse
                            label = 'REV';
                        else
                            label = 'FWD';
                        end
                        fprintf('Run %d [%s] | Phase: %.1f° | R^2: %.3f\n', ...
                            r, label, rad2deg(phase_rad), r2_val);
                    end

                    phases(end+1) = phase_rad; %#ok<AGROW>
                    r2_vals(end+1) = r2_val;   %#ok<AGROW>
                    amps(end+1) = amp_val;     %#ok<AGROW>
                end
            end

            %% save voxelwise averages across valid runs
            if ~isempty(r2_vals)
                r2(x, y, z) = mean(r2_vals);
                amplitude(x, y, z) = mean(amps);

                phase_mean_rad = atan2(mean(sin(phases)), mean(cos(phases)));
                if phase_mean_rad < 0
                    phase_mean_rad = phase_mean_rad + 2*pi;
                end
                phase_deg(x, y, z) = rad2deg(phase_mean_rad);
            end
        end
    end
end

%% =========================
% PREPARE NIFTI HEADER
%% =========================

info = niftiinfo(nii_file);
info.ImageSize = [nx, ny, nz];
info.PixelDimensions = info.PixelDimensions(1:3);
info.Raw.dim = [3, nx, ny, nz, 1, 1, 1, 1];
info.Datatype = 'single';

%% =========================
% SAVE MAIN OUTPUTS
%% =========================

niftiwrite(single(phase_deg), 'phase_map.nii.gz', info, 'Compressed', true);
niftiwrite(single(r2), 'r2_map.nii.gz', info, 'Compressed', true);
niftiwrite(single(amplitude), 'amplitude_map.nii.gz', info, 'Compressed', true);

fprintf('\nHybrid analysis complete.\n');
fprintf('Saved: phase_map.nii.gz, r2_map.nii.gz, amplitude_map.nii.gz\n');

%% =========================
% SAVE MASKED PHASE MAP
%% =========================

r2_mask = r2 > r2_mask_thresh;

phase_masked = phase_deg;
phase_masked(~r2_mask) = NaN;

niftiwrite(single(phase_masked), 'phase_masked.nii.gz', info, 'Compressed', true);
niftiwrite(single(r2_mask), 'r2_mask.nii.gz', info, 'Compressed', true);

fprintf('Saved: phase_masked.nii.gz, r2_mask.nii.gz\n');
