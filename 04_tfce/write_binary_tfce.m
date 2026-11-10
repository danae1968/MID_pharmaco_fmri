% WRITE_BINARY_TFCE  Binarise TFCE FWE maps into significant-voxel masks.
%
% For every analysis folder, thresholds each TFCE_log_pFWE_000N.nii at
% FWE p < .05 and writes
%   TFCE_significant_voxels_<contrast name>_fwe_p05.nii
% The output is named after the contrast name read from that folder's
% SPM.mat, not the index, so masks are labelled correctly whatever the
% contrast order in each model.
%
% Loops over the analysis directories in includeDirs (no-covariate models
% plus every covariate model). Run after wrapper_tfce.m has finished.
%
% -------------------------------------------------------------------------
% Danae Papadopetraki, April 2026

clear; clc;

% Paths
% =========================================================================
if ispc
    projectDir = 'P:\3017048.01\';
    if exist('H:\common\matlab\spm12','dir')
        addpath H:\common\matlab\spm12
    else
        addpath 'C:\Program Files\MATLAB\spm12-main'
    end
elseif isunix
    projectDir = '/project/3017048.01/';
    addpath /home/common/matlab/spm12
end

hrf_type = 'canonical_hrf';
groupLevelBase = fullfile(projectDir, 'bids', 'derivatives', 'mri', 'groupLevel', 'MID', 'regular', hrf_type);

% Analysis directories to process. If empty, every covariate_* directory
% under groupLevelBase is found automatically (minus excludeDirs).
includeDirs = {'no_covariate_with_exc', ...
               'covariate_ki_VS_new', 'covariate_ki_wholeCaudate', ...
               'covariate_ki_wholePutamen', 'covariate_ki_wholeStriatum', ...
               'covariate_mphBehEffect_early', 'covariate_mphBehEffect_hit', ...
               'covariate_mphBehEffect_rt', 'covariate_mphBehEffect_late', ...
               'covariate_rtEffect_MPHvsSUL', 'covariate_hitEffect_MPHvsSUL', ...
               'covariate_lateEffect_MPHvsSUL', 'covariate_SULvsPBO_late', ...
               'covariate_TauEffect_MPHvsSUL', 'covariate_tauEffect_MPHvsPBO'};

% Directories to skip in auto-discovery mode
excludeDirs = {};

% Overwrite existing masks? Set true after re-running TFCE; false only
% fills in what is missing.
force = false;

% ---- Analysis directories to process --------------------------------------
topDirs = {};
if ~isempty(includeDirs)
    for i = 1:numel(includeDirs)
        if exist(fullfile(groupLevelBase, includeDirs{i}), 'dir')
            topDirs{end+1} = includeDirs{i}; %#ok<SAGROW>
        else
            warning('includeDirs entry not found, skipping: %s', includeDirs{i});
        end
    end
else
    d = dir(groupLevelBase);
    d = d([d.isdir]);
    for i = 1:numel(d)
        if strcmp(d(i).name, 'no_covariate_with_exc') || startsWith(d(i).name, 'covariate_')
            topDirs{end+1} = d(i).name; %#ok<SAGROW>
        end
    end
    topDirs = setdiff(topDirs, excludeDirs, 'stable');
end
if isempty(topDirs)
    error('No analysis directories found under %s', groupLevelBase);
end
fprintf('Found %d analysis directories to process:\n', numel(topDirs));
fprintf('  %s\n', strjoin(topDirs, '\n  '));

% Build analysis list: drug comparisons AND session averages
% =========================================================================
contrasts = {'Target', 'HCvsLC', 'HitvsMiss', 'Cue'};
drugComps = {'MPHvsPBO', 'SULvsPBO', 'MPHvsSUL'};
drugs     = {'MPH', 'PBO', 'SUL'};

% Also process the single-drug models (MPH_<con>, PBO_<con>, SUL_<con>)
% from wrapper_tfce.m's 'singleDrug' batches. Keep in step with
% include_single_drug in wrapper_extract_clusters_MID.m.
include_single_drug = false;

analyses = {};
for iC = 1:numel(contrasts)
    for iD = 1:numel(drugComps)
        analyses{end+1} = sprintf('%s_%s', drugComps{iD}, contrasts{iC}); %#ok<SAGROW>
    end
    analyses{end+1} = sprintf('sessionAverage_%s', contrasts{iC}); %#ok<SAGROW>
    if include_single_drug
        for iD = 1:numel(drugs)
            analyses{end+1} = sprintf('%s_%s', drugs{iD}, contrasts{iC}); %#ok<SAGROW>
        end
    end
end

thr = -log10(0.05);   % FWE p < .05 on the -log10(p) map

% Loop over every top-level analysis dir, then every analysis folder in it
% =========================================================================
for iT = 1:numel(topDirs)
    mainOutputDir = fullfile(groupLevelBase, topDirs{iT});
    fprintf('\n=== %s ===\n', topDirs{iT});

    for iA = 1:numel(analyses)
        dir_ = fullfile(mainOutputDir, analyses{iA});
        if ~exist(dir_, 'dir')
            warning('Missing folder %s, skipping.', analyses{iA});
            continue;
        end

        % Find whichever TFCE FWE maps were actually written here
        pFiles = dir(fullfile(dir_, 'TFCE_log_pFWE_*.nii'));
        if isempty(pFiles)
            warning('No TFCE FWE maps in %s, skipping.', analyses{iA});
            continue;
        end

        % Contrast names for this model (to label outputs by name, not index)
        S = load(fullfile(dir_, 'SPM.mat'), 'SPM');
        conNames = {S.SPM.xCon.name};

        for iF = 1:numel(pFiles)
            pName = pFiles(iF).name;
            cIx   = sscanf(pName, 'TFCE_log_pFWE_%04d.nii');

            % Label from the real contrast name (spaces -> underscores)
            label   = regexprep(conNames{cIx}, '\s+', '_');
            outName = sprintf('TFCE_significant_voxels_%s_fwe_p05.nii', label);
            outPath = fullfile(dir_, outName);

            if ~force && exist(outPath, 'file')
                fprintf('%s | %s | con %d (%s): %s already exists, skipping (force=false)\n', ...
                    topDirs{iT}, analyses{iA}, cIx, conNames{cIx}, outName);
                continue;
            end

            Vp    = spm_vol(fullfile(dir_, pName));
            p_img = spm_read_vols(Vp);
            sig_mask = double(p_img >= thr);

            Vout       = Vp;
            Vout.fname = outPath;
            Vout.dt    = [spm_type('uint8') 0];   % binary mask
            spm_write_vol(Vout, sig_mask);

            fprintf('%s | %s | con %d (%s): %d significant voxels -> %s\n', ...
                topDirs{iT}, analyses{iA}, cIx, conNames{cIx}, nnz(sig_mask), outName);
        end
    end
end

fprintf('\nDone. Wrote TFCE significant-voxel masks for all analyses in %d directories.\n', ...
    numel(topDirs));
