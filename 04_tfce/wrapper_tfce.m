% WRAPPER_TFCE  Threshold-free cluster enhancement for the group models.
%
% Runs permutation-based TFCE (TFCE toolbox, Gaser) on the group-level
% models estimated in 03_group_level. What to run is defined as a list of
% BATCHES below -- one row per analysis family (drug comparisons, session
% average or single drugs; with or without a covariate) -- and every
% enabled row runs in one call.
%
% Contrast indices are resolved BY NAME in each model's SPM.mat (see
% run_tfce.m), because their order is not guaranteed to be identical
% across covariate models.
%
% Writes TFCE_log_pFWE_<idx>.nii (-log10 FWE-corrected p) into each
% analysis folder. Run write_binary_tfce.m afterwards.
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

addpath(fullfile(projectDir, 'code', 'analysis', 'mri', '2ndLevel', 'MID'));
spm('fmri','defaults');

% =========================================================================
% USER CONFIG
% =========================================================================
hrf_type    = 'canonical_hrf';                              % 'canonical_hrf' | 'fir'
contrasts   = {'HCvsLC', 'HitvsMiss', 'Cue', 'Target'};      % task contrasts
drugComps   = {'MPHvsPBO', 'SULvsPBO', 'MPHvsSUL'};          % drug comparisons
drugs       = {'MPH', 'PBO', 'SUL'};                         % single-drug models

% Families of batches to run. run_tfce.m recomputes everything it is
% given, so switch off families that are already done when adding new ones.
runDrugComparisons = true;    % MPHvsPBO_/SULvsPBO_/MPHvsSUL_<contrast>
runAverages        = true;    % sessionAverage_<contrast>
runSingleDrugs     = true;    % MPH_/PBO_/SUL_<contrast>

% TFCE parameters shared by every batch
baseTfce.nperm    = 5000;     % permutations
baseTfce.E_weight = 0.5;      % TFCE extent weighting (E)
baseTfce.nproc    = 4;        % cores per job
baseTfce.parallel = true;     % submit one Slurm job per analysis

% ---- Batches --------------------------------------------------------------
% type          : 'drugComparisons' | 'average' | 'singleDrug'
% covariateName : '' for no covariate, else the covariate stem used in the
%                 group-level folder and contrast names, e.g. 'ki_VS_new'
% mask          : 'auto'       -> 'none' for the plain session average and
%                                 plain single drugs (whole brain),
%                                 'functional' for everything else
%                 'none' | 'functional' | path to a mask image
% enabled       : switch a row on/off
%
% 'functional' = the session-average p < .001 mask of the same task
% contrast: drug and covariate effects are only tested where the task
% itself activates.
batches = struct('enabled', {}, 'type', {}, 'covariateName', {}, 'mask', {});

batches(end+1) = struct('enabled', runDrugComparisons, 'type', 'drugComparisons', 'covariateName', '', 'mask', 'auto');
batches(end+1) = struct('enabled', runAverages,        'type', 'average',         'covariateName', '', 'mask', 'auto');
batches(end+1) = struct('enabled', runSingleDrugs,     'type', 'singleDrug',      'covariateName', '', 'mask', 'auto');

% Baseline dopamine synthesis capacity (ki) in four striatal ROIs
kiCovariates = {'ki_VS_new', 'ki_wholeCaudate', 'ki_wholePutamen', 'ki_wholeStriatum'};
for c = 1:numel(kiCovariates)
    batches(end+1) = struct('enabled', runDrugComparisons, 'type', 'drugComparisons', 'covariateName', kiCovariates{c}, 'mask', 'auto'); %#ok<SAGROW>
    batches(end+1) = struct('enabled', runAverages,        'type', 'average',         'covariateName', kiCovariates{c}, 'mask', 'auto'); %#ok<SAGROW>
    batches(end+1) = struct('enabled', runSingleDrugs,     'type', 'singleDrug',      'covariateName', kiCovariates{c}, 'mask', 'auto'); %#ok<SAGROW>
end

% Individual behavioural drug effects. Stems must match the conName of the
% todo.covariate.beh* blocks in run_groupLevel_vici_MID.m. Analyses whose
% SPM.mat does not exist yet are skipped with a warning.
behCovariates = {'mphBehEffect_early', 'mphBehEffect_rt', 'mphBehEffect_hit', ...
                  'mphBehEffect_late', 'rtEffect_MPHvsSUL', 'hitEffect_MPHvsSUL', ...
                  'lateEffect_MPHvsSUL', 'SULvsPBO_late', 'TauEffect_MPHvsSUL', ...
                  'tauEffect_MPHvsPBO'};
for c = 1:numel(behCovariates)
    batches(end+1) = struct('enabled', runDrugComparisons, 'type', 'drugComparisons', 'covariateName', behCovariates{c}, 'mask', 'auto'); %#ok<SAGROW>
end

% Example of an explicit mask override:
% batches(end+1) = struct('enabled', true, 'type', 'average', 'covariateName', 'ki_VS_new', 'mask', 'none');

% =========================================================================
% Run enabled batches
% =========================================================================
groupLevelBase = fullfile(projectDir, 'bids', 'derivatives', 'mri', 'groupLevel', 'MID', 'regular', hrf_type);

nEnabled = sum([batches.enabled]);
fprintf('\n%d/%d batches enabled.\n', nEnabled, numel(batches));

for b = 1:numel(batches)
    B = batches(b);
    if ~B.enabled
        continue
    end
    hasCov = ~isempty(B.covariateName);

    settings = struct();
    settings.io.jobDir  = fullfile(projectDir, 'bids', 'derivatives', 'mri', 'groupLevel', 'MID', 'jobs');
    settings.io.codeDir = fullfile(projectDir, 'code', 'analysis', 'mri', '2ndLevel', 'MID');
    settings.tfce = baseTfce;
    settings.tfce.hrfType = hrf_type;
    settings.tfce.covariateName = B.covariateName;

    if hasCov
        settings.io.mainOutputDir = fullfile(groupLevelBase, sprintf('covariate_%s', B.covariateName));
    else
        settings.io.mainOutputDir = fullfile(groupLevelBase, 'no_covariate_with_exc');
    end

    switch B.type
        case 'drugComparisons'
            settings.tfce.analyses = {};
            for iC = 1:numel(contrasts)
                for iD = 1:numel(drugComps)
                    settings.tfce.analyses{end+1} = sprintf('%s_%s', drugComps{iD}, contrasts{iC}); %#ok<SAGROW>
                end
            end
        case 'average'
            settings.tfce.analyses = {};
            for iC = 1:numel(contrasts)
                settings.tfce.analyses{end+1} = sprintf('sessionAverage_%s', contrasts{iC}); %#ok<SAGROW>
            end
        case 'singleDrug'
            % One model per drug session. With a covariate this asks how
            % the covariate relates to task activation WITHIN each drug
            % state. It does not test whether that relation differs
            % between drugs (that is the drugComparisons covariate
            % contrast): "significant under MPH but not under PBO" is not
            % evidence of a difference (Nieuwenhuis et al., 2011).
            settings.tfce.analyses = {};
            for iC = 1:numel(contrasts)
                for iD = 1:numel(drugs)
                    settings.tfce.analyses{end+1} = sprintf('%s_%s', drugs{iD}, contrasts{iC}); %#ok<SAGROW>
                end
            end
        otherwise
            error('Unknown batch type: %s', B.type);
    end

    % Resolve the 'auto' mask rule
    if ~isfield(B, 'mask') || isempty(B.mask) || strcmp(B.mask, 'auto')
        if strcmp(B.type, 'average') && ~hasCov
            settings.tfce.mask = 'none';
        elseif strcmp(B.type, 'singleDrug') && ~hasCov
            settings.tfce.mask = 'none';
        else
            settings.tfce.mask = 'functional';
        end
    else
        settings.tfce.mask = B.mask;
    end

    % Permutation scheme: Draper-Stoneman (0). It matters only for designs
    % with a nuisance regressor (covariate models); see estimate_tfce_mid.m.
    settings.tfce.nuisance_method = 0;

    if hasCov
        label = sprintf('%s x %s', B.type, B.covariateName);
    else
        label = B.type;
    end
    fprintf('\n=== Batch %d/%d: %s (%s, mask=%s) ===\n', ...
        b, numel(batches), label, hrf_type, settings.tfce.mask);

    run_tfce(settings);
end
