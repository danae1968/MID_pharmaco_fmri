% =========================================================================
% WRAPPER_EXTRACT_CLUSTERS_MID  Cluster tables from significant-voxel masks.
%
% For every analysis folder and every significant-voxel mask in it, calls
% extract_clusters_with_labels.m to write two tables to
% <analysis dir>/cluster_tables/:
%   clusters_<analysis>_<contrast>.csv  one row per peak (up to 3 subpeaks
%                                       per cluster, >= 8 mm apart), with
%                                       Neuromorphometrics label, MNI
%                                       coordinates, peak t, k and pFWE
%   regions_<analysis>_<contrast>.csv   one row per cluster x atlas region:
%                                       voxels in the region, share of the
%                                       cluster, strongest peak in it
%
%   corr_method 'tfce' -> TFCE_significant_voxels_<name>_fwe_p05.nii,
%                         pFWE from the TFCE permutation map
%               'fwe'  -> spmT_NNNN_significant_voxels_<name>_fwe_p05.nii,
%                         pFWE from random field theory
%
% Contrasts are resolved BY NAME in each SPM.mat. Loops over every
% directory in includeDirs (no-covariate + covariate models). Afterwards,
% run combine_cluster_tables.m.
% =========================================================================

% ---- USER CONFIG --------------------------------------------------------
corr_method    = 'tfce';              % 'fwe' | 'tfce'
hrf_type       = 'canonical_hrf';     % 'canonical_hrf' | 'fir'

% Minimum cluster size (voxels) to report. 0 reports every FWE-significant
% cluster; a higher value silently drops small but significant clusters.
opts.cluster_extent   = 0;
opts.min_peak_dist_mm = 8;   % subpeak separation -- SPM's own cluster-report default
opts.max_peaks        = 3;   % subpeaks per cluster -- SPM's own cluster-report default
opts.max_regions      = 12;  % regions per cluster in the region breakdown
opts.min_region_k     = 5;   % ignore regions holding fewer voxels than this

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

% Overwrite existing tables?
% false : skip tables that already exist
% true  : DELETE each cluster_tables folder's CSVs first, then recompute.
%
% Use true whenever corr_method or hrf_type changes. File names do not
% encode corr_method, and a contrast with no surviving voxels writes no
% file, so without the delete a table from an earlier run with other
% settings would survive and end up in the combined table.
force = true;

% Audit mode: extract nothing; instead print one line per contrast (and
% write tfce_audit_<hrf_type>.csv) with the search volume, peak -log10(pFWE),
% number of voxels at p < .05 and number of voxels in the binary mask. A
% quick check that every TFCE run did what it should.
audit_only = false;
% TFCE_log_pFWE_*.nii is SIGNED: positive values are significant in the
% contrast's own direction, negative values in the opposite direction. The
% binary mask uses the positive side only, so always threshold the signed
% map, never abs(). The toolbox interpolates the permutation null, so
% -log10(p) can exceed log10(nperm+1).

% Also table the single-drug models (MPH_<con>, PBO_<con>, SUL_<con>).
% With a covariate these describe the covariate relation WITHIN one drug
% state; they do not test a difference between drugs.
include_single_drug = false;

% Optional local copy of the group-level results: the folder that directly
% contains the analysis directories. Leave empty for the project layout.
local.groupLevel = '';
% -------------------------------------------------------------------------

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
spm('defaults', 'fmri');

if ~isempty(local.groupLevel)
    groupLevelBase = local.groupLevel;
else
    groupLevelBase = fullfile(projectDir, 'bids', 'derivatives', 'mri', 'groupLevel', 'MID', 'regular', hrf_type);
end

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

% ---- Mask naming per correction method ----------------------------------
switch corr_method
    case 'fwe'
        % One file per direction, prefixed by the stat image:
        %   spm[TF]_<idx>_significant_voxels_<name>_fwe_p05.nii
        % The unprefixed significant_voxels_<name>_fwe_p05.nii is the union
        % of both directions (for figure contours) and is not tabled.
        maskGlob = 'spm*_significant_voxels_*_fwe_p05.nii';
        maskPat  = '^spm[TF]_\d{4}_significant_voxels_(.+)_fwe_p05\.nii$';
    case 'tfce'
        maskGlob = 'TFCE_significant_voxels_*_fwe_p05.nii';
        maskPat  = '^TFCE_significant_voxels_(.+)_fwe_p05\.nii$';
    otherwise
        error('Unknown corr_method: %s', corr_method);
end

% ---- Analyses to scan ---------------------------------------------------
contrasts = {'Target', 'HCvsLC', 'HitvsMiss', 'Cue'};
drugComps = {'MPHvsPBO', 'SULvsPBO', 'MPHvsSUL'};
drugs     = {'MPH', 'PBO', 'SUL'};

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

auditRows  = {};
logpThresh = -log10(0.05);             % 1.301, the FWE p<.05 line

% ---- Loop over every top-level analysis dir ------------------------------
for iT = 1:numel(topDirs)
    topDir = topDirs{iT};
    use_covariate = ~strcmp(topDir, 'no_covariate_with_exc');
    if use_covariate
        covariate_name = erase(topDir, 'covariate_');
    else
        covariate_name = '';
    end

    data_dir = fullfile(groupLevelBase, topDir);
    outDir = fullfile(data_dir, 'cluster_tables');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    % force: remove old tables first (see the force setting above)
    if force
        stale = [dir(fullfile(outDir, 'clusters_*.csv')); ...
                 dir(fullfile(outDir, 'regions_*.csv'))];
        for iS = 1:numel(stale)
            delete(fullfile(outDir, stale(iS).name));
        end
        if ~isempty(stale)
            fprintf('force=true: removed %d stale table(s) from %s\n', ...
                numel(stale), outDir);
        end
    end

    if use_covariate
        fprintf('\n=== %s (covariate: %s) ===\n', topDir, covariate_name);
    else
        fprintf('\n=== %s (no covariate) ===\n', topDir);
    end

    for iA = 1:numel(analyses)
        dir_ = fullfile(data_dir, analyses{iA});
        if ~exist(dir_, 'dir')
            warning('Missing folder %s, skipping.', analyses{iA}); continue;
        end

        if audit_only
            lp = dir(fullfile(dir_, 'TFCE_log_pFWE_*.nii'));
            if isempty(lp)
                fprintf('%-24s | %-24s | %s\n', topDir, analyses{iA}, ...
                    'NO TFCE OUTPUT (TFCE never run on this folder)');
                auditRows(end+1, :) = {topDir, analyses{iA}, '', NaN, NaN, ...
                    NaN, NaN, 'no TFCE output'}; %#ok<SAGROW>
                continue
            end
            Sa      = load(fullfile(dir_, 'SPM.mat'), 'SPM');
            conNmsA = strrep({Sa.SPM.xCon.name}, ' ', '_');
            for iL = 1:numel(lp)
                ix = str2double(regexp(lp(iL).name, '\d{4}', 'match', 'once'));
                if ~isnan(ix) && ix >= 1 && ix <= numel(conNmsA)
                    conLab = conNmsA{ix};
                else
                    conLab = sprintf('con%04d', ix);
                end

                Y   = spm_read_vols(spm_vol(fullfile(dir_, lp(iL).name)));
                ok  = isfinite(Y);
                vol = sum(ok(:));
                mx  = max(Y(ok));           % signed: this contrast's own tail
                if isempty(mx), mx = NaN; end
                kSurv = sum(Y(ok) >= logpThresh);

                binFile = fullfile(dir_, ...
                    sprintf('TFCE_significant_voxels_%s_fwe_p05.nii', conLab));
                if exist(binFile, 'file')
                    B    = spm_read_vols(spm_vol(binFile));
                    kBin = sum(B(:) > 0);
                else
                    kBin = NaN;
                end

                flag = '';
                if vol == 0
                    flag = 'EMPTY MAP -- all NaN, TFCE produced nothing';
                elseif isnan(kBin) && kSurv > 0
                    flag = 'survives but no binary mask written';
                elseif ~isnan(kBin) && kBin ~= kSurv
                    flag = sprintf('binary/map mismatch (%d vs %d)', kBin, kSurv);
                end

                fprintf(['%-24s | %-24s | %-38s | vol %6d | max %6.3f | ' ...
                         'k %6d | bin %6d | %s\n'], ...
                    topDir, analyses{iA}, conLab, vol, mx, kSurv, kBin, flag);
                auditRows(end+1, :) = {topDir, analyses{iA}, conLab, vol, mx, ...
                    kSurv, kBin, flag}; %#ok<SAGROW>
            end
            continue
        end

        masks = dir(fullfile(dir_, maskGlob));
        if isempty(masks)
            warning('No %s masks in %s, skipping.', corr_method, analyses{iA}); continue;
        end

        S = load(fullfile(dir_, 'SPM.mat'), 'SPM');
        conNms = strrep({S.SPM.xCon.name}, ' ', '_');

        for iM = 1:numel(masks)
            maskName = masks(iM).name;
            tok = regexp(maskName, maskPat, 'tokens', 'once');
            if isempty(tok)
                continue    % union file, SVC_ exports, or anything else non-matching
            end
            label    = tok{1};

            % Skip the both-directions union written by
            % wrapper_slice_display.m: a figure byproduct, not a contrast
            if endsWith(label, '_bothdir')
                continue
            end

            maskPath = fullfile(dir_, maskName);

            % A covariate model also contains the task main effect, which
            % is reported from the no-covariate model; keep only the
            % covariate contrasts here.
            isCovCon = use_covariate && ~isempty(regexp(label, ...
                sprintf('_x_%s$', regexptranslate('escape', covariate_name)), 'once'));
            if use_covariate ~= isCovCon
                continue
            end

            outCsv = fullfile(outDir, sprintf('clusters_%s_%s.csv', analyses{iA}, label));
            if ~force && exist(outCsv, 'file')
                fprintf('%s | %s: %s already exists, skipping (force=false)\n', ...
                    analyses{iA}, label, outCsv);
                continue;
            end

            cIx = find(strcmp(conNms, label), 1);
            if isempty(cIx)
                warning('No contrast matches mask %s in %s, skipping.', maskName, analyses{iA});
                continue;
            end

            % Stat image: spmF for F-contrasts (e.g. TimeBin_effect), spmT otherwise
            if S.SPM.xCon(cIx).STAT == 'F'
                statPath = fullfile(dir_, sprintf('spmF_%04d.nii', cIx));
            else
                statPath = fullfile(dir_, sprintf('spmT_%04d.nii', cIx));
            end

            thisOpts               = opts;
            thisOpts.spmT_path     = statPath;
            thisOpts.contrast_name = analyses{iA};
            % Per-region breakdown: a few peaks cannot describe a TFCE
            % cluster of thousands of voxels spanning many regions
            thisOpts.region_csv    = fullfile(outDir, ...
                sprintf('regions_%s_%s.csv', analyses{iA}, label));
            if ~isempty(regexp(label, '^negative_', 'once'))
                thisOpts.direction = 'negative';
            else
                thisOpts.direction = 'positive';
            end

            % Peak-level FWE p source
            switch corr_method
                case 'tfce'
                    thisOpts.tfce_p_path  = fullfile(dir_, sprintf('TFCE_log_pFWE_%04d.nii', cIx));
                case 'fwe'
                    thisOpts.spm          = S.SPM;   % RFT peak-level FWE p
                    thisOpts.contrast_idx = cIx;
            end

            fprintf('\n=== %s | %s (con %d, %s, %s) ===\n', ...
                analyses{iA}, label, cIx, thisOpts.direction, corr_method);
            extract_clusters_with_labels(maskPath, outCsv, thisOpts);
        end
    end
end

if audit_only
    if isempty(auditRows)
        fprintf('\nAudit found no analysis folders at all under %s\n', groupLevelBase);
    else
        A = cell2table(auditRows, 'VariableNames', {'Analysis_dir', 'Folder', ...
            'Contrast', 'SearchVolume', 'Max_logP', ...
            'k_above_threshold', 'k_in_binary', 'Flag'});
        auditCsv = fullfile(groupLevelBase, sprintf('tfce_audit_%s.csv', hrf_type));
        writetable(A, auditCsv);
        nSurv = sum(cell2mat(auditRows(:, 6)) > 0);
        nFlag = sum(~cellfun(@isempty, auditRows(:, 8)));
        fprintf(['\nAudit: %d contrasts scanned, %d with surviving voxels, ' ...
                 '%d flagged.\n  %s\n'], size(auditRows, 1), nSurv, nFlag, auditCsv);
    end
else
    fprintf('\nDone (%s, %s). Cluster tables written under:\n  %s\n', ...
        hrf_type, corr_method, groupLevelBase);
end
