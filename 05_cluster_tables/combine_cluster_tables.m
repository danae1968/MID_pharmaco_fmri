% COMBINE_CLUSTER_TABLES  Merge all cluster tables into one file.
%
% Merges every clusters_*.csv and regions_*.csv written by
% wrapper_extract_clusters_MID.m, across all analysis directories, into
%   all_clusters_<hrf_type>.csv  and  all_regions_<hrf_type>.csv
% in the group-level folder. These are read by the "Neural cluster
% tables" section of 02_behaviour/MID_analysis.Rmd.
%
% Adds three columns so the combined table can be filtered directly:
%   Covariate      - '' for no_covariate_with_exc, else the covariate stem
%                    (e.g. 'ki_VS_new'), taken from the directory name.
%   DrugComparison - 'MPHvsPBO' | 'SULvsPBO' | 'MPHvsSUL' | 'sessionAverage',
%                    parsed from Contrast.
%   TaskContrast   - 'Cue' | 'HCvsLC' | 'HitvsMiss' | 'Target', parsed from
%                    Contrast.
%
% -------------------------------------------------------------------------
% Danae Papadopetraki, August 2026

% ---- USER CONFIG ----------------------------------------------------------
hrf_type = 'canonical_hrf';   % 'canonical_hrf' | 'fir' -- must match what
                               % wrapper_extract_clusters_MID.m was run with
% Analysis directories to combine. If empty, every covariate_* directory
% under groupLevelBase is found automatically (minus excludeDirs).
includeDirs = {'no_covariate_with_exc', ...
               'covariate_ki_VS_new', 'covariate_ki_wholeCaudate', ...
               'covariate_ki_wholePutamen', 'covariate_ki_wholeStriatum', ...
               'covariate_mphBehEffect_early', 'covariate_mphBehEffect_hit', ...
               'covariate_mphBehEffect_rt', 'covariate_mphBehEffect_late', ...
               'covariate_rtEffect_MPHvsSUL', 'covariate_hitEffect_MPHvsSUL', ...
               'covariate_lateEffect_MPHvsSUL', 'covariate_SULvsPBO_late', ...
               'covariate_TauEffect_MPHvsSUL', 'covariate_tauEffect_MPHvsPBO'};

excludeDirs = {};              % only used when includeDirs is empty
outCsv = [];                   % leave empty for the default path below

% Optional local copy of the group-level results: the folder that directly
% contains the analysis directories. Leave empty for the project layout.
local.groupLevel = '';
% -------------------------------------------------------------------------

if ispc
    projectDir = 'P:\3017048.01\';
elseif isunix
    projectDir = '/project/3017048.01/';
else
    error('Unknown OS')
end

if ~isempty(local.groupLevel)
    groupLevelBase = local.groupLevel;
else
    groupLevelBase = fullfile(projectDir, 'bids', 'derivatives', 'mri', 'groupLevel', 'MID', 'regular', hrf_type);
end

if isempty(outCsv)
    outCsv = fullfile(groupLevelBase, sprintf('all_clusters_%s.csv', hrf_type));
end

% ---- Analysis directories to combine --------------------------------------
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
fprintf('Combining cluster tables from %d directories:\n', numel(topDirs));
fprintf('  %s\n', strjoin(topDirs, '\n  '));

% ---- Read and tag every clusters_*.csv -------------------------------------
allTables = {};
allRegionTables = {};
for iT = 1:numel(topDirs)
    topDir = topDirs{iT};
    if strcmp(topDir, 'no_covariate_with_exc')
        covariateName = '';
    else
        covariateName = erase(topDir, 'covariate_');
    end

    ctDir = fullfile(groupLevelBase, topDir, 'cluster_tables');
    if ~exist(ctDir, 'dir')
        warning('No cluster_tables folder in %s, skipping.', topDir);
        continue;
    end

    % Peak tables and region tables are tagged identically and written to
    % separate combined files.
    files = [dir(fullfile(ctDir, 'clusters_*.csv')); ...
             dir(fullfile(ctDir, 'regions_*.csv'))];
    if isempty(files)
        warning('No cluster CSVs in %s, skipping.', ctDir);
        continue;
    end

    for iF = 1:numel(files)
        fPath = fullfile(ctDir, files(iF).name);
        try
            % Read text as string so column types match across files
            T = readtable(fPath, 'TextType', 'string');
        catch ME
            warning('Could not read %s (%s), skipping.', fPath, ME.message);
            continue;
        end
        if height(T) == 0
            continue;   % file exists but no surviving clusters
        end

        % Harmonise the p-value column name (e.g. 'TFCE_pFWE') to 'pFWE' so
        % the tables can be concatenated
        vn = T.Properties.VariableNames;
        pCol = vn(~cellfun('isempty', regexpi(vn, 'pFWE')));
        if numel(pCol) == 1 && ~strcmp(pCol{1}, 'pFWE')
            T = renamevars(T, pCol{1}, 'pFWE');
        elseif numel(pCol) > 1
            warning('%s has multiple pFWE-like columns (%s), leaving as-is -- vertcat may fail.', ...
                fPath, strjoin(pCol, ', '));
        end

        % Contrast is always "<DrugComparison-or-sessionAverage>_<TaskContrast>"
        % with exactly one underscore (neither half ever contains one).
        parts = split(T.Contrast, '_');
        if size(parts, 1) ~= height(T)
            % For a single-row table split() returns a 2x1 column instead
            % of 1x2: transpose
            parts = parts';
        end
        T.DrugComparison = parts(:, 1);
        T.TaskContrast   = parts(:, 2);
        T.Covariate      = repmat(string(covariateName), height(T), 1);

        if startsWith(files(iF).name, 'regions_')
            allRegionTables{end+1} = T; %#ok<SAGROW>
        else
            allTables{end+1} = T; %#ok<SAGROW>
        end
    end
end

if isempty(allTables) && isempty(allRegionTables)
    error('No cluster tables found to combine under %s', groupLevelBase);
end

% Identifying columns first, then the original per-peak/per-region columns.
firstCols = {'Covariate', 'DrugComparison', 'TaskContrast'};

if ~isempty(allTables)
    combined  = vertcat(allTables{:});
    otherCols = setdiff(combined.Properties.VariableNames, firstCols, 'stable');
    combined  = combined(:, [firstCols, otherCols]);
    writetable(combined, outCsv);
    fprintf('\nCombined %d cluster table files (%d peak rows total) -> %s\n', ...
        numel(allTables), height(combined), outCsv);
end

if ~isempty(allRegionTables)
    regionOut = regexprep(outCsv, 'all_clusters_', 'all_regions_');
    combinedR = vertcat(allRegionTables{:});
    otherCols = setdiff(combinedR.Properties.VariableNames, firstCols, 'stable');
    combinedR = combinedR(:, [firstCols, otherCols]);
    writetable(combinedR, regionOut);
    fprintf('Combined %d region table files (%d cluster-region rows total) -> %s\n', ...
        numel(allRegionTables), height(combinedR), regionOut);
else
    fprintf(['No regions_*.csv found -- re-run wrapper_extract_clusters_MID.m ' ...
             '(force = true) to write the per-region breakdown.\n']);
end
