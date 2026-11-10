function run_tfce(settings)
% RUN_TFCE  Build and run (or submit) TFCE jobs for a list of group models.
%
% USAGE
%   run_tfce(settings)      % normally called from wrapper_tfce.m
%
% INPUT
%   settings.io.mainOutputDir   - folder holding the analysis folders
%   settings.io.jobDir          - folder for job files and logs
%   settings.io.codeDir         - added to the path inside cluster jobs
%   settings.tfce.analyses      - cellstr of analysis folder names,
%                                 e.g. {'MPHvsPBO_Target'}
%   settings.tfce.mask          - 'functional' : session-average p < .001
%                                   mask of the same task contrast, always
%                                   taken from the no-covariate models
%                                 'none'       : the model's own brain mask
%                                 or a path to a mask image
%   settings.tfce.hrfType       - 'canonical_hrf' | 'fir'; sets the
%                                 contrast name to look up ('<con>' or
%                                 'Mean_<con>')
%   settings.tfce.covariateName - covariate stem; if set, the covariate
%                                 contrast '<con>_x_<cov>' is tested
%                                 instead of the task effect
%   settings.tfce.contrasts     - optional explicit contrast indices; if
%                                 omitted, the positive and 'negative_'
%                                 contrasts are found BY NAME in each
%                                 analysis's own SPM.mat
%   settings.tfce.nperm         - permutations (default 5000)
%   settings.tfce.E_weight      - extent weighting E (default 0.5)
%   settings.tfce.nproc         - cores per job (default 4)
%   settings.tfce.parallel      - submit one Slurm job per analysis
%                                 (default true), else run serially
%   settings.tfce.nuisance_method - permutation scheme (0 = Draper-
%                                 Stoneman, default); see estimate_tfce_mid.m
%
% EXAMPLE
%   settings.tfce.analyses = {'MPHvsPBO_Target', 'SULvsPBO_Target'};
%   settings.tfce.mask     = 'functional';
%   settings.tfce.hrfType  = 'canonical_hrf';
%   run_tfce(settings);
%
% -------------------------------------------------------------------------
% Danae Papadopetraki, April 2026

% Defaults
% =========================================================================
if ~isfield(settings.tfce, 'nperm'), settings.tfce.nperm = 5000; end
if ~isfield(settings.tfce, 'nproc'), settings.tfce.nproc = 4; end
if ~isfield(settings.tfce, 'E_weight'), settings.tfce.E_weight = 0.5; end
if ~isfield(settings.tfce, 'parallel'), settings.tfce.parallel = true; end
if ~isfield(settings.tfce, 'nuisance_method'), settings.tfce.nuisance_method = 0; end
if ~isfield(settings.tfce, 'covariateName'), settings.tfce.covariateName = ''; end
if ~isfield(settings.tfce, 'hrfType'), settings.tfce.hrfType = ''; end
% settings.tfce.contrasts stays unset unless given: contrasts are then
% resolved by name per analysis folder below.

% Header
% =========================================================================
fprintf('\n=================================================================\n');
fprintf('Threshold-free cluster enhancement (TFCE)\n');
fprintf('=================================================================\n');
fprintf('Analyses: %d directories\n', numel(settings.tfce.analyses));
fprintf('Permutations: %d\n', settings.tfce.nperm);
fprintf('Cores/job: %d\n', settings.tfce.nproc);
fprintf('=================================================================\n\n');

% Build job list
% =========================================================================
job_list = {};

for i = 1:numel(settings.tfce.analyses)
    analysisName = settings.tfce.analyses{i};
    
    % Get SPM.mat
    spmmat = fullfile(settings.io.mainOutputDir, analysisName, 'SPM.mat');
    if ~exist(spmmat, 'file')
        warning('SPM.mat not found for %s, skipping', analysisName);
        continue;
    end
    
    % Task contrast from the folder name (MPHvsPBO_Target -> Target)
    parts = strsplit(analysisName, '_');
    if numel(parts) >= 2
        conNm = parts{end};  % Last part is contrast name
    else
        conNm = analysisName;
    end

    % Get mask
    if strcmp(settings.tfce.mask, 'functional')
        % Session-average p < .001 mask from the NO-covariate models, for
        % every batch: the mask defines where the task activates at all
        % and should not depend on which covariate is being tested.
        % (Same mask as the SVC step in groupLevel_vici_MID_s.m.)
        groupLevelBase = fileparts(settings.io.mainOutputDir);
        mask_file = fullfile(groupLevelBase, 'no_covariate_with_exc', ...
            sprintf('sessionAverage_%s', conNm), ...
            sprintf('significant_voxels_%s_uncorrected_p001.nii', conNm));
    elseif strcmp(settings.tfce.mask, 'none')
        % Whole brain: the model's own estimation mask (in-brain voxels)
        mask_file = fullfile(settings.io.mainOutputDir, analysisName, 'mask.nii');
    else
        % Use provided mask path
        mask_file = settings.tfce.mask;
    end
    
    if ~exist(mask_file, 'file')
        warning('Mask not found for %s: %s\nSkipping', analysisName, mask_file);
        continue;
    end

    % Contrast indices: explicit override, or by name from this
    % analysis's own SPM.mat
    if isfield(settings.tfce, 'contrasts') && ~isempty(settings.tfce.contrasts)
        theseContrasts = settings.tfce.contrasts;
    else
        S = load(spmmat, 'SPM');
        conNms = strrep({S.SPM.xCon.name}, ' ', '_');

        switch settings.tfce.hrfType
            case 'fir',           posBase = sprintf('Mean_%s', conNm);
            case 'canonical_hrf', posBase = conNm;
            otherwise
                error(['settings.tfce.hrfType must be ''fir'' or ' ...
                    '''canonical_hrf'' when settings.tfce.contrasts is not given']);
        end

        if ~isempty(settings.tfce.covariateName)
            % Covariate contrasts are named '<COND>_x_<cov>' for both HRF
            % models (see groupLevel_vici_MID_s.m)
            posName = sprintf('%s_x_%s', conNm, settings.tfce.covariateName);
        else
            posName = posBase;
        end
        negName = ['negative_' posName];

        posIx = find(strcmp(conNms, posName), 1);
        negIx = find(strcmp(conNms, negName), 1);
        theseContrasts = [posIx negIx];

        if isempty(theseContrasts)
            warning('No contrast matching "%s"/"%s" in %s, skipping', ...
                posName, negName, spmmat);
            continue;
        end
        fprintf('  %s: resolved contrasts %s -> [%s]\n', analysisName, ...
            strjoin({posName, negName}, ' / '), num2str(theseContrasts));
    end

    % Add to job list
    job.spmmat = spmmat;
    job.mask = mask_file;
    job.name = analysisName;
    job.contrasts = theseContrasts;
    job.settings = settings;
    job_list{end+1} = job;

    fprintf('Job %d: %s\n', numel(job_list), analysisName);
end

if isempty(job_list)
    warning('No valid jobs to run!');
    return;
end

fprintf('\nTotal jobs: %d\n\n', numel(job_list));

% Execute jobs
% =========================================================================
if settings.tfce.parallel
    fprintf('Submitting to cluster...\n');

    jobDir = fullfile(settings.io.jobDir, 'slurm_jobs');
    if ~exist(jobDir, 'dir'), mkdir(jobDir); end

    % MATLAB launcher inside each job; must support -batch (R2019a+).
    % Override with settings.tfce.matlabCmd, e.g. 'matlab2023b'.
    if isfield(settings.tfce, 'matlabCmd')
        matlabCmd = settings.tfce.matlabCmd;
    else
        matlabCmd = 'matlab';
    end

    % Each job starts a fresh MATLAB, so it must add and initialise SPM
    % itself; use the SPM of the submitting session.
    spmDir = spm('Dir');

    for i = 1:numel(job_list)
        job = job_list{i};
        tag = sprintf('tfce_%d_%s', i, datestr(now, 'HHMMSSFFF'));
        jobFile = fullfile(jobDir, [tag '.mat']);
        save(jobFile, 'job');

        % (1) MATLAB script run by the job (a file avoids nested
        % shell/MATLAB quoting)
        scriptFile = fullfile(jobDir, [tag '.m']);
        fid = fopen(scriptFile, 'w');
        fprintf(fid, 'addpath(''%s'');\n', settings.io.codeDir);
        fprintf(fid, 'addpath(''%s'');\n', spmDir);
        fprintf(fid, 'spm(''defaults'', ''fmri'');\n');
        fprintf(fid, 'spm_jobman(''initcfg'');\n');   % registers cfg_util / matlabbatch
        fprintf(fid, 'S = load(''%s''); job = S.job;\n', jobFile);
        fprintf(fid, 'tfce.mask = {job.mask};\n');
        fprintf(fid, 'tfce.results_title_prefix = sprintf(''TFCE_%%s'', job.name);\n');
        fprintf(fid, 'tfce.nperm = job.settings.tfce.nperm;\n');
        fprintf(fid, 'tfce.weight_clustsize = job.settings.tfce.E_weight;\n');
        fprintf(fid, 'tfce.modelConIxs = job.contrasts;\n');
        fprintf(fid, 'tfce.nproc = job.settings.tfce.nproc;\n');
        fprintf(fid, 'tfce.nuisance_method = job.settings.tfce.nuisance_method;\n');
        fprintf(fid, 'jd = fullfile(job.settings.io.jobDir, ''tfce_jobs'');\n');
        fprintf(fid, 'if ~exist(jd, ''dir''), mkdir(jd); end\n');
        fprintf(fid, 'estimate_tfce_mid({job.spmmat}, jd, tfce);\n');
        fclose(fid);

        % (2) Slurm batch script
        [~, scriptName] = fileparts(scriptFile);
        shFile = fullfile(jobDir, [tag '.sh']);
        fid = fopen(shFile, 'w');
        fprintf(fid, '#!/bin/bash\n');
        fprintf(fid, '#SBATCH --job-name=%s\n', tag);
        fprintf(fid, '#SBATCH --partition=batch\n');
        fprintf(fid, '#SBATCH --time=04:00:00\n');
        fprintf(fid, '#SBATCH --mem=12GB\n');
        fprintf(fid, '#SBATCH --cpus-per-task=%d\n', settings.tfce.nproc);
        fprintf(fid, '#SBATCH --output=%s/%s_%%j.out\n', jobDir, tag);
        fprintf(fid, 'mkdir -p %s\n', jobDir);
        fprintf(fid, 'cd %s || exit 1\n', jobDir);  % so -batch finds the script
        fprintf(fid, '%s -batch %s\n', matlabCmd, scriptName);
        fclose(fid);

        [status, result] = system(sprintf('sbatch %s', shFile));
        if status ~= 0
            warning('Failed to submit job %d (%s): %s', i, job.name, strtrim(result));
        else
            fprintf('  job %d (%s): %s', i, job.name, result);  % prints "Submitted batch job <id>"
        end
    end
    fprintf(['Jobs submitted. Monitor: squeue --me   (or turm --me)\n' ...
             'Failed-fast jobs will not show in squeue; check sacct --me and the\n' ...
             '.out logs in:\n  %s\n\n'], jobDir);
else
    fprintf('Running serially...\n\n');
    for i = 1:numel(job_list)
        try
            run_one_job(job_list{i});
        catch ME
            warning(sprintf('Job %d failed: %s', i, ME.message));
        end
    end
end

fprintf('TFCE complete!\n\n');

end

%==========================================================================
% Helper function
%==========================================================================
function run_one_job(job)
fprintf('\n=== TFCE: %s ===\n', job.name);

% Setup TFCE
tfce.mask = {job.mask};
tfce.results_title_prefix = sprintf('TFCE_%s', job.name);
tfce.nperm = job.settings.tfce.nperm;
tfce.weight_clustsize = job.settings.tfce.E_weight;
tfce.modelConIxs = job.contrasts;
tfce.nproc = job.settings.tfce.nproc;
tfce.nuisance_method = job.settings.tfce.nuisance_method;

% Job directory
jobDir = fullfile(job.settings.io.jobDir, 'tfce_jobs');
if ~exist(jobDir, 'dir'), mkdir(jobDir); end

% Run TFCE
estimate_tfce_mid({job.spmmat}, jobDir, tfce);

fprintf('=== COMPLETE ===\n\n');
end