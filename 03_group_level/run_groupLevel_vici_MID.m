% RUN_GROUPLEVEL_VICI_MID  Driver for the group-level (second-level) models.
%
% For every selected task contrast this calls groupLevel_vici_MID_s.m,
% which specifies and estimates SPM one-sample t-tests on the first-level
% images of:
%   - the session average                 (task effects across drugs)
%   - each drug comparison                (MPHvsPBO, SULvsPBO, MPHvsSUL)
%   - each single drug session            (MPH, PBO, SUL)
% optionally with ONE between-subject covariate (baseline dopamine
% synthesis capacity, ki, from [18F]-FDOPA PET, or an individual
% behavioural drug effect). It then exports binary masks of significant
% voxels (p < .001 uncorrected and FWE p < .05).
%
% Run the no-covariate session-average models first: their p < .001 masks
% are the functional mask used by the covariate and TFCE analyses.
%
% Output: groupLevel/MID/<denoise>/<hrfType>/<no_covariate_with_exc |
%         covariate_<name>>/<analysis>_<contrast>/

clear

% Paths
% -------------------------------------------------------------------------
if isunix
    projectDir = '/project/3017048.01/';
    addpath /home/common/matlab/spm12/
    addpath /home/common/matlab/fieldtrip/qsub
elseif ispc
    projectDir = 'P:\3017048.01';
    if exist('\\fileserver.dccn.nl\home\common\matlab\spm12','dir')
        addpath \\fileserver.dccn.nl\home\common\matlab\spm12
        addpath \\fileserver.dccn.nl\home\common\matlab\fieldtrip\qsub
    else
        addpath 'C:\Program Files\MATLAB\spm12-main'
    end
else
    error('Unknown OS');
end

spm('fmri','defaults');

iostruct = struct('projectDir', projectDir, ...
    'derivMRIdir',fullfile(projectDir,'bids','derivatives','mri'), ...
    'derivPETdir', fullfile(projectDir, 'bids','derivatives','pet'), ...
    'functionDir',fullfile(projectDir,'code', 'analysis', 'mri','functions'),...
    'codeDir', fullfile(projectDir,'code','analysis','mri','2ndLevel','MID'), ...
    'functionsDir', fullfile(projectDir,'code','analysis','mri'), ...
    'groupLevelMIDdir', fullfile(projectDir,'bids','derivatives','mri','groupLevel','MID'));

addpath(iostruct.functionDir)

% Unblinding table: session number of each drug per subject
% -------------------------------------------------------------------------
todo.drugs.drugcodes = fullfile(projectDir,'admin','9.Deblinding','drugcodes_sessByDrug.csv');

% Subjects
% -------------------------------------------------------------------------
subjectIx = 1:100;

% Session-level exclusions (imaging quality control + behavioural RT
% outliers, see 02_behaviour/MID_analysis.Rmd)
% .........................................................................
% exclude.subject : subject numbers
% exclude.session : one cell per subject, listing the sessions to exclude.
%                   A drug comparison drops the subject if EITHER of its
%                   two sessions is excluded.
todo.excludeSubjects = true;
exclude.subject = [3,6,9,12,16,20,21,24,26,27,30,36,59,65,76,77,90,94,95,99];
exclude.session = {[2],[3],[3],[2,3],[2],[1,3],[3],[3],[2],[1],[1],[2,3],[3],[1],[1],[1],[2],[1],[3],[1,3]};

% First-level GLMs: one per session ('separate') or sessions concatenated
% -------------------------------------------------------------------------
firstLevelType = 'separate';

% Denoising level of the first-level data
% -------------------------------------------------------------------------
todo.denoise.regular        = true;
todo.denoise.aggressive     = false;
todo.denoise.nonaggressive  = false;
todo.denoise.semiaggressive = false;

% Regular task analysis (or 'ppi')
% -------------------------------------------------------------------------
todo.analysisType = 'regular';

% Model
% -------------------------------------------------------------------------
todo.designEstimation = true;            % false: only (re)compute contrasts
todo.hrfType = 'canonical_hrf';          % 'canonical_hrf' | 'fir'
todo.fir.Nbins = 11;
todo.parameterModulator = false;

% Which session averages / drug comparisons / single drugs to model
% -------------------------------------------------------------------------
% drugs.do       : 'average' = session average; 'XvsY' = drug comparison;
%                  'X' = single drug session
% drugs.comptest : 'one sample' = one-sample t-test on the first-level
%                  difference images (canonical HRF); 'paired' = paired
%                  t-test on the two sessions' images; 'factorial' =
%                  Subject x TimeBin design (FIR)
todo.drugs.do       = {'average','MPHvsPBO','SULvsPBO','MPHvsSUL','PBO', 'MPH', 'SUL'};
todo.drugs.comptest     = 'one sample';

% Small-volume correction (optional, independent of TFCE)
% -------------------------------------------------------------------------
todo.svc.do = true;
todo.svc.mask_type = 'functional';     % 'functional' or 'anatomical'
todo.svc.mask_file = fullfile(projectDir, "code", 'analysis','pet', 'roi', 'striatal_subregions', 'K5wholeStriatum.nii');               % used if anatomical
todo.svc.contrasts = [1 2];
todo.svc.threshold = 0.05;             % FWE
todo.svc.extent = 0;

% Task contrasts to analyse
% -------------------------------------------------------------------------
todo.contrast.Target     = true;
todo.contrast.Target_bins = false;
todo.contrast.Target_peak_window = false;
todo.contrast.HCvsLC     = true;
todo.contrast.Cue        = true;
todo.contrast.CuevsOut   = false;
todo.contrast.HitvsMiss  = true;
todo.contrast.MissvsHit  = false;
todo.contrast.TvsOut     = false;
todo.contrast.HHvsLH     = false;
todo.contrast.HMvsLM     = false;

if todo.parameterModulator
    todo.contrast.HitvsMiss_surprise = false;
    todo.contrast.MissvsHit_surprise = false;
    todo.contrast.Target_p_hit = false;
end

% Covariates -- set AT MOST ONE .use = true per run
% -------------------------------------------------------------------------
% Each covariate is mean-centred and enters the model as a between-subject
% regressor; the output folder is named covariate_<conName>.

% Baseline dopamine synthesis capacity (ki, [18F]-FDOPA PET), per striatal ROI
todo.covariate.ki.use           = false;
todo.covariate.ki.datafile      = fullfile(iostruct.derivPETdir,'group','roi_results_native.csv');
todo.covariate.ki.roi           = 'wholeStriatum';
todo.covariate.ki.conName       = sprintf('ki_%s',todo.covariate.ki.roi);

todo.covariate.ki_cau.use       = false;
todo.covariate.ki_cau.datafile  = fullfile(iostruct.derivPETdir,'group','roi_results_native.csv');
todo.covariate.ki_cau.roi       = 'wholeCaudate';
todo.covariate.ki_cau.conName   = sprintf('ki_%s',todo.covariate.ki_cau.roi);

todo.covariate.ki_put.use       = false;
todo.covariate.ki_put.datafile  = fullfile(iostruct.derivPETdir,'group','roi_results_native.csv');
todo.covariate.ki_put.roi       = 'wholePutamen';
todo.covariate.ki_put.conName   = sprintf('ki_%s',todo.covariate.ki_put.roi);

todo.covariate.ki_vs.use        = false;
todo.covariate.ki_vs.datafile   = fullfile(iostruct.derivPETdir,'group','roi_results_native.csv');
todo.covariate.ki_vs.roi        = 'VS';
todo.covariate.ki_vs.conName    = sprintf('ki_%s_new',todo.covariate.ki_vs.roi);

% Individual behavioural drug effects (per-subject drug differences
% written by 02_behaviour/MID_analysis.Rmd). Each is paired with the drug
% comparison it was computed from.
behCovFile = fullfile(projectDir,'bids','derivatives','beh','MID','stats','MID_output','MID_beh_covariates.csv');

% MPH vs PBO, premature responses
todo.covariate.behEarlyMPHvsPBO.use        = false;
todo.covariate.behEarlyMPHvsPBO.datafile   = behCovFile;
todo.covariate.behEarlyMPHvsPBO.variable   = 'MPH_PBO_early';
todo.covariate.behEarlyMPHvsPBO.conName    = 'mphBehEffect_early';

% MPH vs PBO, reaction time
todo.covariate.behRTMPHvsPBO.use        = false;
todo.covariate.behRTMPHvsPBO.datafile   = behCovFile;
todo.covariate.behRTMPHvsPBO.variable   = 'MPH_PBO_rt';
todo.covariate.behRTMPHvsPBO.conName    = 'mphBehEffect_rt';

% MPH vs SUL, reaction time
todo.covariate.behRTMPHvsSUL.use        = false;
todo.covariate.behRTMPHvsSUL.datafile   = behCovFile;
todo.covariate.behRTMPHvsSUL.variable   = 'MPH_SUL_rt';
todo.covariate.behRTMPHvsSUL.conName    = 'rtEffect_MPHvsSUL';

% MPH vs PBO, hit rate
todo.covariate.behHitMPHvsPBO.use        = false;
todo.covariate.behHitMPHvsPBO.datafile   = behCovFile;
todo.covariate.behHitMPHvsPBO.variable   = 'MPH_PBO_hit';
todo.covariate.behHitMPHvsPBO.conName    = 'mphBehEffect_hit';

% MPH vs SUL, hit rate
todo.covariate.behHitMPHvsSUL.use        = false;
todo.covariate.behHitMPHvsSUL.datafile   = behCovFile;
todo.covariate.behHitMPHvsSUL.variable   = 'MPH_SUL_hit';
todo.covariate.behHitMPHvsSUL.conName    = 'hitEffect_MPHvsSUL';

% SUL vs PBO, too-late responses
todo.covariate.behLATEeffects.use        = false;
todo.covariate.behLATEeffects.datafile   = behCovFile;
todo.covariate.behLATEeffects.variable   = 'SUL_PBO_late';
todo.covariate.behLATEeffects.conName    = 'SULvsPBO_late';

% MPH vs PBO, too-late responses
todo.covariate.behLateMPHvsPBO.use        = false;
todo.covariate.behLateMPHvsPBO.datafile   = behCovFile;
todo.covariate.behLateMPHvsPBO.variable   = 'MPH_PBO_late';
todo.covariate.behLateMPHvsPBO.conName    = 'mphBehEffect_late';

% MPH vs SUL, too-late responses
todo.covariate.behLateMPHvsSUL.use        = false;
todo.covariate.behLateMPHvsSUL.datafile   = behCovFile;
todo.covariate.behLateMPHvsSUL.variable   = 'MPH_SUL_late';
todo.covariate.behLateMPHvsSUL.conName    = 'lateEffect_MPHvsSUL';

% MPH vs SUL, ex-Gaussian tau (slow-response tail)
todo.covariate.behTauEffects.use         = false;
todo.covariate.behTauEffects.datafile    = behCovFile;
todo.covariate.behTauEffects.variable    = 'MPH_SUL_tau';
todo.covariate.behTauEffects.conName     = 'TauEffect_MPHvsSUL';

% MPH vs PBO, ex-Gaussian tau
todo.covariate.behTauMPHvsPBO.use         = false;
todo.covariate.behTauMPHvsPBO.datafile    = behCovFile;
todo.covariate.behTauMPHvsPBO.variable    = 'MPH_PBO_tau';
todo.covariate.behTauMPHvsPBO.conName     = 'tauEffect_MPHvsPBO';

% Results
% -------------------------------------------------------------------------
% groupContrast    : create the group-level t-contrasts (task effect and,
%                    if applicable, the covariate effect)
% export.do        : write binary masks of significant voxels per contrast
% export.threshold / thresholdType : one mask per (threshold, type) pair
% export.extent    : minimum cluster size
todo.results.groupContrast         = true;
todo.results.export.do             = true;
todo.results.export.threshold      = {0.001, 0.05};
todo.results.export.thresholdType  = {'uncorrected', 'fwe'};
todo.results.export.extent         = 0;


% Parallel processing: one Slurm job per contrast
% -------------------------------------------------------------------------
parallel = false;

% Delete and rebuild existing SPM.mat files
force = 1;

% SPM GUI (serial processing only)
% -------------------------------------------------------------------------
gui = false;

if parallel
    gui = false;
end

% Collect everything in one settings structure
% -------------------------------------------------------------------------
settings = struct('io',iostruct, ...
    'todo',todo, ...
    'subjectIx',subjectIx, ...
    'exclude',exclude, ...
    'firstLevelType',firstLevelType, ...
    'parallel',parallel, ...
    'gui',gui, ...
    'force', force);

% Start group level analysis
% =========================================================================

if (strcmpi(settings.firstLevelType,'concatenate') || strcmpi(settings.firstLevelType,'concatenate_dcm')) && ...
        strcmpi(settings.todo.drugs.comptest,'paired')
    error('Cannot do paired t-tests when sessions were concatenated in firstLevel GLMs.')
end

% Contrasts switched on above
% -------------------------------------------------------------------------
conNms = fieldnames(todo.contrast);

clear cons
j = 1;
for ix = 1:numel(conNms)
    if todo.contrast.(conNms{ix})
        cons{j} = conNms(ix);
        j = j+1;
    end
end
assert(numel(cons)>0,'No contrasts were set to true for processing');

if ~settings.parallel
    if settings.gui
        % SPM GUI would be opened here
    else
        spm('defaults', 'FMRI');
        spm_jobman('initcfg');
    end
end

% Loop over denoising levels
denoise_lvls = fieldnames(settings.todo.denoise);

for iLvl = 1:numel(denoise_lvls)

    denoise = denoise_lvls{iLvl};

    if ~settings.todo.denoise.(denoise)
        continue
    end

    % Output directory: first-level type / denoising / HRF model / covariate
    % ---------------------------------------------------------------------
    switch settings.firstLevelType
        case 'separate'
            settings.io.mainOutputDir = settings.io.groupLevelMIDdir;
        case 'concatenate'
            settings.io.mainOutputDir = fullfile(settings.io.groupLevelMIDdir,'sessionsConcatenated');
        case 'concatenate_dcm'
            settings.io.mainOutputDir = fullfile(settings.io.groupLevelMIDdir,'dcm','concatGLMs_groupResults');
    end

    settings.io.mainOutputDir = fullfile(settings.io.mainOutputDir,denoise,todo.hrfType);

    covs = fieldnames(settings.todo.covariate);
    dos = false(1, numel(covs));

    for icov = 1:numel(covs)
        dos(icov) = settings.todo.covariate.(covs{icov}).use;
    end

    if any(dos)
        dirNm = 'covariate';
        for icov = 1:numel(covs)
            if dos(icov)
                dirNm = [dirNm '_' settings.todo.covariate.(covs{icov}).conName];
            end
        end
        settings.io.mainOutputDir = fullfile(settings.io.mainOutputDir, dirNm);
    else
        settings.io.mainOutputDir = fullfile(settings.io.mainOutputDir, 'no_covariate_with_exc');

    end

    % Run analyses
    % ---------------------------------------------------------------------
    if settings.parallel

        % One Slurm job per contrast
        % -----------------------------------------------------------------
        tmpDir = fullfile(settings.io.groupLevelMIDdir, denoise_lvls{iLvl}, 'slurm_jobs');
        if ~exist(tmpDir, 'dir')
            mkdir(tmpDir);
        end

        logDir = fullfile(tmpDir, 'logs');
        if ~exist(logDir, 'dir')
            mkdir(logDir);
        end

        jobIDs = [];
        for icon = 1:numel(cons)
            con = char(cons{icon});

            timestamp = datestr(now, 'yyyymmdd_HHMMSS');
            jobName = sprintf('%s_%s', con, denoise);
            jobFile = fullfile(tmpDir, sprintf('job_%s_%s.sh', jobName, timestamp));
            settingsFile = fullfile(tmpDir, sprintf('settings_%s_%s.mat', jobName, timestamp));
            matlabScript = fullfile(tmpDir, sprintf('run_%s_%s.m', jobName, timestamp));

            save(settingsFile, 'settings', 'con', 'denoise');

            % MATLAB script executed by the job
            mfid = fopen(matlabScript, 'w');
            fprintf(mfid, '%% Auto-generated MATLAB script for Slurm job\n');
            fprintf(mfid, 'load(''%s'');\n', settingsFile);
            fprintf(mfid, 'addpath(''%s'');\n', settings.io.codeDir);
            fprintf(mfid, 'addpath(''%s'');\n', settings.io.functionsDir);
            fprintf(mfid, 'addpath(''%s'');\n', settings.io.functionDir);
            fprintf(mfid, 'addpath(''/home/common/matlab/spm12'');\n');
            fprintf(mfid, 'try\n');
            fprintf(mfid, '    groupLevel_vici_MID_s(settings, con, denoise);\n');
            fprintf(mfid, '    fprintf(''SUCCESS: Completed analysis for contrast: %%s\\n'', con);\n');
            fprintf(mfid, 'catch ME\n');
            fprintf(mfid, '    fprintf(''ERROR: %%s\\n'', ME.message);\n');
            fprintf(mfid, '    for k = 1:numel(ME.stack)\n');
            fprintf(mfid, '        fprintf(''  %%s (line %%d)\\n'', ME.stack(k).file, ME.stack(k).line);\n');
            fprintf(mfid, '    end\n');
            fprintf(mfid, '    exit(1);\n');
            fprintf(mfid, 'end\n');
            fprintf(mfid, 'exit(0);\n');
            fclose(mfid);

            % Slurm batch script
            fid = fopen(jobFile, 'w');
            fprintf(fid, '#!/bin/bash\n');
            fprintf(fid, '#SBATCH --job-name=%s\n', jobName);
            fprintf(fid, '#SBATCH --time=24:00:00\n');
            fprintf(fid, '#SBATCH --mem=64GB\n');
            fprintf(fid, '#SBATCH --cpus-per-task=1\n');
            fprintf(fid, '#SBATCH --partition=batch\n');
            fprintf(fid, '#SBATCH --output=%s/%s_%%j.out\n', logDir, jobName);
            fprintf(fid, '#SBATCH --error=%s/%s_%%j.err\n', logDir, jobName);
            fprintf(fid, '#SBATCH --mail-type=FAIL\n');
            fprintf(fid, '\n');
            fprintf(fid, 'echo "Starting job for contrast: %s"\n', con);
            fprintf(fid, 'echo "MATLAB script: %s"\n', matlabScript);
            fprintf(fid, 'echo "Settings file: %s"\n', settingsFile);
            fprintf(fid, '\n');
            fprintf(fid, '# Load modules\n');
            fprintf(fid, 'module load matlab\n');
            fprintf(fid, '\n');
            fprintf(fid, '# Run MATLAB script\n');
            fprintf(fid, 'cd %s\n', settings.io.codeDir);
            fprintf(fid, 'matlab -batch "run(''%s'')"\n', matlabScript);
            fprintf(fid, '\n');
            fprintf(fid, 'echo "Job completed with exit code: $?"\n');
            fclose(fid);

            system(sprintf('chmod +x %s', jobFile));

            % Submit and record the job ID ("Submitted batch job 12345")
            [status, result] = system(sprintf('sbatch %s', jobFile));
            if status == 0
                tokens = regexp(result, 'Submitted batch job (\d+)', 'tokens');
                if ~isempty(tokens)
                    jobID = str2double(tokens{1}{1});
                    jobIDs = [jobIDs; jobID];
                    fprintf('Submitted job %d for contrast: %s (denoise: %s)\n', jobID, con, denoise);
                end
            else
                warning('Failed to submit job for contrast: %s', con);
            end
        end

        fprintf('\n========================================\n');
        fprintf('Submitted %d jobs to Slurm cluster\n', numel(jobIDs));
        fprintf('Job IDs: %s\n', num2str(jobIDs'));
        fprintf('========================================\n');
        fprintf('Monitor jobs with: squeue --me\n');
        fprintf('Check logs in: %s\n', logDir);
        fprintf('Check MATLAB scripts in: %s/run_*.m\n', tmpDir);
        fprintf('Cancel all jobs: scancel %s\n', num2str(jobIDs'));
        fprintf('========================================\n');

    else

        % Serial processing, one contrast at a time
        % -----------------------------------------------------------------
        for icon = 1:numel(cons)

            con = char(cons{icon});

            % SPM batch jobs can change the working directory
            cd(settings.io.codeDir);

            try
                groupLevel_vici_MID_s(settings,con, denoise)
            catch ME
                fprintf(1,'Error for contrast: %s\n',con);
                disp(ME)
                for k = 1:numel(ME.stack)
                    fprintf('%s (line %d)\n', ME.stack(k).file, ME.stack(k).line);
                end
                rethrow(ME)
            end
        end
    end
end
