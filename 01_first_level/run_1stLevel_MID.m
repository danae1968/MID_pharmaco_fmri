% RUN_1STLEVEL_MID  Driver for the subject-level (first-level) fMRI GLMs.
%
% Sets all options for the first-level analysis of the MID task and calls
% firstLevel_MID_fir.m for every subject, either serially or as one cluster
% job per subject. Per subject and drug session it:
%   1. builds the SPM multiple-conditions file from the task onsets
%   2. specifies the GLM (canonical HRF or FIR basis set)
%   3. estimates the GLM
%   4. estimates the task contrasts
% and then, across the three sessions of a subject,
%   5. computes drug-difference contrast images (MPH-PBO, SUL-PBO, MPH-SUL)
%   6. computes the session-average contrast image
%
% Steps 5 and 6 produce the per-subject images that enter the group-level
% one-sample t-tests (03_group_level).

% Paths
% -------------------------------------------------------------------------
starttime = tic;

if ispc
    projectDir = 'P:\3017048.01\';
    if exist('H:\common\matlab\spm12','dir')
        addpath H:\common\matlab\spm12
        addpath H:\common\matlab\fieldtrip\qsub
    else
        addpath 'C:\Program Files\MATLAB\spm12-main'
    end

elseif isunix
    projectDir = '/project/3017048.01/';
    addpath /home/common/matlab/spm12
    addpath /home/common/matlab/fieldtrip/qsub
else
    error('Unknown OS')
end

addpath(fullfile(projectDir,'code'))

% Input/output structure (BIDS layout, fMRIPrep derivatives)
% -------------------------------------------------------------------------
iostruct = struct('projectDir', projectDir, ...
                  'bidsDir', fullfile(projectDir,'bids'), ...
                  'derivDir', fullfile(projectDir,'bids','derivatives'), ...
                  'derivMRIdir',fullfile(projectDir,'bids','derivatives','mri'), ...
                  'codeDir', fullfile(projectDir,'code','analysis','mri','1stLevel','MID'));

% Unblinding table: which session number each subject received each drug
% -------------------------------------------------------------------------
drugcodes = fullfile(projectDir,'admin','9.Deblinding','drugcodes_sessByDrug.csv');
addpath(iostruct.codeDir)

% Number of fMRI runs per session
% -------------------------------------------------------------------------
nRuns = 1;

% Haemodynamic response model
% -------------------------------------------------------------------------
% 'canonical_hrf' : SPM canonical HRF (used for the reported results)
% 'fir'           : finite impulse response basis set (exploratory)
hrf_response                = 'canonical_hrf';
firCondition                = 'Target';
window                      = 25.52; % FIR window length (s) = 11 TRs
order                       = 11;    % number of FIR bins

% Steps to run
% -------------------------------------------------------------------------
todo.multipleConditions     = true;
todo.modelSpecification     = true;
todo.modelEstimation        = true;
todo.contrastEstimation     = true;

% Regressors / contrasts to leave out (by name). Button presses are not
% modelled as a separate regressor in the final model.
todo.exclude.regressors      = {'Resp'};
todo.exclude.contrasts       = {'buttonPress'};

% Drug-difference contrast images per subject (drug1 con - drug2 con):
%   MPHvsPBO, SULvsPBO, MPHvsSUL
todo.drugComparisons        = true;

% Average contrast image over the three sessions per subject
todo.averageContrasts       = true;

% Spatial smoothing (fMRIPrep output is unsmoothed)
todo.smooth.use             = true;
todo.smooth.kernel          = [6 6 6];   % FWHM in mm

% Denoising strategy: 'regular' = motion, FD, CSF/WM and aCompCor
% confounds from fMRIPrep (see firstLevel_MID_fir.m)
todo.denoise.regular        = true;
todo.denoise.aggressive     = false;
todo.denoise.nonaggressive  = false;
todo.denoise.semiaggressive = false;

% Optional trial-wise parametric modulator (not used for reported results)
todo.parameterModulator = false;
todo.parameterModulatorAt.target = false;
todo.parameterModulatorAt.outcome = false;
todo.parameterModulatorName = 'P_hit';

% Subjects and sessions
% -------------------------------------------------------------------------
subjectIx = 1:100;
sessions  = [1 2 3];

% Subjects with incomplete data: only these sessions are analysed
deviantSessions.sub  = [42 90 95];
deviantSessions.sess = {[1] [2 3] [1 2]};

% Overwrite existing SPM.mat files?
% -------------------------------------------------------------------------
force = true;

% Parallel processing: one cluster job per subject (FieldTrip qsub)
% -------------------------------------------------------------------------
parallel = false;

% Collect everything in one settings structure
% -------------------------------------------------------------------------
settings.io                 = iostruct;
settings.drugcodes          = drugcodes;
settings.nRuns              = nRuns;
settings.todo               = todo;
settings.sessions           = sessions;
settings.deviantSessions    = deviantSessions;
settings.force              = force;
settings.parallel           = parallel;
settings.window             = window;
settings.order              = order;
settings.hrf_response       = hrf_response;
settings.firCondition       = firCondition;

if settings.parallel

    % One job per subject
    % ---------------------------------------------------------------------
    cellarrSettings = repmat({settings},1,numel(subjectIx));

    clear subjects
    for s = 1:numel(subjectIx)
        subjects{s} = subjectIx(s);

        % Subject-specific session list, if applicable
        if ismember(subjects{s},deviantSessions.sub)
            deviantIx = find(subjects{s} == deviantSessions.sub);
            cellarrSettings{s}.sessions = deviantSessions.sess{deviantIx};
        end
    end

    % memreq = 6 GB, timreq = 600 min per job
    qsubcellfun(@firstLevel_MID_fir, cellarrSettings, subjects, 'memreq', 6*1024^3, 'timreq', 600*60, 'StopOnError', false);

else

    % Serial processing
    % ---------------------------------------------------------------------
    spm('fmri','defaults')

    for iSub = 1:numel(subjectIx)

        subject = subjectIx(iSub);

        % SPM batch jobs can change the working directory
        cd(settings.io.codeDir);

        try
            firstLevel_MID_fir(settings,subject)
        catch ME
            fprintf(1,'Error for sub-%.3d:\n',subject);
            disp(ME)
        end
    end

end

elapsedTime_minutes = toc(starttime)/60
