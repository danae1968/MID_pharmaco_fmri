function firstLevel_MID_fir(settings,iSubject)
% FIRSTLEVEL_MID_FIR  First-level fMRI analysis of the MID task for one subject.
%
% Despite the name, handles both the canonical-HRF and the FIR model
% (settings.hrf_response). For each drug session of the subject:
%   - build the SPM multiple-conditions .mat file (names/onsets/durations)
%   - smooth the fMRIPrep output if needed
%   - specify the GLM: task regressors + fMRIPrep confounds
%   - estimate the GLM
%   - define and estimate the task contrasts
%   - re-compress the functional data
% Then, across the subject's sessions:
%   - drug-difference images (MPHvsPBO, SULvsPBO, MPHvsSUL) per contrast
%   - session-average image per contrast
%
% INPUTS
%   settings : structure built by run_1stLevel_MID.m (paths, steps to run,
%              sessions, model options)
%   iSubject : subject number
%
% DEPENDENCIES
%   define_contrasts.m (canonical HRF), define_contrasts_fir.m (FIR),
%   run_spm_jobs_MID.m, smoothing.m (lab utility, not included)
% .........................................................................
% Danae Papadopetraki, 2019-2026
% Initial structure based on a lab template by Ruben van den Bosch (2017)

if ispc
    if exist('H:\common\matlab\spm12','dir')
        addpath H:\common\matlab\spm12
    else
        addpath 'C:\Program Files\MATLAB\spm12-main'
    end
elseif isunix
    addpath /home/common/matlab/spm12
else
    error('Unknown OS')
end

% A parallel job starts from a fresh MATLAB session, so SPM must be opened
if settings.parallel
    spm fmri
end

% Loop over the requested denoising levels
% =========================================================================
denoise_lvls = fieldnames(settings.todo.denoise);

% The multiple-conditions files do not depend on the denoising level, so
% they are only created on the first pass through this loop.
iter1 = true;

for iLvl = 1:numel(denoise_lvls)

    if ~settings.todo.denoise.(denoise_lvls{iLvl})
        continue
    end

    denoise = denoise_lvls{iLvl};

    % Sessions to analyse for this subject
    % ---------------------------------------------------------------------
    deviantIx = find(settings.deviantSessions.sub == iSubject);
    if ~isempty(deviantIx)
        sessions = settings.deviantSessions.sess{deviantIx};
    else
        sessions = settings.sessions;
    end

    for sesIx = 1:numel(sessions)

        iSession = sessions(sesIx);

        if iSession == 0
            warning('Skipping ses-intake for sub-%.3d. No fMRI on the intake session, only an anatomical scan.',iSubject);
            continue
        end

        fprintf(1,'Running first level analysis for sub-%.3d ses-drug%d. Denoising level: %s\n', iSubject, iSession, denoise);

        % Subject directories
        % -----------------------------------------------------------------
        dirs.bidsFunc    = fullfile(settings.io.bidsDir,sprintf('sub-%.3d',iSubject),sprintf('ses-drug%d',iSession),'func');
        dirs.preproc     = fullfile(settings.io.derivMRIdir,'fmriprep', sprintf('sub-%.3d',iSubject),sprintf('ses-drug%d',iSession),'func');
        dirs.derivMID     = fullfile(settings.io.derivMRIdir,'firstLevel','MID',sprintf('sub-%.3d',iSubject),sprintf('ses-drug%d',iSession));
        dirs.jobs        = fullfile(settings.io.derivMRIdir,'firstLevel','MID',sprintf('sub-%.3d',iSubject),sprintf('ses-drug%d',iSession),'jobs');
        dirs.reports     = fullfile(settings.io.derivMRIdir,'firstLevel','MID',sprintf('sub-%.3d',iSubject),sprintf('ses-drug%d',iSession),'reports');

        % Output directory for this subject, session and denoising level
        if settings.todo.parameterModulator
            dirs.firstLevel = fullfile(dirs.derivMID, denoise, settings.hrf_response, 'surprise_modulation');
        else
            dirs.firstLevel = fullfile(dirs.derivMID, denoise, settings.hrf_response);
        end

        % Skip sessions that are already specified, unless force is true
        % -----------------------------------------------------------------
        SPMmat = fullfile(dirs.firstLevel,'SPM.mat');
        if exist(SPMmat,'file') && settings.todo.modelSpecification && ~settings.force
            continue
        elseif exist(SPMmat,'file') && settings.todo.modelSpecification && settings.force
            delete(SPMmat);
        end

        % Determine the MID run numbers. Aborted runs mean the runs to
        % analyse are not necessarily numbered 1..n or consecutive.
        % -----------------------------------------------------------------
        MIDbold_files = spm_select('FPList',dirs.preproc,sprintf('^sub-%.3d_ses-drug%d_task-MID_run-._.*_bold.nii',iSubject,iSession));

        if isempty(MIDbold_files)
            warning('No MID bold data found for sub-%.3d, ses-drug%d. Skipping this session', iSubject,iSession)
            continue
        end

        runs = [];
        for iFile = 1:size(MIDbold_files,1)
            ix = strfind(MIDbold_files(iFile,:),'run-');
            runs(1,iFile) = str2num(MIDbold_files(iFile,ix+length('run-')));
        end
        runNrs = unique(runs);

        % Make output dirs if they don't exist
        % -----------------------------------------------------------------
        dirNames = fieldnames(dirs);
        for iDir = 1:numel(dirNames)
            if ~ismember(dirNames{iDir},{'rawPerform','func'}) && ~exist(dirs.(dirNames{iDir}),'dir')
                mkdir(dirs.(dirNames{iDir}));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CREATE MULTIPLE CONDITIONS .MAT FILE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Task regressors (from onsetsMID.m):
        %   LowCue, HighCue, Target, LowHit, HighHit, LowMiss, HighMiss

        if settings.todo.multipleConditions && iter1

            clear MIDonsets
            for iRun = 1:numel(runNrs)
                run = sprintf('run%d',iRun);

                % Load onsets written by onsetsMID.m
                % ---------------------------------------------------------
                eventsFname = fullfile(dirs.bidsFunc,sprintf('sub-%.3d_ses-drug%d_task-MID_onsets.mat',iSubject,iSession));
                MIDonsets      = load(eventsFname);

                % Optional trial-wise parametric modulator values
                if settings.todo.parameterModulator
                    clear modulParams
                    modulFname = fullfile(dirs.bidsFunc,sprintf('sub-%.3d_ses-drug%d_task-MID_surprise_parameters.mat',iSubject,iSession));

                    if ~exist(modulFname, 'file')
                        warning('Surprise parameter file not found: %s', modulFname);
                    else
                        modulParams = load(modulFname);
                    end
                end

                % Condition names, minus any excluded regressors
                % =========================================================
                conditionNames = MIDonsets.names;
                conditionNames(ismember(conditionNames,settings.todo.exclude.regressors)) = [];

                % Drop conditions without onsets (SPM cannot model an
                % empty regressor)
                % .........................................................
                nonempty_cond = conditionNames;
                for i = 1:numel(conditionNames)

                    condNm = char(conditionNames(i));

                    if isempty(MIDonsets.onsets(i))
                        for j = 1:numel(nonempty_cond)
                            if strcmp(nonempty_cond(j),condNm)
                                nonempty_cond(j) = [];
                                break
                            end
                        end
                    end
                end

                % Keep the selected conditions in their original order
                [~,~,index]=intersect(conditionNames,MIDonsets.names);
                sorted_names = MIDonsets.names(sort(index));
                sorted_onsets = MIDonsets.onsets(sort(index));
                sorted_durations = MIDonsets.durations(sort(index));

                if settings.todo.parameterModulator
                    % Build SPM-compatible names/onsets/durations/pmod
                    names = {};
                    onsets = {};
                    durations = {};
                    pmod = struct('name', {}, 'param', {}, 'poly', {});
                    orth = {};

                    for i = 1:numel(sorted_names)
                        cond_name = sorted_names{i};
                        cond_onsets = sorted_onsets{i};
                        cond_durations = sorted_durations{i};

                        % Modulator values for this condition, if any
                        if settings.todo.parameterModulatorAt.outcome
                        switch cond_name
                            case 'LowHit'
                                pmod_values = modulParams.LowHit_surprise;
                            case 'HighHit'
                                pmod_values = modulParams.HighHit_surprise;
                            case 'LowMiss'
                                pmod_values = modulParams.LowMiss_surprise;
                            case 'HighMiss'
                                pmod_values = modulParams.HighMiss_surprise;
                            otherwise
                                pmod_values = [];
                        end
                        end
                        if settings.todo.parameterModulatorAt.target
                              switch cond_name
                            case 'Target'
                                pmod_values = modulParams.Target_p_hit;
                            otherwise
                                pmod_values = [];
                              end
                        end

                        % Drop trials with a missing modulator value
                        if ~isempty(pmod_values)
                            valid_indices = ~isnan(pmod_values);
                            cond_onsets = cond_onsets(valid_indices);
                            pmod_values = pmod_values(valid_indices);
                        end

                        if ~isempty(cond_onsets)
                            names{end+1} = cond_name;
                            onsets{end+1} = cond_onsets;
                            durations{end+1} = cond_durations;

                            if ~isempty(pmod_values)
                                pmod(end+1).name = {settings.todo.parameterModulatorName};
                                pmod(end).param = {pmod_values};
                                pmod(end).poly = {1};      % linear modulation
                                orth{end+1} = 0;           % do not orthogonalise
                            else
                                pmod(end+1).name = {''};
                                pmod(end).param = {[]};
                                pmod(end).poly = {[]};
                                orth{end+1} = 1;
                            end
                        end
                    end
                else
                    names = sorted_names;
                    durations = sorted_durations;
                    onsets = sorted_onsets;
                end

                % Save the multiple-conditions file
                % =========================================================
                if settings.todo.parameterModulator
                    filename = fullfile(dirs.derivMID, sprintf('sub-%.3d_ses-drug%d_task-MID_run-%d_conditions_surprise.mat', iSubject, iSession, iRun));
                    save(filename, 'names', 'onsets', 'durations', 'pmod', 'orth');
                else
                    filename = fullfile(dirs.derivMID, sprintf('sub-%.3d_ses-drug%d_task-MID_run-%d_conditions.mat', iSubject, iSession, iRun));
                    save(filename,'names','onsets','durations');
                end

            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FMRI MODEL SPECIFICATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if settings.todo.modelSpecification

            outDir = dirs.firstLevel;

            % Inputs per run
            % -------------------------------------------------------------
            % runNrs selects the fMRIPrep files; the files written here are
            % numbered 1..number of runs.
            for iRun = 1:numel(runNrs)
                run = sprintf('run%d',iRun);

                % Functional data: unzip, or smooth if no smoothed copy
                % exists yet
                % ---------------------------------------------------------
                switch denoise
                    case {'regular', 'semiaggressive'}
                        if settings.todo.smooth.use
                            func4D.(run) = fullfile(dirs.preproc,sprintf('sub-%.3d_ses-drug%d_task-MID_run-%d_space-MNI152NLin2009cAsym_desc-preprocSmooth_bold.nii',iSubject,iSession,runNrs(iRun)));
                        else
                            func4D.(run) = fullfile(dirs.preproc,sprintf('sub-%.3d_ses-drug%d_task-MID_run-%d_space-MNI152NLin2009cAsym_desc-preproc_bold.nii',iSubject,iSession,runNrs(iRun)));
                        end
                    case 'nonaggressive'
                        func4D.(run) = fullfile(dirs.preproc,sprintf('sub-%.3d_ses-drug%d_task-MID_run-%d_space-MNI152NLin2009cAsym_desc-smoothAROMAnonaggr_bold.nii',iSubject,iSession,runNrs(iRun)));
                    case 'aggressive'
                        error('Aggressive denoising in preprocessing or through intermediate regression not implemented')
                end

                if ~exist(func4D.(run),'file')

                    gzfile = [func4D.(run) '.gz'];

                    % ICA-AROMA output is already smoothed, so smoothing is
                    % only needed for the regular/semiaggressive data
                    if settings.todo.smooth.use && ismember(denoise,{'regular','semiaggressive'}) && ~exist(gzfile,'file')

                        unsmoothedFunc    = strrep(func4D.(run),'preprocSmooth','preproc');
                        unsmoothedFunc_gz = [unsmoothedFunc '.gz'];

                        if exist(unsmoothedFunc,'file')
                            smoothing(unsmoothedFunc, settings.todo.smooth.kernel, dirs.jobs);
                        elseif exist(unsmoothedFunc_gz,'file')
                            smoothing(unsmoothedFunc_gz, settings.todo.smooth.kernel, dirs.jobs);
                        else
                            error('Func data file %s not found', unsmoothedFunc)
                        end

                        assert(exist(func4D.(run),'file') > 0, 'Smoothed data %s not found', func4D.(run))

                    else
                        assert(exist(gzfile,'file') > 0, ...
                            'Func data file %s not found',gzfile)

                        if isunix
                            cmd = ['gunzip ' gzfile];
                            system(cmd);
                        elseif ispc
                            warning('gunzip in Matlab is EXTREMELY slow!')
                            gunzip(gzfile);
                        end
                    end
                end

                % All volumes of the 4D file
                funcScans.(run) = cellstr(spm_select('expand', func4D.(run)));

                % =========================================================
                % REGRESSORS
                % =========================================================
                % Task conditions
                % ---------------------------------------------------------
                if settings.todo.parameterModulator
                    multiCond.(run) = fullfile(dirs.derivMID, sprintf('sub-%.3d_ses-drug%d_task-MID_run-%d_conditions_surprise.mat', iSubject, iSession, iRun));
                else
                    multiCond.(run) = fullfile(dirs.derivMID, sprintf('sub-%.3d_ses-drug%d_task-MID_run-%d_conditions.mat', iSubject, iSession, iRun));
                end

                % Nuisance regressors from the fMRIPrep confounds file
                % ---------------------------------------------------------
                confreg = dataset('file', ...
                    fullfile(dirs.preproc,sprintf('sub-%.3d_ses-drug%d_task-MID_run-%d_desc-confounds_regressors.tsv',iSubject,iSession,runNrs(iRun))), ...
                    'delimiter', '\t', 'TreatAsEmpty','n/a');

                confregNms = confreg.Properties.VarNames;

                % Confound set depends on the denoising level
                % .........................................................
                regressors = [];
                switch denoise
                    case {'regular', 'nonaggressive'}

                        % 6 realignment parameters (translation, rotation)
                        regressors = strvcat(confregNms{contains(confregNms,'trans')}, ...
                            confregNms{contains(confregNms,'rot')});

                        % Framewise displacement
                        regressors = strvcat(regressors, ...
                            'framewise_displacement');

                        % Mean CSF and white-matter signal
                        regressors = strvcat(regressors, ...
                            'csf', ...
                            'white_matter');

                        % aCompCor components
                        regressors = strvcat(regressors, ...
                            confregNms{contains(confregNms,'a_comp_cor')});

                    case 'aggressive'
                        error('Aggressive denoising in preprocessing or through intermediate regression not implemented at the moment')

                    case 'semiaggressive'

                        % As 'regular', plus the ICA-AROMA motion components
                        regressors = strvcat(confregNms{contains(confregNms,'trans')}, ...
                            confregNms{contains(confregNms,'rot')});

                        regressors = strvcat(regressors, ...
                            'framewise_displacement');

                        regressors = strvcat(regressors, ...
                            'csf', ...
                            'white_matter');

                        regressors = strvcat(regressors, ...
                            confregNms{contains(confregNms,'a_comp_cor')});

                        regressors = strvcat(regressors, ...
                            confregNms{contains(confregNms,'aroma_motion')});
                end

                if ~isempty(regressors)
                    regressors = cellstr(regressors);
                    regressors(ismember(regressors,settings.todo.exclude.regressors)) = [];

                    % Store names and vectors (matrix R) in a .mat file
                    % .....................................................
                    names = regressors;
                    R = [];
                    for iNm = 1:numel(names)

                        regVector = confreg.(names{iNm});

                        % SPM cannot handle NaN in nuisance regressors (e.g.
                        % FD in the first volume): replace with the mean
                        if any(isnan(regVector))
                            regVector(isnan(regVector)) = nanmean(regVector);
                        end

                        R = [R regVector];
                    end

                    filename = fullfile(dirs.firstLevel,sprintf('sub-%.3d_ses-drug%d_task-MID_run-%d_regressors.mat',iSubject,iSession,iRun));
                    save(filename,'names','R');

                    multiReg.(run) = filename;
                else
                    warning('No confound regressors selected. sub-%.3d ses-drug%d',iSubject,iSession)
                    multiReg.(run) = '';
                end

                % High-pass filter cutoff: SPM default 128 s. ICA-AROMA
                % data are already high-pass filtered, so a very long
                % cutoff is used there to effectively switch it off.
                % ---------------------------------------------------------
                if strcmp(denoise,'regular') || strcmp(denoise,'semiaggressive')
                    HPcutoff.(run) = 128;
                else
                    HPcutoff.(run) = 1e6;
                end
            end

            % Fill in the SPM model specification batch
            % -------------------------------------------------------------
            clear jobs
            jobs{1}.spm.stats.fmri_spec.dir             = cellstr(outDir);
            jobs{1}.spm.stats.fmri_spec.timing.units    = 'secs';
            jobs{1}.spm.stats.fmri_spec.timing.RT       = 2.32;
            jobs{1}.spm.stats.fmri_spec.timing.fmri_t   = 38;
            jobs{1}.spm.stats.fmri_spec.timing.fmri_t0  = 19;
            for iRun = 1:numel(runNrs)
                run = sprintf('run%d',iRun);
                jobs{1}.spm.stats.fmri_spec.sess(iRun).scans        = funcScans.(run);
                jobs{1}.spm.stats.fmri_spec.sess(iRun).cond         = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
                jobs{1}.spm.stats.fmri_spec.sess(iRun).multi        = cellstr(multiCond.(run));
                jobs{1}.spm.stats.fmri_spec.sess(iRun).regress      = struct('name', {}, 'val', {});
                jobs{1}.spm.stats.fmri_spec.sess(iRun).multi_reg    = cellstr(multiReg.(run));
                jobs{1}.spm.stats.fmri_spec.sess(iRun).hpf          = HPcutoff.(run);
            end
            jobs{1}.spm.stats.fmri_spec.fact                = struct('name', {}, 'levels', {});
            jobs{1}.spm.stats.fmri_spec.volt                = 1;
            jobs{1}.spm.stats.fmri_spec.global              = 'None';
            jobs{1}.spm.stats.fmri_spec.mthresh             = 0.8;
            jobs{1}.spm.stats.fmri_spec.mask                = {''};
            jobs{1}.spm.stats.fmri_spec.cvi                 = 'AR(1)';
            if strcmp(settings.hrf_response,'canonical_hrf')
                jobs{1}.spm.stats.fmri_spec.bases.hrf.derivs    = [0 0];
            elseif strcmp(settings.hrf_response,'fir')
                jobs{1}.spm.stats.fmri_spec.bases.fir.order     = settings.order;
                jobs{1}.spm.stats.fmri_spec.bases.fir.length    = settings.window;
            else
                error('unknown hrf response')
            end

            jobName = sprintf('fmri_model_spec_%s_%s',denoise,settings.hrf_response);
            run_spm_jobs_MID(jobName,dirs.jobs,jobs);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MODEL ESTIMATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if settings.todo.modelEstimation

            SPMmat = fullfile(dirs.firstLevel, 'SPM.mat');

            clear jobs
            jobs{1}.spm.stats.fmri_est.spmmat           = cellstr(SPMmat);
            jobs{1}.spm.stats.fmri_est.write_residuals  = 0;
            jobs{1}.spm.stats.fmri_est.method.Classical = 1;

            jobName = sprintf('fmri_model_estimation_%s_%s',denoise,settings.hrf_response);
            run_spm_jobs_MID(jobName,dirs.jobs,jobs);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CONTRAST ESTIMATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if settings.todo.contrastEstimation

            SPMmat = fullfile(dirs.firstLevel, 'SPM.mat');

            % Contrast vectors, found by regressor name in SPM.mat:
            % one contrast per condition (canonical HRF) or one per
            % condition and time bin (FIR)
            if strcmp(settings.hrf_response,'fir')
                contrasts = define_contrasts_fir(dirs.firstLevel,SPMmat, settings);
            else
                contrasts = define_contrasts(dirs.firstLevel,SPMmat, settings);
            end

            clear jobs
            jobs{1}.spm.stats.con.spmmat = cellstr(SPMmat);

            conNames = fieldnames(contrasts);
            conNames(ismember(conNames,settings.todo.exclude.contrasts)) = [];

            for icon = 1:numel(conNames)
                jobs{1}.spm.stats.con.consess{icon}.tcon.name    = conNames{icon};
                jobs{1}.spm.stats.con.consess{icon}.tcon.weights = contrasts.(conNames{icon});
                jobs{1}.spm.stats.con.consess{icon}.tcon.sessrep = 'none';
            end
            jobs{1}.spm.stats.con.delete = 1;   % replace existing contrasts

            jobName = sprintf('contrasts_%s_%s',denoise,settings.hrf_response);
            run_spm_jobs_MID(jobName,dirs.jobs,jobs);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % RE-COMPRESS FUNC DATA
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        boldNii = spm_select('FPList',dirs.preproc,sprintf('^sub-%.3d_ses-drug%d_task-MID_.*_bold.nii$',iSubject,iSession));

        if ~isempty(boldNii)
            boldNii = cellstr(boldNii);
            for iNii = 1:numel(boldNii)
                if isunix
                    cmd = ['gzip ' boldNii{iNii}];
                    system(cmd);
                elseif ispc
                    warning('gzip in Matlab is EXTREMELY slow!')
                    gzip(boldNii{iNii});
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DRUG COMPARISONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Difference images (drug1 - drug2) of every contrast. These enter the
    % group-level one-sample t-tests of the drug effects.

    if settings.todo.drugComparisons

        dirs.firstLevel = fullfile(settings.io.derivMRIdir,'firstLevel','MID',sprintf('sub-%.3d',iSubject));

        if settings.todo.parameterModulator
            dirs.output     = fullfile(dirs.firstLevel,'drugComparisons',denoise,settings.hrf_response, 'surprise_modulation');
        else
            dirs.output     = fullfile(dirs.firstLevel,'drugComparisons',denoise,settings.hrf_response);
        end
        dirs.jobs       = fullfile(dirs.firstLevel,'drugComparisons','jobs');

        dirNames = fieldnames(dirs);
        for iDir = 1:numel(dirNames)
            if ~exist(dirs.(dirNames{iDir}),'dir')
                mkdir(dirs.(dirNames{iDir}));
            end
        end

        % Session number of each drug for this subject
        drugcodes = readtable(settings.drugcodes);

        drugComps = {'MPHvsPBO','SULvsPBO','MPHvsSUL'};

        for iComp = 1:numel(drugComps)

            comparison = drugComps{iComp};
            drugs      = regexp(comparison,'vs','split');

            sessdrug1 = drugcodes.(sprintf('%ssession',drugs{1}))(drugcodes.subject==iSubject);
            sessdrug2 = drugcodes.(sprintf('%ssession',drugs{2}))(drugcodes.subject==iSubject);

            if settings.todo.parameterModulator
                spmDir1  = fullfile(dirs.firstLevel,sprintf('ses-drug%d',sessdrug1),denoise,settings.hrf_response, 'surprise_modulation');
            else
                spmDir1  = fullfile(dirs.firstLevel,sprintf('ses-drug%d',sessdrug1),denoise,settings.hrf_response);
            end

            SPMmat1  = fullfile(spmDir1,'SPM.mat');

            if settings.todo.parameterModulator
                spmDir2  = fullfile(dirs.firstLevel,sprintf('ses-drug%d',sessdrug2),denoise,settings.hrf_response, 'surprise_modulation');
            else
                spmDir2  = fullfile(dirs.firstLevel,sprintf('ses-drug%d',sessdrug2),denoise,settings.hrf_response);
            end
            SPMmat2  = fullfile(spmDir2,'SPM.mat');

            if exist(SPMmat1,'file') && exist(SPMmat2,'file')
                SPMdrug1 = load(SPMmat1);
                SPMdrug2 = load(SPMmat2);
            else
                warning('One or both SPMs not found for drug comparison %s. Skipping.',comparison)
                continue
            end

            assert(numel(SPMdrug1.SPM.xCon) == numel(SPMdrug2.SPM.xCon), ...
                'Different number of contrasts found for the drugs in comparison %s', comparison)

            % One imcalc job per contrast: drug1 - drug2
            % -------------------------------------------------------------
            clear jobs
            for icon = 1:numel(SPMdrug1.SPM.xCon)

                conImg_drug1 = fullfile(spmDir1, SPMdrug1.SPM.xCon(icon).Vcon.fname);
                conImg_drug2 = fullfile(spmDir2, SPMdrug2.SPM.xCon(icon).Vcon.fname);

                jobs{icon}.spm.util.imcalc.input = {conImg_drug1;conImg_drug2};

                conNm = SPMdrug1.SPM.xCon(icon).name;

                jobs{icon}.spm.util.imcalc.output = [comparison '_' conNm];
                jobs{icon}.spm.util.imcalc.outdir = cellstr(dirs.output);
                jobs{icon}.spm.util.imcalc.expression = 'i1 - i2';

                jobs{icon}.spm.util.imcalc.var              = struct('name', {}, 'value', {});
                jobs{icon}.spm.util.imcalc.options.dmtx     = 0;
                jobs{icon}.spm.util.imcalc.options.mask     = 0;
                jobs{icon}.spm.util.imcalc.options.interp   = 0;
                jobs{icon}.spm.util.imcalc.options.dtype    = 4;
            end

            jobName = sprintf('drugComp_imCalc_%s_%s',denoise,settings.hrf_response);
            run_spm_jobs_MID(jobName,dirs.jobs,jobs);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AVERAGE CONTRASTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mean of every contrast over the available sessions. These enter the
    % group-level test of the task effects across drug sessions.

    if settings.todo.averageContrasts

        dirs.firstLevel = fullfile(settings.io.derivMRIdir,'firstLevel','MID',sprintf('sub-%.3d',iSubject));

        if settings.todo.parameterModulator
            dirs.output     = fullfile(dirs.firstLevel,'average',denoise,settings.hrf_response, 'surprise_modulation');
        else
            dirs.output     = fullfile(dirs.firstLevel,'average',denoise,settings.hrf_response);
        end

        dirs.jobs       = fullfile(dirs.firstLevel,'average','jobs');

        dirNames = fieldnames(dirs);
        for iDir = 1:numel(dirNames)
            if ~exist(dirs.(dirNames{iDir}),'dir')
                mkdir(dirs.(dirNames{iDir}));
            end
        end

        % SPMs of the sessions that exist for this subject
        % -----------------------------------------------------------------
        SPMmats = {};
        for ises = 1:numel(settings.sessions)
            iSession = settings.sessions(ises);
            SPMmat   = fullfile(dirs.firstLevel,sprintf('ses-drug%d',iSession),denoise,settings.hrf_response,'SPM.mat');
            if exist(SPMmat,'file')
                SPMmats{end+1,1} = SPMmat;
            end
        end

        if numel(SPMmats) < 2
            warning('Not enough sessions available for sub-%.3d to create average contrasts. Skipping',iSubject)
            continue
        end

        for ispm = 1:numel(SPMmats)
            SPMs.(sprintf('spm%d',ispm)) = load(SPMmats{ispm});
        end

        for ispm = 1:numel(SPMmats)
            ncons(ispm) = numel(SPMs.(sprintf('spm%d',ispm)).SPM.xCon);
        end
        assert(numel(unique(ncons)) == 1, 'Different number of contrasts found in SPMs to average')

        % imcalc expression: (i1+i2+...)/number of sessions
        % =================================================================
        expression = '(i1';
        for i = 2:numel(SPMmats)
            expression = sprintf('%s+i%d',expression,i);
        end
        expression = sprintf('%s)/%d',expression,numel(SPMmats));

        clear jobs
        for icon = 1:unique(ncons)

            for iImg = 1:numel(SPMmats)
                [p,~,~]         = fileparts(SPMmats{iImg});
                conImgs{iImg,1} = fullfile(p,sprintf('con_%.4d.nii',icon));
            end

            outputname = ['average_' SPMs.spm1.SPM.xCon(icon).name];

            jobs{icon}.spm.util.imcalc.input            = conImgs;
            jobs{icon}.spm.util.imcalc.output           = outputname;
            jobs{icon}.spm.util.imcalc.outdir           = cellstr(dirs.output);
            jobs{icon}.spm.util.imcalc.expression       = expression;
            jobs{icon}.spm.util.imcalc.var              = struct('name', {}, 'value', {});
            jobs{icon}.spm.util.imcalc.options.dmtx     = 0;
            jobs{icon}.spm.util.imcalc.options.mask     = 0;
            jobs{icon}.spm.util.imcalc.options.interp   = 0;
            jobs{icon}.spm.util.imcalc.options.dtype    = 16;
        end

        jobName = 'averageContrasts';
        run_spm_jobs_MID(jobName,dirs.jobs,jobs);
    end
    iter1 = false;
end
end
