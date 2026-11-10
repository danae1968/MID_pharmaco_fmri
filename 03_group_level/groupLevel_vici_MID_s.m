function groupLevel_vici_MID_s(settings,contrast,denoise)
%GROUPLEVEL_VICI_MID_S  Group-level SPM models for one MID task contrast.
%
% For one task contrast, loops over the requested session average / drug
% comparisons / single drugs (settings.todo.drugs.do) and for each:
%   1. collects the first-level contrast images of all included subjects,
%      applying the session-level exclusions
%   2. adds the (mean-centred) covariate, if one is selected
%   3. specifies and estimates the group model:
%        canonical HRF : one-sample t-test (or paired t-test)
%        FIR           : Subject x TimeBin flexible factorial design
%   4. defines the group contrasts (positive and negative task effect, and
%      covariate effect if applicable)
%   5. exports binary significant-voxel masks (create_significant_voxels_dp.m)
%   6. optionally runs small-volume correction (extract_svc_for_analysis.m)
%
% INPUTS
% - settings : structure from run_groupLevel_vici_MID.m (paths, subjects,
%              exclusions, model options)
% - contrast : char; name of the first-level contrast, e.g. 'Target'
% - denoise  : char; denoising level of the first-level data
% .........................................................................
% Danae Papadopetraki, 2019-2026
% Initial structure based on a lab template by Ruben van den Bosch
%

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

addpath(settings.io.functionsDir)
addpath(settings.io.codeDir)


% A parallel job starts from a fresh MATLAB session
if settings.parallel
    addpath /home/common/matlab/spm12/
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');
end

% First-level derivatives live in firstLevel/MID/sub-XXX/ on the cluster.
% A local copy may instead use firstLevel/<hrfType>/sub-XXX/. Detect which
% layout is present so the same code runs against either.
if exist(fullfile(settings.io.derivMRIdir,'firstLevel','MID'),'dir')
    firstLevelTaskDir = 'MID';
elseif exist(fullfile(settings.io.derivMRIdir,'firstLevel',settings.todo.hrfType),'dir')
    firstLevelTaskDir = settings.todo.hrfType;
else
    error(['Could not find firstLevel' filesep 'MID or firstLevel' filesep '%s under %s'], ...
        settings.todo.hrfType, settings.io.derivMRIdir);
end

set(0, 'DefaultFigureVisible', 'off');

% Directories
% -------------------------------------------------------------------------
dirs.groupLevel = settings.io.groupLevelMIDdir;
dirs.jobs       = fullfile(dirs.groupLevel,'jobs');

dirNames = fieldnames(dirs);
for iDir = 1:numel(dirNames)
    if ~exist(dirs.(dirNames{iDir}),'dir')
        mkdir(dirs.(dirNames{iDir}));
    end
end

% Session number of each drug per subject
% -------------------------------------------------------------------------
drugcodes = readtable(settings.todo.drugs.drugcodes);

% FIR with a factorial design?
% -------------------------------------------------------------------------
if strcmpi(settings.todo.drugs.comptest, 'factorial') & strcmpi(settings.todo.hrfType, 'fir')
    isFactorialFir = 1;
else
    isFactorialFir = 0;
end

conNm = contrast;


% Loop over session average / drug comparisons / single drugs
% =========================================================================
for idrug = 1:numel(settings.todo.drugs.do)

    % 'average'  -> session average
    % 'XvsY'     -> drug comparison
    % otherwise  -> single drug session
    % ---------------------------------------------------------------------
    if idrug <= numel(settings.todo.drugs.do)

        drug = settings.todo.drugs.do{idrug};
        if strcmpi(drug, 'average')
            average  = true;
            drugComp = false;
        else
            average = false;
            drugComp = contains(drug,'vs');
        end
    end

    fprintf('Running group level analysis for %s contrasts (denoise: %s)\n',drug,denoise)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Design specification and Estimation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if settings.todo.designEstimation

        % Output directory
        % -----------------------------------------------------------------
        if average
            outDir = fullfile(settings.io.mainOutputDir,['sessionAverage_' conNm]);
        elseif drugComp && strcmpi(settings.todo.drugs.comptest,'paired')
            outDir = fullfile(settings.io.mainOutputDir,'paired_ttest',[drug '_' conNm]);
        else
            outDir = fullfile(settings.io.mainOutputDir,[drug '_' conNm]);
        end


        if ~exist(outDir,'dir')
            mkdir(outDir)
        end

        % Index of this contrast in the first-level SPM.mat
        % -----------------------------------------------------------------
        % Read from the first subject's first available session
        if settings.todo.parameterModulator
            firstLevelDir = fullfile(settings.io.derivMRIdir,'firstLevel',firstLevelTaskDir, ...
                sprintf('sub-%.3d',settings.subjectIx(1)),'ses-drug1',denoise,settings.todo.hrfType, 'surprise_modulation');
        else
            firstLevelDir = fullfile(settings.io.derivMRIdir,'firstLevel',firstLevelTaskDir, ...
                sprintf('sub-%.3d',settings.subjectIx(1)),'ses-drug1',denoise,settings.todo.hrfType);
        end


        SPMfile = fullfile(firstLevelDir,'SPM.mat');
        if exist(SPMfile,'file')
            load(SPMfile,'SPM')
        elseif exist(strrep(SPMfile, 'ses-drug1', 'ses-drug2'),'file')
            SPMfile = strrep(SPMfile, 'ses-drug1', 'ses-drug2');
            load(SPMfile,'SPM')
        elseif exist(strrep(SPMfile, 'ses-drug1', 'ses-drug3'),'file')
            SPMfile = strrep(SPMfile, 'ses-drug1', 'ses-drug3');
            load(SPMfile,'SPM')
        else
            error('No SPM.mat file found for first subject (sub-%.3d). Needed to find contrast index numbers',settings.subjectIx(1))
        end

        if isFactorialFir
            % FIR: one first-level contrast per time bin, '<con>_binNN'
            Nbins = settings.todo.fir.Nbins;
            conNames = {SPM.xCon.name};
            bin_indices = find(contains(conNames,sprintf('%s_bin', conNm)));
            bin_cons = conNames(bin_indices);

        else
            for ix = 1:numel(SPM.xCon)
                if strcmpi(SPM.xCon(ix).name,conNm)
                    conIx = ix;
                    break
                end
            end
        end


        if ~isFactorialFir

            % Collect one image per subject
            % -----------------------------------------------------------------
            con_imgs.subs = [];
            con_imgs.imgs = '';
            for ix = 1:numel(settings.subjectIx)

                iSubject = settings.subjectIx(ix);

                % Single drug: contrast image of that drug's session
                if ~drugComp && ~average

                    sessNr = drugcodes.(sprintf('%ssession',upper(drug)))(drugcodes.subject == iSubject);

                    % Skip excluded sessions
                    if any(settings.exclude.subject==iSubject)
                        if ismember(sessNr,settings.exclude.session{settings.exclude.subject==iSubject})
                            continue
                        end
                    end

                    if settings.todo.parameterModulator
                        sub_firstLevel_dir = fullfile(settings.io.derivMRIdir,'firstLevel',firstLevelTaskDir,sprintf('sub-%.3d',iSubject), ...
                            sprintf('ses-drug%d',sessNr),denoise, settings.todo.hrfType, 'surprise_modulation');
                    else

                        sub_firstLevel_dir = fullfile(settings.io.derivMRIdir,'firstLevel',firstLevelTaskDir,sprintf('sub-%.3d',iSubject), ...
                            sprintf('ses-drug%d',sessNr),denoise, settings.todo.hrfType);
                    end

                    conImg = fullfile(sub_firstLevel_dir,sprintf('con_%.4d.nii',conIx));


                % Drug comparison: first-level difference image
                elseif drugComp && ~average

                    % Skip the subject if either session is excluded
                    % .................................................
                    drugs = strsplit(drug,'vs');

                    sessNrs(1,1) = drugcodes.(sprintf('%ssession',upper(drugs{1})))(drugcodes.subject == iSubject);
                    sessNrs(1,2) = drugcodes.(sprintf('%ssession',upper(drugs{2})))(drugcodes.subject == iSubject);

                    if any(settings.exclude.subject==iSubject)
                        if ismember(sessNrs(1,1),settings.exclude.session{settings.exclude.subject==iSubject}) || ...
                                ismember(sessNrs(1,2),settings.exclude.session{settings.exclude.subject==iSubject})
                            continue
                        end
                    end

                    if settings.todo.parameterModulator
                        drugCompDir = fullfile(settings.io.derivMRIdir,'firstLevel',firstLevelTaskDir, ...
                            sprintf('sub-%.3d',iSubject),'drugComparisons',denoise,settings.todo.hrfType, 'surprise_modulation');

                    else
                        drugCompDir = fullfile(settings.io.derivMRIdir,'firstLevel',firstLevelTaskDir, ...
                            sprintf('sub-%.3d',iSubject),'drugComparisons',denoise,settings.todo.hrfType);
                    end
                    if strcmpi(settings.todo.analysisType,'ppi')
                        drugCompDir = fullfile(drugCompDir,'ppi',sprintf('PPI_%s',settings.todo.ppi.seed));
                    end

                    conImg = fullfile(drugCompDir,[drug '_' conNm '.nii']);

                % Session average: first-level average image
                elseif average

                    averageDir = fullfile(settings.io.derivMRIdir,'firstLevel',firstLevelTaskDir,sprintf('sub-%.3d',iSubject),'average',denoise,settings.todo.hrfType);

                    conImg = fullfile(averageDir,['average_' conNm '.nii']);
                end


                % Keep the image if it exists
                % .............................................................
                if ~(exist(conImg,'file') == 2)
                    warning('Cannot find contrast image for contrast %s of subject %d. Skipping.',conNm,iSubject)
                    continue
                end
                con_imgs.subs = [con_imgs.subs;iSubject];
                con_imgs.imgs = strvcat(con_imgs.imgs,conImg);
            end

            con_imgs.imgs = cellstr(con_imgs.imgs);

        elseif isFactorialFir

            % Collect one image per subject and time bin (Nsubjects x Nbins)
            con_imgs.subs = [];
            con_imgs.imgs = {};

            Nbins = settings.todo.fir.Nbins;
            for ix = 1:numel(settings.subjectIx)

                iSubject = settings.subjectIx(ix);

                if  ~drugComp && ~average % single drug session

                    sessNr = drugcodes.(sprintf('%ssession',upper(drug)))(drugcodes.subject == iSubject);

                    if any(settings.exclude.subject==iSubject)
                        if ismember(sessNr,settings.exclude.session{settings.exclude.subject==iSubject})
                            continue
                        end
                    end

                    sub_firstLevel_dir = fullfile(settings.io.derivMRIdir,'firstLevel',firstLevelTaskDir,sprintf('sub-%.3d',iSubject), ...
                        sprintf('ses-drug%d',sessNr),denoise, settings.todo.hrfType);

                    bin_paths = cell(1, Nbins);
                    for bin = 1:Nbins
                        ind = bin_indices(bin);
                        bin_paths{bin} = fullfile(sub_firstLevel_dir,  sprintf('con_%04d.nii', ind));
                    end


                elseif drugComp && ~average % drug comparison

                    % Per-bin first-level difference images (drug1 - drug2)
                    drugComparison = settings.todo.drugs.do{idrug};
                    drugsInComparison = strsplit(drugComparison, 'vs');

                    if settings.todo.excludeSubjects
                        sessNrs(1,1) = drugcodes.(sprintf('%ssession',upper(drugsInComparison{1})))(drugcodes.subject == iSubject);
                        sessNrs(1,2) = drugcodes.(sprintf('%ssession',upper(drugsInComparison{2})))(drugcodes.subject == iSubject);

                        if any(settings.exclude.subject==iSubject)
                            if ismember(sessNrs(1,1),settings.exclude.session{settings.exclude.subject==iSubject}) || ...
                                    ismember(sessNrs(1,2),settings.exclude.session{settings.exclude.subject==iSubject})
                                continue
                            end
                        end
                    end

                    drugCompDir = fullfile(settings.io.derivMRIdir,'firstLevel',firstLevelTaskDir, ...
                        sprintf('sub-%.3d',iSubject),'drugComparisons',denoise,settings.todo.hrfType);

                    bin_paths = cell(1, Nbins);
                    skip_subject = 0;

                    for bin = 1:Nbins
                        % e.g. MPHvsPBO_Target_bin01.nii
                        bin_file = fullfile(drugCompDir, sprintf('%s_%s_bin%02d.nii', drugComparison, conNm, bin));

                        if ~exist(bin_file, 'file')
                            warning('Cannot find bin file for %s subject %d bin %d. Skipping subject.', ...
                                drugComparison, iSubject, bin);
                            skip_subject = 1;
                            break;
                        end

                        bin_paths{bin} = bin_file;
                    end

                    if skip_subject
                        continue
                    end

                elseif average

                    averageDir = fullfile(settings.io.derivMRIdir,'firstLevel',firstLevelTaskDir,sprintf('sub-%.3d',iSubject),'average',denoise,settings.todo.hrfType);

                    bin_paths = cell(1, Nbins);

                    % Average images are named by bin, not contrast index
                    for bin = 1:Nbins
                        bin_paths{bin} = fullfile(averageDir, sprintf('average_%s_bin%02d.nii', conNm, bin));
                    end

                end

                % Keep the subject if its images exist
                % .............................................................
                if ~(exist(bin_paths{1},'file') == 2)
                    warning('Cannot find contrast image for contrast %s of subject %d. Skipping.',conNm,iSubject)
                    continue
                end

                con_imgs.subs = [con_imgs.subs; iSubject];
                con_imgs.imgs = [con_imgs.imgs; bin_paths];
            end
        end

        % Covariate vector(s), if any
        % -----------------------------------------------------------------
        clear covariates

        covnms = fieldnames(settings.todo.covariate);

        for icov = 1:numel(covnms)

            cov = covnms{icov};

            if ~settings.todo.covariate.(cov).use
                continue
            end

            covdata = readtable(settings.todo.covariate.(cov).datafile, 'TreatAsEmpty','NA');

            % ki file has one row per subject x ROI: keep the requested ROI
            if startsWith(cov,'ki')
                covdata = covdata(strcmp(covdata.ROI,settings.todo.covariate.(cov).roi),:);
            end

            % First column must hold the subject numbers
            if strcmp(covdata.Properties.VariableNames{1},'vID')
                covdata.Properties.VariableNames{1} = 'subject';
            elseif ~strcmp(covdata.Properties.VariableNames{1},'subject')
                error('First column of covariate data (%s) is not named "subject" or "vID". Check that it contains subject numbers',cov)
            end

            % Keep only subjects with BOTH an image and a valid covariate
            % value, so images and covariate values stay aligned
            subjects_to_keep = [];
            indices_to_keep = [];

            for ix = 1:numel(con_imgs.subs)
                iSubject = con_imgs.subs(ix);

                if ~any(covdata.subject == iSubject)
                    fprintf('Removing subject %d: no covariate data\n', iSubject);
                    continue
                end

                if startsWith(cov,'ki')
                    covvalue = covdata.meanRoiContrastEstimate(covdata.subject == iSubject);
                else
                    covvalue = covdata.(settings.todo.covariate.(cov).variable)(covdata.subject == iSubject);
                end

                if isnan(covvalue)
                    fprintf('Removing subject %d: NaN covariate value\n', iSubject);
                    continue
                end

                subjects_to_keep = [subjects_to_keep; iSubject];
                indices_to_keep = [indices_to_keep; ix];
            end

            con_imgs.subs = con_imgs.subs(indices_to_keep);
            if isFactorialFir
                con_imgs.imgs = con_imgs.imgs(indices_to_keep, :);
            else
                con_imgs.imgs = con_imgs.imgs(indices_to_keep);
            end

            % Covariate vector for the kept subjects
            cov_vector = [];
            for ix = 1:numel(subjects_to_keep)
                iSubject = subjects_to_keep(ix);

                if startsWith(cov,'ki')
                    covvalue = covdata.meanRoiContrastEstimate(covdata.subject == iSubject);
                else
                    covvalue = covdata.(settings.todo.covariate.(cov).variable)(covdata.subject == iSubject);
                end

                if isnumeric(covvalue)
                    cov_vector = [cov_vector; covvalue];
                else
                    error('Non-numeric covariate value for sub-%.3d, covariate: %s', iSubject, cov)
                end
            end

            % Mean-centre here and pass to SPM with iCC = 5 (no further
            % centring), so the task effect is estimated at the mean
            % covariate value
            cov_mean = mean(cov_vector);
            cov_vector_centered = cov_vector - cov_mean;

            fprintf('Covariate %s: mean = %.4f, SD = %.4f (before centering)\n', ...
                cov, cov_mean, std(cov_vector));

            covstruct.name   = settings.todo.covariate.(cov).conName;
            covstruct.vector = cov_vector_centered;

            if ~exist('covariates','var')
                covariates = covstruct;
            else
                covariates = [covariates covstruct];
            end

            fprintf('Covariate %s: using %d subjects (mean-centered)\n', cov, numel(cov_vector));
        end

        % Paired t-test: the two drug sessions' images as pairs
        % -----------------------------------------------------------------
        if drugComp && strcmpi(settings.todo.drugs.comptest,'paired')

            pairs = [];

            for ix = 1:numel(con_imgs.subs)

                iSubject = con_imgs.subs(ix);

                drugs = strsplit(drug,'vs');

                sessNrs(1,1) = drugcodes.(sprintf('%ssession',upper(drugs{1})))(drugcodes.subject == iSubject);
                sessNrs(1,2) = drugcodes.(sprintf('%ssession',upper(drugs{2})))(drugcodes.subject == iSubject);



                if settings.todo.parameterModulator
                    conImgDirs = cellstr(char(fullfile(settings.io.derivMRIdir,'firstLevel',firstLevelTaskDir,sprintf('sub-%.3d',iSubject), ...
                        sprintf('ses-drug%d',sessNrs(1,1)),denoise,settings.todo.hrfType, 'surprise_modulation'), ...
                        fullfile(settings.io.derivMRIdir,'firstLevel',firstLevelTaskDir,sprintf('sub-%.3d',iSubject), ...
                        sprintf('ses-drug%d',sessNrs(1,2)),denoise,settings.todo.hrfType, 'surprise_modulation')));
                else
                    conImgDirs = cellstr(char(fullfile(settings.io.derivMRIdir,'firstLevel',firstLevelTaskDir,sprintf('sub-%.3d',iSubject), ...
                        sprintf('ses-drug%d',sessNrs(1,1)),denoise,settings.todo.hrfType), ...
                        fullfile(settings.io.derivMRIdir,'firstLevel',firstLevelTaskDir,sprintf('sub-%.3d',iSubject), ...
                        sprintf('ses-drug%d',sessNrs(1,2)),denoise,settings.todo.hrfType)));
                end
                conImgs = cellstr(char(fullfile(conImgDirs{1},sprintf('con_%.4d.nii',conIx)), ...
                    fullfile(conImgDirs{2},sprintf('con_%.4d.nii',conIx))));


                pairs(ix).scans = conImgs;
            end

        end

        % Fill in batch
        % -----------------------------------------------------------------
        clear jobs

        % =================================================================
        % FIR: Subject x TimeBin flexible factorial design
        % =================================================================
        if strcmpi(settings.todo.hrfType, 'fir') && strcmpi(settings.todo.drugs.comptest, 'factorial')

            Nbins = settings.todo.fir.Nbins;

            % A between-subject covariate would have to interact with the
            % TimeBin factor (cov.iCFI); this is not implemented, so stop
            % rather than build a wrong design.
            if exist('covariates', 'var') && ~isempty(covariates)
                error(['FIR factorial design does not yet support covariates. ' ...
                       'Set the relevant todo.covariate.*.use to false, or use ' ...
                       'hrfType=''canonical_hrf'' for covariate analyses.']);
            end

            % -----------------------------------------------------------------
            % Repeated-measures design: Subject x TimeBin.
            %   Subject is a blocking factor (one level per subject) and
            %   TimeBin the repeated factor. Modelling Subject explicitly lets
            %   SPM account for the within-subject correlation across bins and
            %   gives the correct residual degrees of freedom, rather than
            %   treating each subject's bins as independent scans.
            % -----------------------------------------------------------------
            jobs{1}.spm.stats.factorial_design.dir = {outDir};
            if settings.force
                if exist(fullfile(outDir, 'SPM.mat'), 'file')
                    delete(fullfile(outDir, 'SPM.mat'));
                    fprintf('Deleted existing SPM.mat in: %s\n', outDir);
                end
            end

            % Keep subjects with an image for every bin
            Nsubs = size(con_imgs.imgs, 1);
            valid_indices = [];
            for s = 1:Nsubs
                scans = con_imgs.imgs(s, :)';
                if ~any(cellfun(@isempty, scans))
                    valid_indices = [valid_indices; s]; %#ok<AGROW>
                else
                    fprintf('Excluding subject %d: missing scans\n', con_imgs.subs(s));
                end
            end
            con_imgs.subs = con_imgs.subs(valid_indices);
            con_imgs.imgs = con_imgs.imgs(valid_indices, :);
            Nsubs = numel(valid_indices);
            fprintf('Building Subject x TimeBin factorial design with %d subjects, %d bins\n', Nsubs, Nbins);

            % Flexible factorial (des.fblock), not full factorial (des.fd):
            % a full factorial always includes the Subject x TimeBin
            % interaction, which with one scan per cell leaves zero
            % residual degrees of freedom. The flexible factorial lets the
            % model contain main effects only (maininters below). Scans are
            % passed together with a factor matrix (imatrix, SPM.xX.I
            % format): column 1 replication (always 1), columns 2-3 the
            % Subject and TimeBin levels, column 4 an unused third factor.
            jobs{1}.spm.stats.factorial_design.des.fblock.fac(1).name     = 'Subject';
            jobs{1}.spm.stats.factorial_design.des.fblock.fac(1).dept     = 0;
            jobs{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
            jobs{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca    = 0;
            jobs{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova   = 0;

            % TimeBin: dependent levels (dept = 1) with a single shared
            % variance component (variance = 0). Unequal variances would
            % require estimating a full 11 x 11 covariance matrix (66 ReML
            % hyperparameters), which does not converge reliably with
            % this number of subjects.
            jobs{1}.spm.stats.factorial_design.des.fblock.fac(2).name     = 'TimeBin';
            jobs{1}.spm.stats.factorial_design.des.fblock.fac(2).dept     = 1;
            jobs{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
            jobs{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca    = 0;
            jobs{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova   = 0;

            allScans = cell(Nsubs*Nbins, 1);
            imatrix  = zeros(Nsubs*Nbins, 4);
            row = 0;
            for s = 1:Nsubs
                for bin = 1:Nbins
                    row = row + 1;
                    allScans{row}  = con_imgs.imgs{s, bin};
                    imatrix(row,:) = [1, s, bin, 1];
                end
            end
            jobs{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans   = allScans;
            jobs{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.imatrix = imatrix;

            % Additive model: main effects only, no interaction term
            jobs{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 1;  % Subject
            jobs{1}.spm.stats.factorial_design.des.fblock.maininters{2}.fmain.fnum = 2;  % TimeBin

            jobs{1}.spm.stats.factorial_design.cov       = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            jobs{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});

            % Masking / globals
            jobs{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            jobs{1}.spm.stats.factorial_design.masking.im         = 1;
            jobs{1}.spm.stats.factorial_design.masking.em         = {''};
            jobs{1}.spm.stats.factorial_design.globalc.g_omit     = 1;
            jobs{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            jobs{1}.spm.stats.factorial_design.globalm.glonorm    = 1;

            % Estimation
            jobs{2}.spm.stats.fmri_est.spmmat          = {fullfile(outDir, 'SPM.mat')};
            jobs{2}.spm.stats.fmri_est.method.Classical = 1;

            jobName = sprintf('design_estimate_factorial_%s_%s', drug, conNm);
            run_spm_jobs(jobName, dirs.jobs, jobs);
            % =================================================================
            % Canonical HRF: one-sample (or paired) t-test
            % =================================================================
        else

            jobs{1}.spm.stats.factorial_design.dir = cellstr(outDir);

            if settings.force
                if exist(fullfile(outDir, 'SPM.mat'), 'file')
                    delete(fullfile(outDir, 'SPM.mat'));
                    fprintf('Deleted existing SPM.mat in: %s\n', outDir);
                end
            end

            if drugComp && strcmpi(settings.todo.drugs.comptest,'paired')
                jobs{1}.spm.stats.factorial_design.des.pt.pair   = pairs;
                jobs{1}.spm.stats.factorial_design.des.pt.gmsca  = 0;
                jobs{1}.spm.stats.factorial_design.des.pt.ancova = 0;
            else
                jobs{1}.spm.stats.factorial_design.des.t1.scans  = con_imgs.imgs;
            end

            % Covariate(s)
            jobs{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            if exist('covariates','var')
                for icov = 1:numel(covariates)
                    % Paired t-test: one value per scan, i.e. each
                    % subject's value repeated twice
                    if drugComp && strcmpi(settings.todo.drugs.comptest,'paired')
                        jobs{1}.spm.stats.factorial_design.cov(icov).c = repelem(covariates(icov).vector, 2);
                    else
                        jobs{1}.spm.stats.factorial_design.cov(icov).c = covariates(icov).vector(:);
                    end
                    jobs{1}.spm.stats.factorial_design.cov(icov).cname = covariates(icov).name;
                    jobs{1}.spm.stats.factorial_design.cov(icov).iCFI  = 1;
                    jobs{1}.spm.stats.factorial_design.cov(icov).iCC   = 5;   % already centred


                end
            end
            jobs{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});

            jobs{1}.spm.stats.factorial_design.masking.tm.tm_none     = 1;
            jobs{1}.spm.stats.factorial_design.masking.im             = 1;
            jobs{1}.spm.stats.factorial_design.masking.em             = {''};
            jobs{1}.spm.stats.factorial_design.globalc.g_omit         = 1;
            jobs{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            jobs{1}.spm.stats.factorial_design.globalm.glonorm        = 1;

            % Estimation
            jobs{2}.spm.stats.fmri_est.spmmat         = cellstr(fullfile(outDir,'SPM.mat'));
            jobs{2}.spm.stats.fmri_est.write_residuals = 0;
            jobs{2}.spm.stats.fmri_est.method.Classical = 1;

            jobName = sprintf('design_estimate_%s_%s_%s', drug, conNm, denoise);
            run_spm_jobs(jobName, dirs.jobs, jobs);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RESULTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if settings.todo.results.groupContrast

        conNm = contrast;

        % Covariates switched on in settings (one extra contrast each)
        % -----------------------------------------------------------------
        covariate_names = fieldnames(settings.todo.covariate);

        covs = {};
        for i = 1:numel(covariate_names)
            if settings.todo.covariate.(covariate_names{i}).use
                covs = [covs covariate_names{i}];
            end
        end

        clear jobs

        % Group-level SPM.mat
        % -----------------------------------------------------------------
        if average

            SPMmat = fullfile(settings.io.mainOutputDir,['sessionAverage_' conNm],'SPM.mat');

        elseif drugComp && strcmpi(settings.todo.drugs.comptest,'paired')
            SPMmat = fullfile(settings.io.mainOutputDir,'paired_ttest', [drug '_' conNm], 'SPM.mat');
        else
            SPMmat = fullfile(settings.io.mainOutputDir,[drug '_' conNm],'SPM.mat');
        end


        % Define contrasts
        % -----------------------------------------------------------------
        if strcmpi(settings.todo.drugs.comptest, 'factorial') & strcmpi(settings.todo.hrfType, 'fir')

            Nbins = settings.todo.fir.Nbins;

            if exist('covariates', 'var') && ~isempty(covariates)
                error('FIR factorial covariate contrasts are not yet supported.');
            end

            load(SPMmat,'SPM');
            nCols = size(SPM.xX.X,2);

            % TimeBin columns, located by name (SPM names them
            % 'TimeBin_{1}' ... 'TimeBin_{11}'); the Subject columns come
            % first, so their position is not fixed.
            binCol = nan(1, Nbins);
            for b = 1:Nbins
                pattern = sprintf('TimeBin_\\{%d\\}', b);
                idx = find(~cellfun('isempty', regexp(SPM.xX.name, pattern, 'once')));
                if numel(idx) ~= 1
                    fprintf('SPM.xX.name contents:\n'); disp(SPM.xX.name');
                    error(['Expected exactly one TimeBin column for bin %d in %s, found %d. ' ...
                           'Check the TimeBin naming convention in SPM.xX.name above and update ' ...
                           'the regex pattern in this function accordingly.'], b, SPMmat, numel(idx));
                end
                binCol(b) = idx;
            end

            % Subject columns. The number of subjects is read from the
            % design, so contrasts can be (re)computed without re-running
            % the estimation step.
            subjCol = find(~cellfun('isempty', regexp(SPM.xX.name, 'Subject_\{\d+\}', 'once')));
            assert(~isempty(subjCol), ...
                'Expected Subject_{N} columns in %s, found none. Check SPM.xX.name naming.', SPMmat);
            NsubsLocal = numel(subjCol);

            clear matlabbatch
            matlabbatch{1}.spm.stats.con.spmmat = cellstr(SPMmat);
            conIdx = 1;

            % ---- (1) Mean over bins (T-contrast)
            %
            % Subject and TimeBin are both fully dummy-coded without a
            % separate intercept, so the sum of the TimeBin columns equals
            % the sum of the Subject columns. A contrast on the TimeBin
            % columns alone is therefore not estimable. Weighting the
            % Subject columns by 1/Nsubs as well makes it estimable and
            % gives the average response over the whole peristimulus
            % window -- the FIR analogue of the canonical-HRF contrast.
            mean_bins = zeros(1, nCols);
            mean_bins(binCol)  = 1/Nbins;
            mean_bins(subjCol) = 1/NsubsLocal;

            matlabbatch{1}.spm.stats.con.consess{conIdx}.tcon.name     = ...
                sprintf('Mean_%s', conNm);
            matlabbatch{1}.spm.stats.con.consess{conIdx}.tcon.weights  = mean_bins;
            matlabbatch{1}.spm.stats.con.consess{conIdx}.tcon.sessrep  = 'none';
            conIdx = conIdx + 1;

            % ---- (2) Negative mean (deactivations)
            matlabbatch{1}.spm.stats.con.consess{conIdx}.tcon.name     = ...
                sprintf('negative_Mean_%s', conNm);
            matlabbatch{1}.spm.stats.con.consess{conIdx}.tcon.weights  = -mean_bins;
            matlabbatch{1}.spm.stats.con.consess{conIdx}.tcon.sessrep  = 'none';
            conIdx = conIdx + 1;

            % ---- (3) Omnibus TimeBin effect (F-contrast): does the
            % response change anywhere over the time course? Rows are
            % differences between consecutive bins; they sum to zero
            % within the TimeBin block and are estimable as they are.
            timebin_F = zeros(Nbins-1, nCols);
            for b = 1:(Nbins-1)
                timebin_F(b, binCol(b))   =  1;
                timebin_F(b, binCol(b+1)) = -1;
            end

            matlabbatch{1}.spm.stats.con.consess{conIdx}.fcon.name     = ...
                sprintf('TimeBin_effect_%s', conNm);
            matlabbatch{1}.spm.stats.con.consess{conIdx}.fcon.weights  = timebin_F;
            matlabbatch{1}.spm.stats.con.consess{conIdx}.fcon.sessrep  = 'none';
            conIdx = conIdx + 1;

            matlabbatch{1}.spm.stats.con.delete = 1;

            spm_jobman('initcfg');
            spm_jobman('run', matlabbatch);
        elseif ~isFactorialFir
            % Task effect (intercept)
            jobs{1}.spm.stats.con.spmmat                     = cellstr(SPMmat);
            jobs{1}.spm.stats.con.consess{1}.tcon.name       = conNm;
            jobs{1}.spm.stats.con.consess{1}.tcon.weights    = zeros(1,numel(covs)+1);
            jobs{1}.spm.stats.con.consess{1}.tcon.weights(1) = 1;
            jobs{1}.spm.stats.con.consess{1}.tcon.sessrep    = 'none';

            % Covariate effect(s): '<con> x <covariate>'
            for ix = 1:numel(covs)
                jobs{1}.spm.stats.con.consess{ix+1}.tcon.name          = [conNm ' x ' settings.todo.covariate.(covs{ix}).conName];
                jobs{1}.spm.stats.con.consess{ix+1}.tcon.weights       = zeros(1,numel(covs)+1);
                jobs{1}.spm.stats.con.consess{ix+1}.tcon.weights(ix+1) = 1;
                jobs{1}.spm.stats.con.consess{ix+1}.tcon.sessrep       = 'none';
            end

            jobs{1}.spm.stats.con.delete = 1;


            jobName = 'groupContrast';
            run_spm_jobs(jobName,dirs.jobs,jobs);
        end
    end
    % Export binary masks of significant voxels for each T-contrast.
    % Negative contrasts are added automatically where missing.
    % =====================================================================
    if settings.todo.results.export.do

        if average

            SPMmat = fullfile(settings.io.mainOutputDir,['sessionAverage_' conNm],'SPM.mat');

        elseif drugComp && strcmpi(settings.todo.drugs.comptest,'paired')
            SPMmat = fullfile(settings.io.mainOutputDir,'paired_ttest', [drug '_' conNm], 'SPM.mat');
        else
            SPMmat = fullfile(settings.io.mainOutputDir,[drug '_' conNm],'SPM.mat');
        end


        load(SPMmat,'SPM')
        cons = {SPM.xCon(:).name}';

        % T-contrasts only
        conStats = {SPM.xCon.STAT}';
        isTcontrast = strcmpi(conStats, 'T');
        cons = cons(isTcontrast);
        create_significant_voxels_dp(SPMmat,cons,'mri',settings.todo.results.export,true);

    end
% =========================================================================
% SMALL VOLUME CORRECTION (optional)
% =========================================================================
    % Uses the existing SPM.mat. The functional mask is the session-average
    % p < .001 mask of the same task contrast (no-covariate model).
    % =====================================================================
    if settings.todo.svc.do

        if strcmp(settings.todo.svc.mask_type, 'functional') && ~drugComp
            fprintf('\n=== Skipping SVC: functional mask requires drug comparison ===\n');
            continue
        end

        if strcmp(settings.todo.svc.mask_type, 'functional') && drugComp
             if strcmpi(settings.todo.drugs.comptest, 'factorial') & strcmpi(settings.todo.hrfType, 'fir')
                                    mask_file = fullfile(settings.io.groupLevelMIDdir, denoise, ...
                    settings.todo.hrfType, 'no_covariate_with_exc', ...
                    sprintf('sessionAverage_%s', conNm), ...
                    sprintf('significant_voxels_Mean_%s_uncorrected_p001.nii', conNm));
                else
               mask_file = fullfile(settings.io.groupLevelMIDdir, denoise, ...
    settings.todo.hrfType, 'no_covariate_with_exc', ...
                    sprintf('sessionAverage_%s', conNm), ...
                    sprintf('significant_voxels_%s_uncorrected_p001.nii', conNm));
             end
        elseif  strcmp(settings.todo.svc.mask_type, 'anatomical')
               mask_file = settings.todo.svc.mask_file;
        end


        fprintf('\n=== Extracting SVC Results ===\n');

        analysisDir = fullfile(settings.io.mainOutputDir, [drug '_' conNm]);

        SPMmat = fullfile(analysisDir, 'SPM.mat');

        if ~exist(SPMmat, 'file')
            fprintf('=== Skipping SVC: SPM.mat not found: %s ===\n', SPMmat);
            continue
        end
        if ~exist(mask_file, 'file')
            fprintf('=== Skipping SVC: mask not found: %s ===\n', mask_file);
            continue
        end

        try
            extract_svc_for_analysis(SPMmat, mask_file, settings.todo.svc.contrasts, settings.todo.svc.threshold, settings.todo.svc.extent);
            fprintf('SVC results extracted successfully\n');
        catch ME
            warning('SVC extraction failed: %s', ME.message);
        end
    end
end
end

function run_spm_jobs(jobName,jobDir,jobs)
% RUN_SPM_JOBS(jobName,jobDir,jobs)
%
% Save an SPM batch to jobDir and run it.

jobFile = fullfile(jobDir,[jobName,'_',datestr(now,'yyyymmddTHHMMSS'),'.m']);

jobId = cfg_util('initjob', jobs);

sts = cfg_util('isjob_id', jobId);
if sts
    cfg_util('savejob', jobId, jobFile);
    cfg_util('run', jobId);
else
    error('Error in initialising %s job.',jobName)
end
end
