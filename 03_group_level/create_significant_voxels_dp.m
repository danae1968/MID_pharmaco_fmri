function create_significant_voxels_dp(SPMmat, contrast, modality, sigOptions, varargin)
% CREATE_SIGNIFICANT_VOXELS_DP  Export binary masks of significant voxels.
%
% create_significant_voxels_dp(SPMmat, contrast, modality, sigOptions, [isInitCfg])
%
% For each requested contrast in SPMmat, thresholds the statistic image
% (SPM results, parametric/RFT) and writes binary masks:
%   spmT_NNNN_significant_voxels_<con>_<type>_p<thr>.nii   one direction
%   significant_voxels_<con>_<type>_p<thr>.nii             both directions
% If a T-contrast has no negative counterpart ('negative_<con>') and the
% design is not factorial, the negative contrast is added to SPM.mat first.
%
% INPUTS
%   SPMmat     : path to the group-level SPM.mat
%   contrast   : cellstr of contrast names to export
%   modality   : 'mri' / 'fmri' (sets SPM defaults)
%   sigOptions : .threshold      p-value, or cell array of p-values
%                .thresholdType  'uncorrected' | 'fwe' (or cell array,
%                                one per threshold)
%                .extent         minimum cluster size (voxels)
%                .svc_prefix     optional; true prefixes outputs with 'SVC_'
%   isInitCfg  : true if spm_jobman is already initialised
%
% Handles T- and F-contrasts; F-contrasts are exported in one direction.

if ~isempty(varargin)
    assert(numel(varargin) == 1, 'Too many input arguments');
    isInitCfg = varargin{1};
else
    isInitCfg = false;
end

clear jobs

if strcmpi(modality,'MRI')
    modality = 'FMRI';
end
spm('defaults', upper(modality));

if ~isInitCfg
    spm_jobman('initcfg');
end

load(SPMmat)
conNms = {SPM.xCon(:).name}';
conNms = strrep(conNms,' ','_');
contrast = strrep(contrast,' ','_');

conCount = 1;
jobCount = 1;

% Optional 'SVC_' prefix on the output names
if isfield(sigOptions, 'svc_prefix') && sigOptions.svc_prefix
    filePrefix = 'SVC_';
else
    filePrefix = '';
end

% Allow several thresholds in one call
if ~iscell(sigOptions.threshold)
    sigOptions.threshold = {sigOptions.threshold};
end
if ~iscell(sigOptions.thresholdType)
    sigOptions.thresholdType = repmat({sigOptions.thresholdType}, size(sigOptions.threshold));
end
assert(numel(sigOptions.threshold) == numel(sigOptions.thresholdType), ...
    'threshold and thresholdType must match in length');

% Factorial designs define their negative contrasts explicitly
isFactorial = isfield(SPM.xX, 'iB') && ~isempty(SPM.xX.iB);

for iCon = 1:numel(conNms)
    conName = conNms{iCon};

    if ~ismember(cellstr(conName),contrast)
        continue
    end

    % A negative contrast is handled together with its positive one
    if startsWith(conName,'negative_') && ismember(cellstr(strrep(conName,'negative_','')),contrast)
        continue
    elseif startsWith(conName,'inverse_') && ismember(cellstr(strrep(conName,'inverse_','')),contrast)
        continue
    elseif startsWith(conName,'Negative_') && ismember(cellstr(strrep(conName,'Negative_','')),contrast)
        continue
    end

    % Find, or create, the negative counterpart of a T-contrast
    negPrefix = [];
    negConName = [];

    if ~strcmpi(SPM.xCon(iCon).STAT, 'F')

        if ismember(['negative_' conName], conNms)
            negPrefix = 'negative_';
            negConName = ['negative_' conName];
        elseif ismember(['inverse_' conName], conNms)
            negPrefix = 'inverse_';
            negConName = ['inverse_' conName];
        elseif ismember(['Negative_' conName], conNms)
            negPrefix = 'Negative_';
            negConName = ['Negative_' conName];

        % Not present and not a factorial design: add it to SPM.mat
        elseif ~isFactorial
            clear newConJob
            newConJob{1}.spm.stats.con.spmmat = cellstr(SPMmat);
            negPrefix = 'negative_';
            negConName = ['negative_' conName];
            newConJob{1}.spm.stats.con.consess{1}.tcon.name = [negPrefix SPM.xCon(iCon).name];

            c_in = full(SPM.xCon(iCon).c);
            c_in = c_in(:);

            % Negated weights, projected onto the estimable space
            xNeg = spm_FcUtil('Set', [negPrefix SPM.xCon(iCon).name], ...
                              'T','c', -c_in, SPM.xX.xKXs);
            wNeg = full(xNeg.c)';
            newConJob{1}.spm.stats.con.consess{1}.tcon.weights = wNeg;

            newConJob{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            newConJob{1}.spm.stats.con.delete = 0;

            run_spm_jobs('create_negative_contrast', newConJob);
            load(SPMmat);
        end
    end

    % One results/export job per threshold and direction
    for ip = 1:numel(sigOptions.threshold)
        thisThresh = sigOptions.threshold{ip};
        thisType   = sigOptions.thresholdType{ip};

        pstr = regexp(num2str(thisThresh), '\.', 'split');
        p = pstr{2};

        % F-contrasts and T-contrasts without a negative: one direction;
        % T-contrasts with a negative: both
        loop_range = 1;
        if ~strcmpi(SPM.xCon(iCon).STAT, 'F') && ~isempty(negConName)
            loop_range = 1:2;
        end

        for ix = loop_range
            if ix == 1
                name = conName;
                outname = conName;
            else
                name = negConName;
                outname = negConName;
            end

            jobs{jobCount}.spm.stats.results.spmmat = cellstr(SPMmat);
            jobs{jobCount}.spm.stats.results.conspec(1).titlestr = '';

            for conIx = 1:numel(SPM.xCon)
                if strcmp(name, strrep(SPM.xCon(conIx).name, ' ', '_'))
                    jobs{jobCount}.spm.stats.results.conspec(1).contrasts = conIx;
                    break
                end
            end

            if strcmpi(thisType, 'uncorrected')
                jobs{jobCount}.spm.stats.results.conspec(1).threshdesc = 'none';
            elseif strcmpi(thisType, 'fwe')
                jobs{jobCount}.spm.stats.results.conspec(1).threshdesc = 'FWE';
            end
            jobs{jobCount}.spm.stats.results.conspec(1).thresh = thisThresh;
            jobs{jobCount}.spm.stats.results.conspec(1).extent = sigOptions.extent;
            jobs{jobCount}.spm.stats.results.conspec(1).mask.none = 1;

            basename = sprintf('%ssignificant_voxels_%s_%s_p%s', filePrefix, outname, thisType, p);
            fprintf("Creating significant voxels for %s\n", basename)

            jobs{jobCount}.spm.stats.results.export{1}.binary.basename = basename;

            % SPM prefixes the exported file with the stat image name
            [outDir,~,~] = fileparts(SPMmat);
            if ix == 1
                pos = fullfile(outDir,sprintf('spmT_%.4d_%s.nii', conIx, basename));
            elseif ix == 2
                neg = fullfile(outDir,sprintf('spmT_%.4d_%s.nii', conIx, basename));
            end

            jobs{jobCount}.spm.stats.results.units = 1;
            jobCount = jobCount + 1;
        end

        % Union of both directions (used for figure contours)
        if ~strcmpi(SPM.xCon(iCon).STAT, 'F') && loop_range(end) == 2
            outputname = sprintf('significant_voxels_%s_%s_p%s', conName, thisType, p);
            jobs{jobCount}.spm.util.imcalc.input = cellstr(char(pos, neg));
            jobs{jobCount}.spm.util.imcalc.output = outputname;
            jobs{jobCount}.spm.util.imcalc.outdir = cellstr(outDir);
            jobs{jobCount}.spm.util.imcalc.expression = 'i1 | i2';
            jobs{jobCount}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            jobs{jobCount}.spm.util.imcalc.options.dmtx = 0;
            jobs{jobCount}.spm.util.imcalc.options.mask = 0;
            jobs{jobCount}.spm.util.imcalc.options.interp = 0;
            jobs{jobCount}.spm.util.imcalc.options.dtype = 2;
            jobCount = jobCount + 1;
        end
    end

    conCount = conCount + 1;
    if conCount > numel(contrast)
        break
    end
end

run_spm_jobs('export_contrasts_binary', jobs);

end

function run_spm_jobs(jobName, jobs)
jobId = cfg_util('initjob', jobs);
sts = cfg_util('isjob_id', jobId);
if sts
    cfg_util('run', jobId);
else
    error('Error in initialising %s job.', jobName)
end
end
