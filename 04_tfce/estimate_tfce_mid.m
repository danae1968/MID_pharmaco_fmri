function estimate_tfce_mid(SPMmat,jobDir,tfce)
% ESTIMATE_TFCE_MID  Run TFCE estimation (TFCE toolbox) on group SPM models.
%
% estimate_tfce_mid(SPMmat, jobDir, tfce)
%
% INPUT
% SPMmat    : cell array of full paths to SPM.mat file(s)
% jobDir    : folder where the SPM batch file is saved
% tfce      : structure with fields
%  - mask   : cell with path to the analysis mask ('' = SPM's model mask)
%  - results_title_prefix : title prefix for the results
%  - nperm  : number of permutations (e.g. 5000)
%  - weight_clustsize : TFCE extent weighting E (0.5 here; lower values
%             weight focal effects more than large clusters)
%  - modelConIxs : indices of the contrasts to process, e.g. [1 2]
%  - nproc  : number of cores
%  - nuisance_method : permutation scheme when the design contains a
%             nuisance regressor besides the tested effect
%             (0 = Draper-Stoneman, 1 = Freedman-Lane, 2 = Smith; Winkler
%             et al., 2014, NeuroImage 92:381-397). Irrelevant for plain
%             one-sample designs; for the covariate models (intercept +
%             covariate) Draper-Stoneman is used. Default: 0.
%
% OUTPUT
% TFCE output (e.g. TFCE_log_pFWE_000N.nii) is written next to SPM.mat.
%
% -------------------------------------------------------------------------
% Ruben van den Bosch, August 2023; adapted by Danae Papadopetraki, 2026
%
% Process input arguments
% -------------------------------------------------------------------------
assert(iscell(SPMmat),'SPMmat input argument should be of type cell');
assert(ischar(jobDir),'jobDir input argument should be of type char');
if ~isfield(tfce,'nuisance_method') || isempty(tfce.nuisance_method)
    tfce.nuisance_method = 0;  % Draper-Stoneman default
end
% Create job(s) for TFCE estimation
% -------------------------------------------------------------------------
jobs = cell(1,numel(SPMmat));
for iJob = 1:numel(jobs)
    jobs{iJob}.spm.tools.tfce_estimate.data              = SPMmat(iJob);
    jobs{iJob}.spm.tools.tfce_estimate.nproc             = tfce.nproc;
    jobs{iJob}.spm.tools.tfce_estimate.mask              = tfce.mask;
    jobs{iJob}.spm.tools.tfce_estimate.conspec.titlestr  = tfce.results_title_prefix;
    jobs{iJob}.spm.tools.tfce_estimate.conspec.contrasts = tfce.modelConIxs;
    jobs{iJob}.spm.tools.tfce_estimate.conspec.n_perm    = tfce.nperm;
    jobs{iJob}.spm.tools.tfce_estimate.nuisance_method   = tfce.nuisance_method;
    jobs{iJob}.spm.tools.tfce_estimate.tbss              = 0;
    jobs{iJob}.spm.tools.tfce_estimate.E_weight          = tfce.weight_clustsize;
    jobs{iJob}.spm.tools.tfce_estimate.singlethreaded    = 0;
end
% Save and run job
% -------------------------------------------------------------------------
jobName = 'TFCE_estimation';
run_spm_jobs(jobName, jobDir, jobs);
end