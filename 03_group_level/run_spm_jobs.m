function run_spm_jobs(jobName,jobDir,jobs)
% RUN_SPM_JOBS(jobName,jobDir,jobs)
%
% Save an SPM batch to jobDir (timestamped, for provenance) and run it.
% Used by estimate_tfce_mid.m.

% Filename for the saved job
jobFile = fullfile(jobDir,[jobName,'_',datestr(now,'yyyymmddTHHMMSS'),'.m']);

% Initialise job
jobId = cfg_util('initjob', jobs);

% If successful, save and run
sts = cfg_util('isjob_id', jobId);
if sts
    cfg_util('savejob', jobId, jobFile);
    cfg_util('run', jobId);
else
    error('Error in initialising %s job.',jobName)
end
end
