% SETUP_PATHS  Add all pipeline folders of this repository to the MATLAB path.
%
% Run once per MATLAB session, from any directory. SPM12, the TFCE toolbox
% and (for figures) slice_display and panel must be added separately.

repoDir = fileparts(mfilename('fullpath'));
steps   = dir(fullfile(repoDir, '0*_*'));
for i = 1:numel(steps)
    if steps(i).isdir
        addpath(fullfile(repoDir, steps(i).name));
    end
end
fprintf('Added %d pipeline folders to the MATLAB path.\n', sum([steps.isdir]));
