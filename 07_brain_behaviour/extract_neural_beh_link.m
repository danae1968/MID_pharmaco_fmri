%% extract_neural_beh_link.m
%
% Per-subject values for the brain-behaviour scatter plots (Figure 8B and
% Supplementary Figure S1B). For each of the two Cue covariate clusters
% (MPH-PBO anticipation contrast x MPH-PBO change in hit rate / premature
% responses), extracts every subject's mean value from their own
% first-level MPH-vs-PBO Cue image within the principal (largest
% connected) cluster of the TFCE mask. Plotted in
% neural_beh_scatter.R.
%
% The clusters were defined by their association with the behavioural
% measure, so these values only display the form of that association;
% they are not an independent test of it.
%
% ROIs are the TFCE_significant_voxels_*_fwe_p05.nii masks written by
% write_binary_tfce.m; no new ROI definition is introduced.
%
% Output: one long-format CSV, one row per subject x cluster:
%   subject, drugComparison, contrast, cluster, roiMean, nVoxInMask
%
% -------------------------------------------------------------------------

%% Paths -----------------------------------------------------------------
if ispc
    projectDir = 'P:\3017048.01';
    addpath('\\fileserver.dccn.nl\home\common\matlab\spm12');
elseif isunix
    projectDir = '/project/3017048.01/';
    addpath('/home/common/matlab/spm12');
else
    error('Unknown OS');
end
spm('fmri','defaults');

derivMRIdir   = fullfile(projectDir,'bids','derivatives','mri');
groupLevelBase = fullfile(derivMRIdir,'groupLevel','MID','regular','canonical_hrf');
drugcodesFile = fullfile(projectDir,'admin','9.Deblinding','drugcodes_sessByDrug.csv');
outCsv        = fullfile(projectDir,'bids','derivatives','beh','MID','stats','neural_roi_values.csv');

denoise = 'regular';       % matches todo.denoise.regular = true in run_groupLevel_vici_MID.m
hrfType = 'canonical_hrf';

% First-level folder layout: 'MID' on the cluster, hrfType-named in a
% local copy (as in groupLevel_vici_MID_s.m)
if exist(fullfile(derivMRIdir,'firstLevel','MID'),'dir')
    firstLevelTaskDir = 'MID';
else
    firstLevelTaskDir = hrfType;
end

% Same session-level exclusions as run_groupLevel_vici_MID.m (keep in sync)
exclude.subject = [3,6,9,12,16,20,21,24,26,27,30,36,59,65,76,77,90,94,95,99];
exclude.session  = {[2],[3],[3],[2,3],[2],[1,3],[3],[3],[2],[1],[1],[2,3],[3],[1],[1],[1],[2],[1],[3],[1,3]};

subjectIx = 1:100;

drugcodes = readtable(drugcodesFile);

%% Clusters to extract -----------------------------------------------------
% One row per cluster:
%   1 drug comparison ('MPHvsPBO'; split on 'vs' to find the sessions)
%   2 task contrast
%   3 'pos' | 'neg' : positive or negative_ TFCE mask
%   4 cluster label used in the output
%   5 covariate stem: the mask is '<con>_x_<stem>' from covariate_<stem>/.
%     The per-subject image is the first-level drug-comparison image.
clusters = {
    'MPHvsPBO', 'Cue', 'pos', 'parietal_MPHvsPBO_Cue_x_hit',   'mphBehEffect_hit'    % Figure 8
    'MPHvsPBO', 'Cue', 'neg', 'parietal_MPHvsPBO_Cue_x_early', 'mphBehEffect_early'  % Figure S1
};

% Only the largest connected cluster of each mask is used: a TFCE mask
% holds every surviving cluster of a contrast, and averaging over disjoint
% regions would mix them. Connectivity as in extract_clusters_with_labels.m
% (spm_clusters).

rows = {};  % subject, drugComparison, contrast, cluster, roiMean, nVoxInMask

for c = 1:size(clusters,1)
    drugComparison = clusters{c,1};
    conNm          = clusters{c,2};
    sgn            = clusters{c,3};
    label          = clusters{c,4};
    covStem        = clusters{c,5};

    maskDir = fullfile(groupLevelBase, ['covariate_' covStem]);
    maskCon = [conNm '_x_' covStem];

    if strcmp(sgn,'pos')
        maskFile = fullfile(maskDir, [drugComparison '_' conNm], ...
            sprintf('TFCE_significant_voxels_%s_fwe_p05.nii', maskCon));
    else
        maskFile = fullfile(maskDir, [drugComparison '_' conNm], ...
            sprintf('TFCE_significant_voxels_negative_%s_fwe_p05.nii', maskCon));
    end

    if ~exist(maskFile,'file')
        warning('Mask not found, skipping: %s', maskFile);
        continue
    end

    Vm   = spm_vol(maskFile);
    mask = spm_read_vols(Vm) > 0;
    nVoxAll = sum(mask(:));

    maskIdx = find(mask);
    [mx, my, mz] = ind2sub(size(mask), maskIdx);
    CL = spm_clusters([mx'; my'; mz']);
    counts = accumarray(CL(:), 1);
    [kBig, biggest] = max(counts);
    mask = false(size(mask));
    mask(maskIdx(CL == biggest)) = true;
    fprintf(['\n%s (%s, %s): restricted to the largest of %d clusters, ' ...
             '%d of %d voxels\n'], label, drugComparison, conNm, ...
             numel(counts), kBig, nVoxAll);
    nVox = sum(mask(:));

    drugsInComparison = strsplit(drugComparison,'vs');   % e.g. {'MPH','PBO'}

    for ix = 1:numel(subjectIx)
        iSubject = subjectIx(ix);

        % Skip the subject if either session is excluded (as in
        % groupLevel_vici_MID_s.m)
        excluded = false;
        if any(exclude.subject == iSubject)
            sessNrs = [ ...
                drugcodes.(sprintf('%ssession',upper(drugsInComparison{1})))(drugcodes.subject==iSubject), ...
                drugcodes.(sprintf('%ssession',upper(drugsInComparison{2})))(drugcodes.subject==iSubject) ];
            if any(ismember(sessNrs, exclude.session{exclude.subject==iSubject}))
                excluded = true;
            end
        end
        if excluded, continue; end

        drugCompDir = fullfile(derivMRIdir,'firstLevel',firstLevelTaskDir, ...
            sprintf('sub-%.3d',iSubject),'drugComparisons',denoise,hrfType);
        conImg = fullfile(drugCompDir, [drugComparison '_' conNm '.nii']);

        if ~(exist(conImg,'file') == 2)
            continue   % no image for this subject (incomplete session)
        end

        Vc = spm_vol(conImg);
        Yc = spm_read_vols(Vc);
        if ~isequal(size(Yc), size(mask))
            warning('Dimension mismatch for sub-%.3d, %s -- skipping', iSubject, label);
            continue
        end

        roiVals = Yc(mask);
        roiVals = roiVals(~isnan(roiVals));
        if isempty(roiVals), continue; end

        rows(end+1,:) = {iSubject, drugComparison, conNm, label, mean(roiVals), nVox}; %#ok<AGROW>
    end
end

if isempty(rows)
    error(['No rows extracted at all -- every mask/subject lookup failed. ' ...
        'Check the "Mask not found" warnings above for the paths tried.']);
end

T = cell2table(rows, 'VariableNames', ...
    {'subject','drugComparison','contrast','cluster','roiMean','nVoxInMask'});
writetable(T, outCsv);
fprintf('\nSaved %d rows (%d subjects x %d clusters) to:\n  %s\n', ...
    height(T), numel(unique(T.subject)), size(clusters,1), outCsv);
