function [slices, tbl] = sd_pick_slices(maskFile, orientation, nSlices, opts)
% SD_PICK_SLICES Choose a few representative slices from a significance mask.
%
% DESCRIPTION
% Large effects can contain significant voxels in 30+ slices per plane; a
% figure needs a handful that show the effect. This scores every candidate
% slice position along the display axis and greedily picks the best ones,
% keeping them a minimum distance apart so they do not all cut through the
% same blob.
%
% SYNTAX
%   [slices, tbl] = sd_pick_slices(maskFile, orientation, nSlices, opts)
%
% maskFile    - binary significance mask (.nii); voxels > 0 count as significant
% orientation - 'axial' | 'coronal' | 'sagittal'
% nSlices     - how many slices to return
% opts        - struct, all fields optional:
%       .metric      'peaks'  (default) score each slice by the strongest
%                             |t| it contains -> slices through the peaks of
%                             the biggest effects
%                    'extent' score by the number of significant voxels ->
%                             slices that best show the spatial extent
%                    'spread' ignore scores, space the slices evenly across
%                             the significant range
%       .stat_file   spmT image, REQUIRED for metric 'peaks'
%       .min_sep_mm  minimum spacing between chosen slices (default 10)
%       .round_to    snap candidate positions to this step in mm (default 2)
%
% slices      - ascending vector of slice coordinates in mm
% tbl         - struct array over ALL candidate slices, with fields .mm,
%               .n_voxels, .max_t and .chosen -- useful for sanity-checking
%               or for hand-picking a different set afterwards.
%
% NOTE ON THE MASK
% Pass the mask you are actually going to draw the contour of. For the TFCE
% figures that is the corrected mask, so the chosen slices are the ones that
% show the corrected clusters -- not the (much larger) uncorrected blob.

if nargin < 4, opts = struct(); end
if ~isfield(opts,'metric')     || isempty(opts.metric),     opts.metric     = 'peaks'; end
if ~isfield(opts,'stat_file'),                              opts.stat_file  = '';      end
if ~isfield(opts,'min_sep_mm') || isempty(opts.min_sep_mm), opts.min_sep_mm = 10;      end
if ~isfield(opts,'round_to')   || isempty(opts.round_to),   opts.round_to   = 2;       end

axisRow = struct('sagittal', 1, 'coronal', 2, 'axial', 3);
if ~isfield(axisRow, lower(orientation))
    error('sd_pick_slices:badOrientation','Unknown orientation "%s".', orientation);
end
row = axisRow.(lower(orientation));

% ---- Significant voxels and their world coordinates ----------------------
Vm            = spm_vol(maskFile);
[Ym, XYZmm]   = spm_read_vols(Vm);
sig           = find(Ym > 0);
if isempty(sig)
    slices = []; tbl = struct('mm',{},'n_voxels',{},'max_t',{},'chosen',{});
    return
end

% ---- |t| at those voxels (for the 'peaks' metric) ------------------------
absT = [];
if strcmpi(opts.metric,'peaks')
    if isempty(opts.stat_file) || ~exist(opts.stat_file,'file')
        warning('sd_pick_slices:noStat', ...
            'metric ''peaks'' needs opts.stat_file; falling back to ''extent''.');
        opts.metric = 'extent';
    else
        Vt = spm_vol(opts.stat_file);
        if ~isequal(Vt.dim, Vm.dim) || ~isequal(round(Vt.mat*1e4), round(Vm.mat*1e4))
            warning('sd_pick_slices:gridMismatch', ...
                ['Stat image and mask are not on the same grid ' ...
                 '(%s vs %s); falling back to ''extent''.'], opts.stat_file, maskFile);
            opts.metric = 'extent';
        else
            Yt   = spm_read_vols(Vt);
            absT = abs(Yt(sig));
            absT(~isfinite(absT)) = 0;
        end
    end
end

% ---- Candidate slice positions and their scores --------------------------
pos = round(XYZmm(row, sig) / opts.round_to) * opts.round_to;
mm  = unique(pos);

nVox = zeros(size(mm));
maxT = zeros(size(mm));
for i = 1:numel(mm)
    here    = (pos == mm(i));
    nVox(i) = sum(here);
    if ~isempty(absT)
        maxT(i) = max(absT(here));
    end
end

switch lower(opts.metric)
    case 'peaks',  score = maxT;
    case 'extent', score = nVox;
    case 'spread', score = [];
    otherwise
        error('sd_pick_slices:badMetric','Unknown metric "%s".', opts.metric);
end

% ---- Pick --------------------------------------------------------------
if strcmpi(opts.metric,'spread')
    if numel(mm) <= nSlices
        slices = mm;
    else
        want   = linspace(mm(1), mm(end), nSlices);
        [~,ix] = min(abs(bsxfun(@minus, mm(:), want)), [], 1);
        slices = unique(mm(ix));
    end
else
    % Greedy: best score first, never within min_sep_mm of an already
    % chosen slice. Ties broken towards the slice with more voxels, then
    % towards the middle of the significant range.
    [~, ord] = sortrows([-score(:), -nVox(:), abs(mm(:) - median(mm))]);
    picked   = [];
    for k = 1:numel(ord)
        cand = mm(ord(k));
        if isempty(picked) || all(abs(picked - cand) >= opts.min_sep_mm)
            picked(end+1) = cand; %#ok<AGROW>
        end
        if numel(picked) >= nSlices, break; end
    end
    % If the effect is too compact to yield nSlices at this spacing, relax
    % the spacing rather than returning a single slice.
    if numel(picked) < nSlices && numel(mm) > numel(picked)
        for k = 1:numel(ord)
            cand = mm(ord(k));
            if ~ismember(cand, picked)
                picked(end+1) = cand; %#ok<AGROW>
            end
            if numel(picked) >= nSlices, break; end
        end
    end
    slices = sort(picked);
end

% ---- Report table --------------------------------------------------------
tbl = struct('mm', num2cell(mm), ...
             'n_voxels', num2cell(nVox), ...
             'max_t', num2cell(maxT), ...
             'chosen', num2cell(ismember(mm, slices)));
end
