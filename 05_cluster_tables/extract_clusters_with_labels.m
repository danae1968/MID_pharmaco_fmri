function extract_clusters_with_labels(binary_image, output_csv, opts)
% EXTRACT_CLUSTERS_WITH_LABELS: Cluster table with subpeaks and atlas labels
%
% Builds a results table from a thresholded (binary) significant-voxel image:
% connected clusters, local subpeaks (>min distance apart), Neuromorphometrics
% region label and hemisphere per peak, peak t, cluster size, and the TFCE
% familywise-error p at each peak.
%
% USAGE
%   opts.spmT_path   = '.../spmT_0002.nii';   % stat image for this contrast
%   opts.tfce_p_path = '.../TFCE_log_pFWE_0002.nii';  % -log10(p) FWE map
%   opts.contrast_name = 'Cue';
%   opts.direction     = 'positive';
%   extract_clusters_with_labels(mask, 'clusters.csv', opts)
%
% OPTS FIELDS
%   spmT_path        - stat image; peaks and peak t are read from it (required)
%   tfce_p_path      - TFCE_log_pFWE image; gives TFCE FWE p per peak (optional)
%   contrast_name    - label for the Contrast column (default '')
%   direction        - 'positive' | 'negative' for the Direction column (default '')
%   cluster_extent   - minimum cluster size in voxels (default 0)
%   min_peak_dist_mm - minimum distance between subpeaks (default 8)
%   max_peaks        - max subpeaks reported per cluster (default 3)
%   region_csv       - if set, also write a per-region breakdown of every
%                      cluster here: one row per (cluster x atlas region),
%                      with the number of significant voxels in that region
%                      and the strongest peak within it. Useful for large
%                      TFCE clusters that a peak list describes poorly.
%   max_regions      - max regions reported per cluster (default 12)
%   min_region_k     - ignore regions holding fewer voxels than this
%                      (default 5), so single-voxel atlas noise is dropped
%   label_fallback_mm - see below (default 10)
%
% NOTE: Neuromorphometrics gives gyrus labels only, not Brodmann areas.
%
% -------------------------------------------------------------------------
% Danae Papadopetraki, April 2026

% ---- Options / defaults -------------------------------------------------
if nargin < 3, opts = struct(); end
if ~isfield(opts, 'cluster_extent'),   opts.cluster_extent   = 0;  end
if ~isfield(opts, 'min_peak_dist_mm'), opts.min_peak_dist_mm = 8;  end
if ~isfield(opts, 'max_peaks'),        opts.max_peaks        = 3;  end
if ~isfield(opts, 'contrast_name'),    opts.contrast_name    = ''; end
if ~isfield(opts, 'direction'),        opts.direction        = ''; end
if ~isfield(opts, 'spmT_path'),        opts.spmT_path        = ''; end
if ~isfield(opts, 'tfce_p_path'),      opts.tfce_p_path      = ''; end
if ~isfield(opts, 'region_csv'),       opts.region_csv       = ''; end
if ~isfield(opts, 'max_regions'),      opts.max_regions      = 12; end
if ~isfield(opts, 'min_region_k'),     opts.min_region_k     = 5;  end
% Peaks on an UNLABELLED atlas voxel (0), which happens at the edges of the
% Neuromorphometrics parcellation, get the nearest labelled grey-matter
% region within this radius instead of 'n/a'. Set to 0 to disable.
% Peaks that are labelled white matter or ventricle keep that label:
% renaming them to a neighbouring gyrus would misreport the peak location.
if ~isfield(opts, 'label_fallback_mm'), opts.label_fallback_mm = 10; end

fprintf('Loading image: %s\n', binary_image);
V   = spm_vol(binary_image);
img = spm_read_vols(V);

maskIdx = find(img > 0);
if isempty(maskIdx)
    fprintf('No significant voxels, nothing to extract.\n');
    return;
end

% ---- Statistic image (peaks and peak t) --------------------------------
if isempty(opts.spmT_path) || ~exist(opts.spmT_path, 'file')
    error('opts.spmT_path is required and must exist: %s', opts.spmT_path);
end
Vt       = spm_vol(opts.spmT_path);
stat_img = spm_read_vols(Vt);
% SPM writes a positive-t map for each contrast (a "negative" contrast has
% positive t where its effect is significant), so peaks are max t either way.

% ---- Peak FWE p source: TFCE map, parametric RFT, else NaN -------------
pfwe_source = 'none';
if isfield(opts,'tfce_p_path') && ~isempty(opts.tfce_p_path) && exist(opts.tfce_p_path,'file')
    Vp   = spm_vol(opts.tfce_p_path);
    logp = spm_read_vols(Vp);        % -log10(p_FWE)
    pfwe_source = 'tfce';
elseif isfield(opts,'spm') && ~isempty(opts.spm) && ...
       isfield(opts,'contrast_idx') && ~isempty(opts.contrast_idx)
    SPMp  = opts.spm;  cP = opts.contrast_idx;
    STATp = SPMp.xCon(cP).STAT;      % 'T' or 'F'
    Rp    = SPMp.xVol.R;             % resel counts, stored at estimation
    if STATp == 'F'
        dfp = [SPMp.xCon(cP).eidf, SPMp.xX.erdf];
    else
        dfp = [1, SPMp.xX.erdf];
    end
    pfwe_source = 'rft';
else
    fprintf('Warning: no p source given, pFWE column will be NaN\n');
end

% ---- Local maxima via spm_max ------------------------------------------
[x, y, z] = ind2sub(size(img), maskIdx);
XYZvox    = [x'; y'; z'];
vals      = stat_img(maskIdx)';

[N, Zmax, Mmax, A] = spm_max(vals, XYZvox);   % sizes, peak vals, peak vox, cluster id
% spm_max may return column vectors; force rows so the cluster loop and
% logical indexing below stay orientation-safe.
N = N(:).'; Zmax = Zmax(:).'; A = A(:).';
if size(Mmax, 1) ~= 3, Mmax = Mmax.'; end     % ensure 3 x nMaxima

% ---- Atlas (Neuromorphometrics, shipped with SPM12) --------------------
[atlasV, atlasImg, labelMap] = load_neuromorphometrics();

% ---- Per-voxel cluster membership and atlas labels (region breakdown) --
% spm_max returns one entry per MAXIMUM, so it cannot say which voxels make
% up a cluster. spm_clusters gives a label per voxel using the same
% connectivity, and each spm_max cluster is matched to its spm_clusters
% label through its strongest peak voxel, so the two numberings agree.
wantRegions = ~isempty(opts.region_csv);
if wantRegions
    CL        = spm_clusters(XYZvox);
    allMm     = V.mat(1:3,:) * [XYZvox; ones(1, size(XYZvox,2))];
    voxLabel  = atlas_label_indices(allMm, atlasV, atlasImg);
    regionRows = {};
    regionHeader = {'Contrast','Direction','Cluster','Size_k','Region','Hemisphere', ...
                    'Region_k','Region_pct','Peak_t','x','y','z','pFWE'};
end

% ---- Build rows: prune subpeaks per cluster ----------------------------
header = {'Contrast','Direction','Region','Hemisphere', ...
          'x','y','z','Peak_t','Cluster','Size_k','pFWE'};
rows = {};
clusterCount = 0;

clusterIds = unique(A);
clusterIds = clusterIds(:).';
for a = clusterIds
    sel = find(A == a);
    kSize = N(sel(1));
    if kSize < opts.cluster_extent, continue; end

    % Peaks of this cluster, strongest first
    [~, ord] = sort(Zmax(sel), 'descend');
    sel = sel(ord);

    Mmm = V.mat(1:3,:) * [Mmax(:, sel); ones(1, numel(sel))];  % maxima in mm

    kept = [];  % indices (into sel) of retained peaks
    for p = 1:numel(sel)
        if numel(kept) >= opts.max_peaks, break; end
        if isempty(kept)
            kept(end+1) = p; %#ok<AGROW>
        else
            d = sqrt(sum((Mmm(:, kept) - Mmm(:, p)).^2, 1));
            if all(d > opts.min_peak_dist_mm)
                kept(end+1) = p; %#ok<AGROW>
            end
        end
    end

    clusterCount = clusterCount + 1;
    for kk = 1:numel(kept)
        p       = kept(kk);
        peakVox = Mmax(:, sel(p));
        peakMm  = Mmm(:, p);
        peakT   = Zmax(sel(p));

       % Peak-level FWE p
        switch pfwe_source
            case 'tfce'
                lv   = logp(peakVox(1), peakVox(2), peakVox(3));
                pFWE = 10.^(-lv);
            case 'rft'
                pFWE = spm_P_RF(1, 0, peakT, dfp, STATp, Rp, 1);
            otherwise
                pFWE = NaN;
        end

        % Atlas label + hemisphere
        [region, hemi] = atlas_label(peakMm, atlasV, atlasImg, labelMap, ...
            opts.label_fallback_mm);

        rows(end+1, :) = { opts.contrast_name, opts.direction, region, hemi, ...
            round(peakMm(1)), round(peakMm(2)), round(peakMm(3)), ...
            peakT, clusterCount, kSize, pFWE };  %#ok<AGROW>

        if kk == 1
            fprintf('  Cluster %d (k=%d): %s %s [%d %d %d] t=%.2f p=%.3g\n', ...
                clusterCount, kSize, hemi, region, ...
                round(peakMm(1)), round(peakMm(2)), round(peakMm(3)), peakT, pFWE);
        else
            fprintf('     subpeak: %s %s [%d %d %d] t=%.2f p=%.3g\n', ...
                hemi, region, round(peakMm(1)), round(peakMm(2)), round(peakMm(3)), peakT, pFWE);
        end
    end

    % ---- Region breakdown of this cluster ------------------------------
    if wantRegions
        topVox = Mmax(:, sel(1));
        vi = find(XYZvox(1,:) == topVox(1) & XYZvox(2,:) == topVox(2) & ...
                  XYZvox(3,:) == topVox(3), 1);
        if isempty(vi), continue; end
        vs   = find(CL == CL(vi));            % every voxel of this cluster
        labs = voxLabel(vs);

        uLab = unique(labs(labs > 0));        % 0 = outside the atlas
        cnt  = arrayfun(@(L) sum(labs == L), uLab);
        [cnt, ordL] = sort(cnt, 'descend');
        uLab = uLab(ordL);

        nShown = 0;
        for r = 1:numel(uLab)
            if cnt(r) < opts.min_region_k || nShown >= opts.max_regions, break; end
            nShown = nShown + 1;

            inR       = vs(labs == uLab(r));
            [rT, rIx] = max(vals(inR));
            rVox      = XYZvox(:, inR(rIx));
            rMm       = allMm(:, inR(rIx));

            switch pfwe_source
                case 'tfce', rP = 10.^(-logp(rVox(1), rVox(2), rVox(3)));
                case 'rft',  rP = spm_P_RF(1, 0, rT, dfp, STATp, Rp, 1);
                otherwise,   rP = NaN;
            end

            if isKey(labelMap, uLab(r))
                rName = regexprep(labelMap(uLab(r)), '^(Left|Right)\s+', '');
            else
                rName = 'n/a';
            end
            if rMm(1) < 0, rHemi = 'L'; elseif rMm(1) > 0, rHemi = 'R'; else, rHemi = 'M'; end

            regionRows(end+1, :) = { opts.contrast_name, opts.direction, ...
                clusterCount, kSize, rName, rHemi, cnt(r), ...
                100 * cnt(r) / numel(vs), rT, ...
                round(rMm(1)), round(rMm(2)), round(rMm(3)), rP }; %#ok<AGROW>
        end
    end
end

% ---- Write --------------------------------------------------------------
if isempty(rows)
    fprintf('No clusters survived extent threshold of %d voxels\n', opts.cluster_extent);
    return;
end
T = cell2table(rows, 'VariableNames', header);
writetable(T, output_csv);
fprintf('Saved %d peaks in %d clusters to: %s\n', size(rows,1), clusterCount, output_csv);

if wantRegions
    if isempty(regionRows)
        fprintf('No regions passed min_region_k = %d, no region table written.\n', ...
            opts.min_region_k);
    else
        TR = cell2table(regionRows, 'VariableNames', regionHeader);
        writetable(TR, opts.region_csv);
        fprintf('Saved %d cluster-region rows to: %s\n', size(regionRows,1), opts.region_csv);
    end
end
end

% -------------------------------------------------------------------------
function lbl = atlas_label_indices(mm, atlasV, atlasImg)
% Atlas label INDEX for a 3xN set of mm coordinates (0 where unlabelled).
% Index rather than name so the two hemispheres of a structure stay
% separate -- the names have their Left/Right prefix stripped for display.
lbl = zeros(1, size(mm, 2));
if isempty(atlasV), return; end
vox = round(atlasV.mat \ [mm; ones(1, size(mm,2))]);
d   = atlasV.dim;
ok  = all(vox(1:3,:) >= 1, 1) & vox(1,:) <= d(1) & vox(2,:) <= d(2) & vox(3,:) <= d(3);
if ~any(ok), return; end
lbl(ok) = atlasImg(sub2ind(d, vox(1,ok), vox(2,ok), vox(3,ok)));
end

% =========================================================================
% Helpers
% =========================================================================
function [atlasV, atlasImg, labelMap] = load_neuromorphometrics()
% Load the Neuromorphometrics atlas nii + xml shipped with SPM12.
niiPath = fullfile(spm('Dir'), 'tpm', 'labels_Neuromorphometrics.nii');
xmlPath = fullfile(spm('Dir'), 'tpm', 'labels_Neuromorphometrics.xml');
if ~exist(niiPath, 'file') || ~exist(xmlPath, 'file')
    warning('Neuromorphometrics atlas not found in %s. Labels will be blank.', ...
        fullfile(spm('Dir'), 'tpm'));
    atlasV = []; atlasImg = []; labelMap = containers.Map('KeyType','double','ValueType','char');
    return;
end
atlasV   = spm_vol(niiPath);
atlasImg = spm_read_vols(atlasV);

% Parse <index>..</index><name>..</name> pairs from the xml
raw = fileread(xmlPath);
tok = regexp(raw, '<index>(\d+)</index>\s*<name>(.*?)</name>', 'tokens');
labelMap = containers.Map('KeyType','double','ValueType','char');
for i = 1:numel(tok)
    labelMap(str2double(tok{i}{1})) = tok{i}{2};
end
end

function [region, hemi] = atlas_label(mm, atlasV, atlasImg, labelMap, fallback_mm)
% Look up the atlas region at an mm coordinate; hemisphere from x sign.
% Falls back to the nearest labelled grey-matter voxel within fallback_mm
% when the peak itself is unlabelled (see opts.label_fallback_mm above).
if nargin < 5, fallback_mm = 0; end
if mm(1) < 0, hemi = 'L'; elseif mm(1) > 0, hemi = 'R'; else, hemi = 'M'; end
region = 'n/a';
if isempty(atlasV), return; end

vox = round(atlasV.mat \ [mm(:); 1]);
d   = atlasV.dim;
if any(vox(1:3) < 1) || any(vox(1:3)' > d), return; end
lblIdx = atlasImg(vox(1), vox(2), vox(3));

if lblIdx == 0 && fallback_mm > 0
    lblIdx = nearest_labelled(vox(1:3)', atlasV, atlasImg, labelMap, fallback_mm);
end
if lblIdx == 0, return; end

if isKey(labelMap, lblIdx)
    region = labelMap(lblIdx);
    % Strip a leading Left/Right so it does not duplicate the hemisphere column
    region = regexprep(region, '^(Left|Right)\s+', '');
end
end

% -------------------------------------------------------------------------
function lblIdx = nearest_labelled(vox, atlasV, atlasImg, labelMap, radius_mm)
% Nearest labelled grey-matter atlas index to a voxel, searching outward.
% Returns 0 if nothing qualifying is found inside radius_mm.
lblIdx = 0;
vsize  = sqrt(sum(atlasV.mat(1:3,1:3).^2, 1));       % mm per voxel
rvox   = ceil(radius_mm ./ vsize);
d      = atlasV.dim;

[ox, oy, oz] = ndgrid(-rvox(1):rvox(1), -rvox(2):rvox(2), -rvox(3):rvox(3));
off  = [ox(:), oy(:), oz(:)];
dist = sqrt(sum((off .* vsize).^2, 2));
keep = dist > 0 & dist <= radius_mm;
off  = off(keep, :);  dist = dist(keep);
[~, order] = sort(dist);
off = off(order, :);

for i = 1:size(off, 1)
    v = vox + off(i, :);
    if any(v < 1) || any(v > d), continue; end
    cand = atlasImg(v(1), v(2), v(3));
    if cand == 0 || ~isKey(labelMap, cand), continue; end
    nm = labelMap(cand);
    % Skip non-informative tissue classes -- see the note at the top.
    if ~isempty(regexpi(nm, 'White Matter|Ventricle|CSF|Vessel', 'once'))
        continue
    end
    lblIdx = cand; return
end
end
