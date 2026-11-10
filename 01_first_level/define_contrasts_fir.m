function contrasts = define_contrasts_fir(dirs, SPMmat, settings)
% DEFINE_CONTRASTS_FIR  Per-time-bin contrast vectors for the FIR MID model.
%
% With a finite impulse response basis set every condition has one
% regressor per time bin. This creates one contrast per bin for:
%   HCvsLC_binNN     high - low reward cue
%   Cue_binNN        mean of high and low reward cue
%   HitvsMiss_binNN  hits - misses, averaged over reward level
%   Outcome_binNN    mean of all feedback conditions
%   Target_binNN     target onset
%
% The bin images are analysed at group level in a Subject x TimeBin
% design (see groupLevel_vici_MID_s.m).
%
% Danae Papadopetraki, 2025
% -------------------------------------------------------------

cd(dirs);
load(SPMmat, 'SPM');

Nbins = settings.order;     % number of FIR bins
nReg  = numel(SPM.xX.name); % number of regressors
contrasts = struct();

% All FIR regressors (columns named '... <name>*bf(k)') of a condition
get_idx = @(name) find(contains(SPM.xX.name, name) & contains(SPM.xX.name, 'bf('));

% =============================================================
% Regressor indices per condition
% =============================================================
conditions = struct();

% Cue (anticipation)
conditions.Cue.High = get_idx('HighCue');
conditions.Cue.Low  = get_idx('LowCue');

% Outcome (feedback)
conditions.Outcome.HighHit  = get_idx('HighHit');
conditions.Outcome.LowHit   = get_idx('LowHit');
conditions.Outcome.HighMiss = get_idx('HighMiss');
conditions.Outcome.LowMiss  = get_idx('LowMiss');

% Target (response)
conditions.Target.Bins = get_idx('Target');

% Sort each condition's regressors by bin number
sort_bins = @(idx) sort_by_bin(SPM, idx);

% =============================================================
% Per-bin contrasts
% =============================================================

% High vs low reward cue
if ~isempty(conditions.Cue.High) && ~isempty(conditions.Cue.Low)
    [~, hc_sorted] = sort_bins(conditions.Cue.High);
    [~, lc_sorted] = sort_bins(conditions.Cue.Low);
    for b = 1:Nbins
        cvec = zeros(1, nReg);
        cvec(hc_sorted(b)) =  1;
        cvec(lc_sorted(b)) = -1;
        contrasts.(sprintf('HCvsLC_bin%02d', b)) = cvec;
    end
end

% Cue, averaged over reward level
if ~isempty(conditions.Cue.High) && ~isempty(conditions.Cue.Low)
    [~, hc_sorted] = sort_bins(conditions.Cue.High);
    [~, lc_sorted] = sort_bins(conditions.Cue.Low);
    for b = 1:Nbins
        cvec = zeros(1, nReg);
        cvec(hc_sorted(b)) =  0.5;
        cvec(lc_sorted(b)) = 0.5;
        contrasts.(sprintf('Cue_bin%02d', b)) = cvec;
    end
end

% Hit vs miss
if all(structfun(@(x) ~isempty(x), conditions.Outcome))
    [~, hh_sorted] = sort_bins(conditions.Outcome.HighHit);
    [~, lh_sorted] = sort_bins(conditions.Outcome.LowHit);
    [~, hm_sorted] = sort_bins(conditions.Outcome.HighMiss);
    [~, lm_sorted] = sort_bins(conditions.Outcome.LowMiss);
    for b = 1:Nbins
        cvec = zeros(1, nReg);
        cvec([hh_sorted(b), lh_sorted(b)]) =  0.5;
        cvec([hm_sorted(b), lm_sorted(b)]) = -0.5;
        contrasts.(sprintf('HitvsMiss_bin%02d', b)) = cvec;
    end
end

% Outcome, averaged over all feedback conditions
if all(structfun(@(x) ~isempty(x), conditions.Outcome))
    [~, hh_sorted] = sort_bins(conditions.Outcome.HighHit);
    [~, lh_sorted] = sort_bins(conditions.Outcome.LowHit);
    [~, hm_sorted] = sort_bins(conditions.Outcome.HighMiss);
    [~, lm_sorted] = sort_bins(conditions.Outcome.LowMiss);
    for b = 1:Nbins
        cvec = zeros(1, nReg);
        cvec([hh_sorted(b), lh_sorted(b)]) =  0.5;
        cvec([hm_sorted(b), lm_sorted(b)]) = 0.5;
        contrasts.(sprintf('Outcome_bin%02d', b)) = cvec;
    end
end

% Target
if ~isempty(conditions.Target.Bins)
    [~, tgt_sorted] = sort_bins(conditions.Target.Bins);
    for b = 1:Nbins
        cvec = zeros(1, nReg);
        cvec(tgt_sorted(b)) = 1;
        contrasts.(sprintf('Target_bin%02d', b)) = cvec;
    end
end


end

%=============================================================
% Local helper
% =============================================================
function [binNumbers, sortedIdx] = sort_by_bin(SPM, idx)
% Sort FIR regressors (idx) by bin number bf(1)...bf(N)
    names = SPM.xX.name(idx);
    binTokens = regexp(names, 'bf\((\d+)\)', 'tokens');
    binNumbers = cellfun(@(x) str2double(x{1}), binTokens);
    [~, sortIdx] = sort(binNumbers);
    sortedIdx = idx(sortIdx);
end
