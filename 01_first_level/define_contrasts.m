function contrasts = define_contrasts(dirs,SPMmat,settings)
% DEFINE_CONTRASTS  Task contrast vectors for the canonical-HRF MID model.
%
% contrasts = define_contrasts(dirs, SPMmat, settings)
%
% Regressors are located BY NAME in SPM.xX.name, so the vectors are
% correct regardless of regressor order or the number of nuisance
% regressors. Contrasts:
%
%   Target      target onset (motor response / response preparation)
%   Cue         mean of high and low reward cue (anticipation)
%   HCvsLC      high - low reward cue (reward anticipation)
%   HitvsMiss   hits - misses, averaged over reward level (feedback)
%   TvsCue      target - cue
%   CuevsOut    cue - feedback
%   HHvsLH      high - low reward hit
%   HMvsLM      high - low reward miss
%   buttonPress button-press regressor (sanity check)
% plus the parametric-modulator contrasts when settings.todo.parameterModulator
% is on.
%
% OUTPUT
%   contrasts : struct, one field per contrast holding its weight vector
% -------------------------------------------------------------------------
% Danae Papadopetraki, 2019-2026
% Initial structure based on a lab template by Ruben van den Bosch (2017)
% =========================================================================
cd(dirs)
load(SPMmat)

% Regressor names used in the contrasts
% -------------------------------------------------------------------------
regs{1} = 'LowCue';
regs{2} = 'HighCue';
regs{3} = 'Target';
regs{4} = 'LowHit';
regs{5} = 'HighHit';
regs{6} = 'LowMiss';
regs{7} = 'HighMiss';
regs{8} = 'Resp';

if settings.todo.parameterModulatorAt.outcome
    regs{9} = 'LowHitxSurprise^1';
    regs{10} = 'HighHitxSurprise^1';
    regs{11} = 'LowMissxSurprise^1';
    regs{12} = 'HighMissxSurprise^1';
end

if settings.todo.parameterModulatorAt.target
    regs{9} = 'TargetxP_hit^1';
end

% Zero vector with one entry per regressor
% -------------------------------------------------------------------------
baseCon = zeros(1,numel(find(strncmp('Sn(',SPM.xX.name,3))));

TvsCue      = baseCon;
Target      = baseCon;
Cue         = baseCon;
HCvsLC      = baseCon;
HitvsMiss   = baseCon;
CuevsOut    = baseCon;
HHvsLH      = baseCon;
HMvsLM      = baseCon;
buttonPress = baseCon;
HitvsMiss_surprise = baseCon;
HHvsLH_surprise = baseCon;
HMvsLM_surprise = baseCon;
Target_p_hit = baseCon;

% Column index of each regressor. SPM names them 'Sn(1) <name>*bf(1)'.
% '^' (from the parametric-modulator names) is not allowed in a struct
% field name, so it is replaced by '_'.
% ---------------------------------------------------------------------
for i = 1:numel(regs)
    safeName = strrep(regs{i}, '^', '_');
    idx = find(contains(SPM.xX.name, ['Sn(1) ' regs{i} '*bf(1)']));
    pos.(safeName) = idx;
end

% Contrast weights
% =====================================================================
% Target vs cue
TvsCue(pos.Target) = 1;
TvsCue(pos.HighCue) = -0.5;
TvsCue(pos.LowCue) = -0.5;

% Target
Target(pos.Target) = 1;

% Cue (anticipation across reward levels)
Cue(pos.HighCue) = 0.5;
Cue(pos.LowCue) = 0.5;

% High vs low reward cue
HCvsLC(pos.HighCue) = 1;
HCvsLC(pos.LowCue) = -1;

% Hit vs miss feedback
HitvsMiss(pos.LowHit) = 0.5;
HitvsMiss(pos.LowMiss) = -0.5;
HitvsMiss(pos.HighHit) = 0.5;
HitvsMiss(pos.HighMiss) = -0.5;

% Cue vs feedback
CuevsOut(pos.HighCue) = 0.5;
CuevsOut(pos.LowCue) = 0.5;
CuevsOut(pos.LowHit) = -0.25;
CuevsOut(pos.LowMiss) = -0.25;
CuevsOut(pos.HighHit) = -0.25;
CuevsOut(pos.HighMiss) = -0.25;

% High vs low reward hit
HHvsLH(pos.HighHit) = 1;
HHvsLH(pos.LowHit) = -1;

% High vs low reward miss
HMvsLM(pos.HighMiss) = 1;
HMvsLM(pos.LowMiss) = -1;

% Button presses
buttonPress(pos.Resp) = 1;

if settings.todo.parameterModulatorAt.outcome
    % Hit vs miss on the parametric modulator
    HitvsMiss_surprise(pos.LowHitxSurprise_1) = 0.5;
    HitvsMiss_surprise(pos.LowMissxSurprise_1) = -0.5;
    HitvsMiss_surprise(pos.HighHitxSurprise_1) = 0.5;
    HitvsMiss_surprise(pos.HighMissxSurprise_1) = -0.5;
end

if settings.todo.parameterModulatorAt.target
    Target_p_hit(pos.TargetxP_hit_1) = 1;
end

% Output structure. Field order sets the contrast numbering in SPM.mat.
% =========================================================================
contrasts.TvsCue   = TvsCue;
contrasts.Target   = Target;
contrasts.HCvsLC   = HCvsLC;
contrasts.Cue      = Cue;
contrasts.HitvsMiss = HitvsMiss;
contrasts.CuevsOut = CuevsOut;
contrasts.HHvsLH = HHvsLH;
contrasts.HMvsLM = HMvsLM;
contrasts.buttonPress = buttonPress;

if settings.todo.parameterModulatorAt.outcome
    contrasts.HitvsMiss_surprise = HitvsMiss_surprise;
end

if settings.todo.parameterModulatorAt.target
    contrasts.Target_p_hit = Target_p_hit;
end
end
