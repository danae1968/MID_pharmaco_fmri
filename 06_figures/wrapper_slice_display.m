% =========================================================================
% WRAPPER_SLICE_DISPLAY  Dual-coded slice figures of the group-level results.
%
% Renders each (analysis, task contrast) as a dual-coded map (Allen et al.,
% 2012): the contrast estimate is shown in colour, the |t| statistic as
% opacity, and the significance mask as a contour (black = positive,
% white = negative effect). Unlike a thresholded map this also shows the
% sub-threshold evidence around significant clusters.
%
% Loops over every analysis directory in includeDirs.
%
%   corr_method : 'fwe'  whole-brain FWE (p < .05)
%                 'svc'  small-volume corrected FWE (p < .05)
%                 'tfce' TFCE permutation FWE (p < .05)
%   hrf_type    : 'canonical_hrf' | 'fir'
%
% Contrasts are resolved BY NAME in each folder's SPM.mat, so the con/spmT
% images and mask names follow whatever index the contrast has per model.
%
% Rendering uses sd_display_multi.m (several planes per figure, fixed
% physical size, vector PDF output) on top of the slice_display toolbox
% (Zandbelt, https://github.com/bramzandbelt/slice_display), and
% sd_pick_slices.m to choose representative slices.
% =========================================================================

% ---- USER CONFIG --------------------------------------------------------
corr_method    = 'tfce';              % 'fwe' | 'svc' | 'tfce'
hrf_type       = 'canonical_hrf';     % 'canonical_hrf' | 'fir'
% Default planes, used for any contrast without its own entry in figSpec
% below. Any of 'axial', 'coronal', 'sagittal'.
orientations   = {'coronal','sagittal'};

% true  : one figure per contrast containing that contrast's planes stacked
%         as row-blocks (one shared colour bar). Filename prefix names them,
%         e.g. 'sagittal-axial_'.
% false : one figure per plane.
combine_orientations = true;

% ---- How slices are chosen when figSpec does not name them ---------------
% 'auto' : sd_pick_slices scores every slice containing significant voxels
%          and picks `n_slices` of them, at least `min_slice_sep_mm` apart.
% 'all'  : every slice that contains a significant voxel.
slice_mode       = 'auto';
n_slices         = 4;        % slices per plane, for 'auto'
slice_metric     = 'peaks';  % 'peaks'  -> slices through the strongest |t|
                             % 'extent' -> slices holding the most voxels
                             % 'spread' -> evenly spaced over the range
min_slice_sep_mm = 10;       % keep chosen slices this far apart

% Print the candidate slices for ALL THREE planes (mm, number of significant
% voxels, peak |t|), with the automatic choice marked '*'. Run once with this
% on to see which plane shows a given effect best and which slices to use,
% then write the winners into figSpec.
print_candidates = false;

% ---- Per-contrast figure recipe ------------------------------------------
% Optional. Everything here is an override; any contrast you leave out keeps
% the defaults above, so a partial figSpec is fine.
%
% Key is <contrast>_<condition>, exactly as in the folder names. Under it:
%   .orientations  which plane(s) to show for this effect
%   .<plane>       the exact slices in mm for that plane (omit to let
%                  slice_mode choose them)
%   .thresholds    'both' | 'corrected' | 'unc' for this contrast only,
%                  overriding thresholds_to_plot
%   .same_row      true/false for this contrast only, overriding same_row
%   .title         replaces the auto-generated title, e.g.
%                  'Anticipation: high vs low reward (session average)'.
%                  The FILENAME still uses the auto-generated description,
%                  so output names stay predictable.
%
%   % striatal reward effect: coronal reads best
%   figSpec.sessionAverage_HCvsLC.orientations = {'coronal'};
%   figSpec.sessionAverage_HCvsLC.coronal      = [-2 6 14 22];
%
%   % midline frontal / PCC drug effect: sagittal, plus an axial row
%   figSpec.MPHvsSUL_Target.orientations = {'sagittal','axial'};
%   figSpec.MPHvsSUL_Target.sagittal     = [-8 -4 0 4];
%   % (no .axial line -> those 4 slices still get picked automatically)
%
% An entry applies in every analysis directory. To pin something for one
% covariate model only, nest it under that directory name first -- that is
% checked before the plain key:
%
%   figSpec.covariate_ki_VS_new.MPHvsPBO_Target.orientations = {'coronal'};
%
% The uncorrected and corrected figure of a contrast share an entry, which
% keeps the two directly comparable.
figSpec = struct();

% ---- Figure recipes used for the paper -----------------------------------
% Session averages: p < .001 uncorrected contour. These masks are also the
% functional mask for the drug and covariate TFCE analyses.
figSpec.sessionAverage_Cue.thresholds   = 'unc';
figSpec.sessionAverage_Cue.orientations = {'sagittal','coronal','axial'};
figSpec.sessionAverage_Cue.sagittal     = 3;
figSpec.sessionAverage_Cue.coronal      = -43;
figSpec.sessionAverage_Cue.axial        = 39;
figSpec.sessionAverage_Cue.same_row     = true;
figSpec.sessionAverage_Cue.title        = 'Anticipation across reward (session average)';

figSpec.sessionAverage_HCvsLC.thresholds    = 'unc';
figSpec.sessionAverage_HCvsLC.orientations  = {'sagittal','coronal','axial'};
figSpec.sessionAverage_HCvsLC.sagittal      = 12;
figSpec.sessionAverage_HCvsLC.coronal       = 0;
figSpec.sessionAverage_HCvsLC.axial         = 0;
figSpec.sessionAverage_HCvsLC.same_row      = true;
figSpec.sessionAverage_HCvsLC.title         = 'Anticipation: high vs low reward (session average)';

figSpec.sessionAverage_Target.thresholds    = 'unc';
figSpec.sessionAverage_Target.orientations  = {'sagittal','coronal','axial'};
figSpec.sessionAverage_Target.sagittal      = 0;
figSpec.sessionAverage_Target.coronal       = -24;
figSpec.sessionAverage_Target.axial         = 0;
figSpec.sessionAverage_Target.same_row      = true;
figSpec.sessionAverage_Target.title         = 'Target (session average)';

figSpec.sessionAverage_HitvsMiss.thresholds   = 'unc';
figSpec.sessionAverage_HitvsMiss.orientations = {'sagittal','coronal','axial'};
figSpec.sessionAverage_HitvsMiss.sagittal     = 50;
figSpec.sessionAverage_HitvsMiss.coronal      = 14;
figSpec.sessionAverage_HitvsMiss.axial        = -4;
figSpec.sessionAverage_HitvsMiss.same_row     = true;
figSpec.sessionAverage_HitvsMiss.title        = 'Feedback: hit vs miss (session average)';

% Drug comparisons: TFCE-corrected contour. A contrast can carry separate
% slices per threshold (.corrected / .unc).
figSpec.MPHvsPBO_Target.thresholds            = 'both';
figSpec.MPHvsPBO_Target.orientations          = {'sagittal','coronal','axial'};
figSpec.MPHvsPBO_Target.same_row              = true;
figSpec.MPHvsPBO_Target.corrected.sagittal    = 0;
figSpec.MPHvsPBO_Target.corrected.coronal     = -50;
figSpec.MPHvsPBO_Target.corrected.axial       = 24;
figSpec.MPHvsPBO_Target.corrected.title       = 'MPH vs PBO, Target (TFCE FWE)';
figSpec.MPHvsPBO_Target.unc.sagittal          = 0;
figSpec.MPHvsPBO_Target.unc.coronal           = 16;
figSpec.MPHvsPBO_Target.unc.axial             = 0;
figSpec.MPHvsPBO_Target.unc.title             = 'MPH vs PBO, Target (p < .001 uncorrected)';

figSpec.MPHvsPBO_HitvsMiss.thresholds   = 'corrected';
figSpec.MPHvsPBO_HitvsMiss.orientations = {'sagittal','coronal','axial'};
figSpec.MPHvsPBO_HitvsMiss.sagittal     = 6;
figSpec.MPHvsPBO_HitvsMiss.coronal      = -20;
figSpec.MPHvsPBO_HitvsMiss.axial        = 10;
figSpec.MPHvsPBO_HitvsMiss.same_row     = true;
figSpec.MPHvsPBO_HitvsMiss.title        = 'MPH vs PBO, Feedback (TFCE FWE)';

figSpec.MPHvsSUL_Target.thresholds   = 'corrected';
figSpec.MPHvsSUL_Target.orientations = {'sagittal','coronal','axial'};
figSpec.MPHvsSUL_Target.sagittal     = -4;
figSpec.MPHvsSUL_Target.coronal      = [-66 4];
figSpec.MPHvsSUL_Target.axial        = 38;
figSpec.MPHvsSUL_Target.same_row     = true;
figSpec.MPHvsSUL_Target.title        = 'MPH vs SUL, Target (TFCE FWE)';

figSpec.SULvsPBO_Target.thresholds   = 'corrected';
figSpec.SULvsPBO_Target.orientations = {'sagittal','coronal','axial'};
figSpec.SULvsPBO_Target.sagittal     = 52;
figSpec.SULvsPBO_Target.coronal      = 4;
figSpec.SULvsPBO_Target.axial        = 12;
figSpec.SULvsPBO_Target.same_row     = true;
figSpec.SULvsPBO_Target.title        = 'SUL vs PBO, Target (TFCE FWE)';

figSpec.MPHvsSUL_HitvsMiss.thresholds   = 'corrected';
figSpec.MPHvsSUL_HitvsMiss.orientations = {'sagittal','coronal','axial'};
figSpec.MPHvsSUL_HitvsMiss.sagittal     = 0;
figSpec.MPHvsSUL_HitvsMiss.coronal      = [-52 4];
figSpec.MPHvsSUL_HitvsMiss.axial        = 10;
figSpec.MPHvsSUL_HitvsMiss.same_row     = true;
figSpec.MPHvsSUL_HitvsMiss.title        = 'MPH vs SUL, Feedback (TFCE FWE)';

% ki covariate models (nested under their covariate directory)
figSpec.covariate_ki_wholePutamen.sessionAverage_Target.thresholds   = 'corrected';
figSpec.covariate_ki_wholePutamen.sessionAverage_Target.orientations = {'sagittal','coronal','axial'};
figSpec.covariate_ki_wholePutamen.sessionAverage_Target.sagittal     = -34;
figSpec.covariate_ki_wholePutamen.sessionAverage_Target.coronal      = 16;
figSpec.covariate_ki_wholePutamen.sessionAverage_Target.axial        = 50;
figSpec.covariate_ki_wholePutamen.sessionAverage_Target.same_row     = true;
figSpec.covariate_ki_wholePutamen.sessionAverage_Target.title        = 'Putamen ki x Target (session average, TFCE FWE)';

figSpec.covariate_ki_wholeStriatum.sessionAverage_Target.thresholds   = 'corrected';
figSpec.covariate_ki_wholeStriatum.sessionAverage_Target.orientations = {'sagittal','coronal','axial'};
figSpec.covariate_ki_wholeStriatum.sessionAverage_Target.sagittal     = -6;
figSpec.covariate_ki_wholeStriatum.sessionAverage_Target.coronal      = 50;
figSpec.covariate_ki_wholeStriatum.sessionAverage_Target.axial        = 10;
figSpec.covariate_ki_wholeStriatum.sessionAverage_Target.same_row     = true;
figSpec.covariate_ki_wholeStriatum.sessionAverage_Target.title        = 'Whole striatum ki x Target (session average, TFCE FWE)';

figSpec.covariate_ki_VS_new.sessionAverage_Target.thresholds   = 'corrected';
figSpec.covariate_ki_VS_new.sessionAverage_Target.orientations = {'sagittal','coronal','axial'};
figSpec.covariate_ki_VS_new.sessionAverage_Target.sagittal     = -10;
figSpec.covariate_ki_VS_new.sessionAverage_Target.coronal      = -52;
figSpec.covariate_ki_VS_new.sessionAverage_Target.axial        = 24;
figSpec.covariate_ki_VS_new.sessionAverage_Target.same_row     = true;
figSpec.covariate_ki_VS_new.sessionAverage_Target.title        = 'Ventral striatum ki x Target (session average, TFCE FWE)';

% Behavioural drug-effect covariate models
figSpec.covariate_mphBehEffect_hit.MPHvsPBO_Cue.thresholds   = 'corrected';
figSpec.covariate_mphBehEffect_hit.MPHvsPBO_Cue.orientations = {'sagittal','coronal','axial'};
figSpec.covariate_mphBehEffect_hit.MPHvsPBO_Cue.sagittal     = -24;
figSpec.covariate_mphBehEffect_hit.MPHvsPBO_Cue.coronal      = -49;
figSpec.covariate_mphBehEffect_hit.MPHvsPBO_Cue.axial        = 52;
figSpec.covariate_mphBehEffect_hit.MPHvsPBO_Cue.same_row     = true;
figSpec.covariate_mphBehEffect_hit.MPHvsPBO_Cue.title        = 'MPH hit-rate benefit x Anticipation cue (MPH vs PBO, TFCE FWE)';

figSpec.covariate_mphBehEffect_early.MPHvsPBO_Cue.thresholds   = 'corrected';
figSpec.covariate_mphBehEffect_early.MPHvsPBO_Cue.orientations = {'sagittal','coronal','axial'};
figSpec.covariate_mphBehEffect_early.MPHvsPBO_Cue.sagittal     = -24;
figSpec.covariate_mphBehEffect_early.MPHvsPBO_Cue.coronal      = -49;
figSpec.covariate_mphBehEffect_early.MPHvsPBO_Cue.axial        = 52;
figSpec.covariate_mphBehEffect_early.MPHvsPBO_Cue.same_row     = true;
figSpec.covariate_mphBehEffect_early.MPHvsPBO_Cue.title        = 'MPH premature-response change x Anticipation cue (MPH vs PBO, TFCE FWE)';

figSpec.covariate_mphBehEffect_late.MPHvsPBO_HitvsMiss.thresholds   = 'corrected';
figSpec.covariate_mphBehEffect_late.MPHvsPBO_HitvsMiss.orientations = {'sagittal','coronal','axial'};
figSpec.covariate_mphBehEffect_late.MPHvsPBO_HitvsMiss.sagittal     = -6;
figSpec.covariate_mphBehEffect_late.MPHvsPBO_HitvsMiss.coronal      = 53;
figSpec.covariate_mphBehEffect_late.MPHvsPBO_HitvsMiss.axial        = -5;
figSpec.covariate_mphBehEffect_late.MPHvsPBO_HitvsMiss.same_row     = true;
figSpec.covariate_mphBehEffect_late.MPHvsPBO_HitvsMiss.title        = 'MPH late-response change x Feedback (MPH vs PBO, TFCE FWE)';

figSpec.covariate_lateEffect_MPHvsSUL.MPHvsSUL_HitvsMiss.thresholds   = 'corrected';
figSpec.covariate_lateEffect_MPHvsSUL.MPHvsSUL_HitvsMiss.orientations = {'sagittal','coronal','axial'};
figSpec.covariate_lateEffect_MPHvsSUL.MPHvsSUL_HitvsMiss.sagittal     = -63;
figSpec.covariate_lateEffect_MPHvsSUL.MPHvsSUL_HitvsMiss.coronal      = -43;
figSpec.covariate_lateEffect_MPHvsSUL.MPHvsSUL_HitvsMiss.axial        = 3;
figSpec.covariate_lateEffect_MPHvsSUL.MPHvsSUL_HitvsMiss.same_row     = true;
figSpec.covariate_lateEffect_MPHvsSUL.MPHvsSUL_HitvsMiss.title        = 'MPH vs SUL late-response effect x Feedback (TFCE FWE)';

% ki x single-drug models
figSpec.covariate_ki_VS_new.MPH_Target.thresholds   = 'corrected';
figSpec.covariate_ki_VS_new.MPH_Target.orientations = {'sagittal','coronal','axial'};
figSpec.covariate_ki_VS_new.MPH_Target.sagittal     = -20;
figSpec.covariate_ki_VS_new.MPH_Target.coronal      = 33;
figSpec.covariate_ki_VS_new.MPH_Target.axial        = 39;
figSpec.covariate_ki_VS_new.MPH_Target.same_row     = true;
figSpec.covariate_ki_VS_new.MPH_Target.title        = 'Ventral striatum ki x Target, under MPH (TFCE FWE)';

figSpec.covariate_ki_wholeStriatum.MPH_Target.thresholds   = 'corrected';
figSpec.covariate_ki_wholeStriatum.MPH_Target.orientations = {'sagittal','coronal','axial'};
figSpec.covariate_ki_wholeStriatum.MPH_Target.sagittal     = -23;
figSpec.covariate_ki_wholeStriatum.MPH_Target.coronal      = 36;
figSpec.covariate_ki_wholeStriatum.MPH_Target.axial        = 36;
figSpec.covariate_ki_wholeStriatum.MPH_Target.same_row     = true;
figSpec.covariate_ki_wholeStriatum.MPH_Target.title        = 'Whole striatum ki x Target, under MPH (TFCE FWE)';

figSpec.covariate_ki_wholeCaudate.MPH_Target.thresholds   = 'corrected';
figSpec.covariate_ki_wholeCaudate.MPH_Target.orientations = {'sagittal','coronal','axial'};
figSpec.covariate_ki_wholeCaudate.MPH_Target.sagittal     = -23;
figSpec.covariate_ki_wholeCaudate.MPH_Target.coronal      = 36;
figSpec.covariate_ki_wholeCaudate.MPH_Target.axial        = 36;
figSpec.covariate_ki_wholeCaudate.MPH_Target.same_row     = true;
figSpec.covariate_ki_wholeCaudate.MPH_Target.title        = 'Caudate ki x Target, under MPH (TFCE FWE)';

figSpec.covariate_ki_wholeCaudate.MPH_HCvsLC.thresholds   = 'corrected';
figSpec.covariate_ki_wholeCaudate.MPH_HCvsLC.orientations = {'sagittal','coronal','axial'};
figSpec.covariate_ki_wholeCaudate.MPH_HCvsLC.sagittal     = 6;
figSpec.covariate_ki_wholeCaudate.MPH_HCvsLC.coronal      = 10;
figSpec.covariate_ki_wholeCaudate.MPH_HCvsLC.axial        = 59;
figSpec.covariate_ki_wholeCaudate.MPH_HCvsLC.same_row     = true;
figSpec.covariate_ki_wholeCaudate.MPH_HCvsLC.title        = 'Caudate ki x Anticipation, under MPH (TFCE FWE)';

figSpec.covariate_ki_wholePutamen.SUL_Target.thresholds   = 'corrected';
figSpec.covariate_ki_wholePutamen.SUL_Target.orientations = {'sagittal','coronal','axial'};
figSpec.covariate_ki_wholePutamen.SUL_Target.sagittal     = 6;
figSpec.covariate_ki_wholePutamen.SUL_Target.coronal      = 59;
figSpec.covariate_ki_wholePutamen.SUL_Target.axial        = -2;
figSpec.covariate_ki_wholePutamen.SUL_Target.same_row     = true;
figSpec.covariate_ki_wholePutamen.SUL_Target.title        = 'Putamen ki x Target, under SUL (TFCE FWE)';

figSpec.covariate_ki_wholeStriatum.SUL_Target.thresholds   = 'corrected';
figSpec.covariate_ki_wholeStriatum.SUL_Target.orientations = {'sagittal','coronal','axial'};
figSpec.covariate_ki_wholeStriatum.SUL_Target.sagittal     = 59;
figSpec.covariate_ki_wholeStriatum.SUL_Target.coronal      = -3;
figSpec.covariate_ki_wholeStriatum.SUL_Target.axial        = -5;
figSpec.covariate_ki_wholeStriatum.SUL_Target.same_row     = true;
figSpec.covariate_ki_wholeStriatum.SUL_Target.title        = 'Whole striatum ki x Target, under SUL (TFCE FWE)';

% tile_width_mm is 39 so that the 4-tile figures (both MPHvsSUL panels and
% the Cue session average) stay inside 180 mm; 3-tile figures are
% unaffected beyond being marginally narrower.

% Output location within each analysis's figure directory:
%   slice_mode = 'all' -> the analysis directory itself
%   otherwise          -> figureSubfolder, so representative-slice figures
%                         and all-slice figures do not overwrite each other
% figureSubfolder = '' puts everything in the analysis directory.
figureSubfolder = 'selected_slices';

% Optional: write every figure into this single folder (file names are
% unique). Empty = one folder per analysis directory.
flat_figure_dir = '';

% Anatomical underlay. 'auto' = the sample's group-average T1
% (bids/derivatives/mri/fmriprep/group/group_average_T1scan.nii), falling
% back to SPM's single-subject canonical T1. Or give a full path.
anat_image = 'auto';

% ---- Dual-coded opacity ceiling -----------------------------------------
% The opacity ceiling should be the statistical threshold the contour marks,
% so voxels inside the outline are fully opaque and evidence fades below it
% (Allen et al., 2012).
%
% 'auto'  : per model, the t for p = opacity_p uncorrected at that model's
%           error df, i.e. the p < .001 threshold. On TFCE figures, voxels
%           significant only through TFCE's extent term stay partly faded,
%           which is informative.
% numeric : the same ceiling for every figure.
opacity_ceiling = 'auto';
opacity_p       = 0.001;

% Colour axis of the dual-coded map.
%   'robust' (default) : +/- the color_range_pct-th percentile of |con|
%                        among the contrast's significant voxels. Robust to
%                        single extreme voxels, which otherwise stretch the
%                        colour axis so far that most of the map looks pale.
%   'auto'             : toolbox behaviour (most extreme voxel).
%   scalar or [lo hi]  : fixed, e.g. 1.5 -> [-1.5 1.5], so every figure
%                        shares one colour bar.
% The range is taken from the CORRECTED mask, so a contrast's uncorrected
% and corrected figures stay on the same scale.
color_range     = 'robust';
color_range_pct = 95;

% Trilinear instead of nearest-neighbour resampling onto the 1 mm
% anatomical grid (display only). Applied to maps and contour masks alike
% so the outline stays aligned with the coloured blob; the 0.5 isocontour
% of the interpolated mask follows the original voxel boundary.
smooth_display  = false;

% |t| below this is drawn fully transparent, to reduce clutter
% (van den Bosch et al., 2022, use 1). 0 = no floor.
opacity_floor   = 0;

% Black outline around positive effects, white around negative ones
% (black is hard to see on blue). Falls back to one black outline if the
% per-direction masks are missing.
contour_by_sign = true;

% Which threshold to render per contrast:
%   'both'      uncorrected p<.001 AND the corrected mask (default)
%   'corrected' corrected only
%   'unc'       uncorrected only
% Overridable per contrast with figSpec.<contrast>.thresholds.
thresholds_to_plot = 'both';

% true : all planes of a contrast on ONE row, tiles sharing a height and
%        taking their width from each plane's aspect ratio.
% false: one row-block per plane, stacked (default).
% Overridable per contrast with figSpec.<contrast>.same_row. Note this puts
% EVERY slice of every plane on that row, so it suits pinned one- or
% two-slice figures, not a 4-slice automatic one (8 tiles = ~340 mm wide).
same_row = false;

% ---- Figure geometry -----------------------------------------------------
fig_columns    = 4;    % max slices per row
tile_width_mm  = 39;   % width of one slice tile -> 4 tiles fit inside 180 mm
fig_dpi        = 300;  % print resolution for PNG
% Output formats: {'png'}, {'pdf'}, or both. PDF is vector -- editable text
% and vector contours, with the slices as embedded rasters -- so it drops
% straight into Illustrator.
fig_formats    = {'png'};
colorbar_label = 'beta (a.u.)';   % short label under the colour bar

% Figures are built off-screen; output size does not depend on the screen.
fig_visible    = 'off';

% Analysis directories to process. If empty, every covariate_* directory
% under groupLevelBase is found automatically (minus excludeDirs).
includeDirs = {'no_covariate_with_exc', ...
               'covariate_ki_VS_new', 'covariate_ki_wholeCaudate', ...
               'covariate_ki_wholePutamen', 'covariate_ki_wholeStriatum', ...
               'covariate_mphBehEffect_early', 'covariate_mphBehEffect_hit', ...
               'covariate_mphBehEffect_rt', 'covariate_mphBehEffect_late', ...
               'covariate_rtEffect_MPHvsSUL', 'covariate_hitEffect_MPHvsSUL', ...
               'covariate_lateEffect_MPHvsSUL', 'covariate_SULvsPBO_late', ...
               'covariate_TauEffect_MPHvsSUL', 'covariate_tauEffect_MPHvsPBO'};

% Directories to skip in auto-discovery mode
excludeDirs = {};

% Overwrite existing figures? false only fills in what is missing.
force = true;
% ---- Local mode ----------------------------------------------------------
% Run against a local copy of the results instead of the project storage,
% e.g. to iterate on figure layout. Leave empty for the project layout.
%   local.groupLevel : folder that directly contains the analysis directories
%   local.figureDir  : output folder for the figures
%   local.extraPaths : SPM, slice_display and panel toolbox folders
%   local.anat       : anatomical underlay
local.groupLevel = '';
local.figureDir  = '';
local.extraPaths = {};
local.anat       = '';   % anatomical underlay in the local copy
% -------------------------------------------------------------------------

localMode = ~isempty(local.groupLevel);

if ispc
    projectDir = 'P:\3017048.01\';
    spmDir     = '\\fileserver.dccn.nl\home\common\matlab\spm12';
elseif isunix
    projectDir = '/project/3017048.01/';
    spmDir     = '/home/common/matlab/spm12';
else
    error('Unknown OS')
end

% sd_display_multi / sd_pick_slices live next to this script
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

if localMode
    for iP = 1:numel(local.extraPaths)
        addpath(genpath(local.extraPaths{iP}));
    end
    % colormaps.mat ships in the slice_display folder next to this script
    code_dir = fullfile(scriptDir, 'slice_display');
    addpath(genpath(code_dir));
    groupLevelBase = local.groupLevel;
else
    addpath(spmDir);
    code_dir  = fullfile(projectDir,'code','analysis','mri','2ndLevel','MID','slice_display');
    panelDir  = fullfile(projectDir,'code','analysis','slice_display');
    addpath(genpath(code_dir));
    addpath(genpath(panelDir));
    groupLevelBase = fullfile(projectDir,'bids','derivatives','mri','groupLevel','MID', ...
        'regular', hrf_type);
end

% Load custom colormaps (used by all dual layers)
load(fullfile(code_dir,'colormaps.mat'));

% ---- Resolve the anatomical underlay -------------------------------------
if localMode && ~isempty(local.anat)
    anatImg = local.anat;
elseif strcmpi(anat_image, 'auto')
    anatImg = fullfile(projectDir,'bids','derivatives','mri','fmriprep', ...
                       'group','group_average_T1scan.nii');
else
    anatImg = anat_image;
end
if ~exist(anatImg, 'file')
    warning(['Anatomical underlay not found (%s); falling back to SPM''s ' ...
             'single-subject canonical.'], anatImg);
    anatImg = fullfile(spm('Dir'),'canonical','single_subj_T1.nii');
end
fprintf('Anatomical underlay: %s\n', anatImg);

% ---- Discover which top-level analysis dirs to process --------------------
% Either exactly includeDirs, or (no_covariate_with_exc + every covariate_*
% dir present under groupLevelBase) minus excludeDirs.
topDirs = {};
if ~isempty(includeDirs)
    for i = 1:numel(includeDirs)
        if exist(fullfile(groupLevelBase, includeDirs{i}), 'dir')
            topDirs{end+1} = includeDirs{i}; %#ok<SAGROW>
        else
            warning('includeDirs entry not found, skipping: %s', includeDirs{i});
        end
    end
else
    d = dir(groupLevelBase);
    d = d([d.isdir]);
    for i = 1:numel(d)
        if strcmp(d(i).name, 'no_covariate_with_exc') || startsWith(d(i).name, 'covariate_')
            topDirs{end+1} = d(i).name; %#ok<SAGROW>
        end
    end
    topDirs = setdiff(topDirs, excludeDirs, 'stable');
end
if isempty(topDirs)
    error('No analysis directories found under %s', groupLevelBase);
end
fprintf('Found %d analysis directories to process:\n', numel(topDirs));
fprintf('  %s\n', strjoin(topDirs, '\n  '));

% Short model label for figure titles/filenames (e.g. 'canonical', 'fir')
hrfLabel = regexprep(hrf_type, '_hrf$', '');

% ---- Loop dimensions ------------------------------------------------------
contrasts  = {'sessionAverage','MPHvsPBO','MPHvsSUL','SULvsPBO'};
conditions = {'Cue','HCvsLC','HitvsMiss','Target'};

contrastLabel = containers.Map( ...
    {'sessionAverage','MPHvsPBO','MPHvsSUL','SULvsPBO'}, ...
    {'Session average','MPH vs PBO','MPH vs SUL','SUL vs PBO'} );

% Positive task contrast name depends on the HRF model:
%   canonical: '<COND>'        fir: 'Mean_<COND>'
switch hrf_type
    case 'fir',           con1 = @(cond) sprintf('Mean_%s', cond);
    case 'canonical_hrf', con1 = @(cond) cond;
    otherwise, error('Unknown hrf_type: %s', hrf_type);
end

axisRow = struct('sagittal', 1, 'coronal', 2, 'axial', 3);
allOrientations = {'axial','coronal','sagittal'};   % for print_candidates

% ---- Loop over every top-level analysis dir --------------------------------
for iT = 1:numel(topDirs)
    topDir = topDirs{iT};
    use_covariate = ~strcmp(topDir, 'no_covariate_with_exc');
    if use_covariate
        covariate_name = erase(topDir, 'covariate_');
        effectLabel    = covariate_name;
    else
        covariate_name = '';
        effectLabel    = 'main';
    end

    data_dir = fullfile(groupLevelBase, topDir);

    figSubdir = topDir;
    switch corr_method
        case 'svc',  figSubdir = [figSubdir '_SVC'];
        case 'tfce', figSubdir = [figSubdir '_TFCE'];
    end
    if ~isempty(flat_figure_dir)
        figureDir = flat_figure_dir;
    elseif localMode
        figureDir = fullfile(local.figureDir, hrf_type, figSubdir);
    else
        figureDir = fullfile(projectDir,'bids','derivatives','mri','groupLevel','MID', ...
            'figure_collections','slice_display', hrf_type, figSubdir);
    end
    if isempty(flat_figure_dir) && ~strcmpi(slice_mode, 'all') && ~isempty(figureSubfolder)
        figureDir = fullfile(figureDir, figureSubfolder);
    end
    if ~exist(figureDir,'dir')
        mkdir(figureDir);
    end

    if use_covariate
        fprintf('\n=== %s (covariate: %s) ===\n', topDir, covariate_name);
    else
        fprintf('\n=== %s (no covariate) ===\n', topDir);
    end

    for c = 1:numel(contrasts)
        for k = 1:numel(conditions)
            contrastName = contrasts{c};
            condName     = conditions{k};
            folderName   = sprintf('%s_%s', contrastName, condName);
            contrast_dir = fullfile(data_dir, folderName);

            if ~exist(contrast_dir,'dir')
                warning('Skipping (folder missing): %s', folderName);
                continue
            end

            spmPath = fullfile(contrast_dir, 'SPM.mat');
            if ~exist(spmPath,'file')
                warning('Skipping (SPM.mat missing): %s', folderName);
                continue
            end

            % ---- Resolve the target contrast BY NAME -----------------------
            S = load(spmPath, 'SPM');
            conNms = strrep({S.SPM.xCon.name}, ' ', '_');

            % Opacity ceiling for this model (see USER CONFIG above).
            if ischar(opacity_ceiling) && strcmpi(opacity_ceiling, 'auto')
                if isfield(S.SPM, 'xX') && isfield(S.SPM.xX, 'erdf') && ...
                        isscalar(S.SPM.xX.erdf) && isfinite(S.SPM.xX.erdf) && S.SPM.xX.erdf > 0
                    % SPM's uncorrected T threshold is one-tailed, matching
                    % how the p<.001 masks were written.
                    tCeiling = spm_invTcdf(1 - opacity_p, S.SPM.xX.erdf);
                else
                    warning(['No usable error df in %s; falling back to 5.77 ' ...
                        'for the opacity ceiling.'], spmPath);
                    tCeiling = 5.77;
                end
            else
                tCeiling = opacity_ceiling;
            end

            baseName = con1(condName);                 % e.g. 'Cue' or 'Mean_Cue'
            if use_covariate
                interBase  = condName;
                pat  = sprintf('^%s_x_%s$', regexptranslate('escape', interBase), ...
                    regexptranslate('escape', covariate_name));
                targetName = sprintf('%s_x_%s', interBase, covariate_name);
            else
                pat  = sprintf('^%s$', regexptranslate('escape', baseName));
                targetName = baseName;
            end

            idx = find(~cellfun('isempty', regexp(conNms, pat, 'once')), 1);
            if isempty(idx)
                warning('Skipping (no contrast matching "%s"): %s', targetName, folderName);
                continue
            end
            targetName = conNms{idx};

            % Index of the matching negative contrast, used to name the
            % per-direction uncorrected/FWE masks (which are written with
            % the contrast index in the filename, unlike the TFCE ones).
            patNeg = sprintf('^negative_%s$', regexptranslate('escape', targetName));
            idxNeg = find(~cellfun('isempty', regexp(conNms, patNeg, 'once')), 1);

            conImg = fullfile(contrast_dir, sprintf('con_%04d.nii',  idx));
            tImg   = fullfile(contrast_dir, sprintf('spmT_%04d.nii', idx));
            if ~exist(conImg,'file') || ~exist(tImg,'file')
                warning('Skipping (con/T image missing): %s', folderName);
                continue
            end

            % ---- Threshold masks (resolved from targetName + idx) ----------
            % The uncorrected mask exists as soon as a model is estimated;
            % for 'tfce'/'svc' the corrected mask is required as well, so a
            % model that was never corrected does not silently produce an
            % uncorrected-only figure.
            clear thresholds
            negIdxStr = '';
            if ~isempty(idxNeg), negIdxStr = sprintf('%04d', idxNeg); end
            switch corr_method
                case 'svc'
                    correctedFile = fullfile(contrast_dir, ...
                        sprintf('spmT_%04d_SVC_significant_voxels_%s_FWE_p05.nii', idx, targetName));
                    correctedLabel = 'SVC FWE';
                    corrPos = ''; corrNeg = '';       % SVC writes one file only
                case 'tfce'
                    posMask = fullfile(contrast_dir, ...
                        sprintf('TFCE_significant_voxels_%s_fwe_p05.nii', targetName));
                    negMask = fullfile(contrast_dir, ...
                        sprintf('TFCE_significant_voxels_negative_%s_fwe_p05.nii', targetName));
                    correctedFile = fullfile(contrast_dir, ...
                        sprintf('TFCE_significant_voxels_%s_bothdir_fwe_p05.nii', targetName));
                    make_union_mask({posMask, negMask}, correctedFile);
                    correctedLabel = 'TFCE FWE';
                    corrPos = posMask; corrNeg = negMask;
                otherwise  % 'fwe'
                    correctedFile = fullfile(contrast_dir, ...
                        sprintf('significant_voxels_%s_fwe_p05.nii', targetName));
                    correctedLabel = 'FWE';
                    corrPos = fullfile(contrast_dir, ...
                        sprintf('spmT_%04d_significant_voxels_%s_fwe_p05.nii', idx, targetName));
                    corrNeg = '';
                    if ~isempty(negIdxStr)
                        corrNeg = fullfile(contrast_dir, sprintf( ...
                            'spmT_%s_significant_voxels_negative_%s_fwe_p05.nii', negIdxStr, targetName));
                    end
            end

            if ~strcmp(corr_method, 'fwe') && ~exist(correctedFile, 'file')
                warning('Skipping (no %s output yet, only uncorrected exists): %s', ...
                    corr_method, folderName);
                continue
            end

            uncPos = fullfile(contrast_dir, sprintf( ...
                'spmT_%04d_significant_voxels_%s_uncorrected_p001.nii', idx, targetName));
            uncNeg = '';
            if ~isempty(negIdxStr)
                uncNeg = fullfile(contrast_dir, sprintf( ...
                    'spmT_%s_significant_voxels_negative_%s_uncorrected_p001.nii', negIdxStr, targetName));
            end

            thresholds(1) = struct('label','unc', 'file', fullfile(contrast_dir, ...
                sprintf('significant_voxels_%s_uncorrected_p001.nii', targetName)), ...
                'pos', uncPos, 'neg', uncNeg);
            thresholds(2) = struct('label', correctedLabel, 'file', correctedFile, ...
                'pos', corrPos, 'neg', corrNeg);

            % Keep only the requested threshold(s) for this contrast.
            whichThresh = spec_lookup(figSpec, topDir, folderName, 'thresholds');
            if isempty(whichThresh), whichThresh = thresholds_to_plot; end
            switch lower(whichThresh)
                case 'corrected', thresholds = thresholds(2);
                case 'unc',       thresholds = thresholds(1);
                case 'both'       % keep both
                otherwise
                    error('Unknown thresholds setting "%s".', whichThresh);
            end


            for th = 1:numel(thresholds)
                threshLab = thresholds(th).label;
                maskPath  = thresholds(th).file;

                if ~exist(maskPath,'file')
                    warning('Skipping (mask missing): %s', maskPath);
                    continue
                end

                % figSpec key for this threshold, so a contrast can carry a
                % different slice set / title per threshold.
                if strcmpi(threshLab, 'unc'), threshKey = 'unc';
                else,                          threshKey = 'corrected'; end

                % ---- Colour range, from the mask this figure outlines ------
                % Taken from the drawn mask rather than always the corrected
                % one, so the scale matches the region the contour encloses.
                % A contrast rendered at both thresholds therefore gets two
                % scales -- each honest about what it is showing.
                if ischar(color_range)
                    switch lower(color_range)
                        case 'robust'
                            cRange = robust_con_range(conImg, maskPath, color_range_pct);
                        case 'auto'
                            cRange = [];   % sd_config_layers derives it
                        otherwise
                            error('Unknown color_range "%s".', color_range);
                    end
                elseif isscalar(color_range)
                    cRange = [-abs(color_range) abs(color_range)];
                else
                    cRange = color_range;
                end

                % Descriptive string: figure title + filename stem.
                descr = sprintf('%s %s %s at %s (%s)', ...
                    contrastLabel(contrastName), hrfLabel, effectLabel, condName, threshLab);

                % ---- Contour layers: one per direction, or one union ---------
                % Black outlines around positive effects, white around
                % negative, so the outline stays visible on blue blobs.
                % Whichever per-direction masks exist get used; if neither
                % does, fall back to a single black outline over the
                % both-directions mask.
                contourFiles = {}; contourCols = {};
                if contour_by_sign
                    if ~isempty(thresholds(th).pos) && exist(thresholds(th).pos, 'file')
                        contourFiles{end+1} = thresholds(th).pos; %#ok<SAGROW>
                        contourCols{end+1}  = [0 0 0];            %#ok<SAGROW>
                    end
                    if ~isempty(thresholds(th).neg) && exist(thresholds(th).neg, 'file')
                        contourFiles{end+1} = thresholds(th).neg; %#ok<SAGROW>
                        contourCols{end+1}  = [1 1 1];            %#ok<SAGROW>
                    end
                end
                if isempty(contourFiles)
                    contourFiles = {maskPath};
                    contourCols  = {[0 0 0]};
                end

                % ---- Initialise layers ---------------------------------------
                layerTypes = [{'truecolor','dual'}, repmat({'contour'}, 1, numel(contourFiles))];
                layers   = sd_config_layers('init', layerTypes);

                % Layer 1: anatomical
                layers(1).color.file = anatImg;
                layers(1).color.map  = gray(256);

                % Layer 2: dual-coded map. The colour-bar label is kept short;
                % the long description goes in the figure title.
                layers(2).color.file    = conImg;
                layers(2).color.map     = CyBuGyRdYl;
                layers(2).color.range   = cRange;   % [] -> toolbox default
                layers(2).color.label   = colorbar_label;
                layers(2).opacity.file  = tImg;
                layers(2).opacity.label = '| t |';
                layers(2).opacity.range = [opacity_floor tCeiling];
                if smooth_display
                    layers(2).color.hold   = 1;   % trilinear, display only
                    layers(2).opacity.hold = 1;
                    for iC = 3:numel(layers)      % contour layers must match
                        layers(iC).color.hold = 1;
                    end
                end

                % Layers 3+: significance contours
                for iC = 1:numel(contourFiles)
                    layers(2+iC).color.file = contourFiles{iC};
                    layers(2+iC).color.map  = contourCols{iC};
                end

                % ---- Check there is anything to show -------------------------
                V = spm_vol(maskPath);
                [maskData, XYZmm] = spm_read_vols(V);
                sigIdx = find(maskData > 0);
                if isempty(sigIdx)
                    warning('No significant voxels in %s. Skipping.', maskPath);
                    continue
                end

                % ---- Choose planes and slices --------------------------------
                pickOpts = struct('metric', slice_metric, ...
                                  'stat_file', tImg, ...
                                  'min_sep_mm', min_slice_sep_mm, ...
                                  'round_to', 2);

                % Candidate tables for ALL THREE planes, not just the ones
                % about to be rendered -- the point is to be able to see
                % which plane shows THIS effect best before committing to it.
                if print_candidates
                    for iC = 1:numel(allOrientations)
                        [~, tbl] = sd_pick_slices(maskPath, allOrientations{iC}, ...
                                                  n_slices, pickOpts);
                        print_slice_candidates(folderName, threshLab, ...
                                               allOrientations{iC}, tbl, slice_metric);
                    end
                end

                % Planes for this contrast: its own figSpec entry if there is
                % one, otherwise the global default.
                oriList = spec_lookup(figSpec, topDir, folderName, 'orientations', threshKey);
                if isempty(oriList)
                    oriList = orientations;
                end
                if ischar(oriList), oriList = {oriList}; end

                blocks = struct('orientation', {}, 'slices', {});
                for iO = 1:numel(oriList)
                    ori = oriList{iO};
                    if ~isfield(axisRow, ori)
                        warning('Unknown orientation "%s", skipping.', ori);
                        continue
                    end

                    % Slices named explicitly for this contrast and plane?
                    sel = spec_lookup(figSpec, topDir, folderName, ori, threshKey);

                    if isempty(sel)
                        switch lower(slice_mode)
                            case 'all'
                                coords = XYZmm(axisRow.(ori), sigIdx);
                                sel    = unique(round(coords/2)*2);
                            case 'auto'
                                sel = sd_pick_slices(maskPath, ori, n_slices, pickOpts);
                            otherwise
                                error('Unknown slice_mode "%s".', slice_mode);
                        end
                    end

                    if isempty(sel)
                        warning('No slices selected for %s / %s.', folderName, ori);
                        continue
                    end
                    blocks(end+1) = struct('orientation', ori, 'slices', sel); %#ok<SAGROW>
                end
                if isempty(blocks)
                    continue
                end

                % ---- One combined figure, or one figure per orientation ------
                if combine_orientations
                    renderSets = {blocks};
                    % Sorted, so reordering `orientations` overwrites the same
                    % file instead of quietly creating a parallel set under a
                    % new name. Alphabetical is also the conventional order
                    % (axial, coronal, sagittal). The ROW order in the figure
                    % still follows `orientations` as written, so the filename
                    % names the planes, not their stacking order.
                    prefixes   = {strjoin(sort({blocks.orientation}), '-')};
                else
                    renderSets = arrayfun(@(b) b, blocks, 'UniformOutput', false);
                    prefixes   = {blocks.orientation};
                end

                for iR = 1:numel(renderSets)

                    fName   = sprintf('%s_%s', prefixes{iR}, regexprep(descr, ' ', '_'));
                    outPath = fullfile(figureDir, [fName '.png']);
                    if ~force && exist(outPath, 'file')
                        fprintf('  %s already exists, skipping (force=false)\n', fName);
                        continue
                    end

                    if isempty(cRange), rangeStr = 'auto';
                    else,                rangeStr = sprintf('%.2f', cRange(2)); end
                    fprintf('\n=== %s | con %d (%s) | %s | %s | |t| ceiling %.2f | colour +/-%s ===\n', ...
                        folderName, idx, targetName, threshLab, prefixes{iR}, tCeiling, rangeStr);
                    for ib = 1:numel(renderSets{iR})
                        fprintf('    %-9s slices: %s\n', renderSets{iR}(ib).orientation, ...
                            strjoin(arrayfun(@(v) sprintf('%d', round(v)), ...
                                renderSets{iR}(ib).slices, 'UniformOutput', false), ', '));
                    end

                    thisSameRow = spec_lookup(figSpec, topDir, folderName, 'same_row', threshKey);
                    if isempty(thisSameRow), thisSameRow = same_row; end

                    figTitle = spec_lookup(figSpec, topDir, folderName, 'title', threshKey);
                    if isempty(figTitle), figTitle = descr; end

                    dispOpts = struct('title',         figTitle, ...
                                      'n_columns',     fig_columns, ...
                                      'tile_width_mm', tile_width_mm, ...
                                      'dpi',           fig_dpi, ...
                                      'formats',       {fig_formats}, ...
                                      'visible',       fig_visible, ...
                                      'same_row',      thisSameRow, ...
                                      'out_file',      fullfile(figureDir, fName));

                    h = sd_display_multi(layers, renderSets{iR}, dispOpts);
                    close(h);
                end
            end
        end
    end
end

if ~isempty(flat_figure_dir)
    outLeaf = '';
elseif strcmpi(slice_mode, 'all') || isempty(figureSubfolder)
    outLeaf = '<analysis>';
else
    outLeaf = ['<analysis>' filesep figureSubfolder];
end
if localMode
    doneRoot = fullfile(local.figureDir, hrf_type, outLeaf);
else
    doneRoot = fullfile(projectDir,'bids','derivatives','mri','groupLevel','MID', ...
        'figure_collections','slice_display', hrf_type, outLeaf);
end
fprintf('\nDone. Figures saved under:\n  %s\n', doneRoot);

% =========================================================================
% Helper: symmetric colour range from a percentile of |con| inside a mask.
% Returns [] (meaning "let the toolbox decide") if anything is missing.
% Uses a sorted index rather than prctile to avoid needing the Statistics
% Toolbox.
% =========================================================================
function r = robust_con_range(conImg, maskFile, pct)
r = [];
if ~exist(conImg, 'file') || isempty(maskFile) || ~exist(maskFile, 'file'), return; end
Vc = spm_vol(conImg);  Yc = spm_read_vols(Vc);
Vm = spm_vol(maskFile); Ym = spm_read_vols(Vm);
if ~isequal(Vc.dim, Vm.dim), return; end
v = abs(Yc(Ym > 0));
v = v(isfinite(v) & v > 0);

% A percentile of only a few voxels is just their maximum, which can be so
% small that the rest of the map saturates. With fewer than 10 significant
% voxels, use a high percentile of the whole map instead.
if numel(v) < 10
    v   = abs(Yc(:));
    v   = v(isfinite(v) & v > 0);
    pct = 99;
end

if isempty(v), return; end
v = sort(v);
q = v(min(numel(v), max(1, ceil(pct/100 * numel(v)))));
if ~isfinite(q) || q <= 0, return; end
r = [-q q];
end

% =========================================================================
% Helper: dump every candidate slice for one contrast/orientation, so a set
% can be hand-picked into `figSpec`. '*' marks the automatic choice.
% =========================================================================
function print_slice_candidates(folderName, threshLab, orientation, tbl, metric, maxRows)
if isempty(tbl), return; end
if nargin < 6 || isempty(maxRows), maxRows = 12; end

% A session-average contrast can span 60+ slices per plane, and printing
% every one for each contrast x threshold x plane buries the useful part.
% Show the best `maxRows` by the same score the picker uses (so the chosen
% slices are always among them) and say how many were left out.
switch lower(metric)
    case 'extent', score = [tbl.n_voxels];
    otherwise,     score = [tbl.max_t];
end
if ~any(score), score = [tbl.n_voxels]; end

show = [tbl.chosen];
if numel(tbl) > maxRows
    [~, ord] = sort(score, 'descend');
    show(ord(1:min(maxRows, numel(ord)))) = true;
else
    show(:) = true;
end

axisChar = struct('sagittal','x', 'coronal','y', 'axial','z');
fprintf('  candidate %s slices for %s (%s):\n', orientation, folderName, threshLab);
fprintf('    %3s  %8s  %8s\n', upper(axisChar.(lower(orientation))), 'voxels', 'peak|t|');
for i = find(show)
    if tbl(i).chosen, mark = '*'; else, mark = ' '; end
    fprintf('   %s%4d  %8d  %8.2f\n', mark, round(tbl(i).mm), tbl(i).n_voxels, tbl(i).max_t);
end
nHidden = numel(tbl) - sum(show);
if nHidden > 0
    fprintf('    ... %d more slices with significant voxels not shown\n', nHidden);
end
end

% =========================================================================
% Helper: look up hand-picked slices for one contrast folder + orientation.
% Returns [] when nothing is configured, so the caller can fall back.
% =========================================================================
function val = spec_lookup(figSpec, topDir, folderName, field, threshKey)
% Resolution order, most specific first:
%   figSpec.<analysisDir>.<contrast>.<unc|corrected>.<field>
%   figSpec.<analysisDir>.<contrast>.<field>
%   figSpec.<contrast>.<unc|corrected>.<field>
%   figSpec.<contrast>.<field>
% The threshold level allows different slices per outlined mask, as the
% corrected and uncorrected effects of a contrast can lie far apart.
val = [];
if ~isstruct(figSpec), return; end
if nargin < 5, threshKey = ''; end

nodes = {};
if isfield(figSpec, topDir) && isstruct(figSpec.(topDir))
    nodes{end+1} = figSpec.(topDir);
end
nodes{end+1} = figSpec;

for iN = 1:numel(nodes)
    if ~isfield(nodes{iN}, folderName), continue; end
    entry = nodes{iN}.(folderName);
    if ~isstruct(entry), continue; end
    if ~isempty(threshKey) && isfield(entry, threshKey) && isstruct(entry.(threshKey)) ...
            && isfield(entry.(threshKey), field)
        val = entry.(threshKey).(field);
        return
    end
    if isfield(entry, field)
        val = entry.(field);
        return
    end
end
end

% =========================================================================
% Helper: union of one or more binary masks -> written to outPath.
% Skips files that are missing or empty. Returns the number of non-empty
% inputs combined; writes outPath only if at least one had voxels.
% =========================================================================
function nUsed = make_union_mask(files, outPath)
nUsed = 0;
acc   = [];
Vref  = [];
for i = 1:numel(files)
    if ~exist(files{i}, 'file'), continue; end
    Vi = spm_vol(files{i});
    Yi = spm_read_vols(Vi) > 0;
    if ~any(Yi(:)), continue; end
    if isempty(acc)
        acc  = Yi;
        Vref = Vi;
    else
        acc = acc | Yi;
    end
    nUsed = nUsed + 1;
end
if nUsed == 0
    % Nothing survived in either direction; make sure no stale union file
    % is left behind to be picked up as significant.
    if exist(outPath, 'file'), delete(outPath); end
    return
end
Vout       = Vref;
Vout.fname = outPath;
Vout.dt    = [spm_type('uint8') 0];
spm_write_vol(Vout, double(acc));
end
