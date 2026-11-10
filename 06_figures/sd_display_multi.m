function [h_fig, info] = sd_display_multi(layers, blocks, opts)
% SD_DISPLAY_MULTI Dual-coded slice display with several orientations in one figure.
%
% DESCRIPTION
% Alternative to sd_display (slice_display toolbox) with three additions:
%
%   1. Several orientations in ONE figure. Each entry of `blocks` becomes a
%      row-block of slices (e.g. a row of axial slices above a row of
%      coronal slices), sharing a single dual-coded colour bar.
%   2. The layout is computed in millimetres and the figure is printed at
%      exactly that physical size, so output does not depend on the screen
%      (text stays in proportion to the slices however many there are).
%   3. The long descriptive string goes in a title across the full figure
%      width; title and colour-bar label shrink automatically to fit.
%
% The per-slice rendering itself is deliberately identical to sd_display:
% same sd_get_image_specs / sd_get_slice / sd_slice_to_rgb /
% sd_slice_to_alpha calls, same layer semantics, same slice labels.
%
% SYNTAX
%   [h_fig, info] = sd_display_multi(layers, blocks, opts)
%
% layers  - Nx1 struct, as for sd_display (see sd_config_layers)
% blocks  - 1xM struct array, one element per orientation block:
%             .orientation  'axial' | 'coronal' | 'sagittal'
%             .slices       vector of slice coordinates in mm (world space)
% opts    - struct, all fields optional:
%             .n_columns          max slices per row          (default 4)
%             .tile_width_mm      width of one slice tile     (default 42)
%             .margin_mm          [left bottom right top]     (default [7 5 7 4])
%             .block_gap_mm       vertical gap between blocks (default 4)
%             .same_row           true: put every plane's slices on ONE row,
%                                 tiles sharing a height and taking their
%                                 width from each plane's aspect ratio.
%                                 false (default): one row-block per plane.
%             .colorbar_width_mm                              (default 55)
%             .colorbar_height_mm                             (default 7)
%             .title              char or cellstr, drawn at the top
%             .title_font_pt                                  (default 9)
%             .label_font_pt      colour bar labels/ticks     (default 8)
%             .show_labels        slice coordinate labels     (default true)
%             .label_style        'axis' ("z = 10 mm") | 'plain' ("10 mm")
%                                                             (default 'axis')
%             .label_color        slice label colour          (default 'w')
%             .label_box          draw the slice label on an opaque box, so
%                                 it stays readable even where it falls on
%                                 bright tissue                (default false)
%             .show_orientation   L/R (or P/A) markers        (default true)
%             .out_file           full path (no extension) -> printed to PNG
%             .dpi                print resolution for PNG    (default 300)
%             .formats            cellstr of {'png','pdf'} (default {'png'}).
%                                 'pdf' prints with the painters renderer,
%                                 giving editable text and vector contours
%                                 with the slices as embedded rasters.
%             .visible            'on' | 'off'                (default 'off')
%             .renderer           'painters' (default) | 'opengl'. Painters
%                                 needs no GL context and renders
%                                 identically here.
%
% info    - struct with .width_mm, .height_mm, .out_file and .blocks (the
%           slice coordinates actually rendered per block, after snapping to
%           the sampling grid).
%
% EXAMPLE
%   blocks(1) = struct('orientation','axial',  'slices',[-6 4 14 24]);
%   blocks(2) = struct('orientation','coronal','slices',[-60 -30 0 12]);
%   opts = struct('title','MPH vs PBO, Target (TFCE FWE)', ...
%                 'out_file','/tmp/fig');
%   sd_display_multi(layers, blocks, opts);
%
% .........................................................................
% Builds on the slice_display toolbox by Bram Zandbelt
% (https://github.com/bramzandbelt/slice_display), which must be on the path.

if nargin < 3, opts = struct(); end

% ---- Option defaults ----------------------------------------------------
opts = set_default(opts, 'n_columns',          4);
opts = set_default(opts, 'tile_width_mm',      42);
opts = set_default(opts, 'margin_mm',          [7 5 7 4]);   % L B R T
opts = set_default(opts, 'block_gap_mm',       4);
opts = set_default(opts, 'same_row',           false);
opts = set_default(opts, 'colorbar_width_mm',  55);
opts = set_default(opts, 'colorbar_height_mm', 7);
opts = set_default(opts, 'title',              '');
opts = set_default(opts, 'title_font_pt',      9);
opts = set_default(opts, 'label_font_pt',      8);
opts = set_default(opts, 'show_labels',        true);
opts = set_default(opts, 'label_style',        'axis');
opts = set_default(opts, 'label_color',        'w');
opts = set_default(opts, 'label_box',          false);
opts = set_default(opts, 'show_orientation',   true);
opts = set_default(opts, 'out_file',           '');
opts = set_default(opts, 'dpi',                300);
opts = set_default(opts, 'formats',            {'png'});
opts = set_default(opts, 'visible',            'off');
opts = set_default(opts, 'renderer',           'painters');

% ---- Fill layer defaults once (headers, colour ranges, ...) -------------
layers = sd_config_layers('fill-defaults', layers);

% ---- Per-block image specs ----------------------------------------------
% sd_get_image_specs writes slice bookkeeping into settings, so each block
% needs a fresh settings struct (a reused one gives wrong slice indices).
nBlock = numel(blocks);
S      = cell(1, nBlock);
keep   = true(1, nBlock);
for b = 1:nBlock
    wanted = snap_slices(layers, blocks(b).orientation, blocks(b).slices);
    if isempty(wanted)
        warning('sd_display_multi:emptyBlock', ...
            'No displayable slices for the %s block; dropping it.', blocks(b).orientation);
        keep(b) = false;
        continue
    end
    s = sd_config_settings('init');
    s.slice.orientation = blocks(b).orientation;
    s.slice.disp_slices = wanted;
    s = sd_config_settings('fill-defaults', s);
    s = sd_get_image_specs(layers, s);
    S{b} = s;
end
S      = S(keep);
blocks = blocks(keep);
nBlock = numel(blocks);
if nBlock == 0
    error('sd_display_multi:noSlices','None of the requested blocks had displayable slices.');
end

% ---- Layout, in millimetres ---------------------------------------------
m     = opts.margin_mm;
tileW = opts.tile_width_mm;

nSlice = cellfun(@(s) numel(s.slice.zmm), S);
ratio  = cellfun(@(s) s.fig_specs.height_width_ratio.slice, S);   % height/width

if opts.same_row
    % One row for everything. Tiles share a height so they line up; each
    % plane's width follows its own aspect ratio (a sagittal slice is wider
    % than a coronal one), so nothing is squashed.
    tileHb = repmat(tileW, 1, nBlock);
    tileWb = tileHb ./ ratio;
    % block_gap_mm separates the planes horizontally here, so the L/R marker
    % of one plane does not butt against the next plane's.
    mapW   = sum(tileWb .* nSlice) + (nBlock-1)*opts.block_gap_mm;
    mapH   = max(tileHb);
    nCols  = nSlice;
    nRows  = ones(1, nBlock);
else
    nColsF = max(1, min(opts.n_columns, max(nSlice)));   % figure-wide column count
    mapW   = nColsF * tileW;
    nCols  = min(nColsF, nSlice);
    nRows  = ceil(nSlice ./ nCols);
    tileWb = repmat(tileW, 1, nBlock);
    tileHb = tileW * ratio;
    mapH   = sum(nRows .* tileHb) + (nBlock-1)*opts.block_gap_mm;
end
figW = mapW + m(1) + m(3);

% Legend strip: colour bar + its tick labels + its xlabel.
hasBar = any_dual_with_label(layers);
if hasBar
    legendH = opts.colorbar_height_mm + 13;
else
    legendH = 0;
end

% Title strip
titleLines = cellstr(opts.title);
titleLines = titleLines(~cellfun(@isempty, titleLines));
if isempty(titleLines)
    titleH = 0;
else
    titleH = numel(titleLines) * (opts.title_font_pt * 0.42) + 3;   % pt -> mm, plus padding
end

figH = m(2) + legendH + mapH + titleH + m(4);

% ---- Figure --------------------------------------------------------------
h_fig = figure('Visible', opts.visible, ...
               'Color', 'w', ...
               'Units', 'centimeters', ...
               'Position', [1 1 figW/10 figH/10], ...
               'PaperUnits', 'centimeters', ...
               'PaperSize', [figW/10 figH/10], ...
               'PaperPositionMode', 'manual', ...
               'PaperPosition', [0 0 figW/10 figH/10], ...
               'InvertHardcopy', 'off');
% InvertHardcopy MUST stay 'off': with 'on', print turns the white slice
% labels and L/R markers black, i.e. invisible on the black background.

% Nothing in the figure uses transparency (image layers are composited
% below), so the painters renderer gives the same result as OpenGL and
% works on headless machines without a GL context.
set(h_fig, 'Renderer', opts.renderer, 'RendererMode', 'manual');

% ---- Draw the slice blocks, top to bottom --------------------------------
yTop = figH - m(4) - titleH;                      % mm above the figure bottom
xCur = m(1);                                      % same_row: left-to-right cursor
for b = 1:nBlock

    s        = S{b};
    zmm      = s.slice.zmm;
    x        = s.slice.x;
    y        = s.slice.y;
    vdims    = s.slice.vdims;
    nvox     = prod(vdims(1:2));
    axisChar = orientation_axis_char(s.slice.orientation);

    for i_slice = 1:numel(zmm)

        if opts.same_row
            pos_mm = [xCur, yTop - tileHb(b), tileWb(b), tileHb(b)];
            xCur   = xCur + tileWb(b);
            if i_slice == numel(zmm), xCur = xCur + opts.block_gap_mm; end
        else
            i_row = ceil(i_slice / nCols(b));
            i_col = i_slice - (i_row-1)*nCols(b);

            % A single-row block is centred; multi-row blocks stay
            % left-aligned so the reading order across rows is unambiguous.
            if nRows(b) == 1
                rowLeft = m(1) + (mapW - numel(zmm)*tileW)/2;
            else
                rowLeft = m(1);
            end

            pos_mm = [rowLeft + (i_col-1)*tileW, ...
                      yTop - i_row*tileHb(b), ...
                      tileWb(b), tileHb(b)];
        end
        ax = axes('Parent', h_fig, 'Units', 'normalized', ...
                  'Position', pos_mm ./ [figW figH figW figH]); %#ok<LAXES>
        hold(ax, 'on');

        xyzmm = [x(:)'; y(:)'; ones(1,nvox)*zmm(i_slice); ones(1,nvox)];

        % Image layers are alpha-blended here into one opaque RGB image
        % instead of stacked with per-pixel AlphaData. Without transparency
        % at render time, painters can write a true vector PDF (editable
        % text, vector contours, slices as embedded rasters).
        canvas = [];
        cJobs  = {};
        for i_layer = 1:numel(layers)
            [Y_rgb, Y_alpha, Y_c] = layer_slice(layers(i_layer), xyzmm, s);
            if strcmpi(layers(i_layer).type, 'contour')
                cJobs{end+1} = {Y_c, layers(i_layer)}; %#ok<AGROW>
            else
                if isempty(canvas)
                    canvas = ones(size(Y_rgb));   % figure background
                end
                a3     = repmat(Y_alpha, [1 1 3]);
                canvas = canvas .* (1 - a3) + Y_rgb .* a3;
            end
        end
        if ~isempty(canvas)
            image(canvas, 'Parent', ax);
        end
        for iJ = 1:numel(cJobs)
            draw_contour(ax, cJobs{iJ}{1}, cJobs{iJ}{2});
        end

        % Fill the axes box exactly: the box was sized to this slice's own
        % height/width ratio, and world-space voxels here are isotropic.
        %
        % YDir MUST be 'normal', not the 'reverse' that image() sets: row 1
        % of the slice is the most posterior (axial) or inferior (coronal)
        % row and belongs at the bottom; 'reverse' draws brains upside down.
        set(ax, 'YDir','normal', ...
                'DataAspectRatio',[1 1 1], ...
                'XLim',[0.5 vdims(1)+0.5], ...
                'YLim',[0.5 vdims(2)+0.5], ...
                'Visible','off');

        % Slice coordinate label
        if opts.show_labels
            switch lower(opts.label_style)
                case 'axis', txt = sprintf('%s = %d mm', axisChar, round(zmm(i_slice)));
                otherwise,   txt = sprintf('%d mm', round(zmm(i_slice)));
            end
            h_lab = text(.025, .025, txt, 'Parent', ax, ...
                'Color', opts.label_color, 'HorizontalAlignment','left', ...
                'VerticalAlignment','bottom', 'Units','normalized', ...
                'FontUnits','normalized', 'FontSize',0.075, 'Interpreter','none');
            if opts.label_box
                set(h_lab, 'BackgroundColor','k', 'Margin',1);
            end
        end

        % Left/right (or posterior/anterior) markers on the last tile of the block
        if opts.show_orientation && i_slice == numel(zmm)
            switch lower(s.slice.orientation)
                case {'axial','coronal'}, marks = {'L','R'};
                case 'sagittal',          marks = {'P','A'};
                otherwise,                marks = {'',''};
            end
            text(.025, .975, marks{1}, 'Parent', ax, 'Color', opts.label_color, ...
                'HorizontalAlignment','left', 'VerticalAlignment','top', ...
                'Units','normalized', 'FontUnits','normalized', 'FontSize',0.1);
            text(.975, .975, marks{2}, 'Parent', ax, 'Color', opts.label_color, ...
                'HorizontalAlignment','right', 'VerticalAlignment','top', ...
                'Units','normalized', 'FontUnits','normalized', 'FontSize',0.1);
        end
    end

    if ~opts.same_row
        yTop = yTop - nRows(b)*tileHb(b) - opts.block_gap_mm;
    end
end

% ---- Colour bar ----------------------------------------------------------
if hasBar
    draw_dual_colorbar(h_fig, layers, opts, figW, figH, m);
end

% ---- Title ---------------------------------------------------------------
if ~isempty(titleLines)
    draw_title(h_fig, titleLines, opts, figW, figH, titleH);
end

% ---- Print at the exact physical size ------------------------------------
info = struct('width_mm', figW, 'height_mm', figH, 'out_file', '');
info.blocks = struct('orientation', {blocks.orientation}, ...
                     'slices',      cellfun(@(s) s.slice.zmm, S, 'UniformOutput', false));
if ~isempty(opts.out_file)
    [outDir, outName] = fileparts(opts.out_file);
    if ~isempty(outDir) && ~exist(outDir, 'dir'), mkdir(outDir); end
    fmts = opts.formats;
    if ischar(fmts), fmts = {fmts}; end
    written = {};
    for iF = 1:numel(fmts)
        switch lower(fmts{iF})
            case 'png'
                f = fullfile(outDir, [outName '.png']);
                print(h_fig, '-dpng', sprintf('-r%d', opts.dpi), f);
            case 'pdf'
                % -painters is what makes it vector; safe because nothing in
                % the figure uses transparency.
                f = fullfile(outDir, [outName '.pdf']);
                print(h_fig, '-dpdf', '-painters', f);
            otherwise
                warning('sd_display_multi:badFormat', ...
                    'Unknown output format "%s", skipping.', fmts{iF});
                continue
        end
        written{end+1} = f; %#ok<AGROW>
    end
    if ~isempty(written), info.out_file = written{1}; end
    info.out_files = written;
end

end % sd_display_multi

% =========================================================================
% Subfunctions
% =========================================================================

function [Y_rgb, Y_alpha, Y_c] = layer_slice(layer, xyzmm, s)
% Colour and opacity for one layer of one slice. Same computation as
% sd_display's layer loop, but returned rather than drawn, so the caller can
% composite the image layers into a single opaque image.

transform = s.slice.transform;
vdims     = s.slice.vdims;

Y_c     = sd_get_slice(layer.color.header, xyzmm, transform, vdims, layer.color.hold);
Y_rgb   = [];
Y_alpha = [];

if strcmpi(layer.type, 'contour')
    return
end

Y_rgb = sd_slice_to_rgb(Y_c, layer, s);

switch lower(layer.type)
    case {'truecolor','blob','cluster'}
        Y_alpha = ones(size(Y_c)) .* layer.color.opacity;
    case 'dual'
        Y_o     = sd_get_slice(layer.opacity.header, xyzmm, transform, vdims, layer.opacity.hold);
        Y_alpha = sd_slice_to_alpha(Y_o, layer);
end

% Masking
if ~isempty(layer.mask.file)
    Y_m = sd_get_slice(layer.mask.header, xyzmm, transform, vdims, layer.color.hold);
else
    switch lower(layer.type)
        case 'cluster'
            Y_m = ones(size(Y_c));
            Y_m(Y_c == 0) = 0;
        otherwise
            Y_m = ones(size(Y_c));
    end
end
Y_alpha = Y_alpha .* Y_m;

% Voxels with no data must not be painted. pr_scaletocmap maps NaN to
% colormap index 1 -- the extreme negative colour -- so without this a NaN
% renders as a saturated cyan pixel. It only becomes visible with trilinear
% resampling (hold = 1), which spreads NaN from outside the analysis mask
% into its neighbours and draws a cyan fringe along every blob edge.
Y_alpha(~isfinite(Y_c))     = 0;
Y_alpha(~isfinite(Y_alpha)) = 0;
end

% -------------------------------------------------------------------------
function draw_contour(ax, Y_c, layer)
h_contour = [];
if strncmp(spm_type(layer.color.header.dt), 'int', 3)
    [~, h_contour] = contour(ax, Y_c, 1);
elseif any(Y_c(:) > 0) && any(Y_c(:) <= 0)
    % Binary masks are written as uint8 (imcalc dtype 2 / spm_type('uint8')),
    % which the test above does not match. Contour them at the mask edge
    % explicitly rather than leaving h_contour undefined.
    [~, h_contour] = contour(ax, Y_c, [0.5 0.5]);
end
if ~isempty(h_contour)
    set(h_contour, 'LineColor', layer.color.map, ...
                   'LineStyle', layer.color.line_style, ...
                   'LineWidth', layer.color.line_width);
end
end

% -------------------------------------------------------------------------
function tf = any_dual_with_label(layers)
tf = false;
for i = 1:numel(layers)
    if strcmpi(layers(i).type, 'dual') && ~isempty(layers(i).color.label)
        tf = true; return
    end
end
end

% -------------------------------------------------------------------------
function draw_dual_colorbar(h_fig, layers, opts, figW, figH, m)
% Single dual-coded colour bar, centred in the legend strip at the bottom.
% Same content as sd_display's 'dual' colour bar, but in an axes whose size
% and position are set explicitly so nothing lands outside the figure.

i_layer = find(arrayfun(@(L) strcmpi(L.type,'dual'), layers), 1);
layer   = layers(i_layer);

cbW = min(opts.colorbar_width_mm, figW - m(1) - m(3) - 14);   % 14 mm for the ylabel
cbH = opts.colorbar_height_mm;
cbL = (figW - cbW)/2;
cbB = m(2) + 10;                                              % 10 mm for ticks + xlabel

ax = axes('Parent', h_fig, 'Units','normalized', ...
          'Position', [cbL cbB cbW cbH] ./ [figW figH figW figH]);
color_vector = linspace(layer.color.range(1), layer.color.range(2), size(layer.color.map,1));
alpha_vector = linspace(layer.opacity.range(1), layer.opacity.range(2), 256);

% Composited to opaque RGB for the same reason the slices are: no
% AlphaData means painters can render this as vector. Row i of the bar is
% the colormap blended onto white at alpha (i-1)/255, which is exactly what
% alpha('scaled') produced against the white background before.
nA  = numel(alpha_vector);
a   = linspace(0, 1, nA)';
cmap = layer.color.map;
barRGB = zeros(nA, size(cmap,1), 3);
for ch = 1:3
    barRGB(:,:,ch) = (1 - a) * ones(1, size(cmap,1)) + a * cmap(:,ch)';
end

image(color_vector, alpha_vector, barRGB, 'Parent', ax);
axis(ax, 'tight');

% Colours set explicitly: MATLAB R2025a+ figure themes would otherwise
% render the labels light grey on white.
set(ax, 'Box','off', ...
        'Color','w', 'XColor','k', 'YColor','k', ...
        'FontSize', opts.label_font_pt, ...
        'XLim',  layer.color.range, ...
        'YLim',  layer.opacity.range, ...
        'YTick', layer.opacity.range, ...
        'YTickLabel', {sprintf('%g', layer.opacity.range(1)), ...
                       sprintf('>%.1f', layer.opacity.range(2))}, ...
        'YDir', 'normal');   % opaque (high |t|) at the top, as in sd_display

h_xlabel = xlabel(ax, layer.color.label,   'FontSize', opts.label_font_pt, ...
                  'Interpreter','tex', 'Color','k');
           ylabel(ax, layer.opacity.label, 'FontSize', opts.label_font_pt, ...
                  'Interpreter','tex', 'Color','k');

% Shrink the xlabel if it is wider than the figure
shrink_text_to_width(h_xlabel, ax, figW, figW - 4, opts.label_font_pt);
end

% -------------------------------------------------------------------------
function draw_title(h_fig, titleLines, opts, figW, figH, titleH)
m  = opts.margin_mm;
ax = axes('Parent', h_fig, 'Units','normalized', ...
          'Position', [m(1) (figH - m(4) - titleH) (figW - m(1) - m(3)) titleH] ...
                      ./ [figW figH figW figH]);
set(ax, 'XLim',[0 1], 'YLim',[0 1], 'Visible','off');

h = text(0.5, 0.5, titleLines, 'Parent', ax, ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
    'FontSize', opts.title_font_pt, 'Interpreter','none', 'Color','k');

% The axes spans the full text width and 0..1 in x, so a normalised extent
% width above 1 means the text overflows. Shrink until it fits.
fs = opts.title_font_pt;
while fs > 5
    e = get(h, 'Extent');
    if e(3) <= 1, break; end
    fs = fs - 0.5;
    set(h, 'FontSize', fs);
end
end

% -------------------------------------------------------------------------
function shrink_text_to_width(h_text, ax, figW, maxW_mm, fs0)
% Reduce a label's font size until its rendered width fits maxW_mm.
if isempty(h_text) || ~ishandle(h_text), return; end
axPos  = get(ax, 'Position');            % normalised figure units
axW_mm = axPos(3) * figW;
xr     = diff(get(ax, 'XLim'));
if ~(xr > 0), return; end                % degenerate colour range
fs     = fs0;
while fs > 5
    e = get(h_text, 'Extent');           % data units
    if (e(3)/xr) * axW_mm <= maxW_mm, break; end
    fs = fs - 0.5;
    set(h_text, 'FontSize', fs);
end
end

% -------------------------------------------------------------------------
function c = orientation_axis_char(orientation)
switch lower(orientation)
    case 'axial',    c = 'z';
    case 'coronal',  c = 'y';
    case 'sagittal', c = 'x';
    otherwise,       c = '?';
end
end

% -------------------------------------------------------------------------
function opts = set_default(opts, field, value)
if ~isfield(opts, field) || isempty(opts.(field))
    opts.(field) = value;
end
end

% -------------------------------------------------------------------------
function [snapped, grid_mm] = snap_slices(layers, orientation, wanted)
% Snap requested slice coordinates onto the grid sd_get_image_specs samples.
%
% sd_get_image_specs keeps only the requested slices that are EXACTLY on
% that grid (`slices(ismember(slices, disp_slices))`), so a coordinate off
% by a fraction of a millimetre is silently dropped. Snapping first, using
% the same grid construction, means what you ask for is what you get.
%
% layers      - layer struct array (only layers(1).color is used, as in
%               sd_get_image_specs) or a single spm_vol header struct
% orientation - 'axial' | 'coronal' | 'sagittal'
% wanted      - vector of slice coordinates in mm (world space)
%
% snapped     - unique, ascending, on-grid slice coordinates
% grid_mm     - the full sampling grid for this orientation
%
% Grid construction copied from sd_get_image_specs (Bram Zandbelt).

% ---- Reference header ----------------------------------------------------
if isstruct(layers) && isfield(layers, 'color')
    V = layers(1).color.header;
    if isempty(V)
        V = spm_vol(layers(1).color.file);
    end
elseif isstruct(layers) && isfield(layers, 'mat')
    V = layers(1);
else
    error('snap_slices:badInput','Pass a layers struct or an spm_vol header.');
end

% ---- Same transform / grid as sd_get_image_specs -------------------------
orientation_ix = strcmpi(orientation, {'axial', 'coronal', 'sagittal'});
if ~any(orientation_ix)
    error('snap_slices:badOrientation','Unknown orientation "%s".', orientation);
end
ts = [0 0 0 0     0 0     1  1 1; ...
      0 0 0 pi/2  0 0     1 -1 1; ...
      0 0 0 pi/2  0 -pi/2 -1 1 1];
transform = spm_matrix(ts(orientation_ix,:));

D = V.dim(1:3);
T = transform * V.mat;

voxel_corners = [1 1 1; D(1) 1 1; 1 D(2) 1; D(1:2) 1; ...
                 1 1 D(3); D(1) 1 D(3); 1 D(2:3); D(1:3)]';
corners        = T * [voxel_corners; ones(1,8)];
sorted_corners = sort(corners, 2);
voxel_size     = sqrt(sum(T(1:3,1:3).^2));

grid_mm = sorted_corners(3,1) : voxel_size(3) : sorted_corners(3,8);

% ---- Snap ----------------------------------------------------------------
wanted = wanted(:)';
wanted = wanted(isfinite(wanted));
if isempty(wanted) || isempty(grid_mm)
    snapped = [];
    return
end

[~, ix]  = min(abs(bsxfun(@minus, grid_mm(:), wanted)), [], 1);
snapped  = unique(grid_mm(ix));

% Drop anything that had to travel more than one slice to land on the grid
% (i.e. was outside the field of view altogether).
tooFar = arrayfun(@(g) min(abs(wanted - g)) > voxel_size(3), snapped);
snapped(tooFar) = [];
end
