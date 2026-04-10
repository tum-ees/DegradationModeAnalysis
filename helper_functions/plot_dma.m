function plot_dma(data, varargin)
% plot_dma  (tiledlayout edition)                                  2025-06-18
%--------------------------------------------------------------------------
%> Author: Mathias Rehm (mathias.rehm@tum.de)
%> additional code by Can Korkmaz (can.korkmaz@tum.de)
%> Date: 2025-06-18
%>
%> Features:
%>   * Keeps digit-based CU filtering and RMSE plotting.
%>   * Uses a nested tiledlayout so the RMSE panel gets extra horizontal
%>     space for its y-label.
%>
% DESCRIPTION
%   * Plots the main degradation modes (LAM / LLI) versus EFC.
%   * Works with data structs whose field names contain digits (CU0, CU1, ...);
%>     non-digit housekeeping fields are skipped.
%   * Optional right-most panel shows RMSE from rmseFullRange.
%   * Supports aggregate-only and blend-aware DMA plots.
%   * All options are name-value pairs.
%>
% Minimal usage
%   plot_dma(dataStruct,'EFC',0:5:15);
%
% INPUTS:
%   data      - DMA result struct with CU-like fields created by main_dma
%   varargin  - name-value options for labels, EFC values, blend toggles,
%               RMSE display, and title behavior
%
% OUTPUTS:
%   none      - creates a DMA summary figure in the current figure context
%--------------------------------------------------------------------------

%% 1 - Set defaults -------------------------------------------------------
useAnodeBlendModel    = false;
useCathodeBlendModel  = false;
EFC                   = [];
plotTitles            = {};
plotRMSE              = true;   % show RMSE panel unless explicitly disabled
calendarOrCyclic      = 0;
labelCathode          = 'Cathode';
labelAnode            = 'Anode';
labelAnodeBlend1      = 'An-blend1';
labelAnodeBlend2      = 'An-blend2';
labelCathodeBlend1    = 'Ca-blend1';
labelCathodeBlend2    = 'Ca-blend2';
labelChargeCarrierInv = 'Charge-carrier-inv';
plotCathode           = true;   % Hide only for LFP, where cathode aging is not meaningful

%% 2 - Parse varargin -----------------------------------------------------
for v = 1:2:numel(varargin)
    key = lower(string(varargin{v}));
    val = varargin{v+1};

    switch key
        case {'useanodeblendmodel','anodeblend'}
            useAnodeBlendModel = logical(val);
        case {'usecathodeblendmodel','cathodeblend'}
            useCathodeBlendModel = logical(val);

        case 'efc'
            EFC = val;

        case {'plottitles','plot_title','titles','title'}
            plotTitles = val;

        case 'plotrmse'
            plotRMSE = logical(val);

        case 'calendarorcyclic'
            calendarOrCyclic = val;

        case 'labels'
            if isstruct(val)
                if isfield(val, 'labelCathode'),          labelCathode = val.labelCathode; end
                if isfield(val, 'labelAnode'),            labelAnode = val.labelAnode; end
                if isfield(val, 'labelAnodeBlend1'),      labelAnodeBlend1 = val.labelAnodeBlend1; end
                if isfield(val, 'labelAnodeBlend2'),      labelAnodeBlend2 = val.labelAnodeBlend2; end
                if isfield(val, 'labelCathodeBlend1'),    labelCathodeBlend1 = val.labelCathodeBlend1; end
                if isfield(val, 'labelCathodeBlend2'),    labelCathodeBlend2 = val.labelCathodeBlend2; end
                if isfield(val, 'labelChargeCarrierInv'), labelChargeCarrierInv = val.labelChargeCarrierInv; end
            end

        case {'labelcathode','cathodelabel'}
            labelCathode = val;
        case {'labelanode','anodelabel'}
            labelAnode = val;
        case {'labelanodeblend1','anodeblend1label'}
            labelAnodeBlend1 = val;
        case {'labelanodeblend2','anodeblend2label'}
            labelAnodeBlend2 = val;
        case {'labelcathodeblend1','cathodeblend1label'}
            labelCathodeBlend1 = val;
        case {'labelcathodeblend2','cathodeblend2label'}
            labelCathodeBlend2 = val;
        case {'labelchargecarrierinv','labellithium','lithiumlabel'}
            labelChargeCarrierInv = val;

        case 'plotcathode'
            plotCathode = logical(val);

        otherwise
            warning('plot_dma: unknown option %s ignored.', varargin{v});
    end
end

%% 3 - Identify CU-like fields and validate EFC --------------------------
allFields = fieldnames(data);

% Keep only fields whose names contain at least one digit.
isCU  = cellfun(@(s) ~isempty(regexp(s, '\d', 'once')), allFields);
cuFld = allFields(isCU);
nCU   = numel(cuFld);
if nCU == 0
    error('plot_dma: no CU-like fields (containing digits) found in data.');
end

% Sort CU fields by numeric suffix (CU1, CU2, ... CU10).
cuNums = zeros(nCU, 1);
for k = 1:nCU
    tok = regexp(cuFld{k}, '\d+', 'match', 'once');
    cuNums(k) = str2double(tok);
end
[~, sortIdx] = sort(cuNums);
cuFld = cuFld(sortIdx);

% If user omitted EFC, fall back to 1:nCU.
if isempty(EFC)
    EFC = 1:nCU;
end
if numel(EFC) ~= nCU
    warning('plot_dma: length(EFC) = %d, but CU fields = %d.', numel(EFC), nCU);
end

%% 4 - Harvest LAM / LI and optionally RMSE ------------------------------
varNames = {};
lamExtractors = {};

if plotCathode
    varNames{end+1} = labelCathode;
    lamExtractors{end+1} = @(d) 100 * d.lamCathode;
    if useCathodeBlendModel
        varNames = [varNames, {labelCathodeBlend1, labelCathodeBlend2}];
        lamExtractors = [lamExtractors, ...
            {@(d) 100 * d.lamCathodeBlend1, @(d) 100 * d.lamCathodeBlend2}];
    end
end

% Anode aggregate always plotted; add blends if enabled.
varNames{end+1} = labelAnode;
lamExtractors{end+1} = @(d) 100 * d.lamAnode;
if useAnodeBlendModel
    varNames = [varNames, {labelAnodeBlend1, labelAnodeBlend2}];
    lamExtractors = [lamExtractors, ...
        {@(d) 100 * d.lamAnodeBlend1, @(d) 100 * d.lamAnodeBlend2}];
end

% Charge-carrier inventory always plotted.
varNames{end+1} = labelChargeCarrierInv;
lamExtractors{end+1} = @(d) 100 * d.li;

nPlots = numel(varNames);

LAM = zeros(nCU, nPlots);   % rows = CU, cols = curves
RMSE = NaN(nCU, 1);         % stays NaN if plotRMSE == false

for k = 1:nCU
    f = cuFld{k};

    for idxLam = 1:nPlots
        LAM(k, idxLam) = lamExtractors{idxLam}(data.(f));
    end

    if plotRMSE && isfield(data.(f), 'rmseFullRange')
        RMSE(k) = 1e3 * data.(f).rmseFullRange;
    end
end

%% 5 - Figure and nested tiledlayout to create RMSE gap ------------------
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

fig = figure('Units', 'centimeters', 'Position', [3 3 20 6]);

outerCols = nPlots + double(plotRMSE);   % +1 column only if RMSE requested
outerTL   = tiledlayout(fig, 1, outerCols, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

% Place an inner tiledlayout spanning the first nPlots tiles of outerTL.
innerTL = tiledlayout(outerTL, 1, nPlots, ...
    'TileSpacing', 'none', 'Padding', 'compact');
innerTL.Layout.Tile     = 1;
innerTL.Layout.TileSpan = [1 nPlots];

% Result: LAM axes have zero internal gap, while outerTL's inter-tile
% spacing provides the extra gap before the RMSE axis.

%% 6 - Style constants ----------------------------------------------------
colors = tum_colors();
markerSet = {'-o', '-s', '-d', '-^', '-v', '-*'};
while numel(markerSet) < outerCols
    markerSet = [markerSet markerSet]; %#ok<AGROW>
end

fontSz = 10;
lnW = 2;
mkSz = 5;
yMin = min(LAM, [], 'all');
yMax = max(LAM, [], 'all');

%% 7 - Plot LAM / LLI panels ---------------------------------------------
for p = 1:nPlots
    ax = nexttile(innerTL, p);
    set(ax, 'ColorOrder', colors.colorOrder);
    hold(ax, 'on');
    plot(ax, EFC, LAM(:, p), markerSet{p}, ...
        'Color', colors.tumBlue, 'LineWidth', lnW, 'MarkerSize', mkSz);

    grid(ax, 'on');
    xlim(ax, 'padded');
    ylim(ax, [yMin yMax]);
    box(ax, 'on');

    if isempty(plotTitles)
        titleText = varNames{p};
    else
        titleText = plotTitles{p};
    end
    title(ax, titleText, 'Interpreter', 'latex', 'FontSize', 1.2 * fontSz);

    if p == 1
        ylabel(ax, 'Capacity Loss / \%', 'Interpreter', 'latex', 'FontSize', fontSz);
    else
        ax.YTickLabel = [];
    end

    if calendarOrCyclic == 0
        xlabel(ax, 'RPT Number / -', 'Interpreter', 'latex', 'FontSize', fontSz);
    else
        xlabel(ax, 'EFC', 'Interpreter', 'latex', 'FontSize', fontSz);
    end

    ax.TickLength = [0 0];
    ax.TickDir    = 'out';
    ax.LineWidth  = lnW;
    ax.FontSize   = fontSz;
    ax.FontName   = 'Times New Roman';
end

%% 8 - RMSE panel ---------------------------------------------------------
if plotRMSE
    axR = nexttile(outerTL, outerCols);
    set(axR, 'ColorOrder', colors.colorOrder);
    hold(axR, 'on');

    if all(isnan(RMSE))
        warning('plot_dma: all RMSE values are NaN. RMSE panel will be empty.');
    end

    plot(axR, EFC, RMSE, '-o', 'Color', colors.tumBlue, ...
        'LineWidth', lnW, 'MarkerSize', mkSz);

    grid(axR, 'on');
    xlim(axR, 'padded');
    ylim(axR, 'padded');
    ylim(axR, [0, max(axR.YLim)]);
    box(axR, 'on');

    title(axR, 'RMSE', 'Interpreter', 'latex', 'FontSize', 1.2 * fontSz);
    ylabel(axR, 'RMSE / mV', 'Interpreter', 'latex', 'FontSize', fontSz);
    if calendarOrCyclic == 0
        xlabel(axR, 'RPT Number / -', 'Interpreter', 'latex', 'FontSize', fontSz);
    else
        xlabel(axR, 'EFC', 'Interpreter', 'latex', 'FontSize', fontSz);
    end

    axR.TickLength = [0 0];
    axR.TickDir    = 'out';
    axR.LineWidth  = lnW;
    axR.FontSize   = fontSz;
    axR.FontName   = 'Times New Roman';
end

end
