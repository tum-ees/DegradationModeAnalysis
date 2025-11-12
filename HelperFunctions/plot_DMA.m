function plot_DMA(Data, varargin)
%==========================================================================
% plotDMA  (tiledlayout edition – **with RMSE gap**)              2025-06-18
%--------------------------------------------------------------------------
%> Author(s) :  Can Korkmaz  &  Mathias Rehm
%> Last edit :  2025-06-18  – MR
%>             • kept digit-based CU filter & robust RMSE detection
%>             • switched to nested *tiledlayout* so the RMSE panel gets
%>               extra horizontal space for its y-label (replicates the old
%>               manual “gapRMSE = 0.07” behaviour)
%>             • figure export code permanently removed (keywords ignored)
%
% DESCRIPTION
%   • Plots the main degradation modes (LAM / LLI) versus *EFC*.
%   • Works with Data structs whose field names contain digits (“CU0” …);
%     non-digit housekeeping fields are skipped.
%   • Optional right-most panel shows RMSE; recognises
%       -  RMSE_0_9   -  RMSE_full   -  RMSE      (all V → mV).
%   • 3-line mode  (NCA / Graphite / Lithium)  **or**
%     4-line mode  (NCA / Graphite / Silicon / Lithium) via the
%     *useBlendElectrodeModel* flag.
%   • All options are key–value pairs; list at call site is order-free.
%
% Minimal usage
%   plotDMA(DataStruct,'EFC',0:5:15);
%
% Accepted key–value options
%   'useBlendElectrodeModel'   logical   (default false)
%   'EFC'                      numeric   (default [])
%   'plotTitles'               cellstr   (auto if empty)
%   'plotRMSE'                 logical   (default true)
%   ***  “savePathDMA” & “saveFig” are parsed only for backward-compat;   ***
%   ***  they DO NOTHING in this edition (no file export).               ***
%==========================================================================

%% 1 — Set defaults --------------------------------------------------------
useBlendElectrodeModel = false;
EFC         = [];
plotTitles  = {};
plotRMSE    = true;       % show RMSE panel unless explicitly disabled
CalendarOrCyclic = 0;

%% 2 — Parse varargin ------------------------------------------------------
for v = 1:2:numel(varargin)
    key = lower(string(varargin{v}));
    val = varargin{v+1};

    switch key
        case {'useblendelectrodemodel','blend','useblend'}
            useBlendElectrodeModel = logical(val);

        case 'efc'
            EFC = val;

        case {'plottitles','plot_title','titles','title'}
            plotTitles = val;

        case 'plotrmse'
            plotRMSE = logical(val);
        
        case 'calendarorcyclic'
            CalendarOrCyclic = val;

        % ---------- legacy / ignored keys --------------------------------
        case {'savepathdma','savefig','savefigure','savepath'}
            % These keys are silently accepted so old scripts don’t error,
            % but they no longer trigger any file export.

        otherwise
            warning('plotDMA: unknown option “%s” – ignored.', varargin{v});
    end
end

%% 3 — Identify CU-like fields & validate EFC -----------------------------
allFields = fieldnames(Data);

% keep only fields whose names contain at least one digit
isCU  = cellfun(@(s) ~isempty(regexp(s,'\d','once')), allFields);
cuFld = allFields(isCU);
nCU   = numel(cuFld);
if nCU == 0
    error('plotDMA: no CU-like fields (containing digits) found in Data.');
end

% if user omitted EFC, fall back to 1:nCU
if isempty(EFC),  EFC = 1:nCU; end
if numel(EFC) ~= nCU
    warning('plotDMA: length(EFC) = %d, but CU fields = %d.', numel(EFC), nCU);
end

%% 4 — Harvest LAM / LI & (optionally) RMSE ------------------------------
if useBlendElectrodeModel
    varNames = {'NCA','Graphite','Silicon','Lithium', 'Anode'};   % 4 curves
else
    varNames = {'NCA','Graphite','Lithium'};             % 3 curves
end
nPlots = numel(varNames);

LAM = zeros(nCU,nPlots);          % rows = CU, cols = curves
RMSE = NaN(nCU,1);                % will stay NaN if plotRMSE == false

for k = 1:nCU
    f = cuFld{k};
    % ---- LAM / LI ------------------------------------------------------
    LAM(k,1) = 100*Data.(f).LAM_Cathode;
    LAM(k,2) = 100*Data.(f).LAM_Anode_Blend1;
    if useBlendElectrodeModel
        LAM(k,3) = 100*Data.(f).LAM_Anode_Blend2;
        LAM(k,4) = 100*Data.(f).LI;
        LAM(k,5) = 100*Data.(f).LAM_Anode;
    else
        LAM(k,3) = 100*Data.(f).LI;
    end
    % ---- RMSE -----------------------------------------------------------
    if plotRMSE
        if     isfield(Data.(f),'RMSE_0_9')
            RMSE(k) = 1e3*Data.(f).RMSE_0_9;      % V → mV
        elseif isfield(Data.(f),'RMSE_full')
            RMSE(k) = 1e3*Data.(f).RMSE_full;
        elseif isfield(Data.(f),'RMSE')
            RMSE(k) = 1e3*Data.(f).RMSE;
        end
        % else leave as NaN (warning later if all missing)
    end
end

%% 5 — Figure & **nested tiledlayout** to create RMSE gap -----------------
set(groot,'defaultAxesTickLabelInterpreter','latex');

fig = figure('Units','centimeters','Position',[3 3 20 6]);

outerCols = nPlots + double(plotRMSE);   % +1 column only if RMSE requested
outerTL   = tiledlayout(fig,1,outerCols, ...
              'TileSpacing','compact','Padding','compact');


% place an inner tiledlayout spanning the *first* nPlots tiles of outerTL
innerTL = tiledlayout(outerTL,1,nPlots, ...
              'TileSpacing','none','Padding','compact');
innerTL.Layout.Tile     = 1;
innerTL.Layout.TileSpan = [1 nPlots];    % occupies first nPlots columns

% Result: LAM axes have zero internal gap, while outerTL’s inter-tile
% spacing provides the *extra* gap before the RMSE axis.

%% 6 - Style constants ----------------------------------------------------
TUM_colors;                        % defines TUM color palette variables
markerSet = {'-o','-s','-d','-^','-v','-*'};
while numel(markerSet) < outerCols
    markerSet = [markerSet markerSet];  %#ok<AGROW>
end

fontSz = 10;  lnW = 2;  mkSz = 5;
yMin = min(LAM,[],'all');  yMax = max(LAM,[],'all');

%% 7 — Plot LAM / LLI panels ----------------------------------------------
for p = 1:nPlots
    ax = nexttile(innerTL,p);  hold(ax,'on');
    plot(ax,EFC,LAM(:,p),markerSet{p}, ...
         'Color',tumBlue,'LineWidth',lnW,'MarkerSize',mkSz);

    grid(ax,'on');  xlim(ax,'padded');  ylim(ax,[yMin yMax]);

    box on;

    % --- titles ----------------------------------------------------------
    if isempty(plotTitles)
        ttxt = varNames{p};
    else
        ttxt = plotTitles{p};
    end
    title(ax,ttxt,'Interpreter','latex','FontSize',1.2*fontSz);

    % --- y-label only on first panel ------------------------------------
    if p == 1
        ylabel(ax,'Capacity Loss / \%','Interpreter','latex','FontSize',fontSz);
    else
        ax.YTickLabel = [];
    end
    if CalendarOrCyclic == 0    % Calendar Aging
        xlabel('RPT Number / -', 'Interpreter','latex', 'FontSize', fontSz);
    else                         % Cyclic Aging
        xlabel('EFC', 'Interpreter','latex', 'FontSize', fontSz);
    end

    % --- axes cosmetics --------------------------------------------------
    ax.TickLength = [0 0];
    ax.TickDir    = 'out';
    ax.LineWidth  = lnW;
    ax.FontSize   = fontSz;
    ax.FontName   = 'Times New Roman';
end

%% 8 — RMSE panel ----------------------------------------------------------
if plotRMSE
    axR = nexttile(outerTL,outerCols);   % last column of outer layout
    hold(axR,'on');

    % warn if *all* RMSE values are missing
    if all(isnan(RMSE))
        warning('plotDMA: all RMSE values are NaN – RMSE panel will be empty.');
    end

    plot(axR,EFC,RMSE,'-o','Color',tumBlue, ...
         'LineWidth',lnW,'MarkerSize',mkSz);

    grid(axR,'on');  xlim(axR,'padded');
    ylim(axR,'padded');  ylim(axR,[0, max(axR.YLim)]);  % force baseline at 0

    box on;

    title(axR,'RMSE','Interpreter','latex','FontSize',1.2*fontSz);
    ylabel(axR,'RMSE / mV','Interpreter','latex','FontSize',fontSz);
    if CalendarOrCyclic == 0    % Calendar Aging
       xlabel('RPT Number / -', 'Interpreter','latex', 'FontSize', fontSz);
    else                         % Cyclic Aging
        xlabel('EFC', 'Interpreter','latex', 'FontSize', fontSz);
    end

    axR.TickLength = [0 0];
    axR.TickDir    = 'out';
    axR.LineWidth  = lnW;
    axR.FontSize   = fontSz;
    axR.FontName   = 'Times New Roman';

    
end



end
%==========================================================================

