function plot_DMA(Data, varargin)
% plotDMA  (tiledlayout edition)                                  2025-06-18
%--------------------------------------------------------------------------
%> Author: Mathias Rehm (mathias.rehm@tum.de)
%> additional code by Can Korkmaz (can.korkmaz@tum.de)
%> Date: 2025-06-18
%>
%>             • kept digit-based CU filter and robust RMSE detection
%>             • uses nested tiledlayout so the RMSE panel gets extra
%>               horizontal space for its y-label (replicates the old
%>               manual gapRMSE = 0.07 behaviour)
%>             • figure export code permanently removed (keywords ignored)
%
% DESCRIPTION
%   • Plots the main degradation modes (LAM / LLI) versus *EFC*.
%   • Works with Data structs whose field names contain digits (“CU0” …);
%     non-digit housekeeping fields are skipped.
%   • Optional right-most panel shows RMSE; recognises
%       -  RMSE_0_9   -  RMSE_full   -  RMSE      (all V → mV).
%   • 3-line mode  (Cathode / Anode / ChargeCarrierInv)  **or**
%     4-line mode  (Cathode / An-blend1 / An-blend2 / ChargeCarrierInv) via the
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
%--------------------------------------------------------------------------

%% 1 — Set defaults --------------------------------------------------------
% Blend toggles are independent; legacy useBlendElectrodeModel maps to the anode blend flag.
useAnodeBlendModel    = false;
useCathodeBlendModel  = false;
EFC         = [];
plotTitles  = {};
plotRMSE    = true;       % show RMSE panel unless explicitly disabled
CalendarOrCyclic = 0;
labelCathode     = 'Cathode';
labelAnode       = 'Anode';
labelAnodeBlend1 = 'An-blend1';
labelAnodeBlend2 = 'An-blend2';
labelCathodeBlend1 = 'Ca-blend1';
labelCathodeBlend2 = 'Ca-blend2';
labelChargeCarrierInv = 'Charge-carrier-inv';
plotCathode      = true;  % Use only to hide cathode for LFP (cathode aging not meaningful)

%% 2 - Parse varargin ------------------------------------------------------
for v = 1:2:numel(varargin)
    key = lower(string(varargin{v}));
    val = varargin{v+1};

    switch key
        case {'useblendelectrodemodel','blend','useblend'} % legacy: maps to anode blend
            useAnodeBlendModel = logical(val);
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
            CalendarOrCyclic = val;

        case 'labels'
            if isstruct(val)
                if isfield(val,'labelCathode'),          labelCathode          = val.labelCathode; end
                if isfield(val,'labelAnode'),            labelAnode            = val.labelAnode; end
                if isfield(val,'labelAnodeBlend1'),      labelAnodeBlend1      = val.labelAnodeBlend1; end
                if isfield(val,'labelAnodeBlend2'),      labelAnodeBlend2      = val.labelAnodeBlend2; end
                if isfield(val,'labelCathodeBlend1'),    labelCathodeBlend1    = val.labelCathodeBlend1; end
                if isfield(val,'labelCathodeBlend2'),    labelCathodeBlend2    = val.labelCathodeBlend2; end
                if isfield(val,'labelChargeCarrierInv'), labelChargeCarrierInv = val.labelChargeCarrierInv; end
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

        % ---------- legacy / ignored keys --------------------------------
        case {'savepathdma','savefig','savefigure','savepath'}
            % These keys are silently accepted so old scripts do not error,
            % but they no longer trigger any file export.

        otherwise
            warning('plotDMA: unknown option %s ignored.', varargin{v});
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
  varNames = {};
  lamExtractors = {};

  if plotCathode
      varNames{end+1} = labelCathode;
      lamExtractors{end+1} = @(d)100*d.LAM_cathode;
      if useCathodeBlendModel
          varNames = [varNames, {labelCathodeBlend1, labelCathodeBlend2}];
          lamExtractors = [lamExtractors, {@(d)100*d.LAM_cathode_blend1, @(d)100*d.LAM_cathode_blend2}];
      end
  end

  % Anode aggregate always plotted; add blends if enabled
  varNames{end+1} = labelAnode;
  lamExtractors{end+1} = @(d)100*d.LAM_anode;
  if useAnodeBlendModel
      varNames = [varNames, {labelAnodeBlend1, labelAnodeBlend2}];
      lamExtractors = [lamExtractors, {@(d)100*d.LAM_anode_blend1, @(d)100*d.LAM_anode_blend2}];
  end

  % Charge-carrier inventory always plotted
  varNames{end+1} = labelChargeCarrierInv;
  lamExtractors{end+1} = @(d)100*d.LI;

  nPlots = numel(varNames);
  
  LAM = zeros(nCU,nPlots);          % rows = CU, cols = curves
  RMSE = NaN(nCU,1);                % will stay NaN if plotRMSE == false

for k = 1:nCU
    f = cuFld{k};
    % ---- LAM / LI ------------------------------------------------------
    for idxLam = 1:nPlots
        LAM(k,idxLam) = lamExtractors{idxLam}(Data.(f));
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

%% 5 — Figure & nested tiledlayout to create RMSE gap -----------------
set(groot,'defaultAxesTickLabelInterpreter','latex');

fig = figure('Units','centimeters','Position',[3 3 20 6]);

outerCols = nPlots + double(plotRMSE);   % +1 column only if RMSE requested
outerTL   = tiledlayout(fig,1,outerCols, ...
              'TileSpacing','compact','Padding','compact');


% place an inner tiledlayout spanning the first nPlots tiles of outerTL
innerTL = tiledlayout(outerTL,1,nPlots, ...
              'TileSpacing','none','Padding','compact');
innerTL.Layout.Tile     = 1;
innerTL.Layout.TileSpan = [1 nPlots];    % occupies first nPlots columns

% Result: LAM axes have zero internal gap, while outerTL’s inter-tile
% spacing provides the extra gap before the RMSE axis.

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

    % warn if all RMSE values are missing
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
