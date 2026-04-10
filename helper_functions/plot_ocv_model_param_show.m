function plot_ocv_model_param_show(solution, i, cellName, labelCfg)
%> Authors : Mathias Rehm
%> E-mail : mathias.rehm@tum.de
%> Date    : 2025-10-09
%
% plot_ocv_model_param_show
% Tiled OCV plus DVA plotter with dynamic limits and parameterized styling
%
% INPUTS:
%   solution  - best-solution struct returned by store_solution_struct
%   i         - current CU or RPT index for titles
%   cellName  - optional cell identifier shown in the figure title
%   labelCfg  - struct with electrode and lithium-inventory display labels
%
% OUTPUTS:
%   none      - creates a parameterized reconstruction plot for one CU

%% 1) Parameters
% 1.1 figure and fonts
figPos_cm            = [3 2];          % lower left corner in cm
figWidth_cm          = 20;             % keep width
figHeight_cm         = 16;             % reduced height, bottom stays fixed
titleFont            = 14;
baseFont             = 12;
axisTickFont         = 11;
legendFont           = 11;

% 1.2 line styles and widths
lineStyleModel       = '-.';           % model and electrodes
lineStyleMeas        = '-';            % measurement
lwMain               = 2;              % default thickness for OCV curves
lwDVA                = 2;              % default thickness for DVA curves
lwGuide              = 1;              % guide lines at SOC 0 and 1
axesBoxLineWidth     = 1.5;            % axes box thickness

% 1.3 arrows
arrowLineWidth       = 2;              % twice lwMain by default
arrowHeadWidth       = 4;
arrowHeadLength      = 4;

% 1.4 colors
colMeas              = 0.5*[1 1 1];    % gray
colModel             = [0 0 0];        % black
colCathode           = [162 173  0]/255;
colAnode             = [ 48 112 179]/255;

% 1.5 dynamic limit padding
xPadLeftFrac         = 0.05;           % reduced left pad of data span
xPadRightFrac        = 0.08;           % reduced right pad of data span
xPadZeroLeftFrac     = 0.02;           % extra pad left of zero for beta labels

yPadLowerOCV         = 0.30;           % more space below for U
yPadUpperOCVFrac     = 0.04;           % top fractional pad for U

yPadLowerDVAFrac     = 0.05;           % bottom pad for DVA
yMaxDVAParam         = 3.2;            % cap for DVA y max

% 1.6 label offsets in V
alphaCatTextDown     = 0.12;           % below cathode max
alphaAnTextUp        = 0.20;           % above anode max
betaCatTextUp        = 0.50;           % slightly higher than before
betaAnTextUp         = 0.00;           % near anode min

% 1.7 arrow vertical offsets in V
alphaArrowYOffset    = 0.05;           % above maxima
betaCatArrowYOffset  = 0.15;           % above cathode min
betaAnArrowYOffset   = 0.12;           % a bit higher for comfort

%% 2) Load data
measFullCellSOC      = solution.solverInput.qCell;
measFullCellVoltage  = solution.solverInput.ocvCell;
calcFullCellSOC      = solution.reconSOC;
calcFullCellVoltage  = solution.fullCellUModel;

measQDVA              = solution.qDVAMeas;
measDVA               = solution.dvaSmoothMeas;
calcQDVA              = solution.qDVACalc;
calcDVA               = solution.dvaSmoothCalc;

params                = solution.params;
cathodeSOC               = solution.cathodeSOC;
cathodeURecon           = solution.cathodeURecon;
anodeSOC              = solution.anodeSOC;
anodeURecon          = solution.anodeURecon;

alphaAn  = params(1);
betaAn   = params(2);
alphaCat = params(3);
betaCat  = params(4);

% 2.1 labels
defaultLabels = struct('labelCathode','Cathode', ...
    'labelAnode','Anode');
if nargin < 4 || isempty(labelCfg)
    labelCfg = defaultLabels;
end
if ~isfield(labelCfg,'labelCathode'), labelCfg.labelCathode = defaultLabels.labelCathode; end
if ~isfield(labelCfg,'labelAnode'),   labelCfg.labelAnode   = defaultLabels.labelAnode;   end

%% 3) Defaults for fonts
set(groot, {'defaultAxesFontSize','defaultTextFontSize','defaultAxesTickLabelInterpreter'}, ...
           {baseFont, baseFont, 'latex'});

%% 4) Precompute DVA of electrodes
[qDVACathode, ~, dvaCathode] = calculate_dva(cathodeSOC, cathodeURecon);
[qDVAAnode, ~, dvaAnode] = calculate_dva(anodeSOC, anodeURecon);
dvaCathodeSmooth = smooth(dvaCathode, 30, 'lowess');
dvaAnodeSmooth   = smooth(dvaAnode, 30, 'lowess');

%% 5) Dynamic limits with arrow space
% 5.1 x limits
xAll = [measFullCellSOC(:); calcFullCellSOC(:); ...
        cathodeSOC(:); anodeSOC(:); ...
        measQDVA(:); calcQDVA(:); ...
        qDVACathode(:); qDVAAnode(:)];
if isempty(xAll) || all(~isfinite(xAll))
    xMinData = 0; xMaxData = 1;
else
    xMinData = min(xAll(isfinite(xAll)));
    xMaxData = max(xAll(isfinite(xAll)));
end
rangeX   = max(xMaxData - xMinData, 1);
xMinNeed = min([xMinData, 0]) - (xPadLeftFrac + xPadZeroLeftFrac)*rangeX;
xMaxNeed = xMaxData + xPadRightFrac*rangeX;
xLimDyn  = [xMinNeed, xMaxNeed];

% 5.2 y limits for OCV
yAllOCV = [measFullCellVoltage(:); calcFullCellVoltage(:); cathodeURecon(:); anodeURecon(:)];
if isempty(yAllOCV) || all(~isfinite(yAllOCV))
    yMinData = 0; yMaxData = 4.4;
else
    yMinData = min(yAllOCV(isfinite(yAllOCV)));
    yMaxData = max(yAllOCV(isfinite(yAllOCV)));
end
rangeY   = max(yMaxData - yMinData, 1);
yNeedTop = max([yMaxData, max(cathodeURecon)+alphaArrowYOffset, max(anodeURecon)+alphaAnTextUp]);
yNeedBot = min([yMinData, min(cathodeURecon)+betaCatArrowYOffset, min(anodeURecon)+betaAnArrowYOffset]);
yLimOCV  = [yNeedBot - yPadLowerOCV, yNeedTop + yPadUpperOCVFrac*rangeY];

% 5.3 y limits for DVA with cap
yAllDVA  = [measDVA(:); calcDVA(:); dvaCathodeSmooth(:); abs(dvaAnodeSmooth(:))];
if isempty(yAllDVA) || all(~isfinite(yAllDVA))
    yMinD = 0; yMaxD = yMaxDVAParam;
else
    yMinD = min(yAllDVA(isfinite(yAllDVA)));
    yMaxD = max(yAllDVA(isfinite(yAllDVA)));
end
rangeYD  = max(yMaxD - yMinD, 1);
yLimDVA  = [max(0, yMinD - yPadLowerDVAFrac*rangeYD), yMaxDVAParam];

%% 6) Figure and layout
fig = figure('Units','centimeters', ...
             'Position',[figPos_cm figWidth_cm figHeight_cm], ...
             'Color','w', 'PaperPositionMode','auto');

tl  = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

if ~isempty(cellName)
    sgtitle(sprintf('%s CU%d',cellName,i), 'Interpreter','latex', ...
            'FontWeight','bold','FontSize',titleFont,'Color','k');
else
    sgtitle(sprintf('CU%d',i), 'Interpreter','latex', ...
            'FontWeight','bold','FontSize',titleFont,'Color','k');
end

%% 7) OCV (top two rows)
ax1 = nexttile(tl,[2 1]); hold(ax1,'on'); box(ax1,'on');
ax1.LineWidth = axesBoxLineWidth;

% OCV curves use lwMain
h1 = plot(ax1, measFullCellSOC, measFullCellVoltage, ...
          'Color', colMeas,  'LineWidth', lwMain, 'LineStyle', lineStyleMeas);
h2 = plot(ax1, calcFullCellSOC, calcFullCellVoltage, ...
          'Color', colModel, 'LineWidth', lwMain, 'LineStyle', lineStyleModel);
h3 = plot(ax1, cathodeSOC, cathodeURecon, ...
          'Color', colCathode, 'LineWidth', lwMain, 'LineStyle', lineStyleModel);
h4 = plot(ax1, anodeSOC, anodeURecon, ...
          'Color', colAnode, 'LineWidth', lwMain, 'LineStyle', lineStyleModel);

xlim(ax1, xLimDyn);
ylim(ax1, yLimOCV);
yticks(ax1, 0:1:4);
set(ax1,'XTickLabel',[]);
ylabel(ax1,'$U$ / V','Interpreter','latex');

% SOC guides
xl_now = xlim(ax1);
if xl_now(1) <= 0 && 0 <= xl_now(2)
    plot(ax1,[0 0],[yLimOCV(1) yLimOCV(2)],'--k','LineWidth',lwGuide,'HandleVisibility','off');
end
if xl_now(1) <= 1 && 1 <= xl_now(2)
    plot(ax1,[1 1],[yLimOCV(1) yLimOCV(2)],'--k','LineWidth',lwGuide,'HandleVisibility','off');
end

ax1.XAxis.FontSize = axisTickFont;
ax1.YAxis.FontSize = axisTickFont;
grid(ax1,'on');

% legend text updated
lgd = legend(ax1,[h1 h2 h3 h4], ...
    'FC measured', 'FC reconstructed', ...
    sprintf('%s reconstructed', labelCfg.labelCathode), ...
    sprintf('%s reconstructed', labelCfg.labelAnode), ...
    'Interpreter','latex','FontSize',legendFont,'Orientation','horizontal');
lgd.Layout.Tile = 'south';

drawnow('nocallbacks');

%% 8) Arrows and labels on OCV
arrowInfo = repmat(struct('h',gobjects(1),'dataX',[],'dataY',[],'ax',ax1),1,4);
arrowOpts = {'LineWidth',arrowLineWidth,'Head1Width',arrowHeadWidth,'Head2Width',arrowHeadWidth, ...
             'Head1Length',arrowHeadLength,'Head2Length',arrowHeadLength,'HeadStyle','plain'};

% alpha cathode
arrowInfo(1) = createArrowAndStore([min(cathodeSOC) max(cathodeSOC)], ...
                                   [max(cathodeURecon)+alphaArrowYOffset max(cathodeURecon)+alphaArrowYOffset], ...
                                   colCathode, arrowOpts, ax1, fig);
text(mean([min(cathodeSOC) max(cathodeSOC)]), max(cathodeURecon)-alphaCatTextDown, ...
     sprintf('$\\alpha_{\\mathrm{cat}} = %.2f$', alphaCat), ...
     'Parent',ax1,'Interpreter','latex','HorizontalAlignment','center', ...
     'Color',colCathode,'FontSize',baseFont);

% alpha anode
arrowInfo(2) = createArrowAndStore([min(anodeSOC) max(anodeSOC)], ...
                                   [max(anodeURecon)+alphaArrowYOffset max(anodeURecon)+alphaArrowYOffset], ...
                                   colAnode, arrowOpts, ax1, fig);
text(mean([min(anodeSOC) max(anodeSOC)]), max(anodeURecon)+alphaAnTextUp, ...
     sprintf('$\\alpha_{\\mathrm{an}} = %.2f$', alphaAn), ...
     'Parent',ax1,'Interpreter','latex','HorizontalAlignment','center', ...
     'Color',colAnode,'FontSize',baseFont);

% beta cathode (shifted upwards a bit)
arrowInfo(3) = createArrowAndStore([min(cathodeSOC) 0], ...
                                   [min(cathodeURecon)+betaCatArrowYOffset min(cathodeURecon)+betaCatArrowYOffset], ...
                                   colCathode, arrowOpts, ax1, fig);
text(0.03, min(cathodeURecon)+betaCatTextUp, ...
     sprintf('$\\beta_{\\mathrm{cat}} = %.2f$', betaCat), ...
     'Parent',ax1,'Interpreter','latex','HorizontalAlignment','left', ...
     'Color',colCathode,'FontSize',baseFont);

% beta anode with more headroom below
arrowInfo(4) = createArrowAndStore([min(anodeSOC) 0], ...
                                   [min(anodeURecon)+betaAnArrowYOffset min(anodeURecon)+betaAnArrowYOffset], ...
                                   colAnode, arrowOpts, ax1, fig);
text(0.03, min(anodeURecon)+betaAnTextUp, ...
     sprintf('$\\beta_{\\mathrm{an}} = %.2f$', round(betaAn,2)+0), ...
     'Parent',ax1,'Interpreter','latex','HorizontalAlignment','left', ...
     'Color',colAnode,'FontSize',baseFont);

% bind listeners without anonymous inline code
updateArrows();
l1 = addlistener(ax1,'Position','PostSet',@cbUpdateArrows);
l2 = addlistener(ax1,'XLim','PostSet',   @cbUpdateArrows);
l3 = addlistener(ax1,'YLim','PostSet',   @cbUpdateArrows);
l4 = addlistener(fig,'SizeChanged',      @cbUpdateArrows);
oldL = getappdata(fig,'arrow_listeners'); if isempty(oldL), oldL = []; end
setappdata(fig,'arrow_listeners',[oldL l1 l2 l3 l4]);

%% 9) DVA (bottom row)
ax2 = nexttile(tl,3); hold(ax2,'on'); box(ax2,'on');
ax2.LineWidth = axesBoxLineWidth;

% DVA curves use lwDVA
plot(ax2, measQDVA, measDVA, ...
     'Color', colMeas,  'LineWidth', lwDVA, 'LineStyle', lineStyleMeas);
plot(ax2, calcQDVA, calcDVA, ...
     'Color', colModel, 'LineWidth', lwDVA, 'LineStyle', lineStyleModel);
plot(ax2, qDVACathode, dvaCathodeSmooth, ...
     'Color', colCathode, 'LineWidth', lwDVA, 'LineStyle', lineStyleModel);
plot(ax2, qDVAAnode, abs(dvaAnodeSmooth), ...
     'Color', colAnode,   'LineWidth', lwDVA, 'LineStyle', lineStyleModel);

xlim(ax2, xLimDyn);
ylim(ax2, yLimDVA);
yticks(ax2, 0:1:4);
xlabel(ax2,'SOC / -','Interpreter','latex');
ylabel(ax2,'$dU(dQ)^{-1}\cdot C_{\mathrm{act}}$ / V','Interpreter','latex');
ax2.XAxis.FontSize = axisTickFont;
ax2.YAxis.FontSize = axisTickFont;
grid(ax2,'on');

%% 10) Helpers
    function s = createArrowAndStore(xDatIn,yDatIn,col,opts,ax,figH)
        hAnn    = annotation(figH,'doublearrow',[0 0],[0 0], 'Color',col, opts{:});
        s.h     = hAnn;
        s.dataX = xDatIn;
        s.dataY = yDatIn;
        s.ax    = ax;
    end

    function cbUpdateArrows(~,~)
        updateArrows();
    end

    function updateArrows()
        % recalc figure normalized coords for each stored arrow
        for kk = 1:numel(arrowInfo)
            axCur = arrowInfo(kk).ax;
            [xFig,yFig] = data2norm(axCur, arrowInfo(kk).dataX, arrowInfo(kk).dataY);
            set(arrowInfo(kk).h, 'X', xFig, 'Y', yFig);
        end
    end

    function [xf,yf] = data2norm(axH,x,y)
        axUnits = get(axH,'Units');
        set(axH,'Units','normalized');
        axPos = get(axH,'Position');
        set(axH,'Units',axUnits);

        xl = get(axH,'XLim'); yl = get(axH,'YLim');
        isLogX = strcmp(get(axH,'XScale'),'log');
        isLogY = strcmp(get(axH,'YScale'),'log');

        if isLogX, x = log10(x); xl = log10(xl); end
        if isLogY, y = log10(y); yl = log10(yl); end

        xnorm = (x - xl(1)) ./ diff(xl);
        ynorm = (y - yl(1)) ./ diff(yl);

        xf = axPos(1) + axPos(3) .* xnorm;
        yf = axPos(2) + axPos(4) .* ynorm;
    end

end
