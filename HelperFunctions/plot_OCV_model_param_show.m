function plot_OCV_model_param_show(solution, i, cellName)
%> Authors : Mathias Rehm
%> E-mail : mathias.rehm@tum.de
%> Date    : 2025-10-09
%
% plot_OCV_model_param_show
% Tiled OCV plus DVA plotter with dynamic limits and parameterized styling

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
yMaxDVA_param        = 3.2;            % cap for DVA y max

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
meas_fullCell_SOC     = solution.myData.Q_cell;
meas_fullCell_voltage = solution.myData.OCV_cell;

calc_fullCell_SOC     = solution.reconSOC;
full_cell_voltage     = solution.fcU_model;

meas_Q_DVA            = solution.Q_DVA_meas;
meas_DVA              = solution.DVA_smooth_meas;
calc_Q_DVA            = solution.Q_DVA_calc;
calc_DVA              = solution.DVA_smooth_calc;

params                = solution.params;
cathSOC               = solution.cathSOC;
cathU_recon           = solution.cathU_recon;
anodeSOC              = solution.anodeSOC;
anodeU_recon          = solution.anodeU_recon;

alpha_an  = params(1);
beta_an   = params(2);
alpha_cat = params(3);
beta_cat  = params(4);

%% 3) Defaults for fonts
set(groot, {'defaultAxesFontSize','defaultTextFontSize','defaultAxesTickLabelInterpreter'}, ...
           {baseFont, baseFont, 'latex'});

%% 4) Precompute DVA of electrodes
[Q_DVA_cath, ~, DVA_cath] = calculate_DVA(cathSOC,  cathU_recon);
[Q_DVA_an  , ~, DVA_an  ] = calculate_DVA(anodeSOC, anodeU_recon);
DVA_cath_s = smooth(DVA_cath,30,'lowess');
DVA_an_s   = smooth(DVA_an  ,30,'lowess');

%% 5) Dynamic limits with arrow space
% 5.1 x limits
xAll = [meas_fullCell_SOC(:); calc_fullCell_SOC(:); ...
        cathSOC(:); anodeSOC(:); ...
        meas_Q_DVA(:); calc_Q_DVA(:); ...
        Q_DVA_cath(:); Q_DVA_an(:)];
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
yAllOCV = [meas_fullCell_voltage(:); full_cell_voltage(:); cathU_recon(:); anodeU_recon(:)];
if isempty(yAllOCV) || all(~isfinite(yAllOCV))
    yMinData = 0; yMaxData = 4.4;
else
    yMinData = min(yAllOCV(isfinite(yAllOCV)));
    yMaxData = max(yAllOCV(isfinite(yAllOCV)));
end
rangeY   = max(yMaxData - yMinData, 1);
yNeedTop = max([yMaxData, max(cathU_recon)+alphaArrowYOffset, max(anodeU_recon)+alphaAnTextUp]);
yNeedBot = min([yMinData, min(cathU_recon)+betaCatArrowYOffset, min(anodeU_recon)+betaAnArrowYOffset]);
yLimOCV  = [yNeedBot - yPadLowerOCV, yNeedTop + yPadUpperOCVFrac*rangeY];

% 5.3 y limits for DVA with cap
yAllDVA  = [meas_DVA(:); calc_DVA(:); DVA_cath_s(:); abs(DVA_an_s(:))];
if isempty(yAllDVA) || all(~isfinite(yAllDVA))
    yMinD = 0; yMaxD = yMaxDVA_param;
else
    yMinD = min(yAllDVA(isfinite(yAllDVA)));
    yMaxD = max(yAllDVA(isfinite(yAllDVA)));
end
rangeYD  = max(yMaxD - yMinD, 1);
yLimDVA  = [max(0, yMinD - yPadLowerDVAFrac*rangeYD), yMaxDVA_param];

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
h1 = plot(ax1, meas_fullCell_SOC, meas_fullCell_voltage, ...
          'Color', colMeas,  'LineWidth', lwMain, 'LineStyle', lineStyleMeas);
h2 = plot(ax1, calc_fullCell_SOC, full_cell_voltage, ...
          'Color', colModel, 'LineWidth', lwMain, 'LineStyle', lineStyleModel);
h3 = plot(ax1, cathSOC, cathU_recon, ...
          'Color', colCathode, 'LineWidth', lwMain, 'LineStyle', lineStyleModel);
h4 = plot(ax1, anodeSOC, anodeU_recon, ...
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
    'FC measured', 'FC reconstructed', 'Cathode reconstructed', 'Anode reconstructed', ...
    'Interpreter','latex','FontSize',legendFont,'Orientation','horizontal');
lgd.Layout.Tile = 'south';

drawnow('nocallbacks');

%% 8) Arrows and labels on OCV
arrowInfo = repmat(struct('h',gobjects(1),'dataX',[],'dataY',[],'ax',ax1),1,4);
arrowOpts = {'LineWidth',arrowLineWidth,'Head1Width',arrowHeadWidth,'Head2Width',arrowHeadWidth, ...
             'Head1Length',arrowHeadLength,'Head2Length',arrowHeadLength,'HeadStyle','plain'};

% alpha cathode
arrowInfo(1) = createArrowAndStore([min(cathSOC) max(cathSOC)], ...
                                   [max(cathU_recon)+alphaArrowYOffset max(cathU_recon)+alphaArrowYOffset], ...
                                   colCathode, arrowOpts, ax1, fig);
text(mean([min(cathSOC) max(cathSOC)]), max(cathU_recon)-alphaCatTextDown, ...
     sprintf('$\\alpha_{\\mathrm{cat}} = %.2f$', alpha_cat), ...
     'Parent',ax1,'Interpreter','latex','HorizontalAlignment','center', ...
     'Color',colCathode,'FontSize',baseFont);

% alpha anode
arrowInfo(2) = createArrowAndStore([min(anodeSOC) max(anodeSOC)], ...
                                   [max(anodeU_recon)+alphaArrowYOffset max(anodeU_recon)+alphaArrowYOffset], ...
                                   colAnode, arrowOpts, ax1, fig);
text(mean([min(anodeSOC) max(anodeSOC)]), max(anodeU_recon)+alphaAnTextUp, ...
     sprintf('$\\alpha_{\\mathrm{an}} = %.2f$', alpha_an), ...
     'Parent',ax1,'Interpreter','latex','HorizontalAlignment','center', ...
     'Color',colAnode,'FontSize',baseFont);

% beta cathode (shifted upwards a bit)
arrowInfo(3) = createArrowAndStore([min(cathSOC) 0], ...
                                   [min(cathU_recon)+betaCatArrowYOffset min(cathU_recon)+betaCatArrowYOffset], ...
                                   colCathode, arrowOpts, ax1, fig);
text(0.03, min(cathU_recon)+betaCatTextUp, ...
     sprintf('$\\beta_{\\mathrm{cat}} = %.2f$', beta_cat), ...
     'Parent',ax1,'Interpreter','latex','HorizontalAlignment','left', ...
     'Color',colCathode,'FontSize',baseFont);

% beta anode with more headroom below
arrowInfo(4) = createArrowAndStore([min(anodeSOC) 0], ...
                                   [min(anodeU_recon)+betaAnArrowYOffset min(anodeU_recon)+betaAnArrowYOffset], ...
                                   colAnode, arrowOpts, ax1, fig);
text(0.03, min(anodeU_recon)+betaAnTextUp, ...
     sprintf('$\\beta_{\\mathrm{an}} = %.2f$', round(beta_an,2)+0), ...
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
plot(ax2, meas_Q_DVA, meas_DVA, ...
     'Color', colMeas,  'LineWidth', lwDVA, 'LineStyle', lineStyleMeas);
plot(ax2, calc_Q_DVA, calc_DVA, ...
     'Color', colModel, 'LineWidth', lwDVA, 'LineStyle', lineStyleModel);
plot(ax2, Q_DVA_cath, DVA_cath_s, ...
     'Color', colCathode, 'LineWidth', lwDVA, 'LineStyle', lineStyleModel);
plot(ax2, Q_DVA_an,   abs(DVA_an_s), ...
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
