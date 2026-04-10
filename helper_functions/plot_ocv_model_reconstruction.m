function plot_ocv_model_reconstruction(measFullCellSOC, measFullCellVoltage, calcFullCellSOC, calcFullCellVoltage, qDVAMeas, dvaSmoothMeas, qDVACalc, dvaSmoothCalc, ~, ...
    cathodeSOC, cathodeURecon, anodeSOC, anodeURecon, finalRMSE, i, cellName)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-10-10
%
% plot_ocv_model_reconstruction
% Tiled OCV plus DVA plotter with dynamic limits and parameterized styling
%
% INPUTS:
%   measFullCellSOC      - measured full-cell SOC vector
%   measFullCellVoltage  - measured full-cell voltage vector
%   calcFullCellSOC      - reconstructed full-cell SOC vector
%   calcFullCellVoltage  - reconstructed full-cell voltage vector
%   qDVAMeas             - measured DVA x-axis
%   dvaSmoothMeas        - measured smoothed DVA curve
%   qDVACalc             - calculated DVA x-axis
%   dvaSmoothCalc        - calculated smoothed DVA curve
%   ~                    - unused legacy ICA placeholder
%   cathodeSOC           - cathode SOC axis used for the electrode plot
%   cathodeURecon        - cathode potential curve
%   anodeSOC             - anode SOC axis used for the electrode plot
%   anodeURecon          - anode potential curve
%   finalRMSE            - optional RMSE value shown in the title
%   i                    - current CU or RPT index
%   cellName             - optional cell identifier shown in the title
%
% OUTPUTS:
%   none                 - creates a tiled reconstruction figure

    colors = tum_colors();

    fontSize = 12;
    lineSize = 2;
    titleSize = 16;
    figSize1 = [3 2 25 15];
    fontName    = 'Times New Roman';

    % Hard-coded switch to show RMSE in title
    showRMSE = ~isempty(finalRMSE) && isfinite(finalRMSE);

    rmse_mV = finalRMSE * 1000;
    set(groot, 'defaultAxesTickLabelInterpreter','latex');

    % Full Cell
    [~, ocvICAFullCellMeas, icaFullCellMeas] = calculate_ica(measFullCellSOC, measFullCellVoltage);
    [~, ocvICAFullCellCalc, icaFullCellCalc] = calculate_ica(calcFullCellSOC, calcFullCellVoltage);
    icaSmoothMeas = smooth(icaFullCellMeas, 30, 'lowess');
    icaSmoothCalc = smooth(icaFullCellCalc, 30, 'lowess');

    % Cathode
    [qDVACathode, ~, dvaCathode] = calculate_dva(cathodeSOC, cathodeURecon);
    dvaCathodeSmooth = smooth(dvaCathode, 30, 'lowess');

    % Anode
    [qDVAAnode, ~, dvaAnode] = calculate_dva(anodeSOC, anodeURecon);
    dvaAnodeSmooth = smooth(dvaAnode, 30, 'lowess');

    fig_offset = [-0.04 0 0 0];
    
    fig = figure();
    set(fig, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
            'PaperUnits', 'Centimeters','Position', figSize1,...
            'PaperSize', [figSize1(3) figSize1(4)]);

    if showRMSE
        if ~isempty(cellName)
            sgtitle(sprintf("C%s CU%d (RMSE: %.1f mV)", cellName, i, rmse_mV), 'Interpreter', 'latex', 'FontSize', titleSize, 'FontName', fontName, 'FontWeight','bold');
        else
            sgtitle(sprintf("CU%d (RMSE: %.1f mV)", i, rmse_mV), 'Interpreter', 'latex', 'FontSize', titleSize, 'FontName', fontName, 'FontWeight','bold');
        end
    else
        if ~isempty(cellName)
            sgtitle(sprintf("C%s CU%d", cellName, i), 'Interpreter', 'latex', 'FontSize', titleSize, 'FontName', fontName, 'FontWeight','bold');
        else
            sgtitle(sprintf("CU%d", i), 'Interpreter', 'latex', 'FontSize', titleSize, 'FontName', fontName, 'FontWeight','bold');
        end
    end
   
    ax = subplot(4,4,[2 3 4 6 7 8 10 11 12]);
    set(ax, 'ColorOrder', colors.colorOrder);
    
    % add extra plots for legend
    h1 = plot(NaN, NaN, 'k-', 'DisplayName', 'model','linewidth', lineSize);
    hold on;
    h2 = plot(NaN, NaN, 'k:', 'DisplayName', 'measurement','linewidth', lineSize);
    hold on;

    plot(calcFullCellSOC, calcFullCellVoltage, '-', 'Color', colors.tumBlue, 'linewidth', lineSize);
    hold on;
    plot(measFullCellSOC, measFullCellVoltage, ':', 'Color', colors.darkGray, 'linewidth', lineSize);
    plot(cathodeSOC, cathodeURecon, '-', 'Color', colors.mediumGray, 'linewidth', lineSize);
    maxVoltage = max([max(measFullCellVoltage), max(calcFullCellVoltage), max(cathodeURecon)]) + 0.05;
    minVoltage = min([min(measFullCellVoltage), min(calcFullCellVoltage)]) - 0.05;
    xlim([min(cathodeSOC)-0.1, max(cathodeSOC)+0.1]);
    ylim([minVoltage, maxVoltage]);
    set(gca, 'XTickLabel', [], 'Position', ax.Position + fig_offset);
    
    [ax, h] = addaxis(measFullCellSOC, (measFullCellVoltage - calcFullCellVoltage)*1000, [-20, 20], ':', 'Color', colors.green, 'linewidth', lineSize*0.7);
    addaxislabel(2, 'Full Cell Fitting Error / mV', 'Fontname', fontName, 'Interpreter', 'latex');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    set(ax, 'box', 'off', 'TickDir', 'out', 'linewidth', lineSize, 'fontSize', fontSize, 'Fontname', fontName, 'Ticklabelinterpreter', 'latex', 'Position', ax.Position + [+0.06 0 0 0]);
    
    yyaxis right;
    plot(anodeSOC, anodeURecon, '-', 'Color', colors.orange, 'linewidth', lineSize);
    hold on;
    ylabel('Anode Potential vs Li/Li$^+$ / V', 'fontSize', fontSize, 'FontName', fontName, 'Interpreter', 'latex');
    set(gca, 'box', 'off', 'TickDir', 'out', 'YColor', colors.orange, 'linewidth', lineSize, 'fontSize', fontSize, 'FontName', fontName, 'TickLabelInterpreter', 'latex');
    
    ax = subplot(4,4, [1 5 9]);
    set(ax, 'ColorOrder', colors.colorOrder);
    plot(icaSmoothCalc, ocvICAFullCellCalc, '-', 'Color', colors.tumBlue, 'linewidth', lineSize);
    hold on;
    maxICA = max([max(icaSmoothMeas), max(ocvICAFullCellMeas)]) + 5;
    plot(icaSmoothMeas, ocvICAFullCellMeas, ':', 'Color', colors.darkGray, 'linewidth', lineSize);
    plot(1000, 1000, '-', 'Color', colors.mediumGray, 'linewidth', lineSize);
    plot(1000, 1000, '-', 'Color', colors.orange, 'linewidth', lineSize);
    ylim([minVoltage, maxVoltage]);
    xlim([0, maxICA]);
    ylabel('Full Cell Voltage / V', 'fontSize', fontSize, 'FontName', fontName, 'Interpreter', 'latex');
    xlabel("$\frac{\Delta Q}{\Delta U}$ / $\frac{\mathrm{Ah}}{\mathrm{V}}$", 'fontSize', fontSize, 'FontName', fontName, 'Interpreter', 'latex');
    legend([h1, h2], 'Position', [0.1 0.15 0.08 0.05], 'Units', 'normalized', 'FontSize', fontSize, 'Interpreter', 'latex');
    
    % JE to MRe: yticks ok? Do we need ylabel againg for full-cell voltage?
    % -> ICA scaling is different to improve graphics
    set(gca, 'box', 'off', 'TickDir', 'out', 'YColor', colors.tumBlue, 'linewidth', lineSize, 'fontSize', fontSize, 'FontName', fontName, 'TickLabelInterpreter', 'latex', 'Position', ax.Position + fig_offset + [0.02 0 -0.02 0]);
    
    [ax, ~] = addaxis([], []);
    set(ax, 'visible', false);
    
    [ax, h] = addaxis(NaN, NaN, [2.3, 4.4], '--', 'Color', colors.mediumGray, 'linewidth', lineSize*0.7);
    addaxislabel(3, 'Cathode Potential vs Li/Li^+ / V');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    set(ax, 'box', 'off', 'TickDir', 'out', 'linewidth', lineSize, 'fontSize', fontSize, 'Fontname', fontName, 'Ticklabelinterpreter', 'latex', 'Position', ax.Position + [-0.05 0 0 0]); %'Fontname', fontName,
     
    ax = subplot(4,4, [14 15 16]);
    set(ax, 'ColorOrder', colors.colorOrder);
    plot(qDVACalc, dvaSmoothCalc, '-', 'Color', colors.tumBlue, 'linewidth', lineSize);
    hold on;
    plot(qDVAMeas, dvaSmoothMeas, ':', 'Color', colors.darkGray, 'linewidth', lineSize);
    plot(qDVACathode, dvaCathodeSmooth, '-', 'Color', colors.mediumGray, 'linewidth', lineSize);
    plot(qDVAAnode, abs(dvaAnodeSmooth), '-', 'Color', colors.orange, 'linewidth', lineSize);
    % Determine overall limits from both datasets, then pad by 0.1
    % Robust one-liner using nested mins/maxs
    xlim([min(min(qDVACathode), min(qDVAAnode)) - 0.1, max(max(qDVACathode), max(qDVAAnode)) + 0.1]);
    ylim([-0.10, 3]);
    ylabel("$\frac{\Delta U}{\Delta Q}$ / $\frac{\mathrm{V}}{\mathrm{Ah}}$", 'fontSize', fontSize, 'FontName', fontName, 'Interpreter', 'latex');
    xlabel("Full cell SOC / -", 'fontSize', fontSize, 'FontName', fontName, 'Interpreter', 'latex');
    set(gca, 'box', 'off', 'TickDir', 'out', 'linewidth', lineSize, 'fontSize', fontSize, 'FontName', fontName, 'TickLabelInterpreter', 'latex', 'Position', ax.Position + fig_offset + [0 -0.02 0 0]);


end
