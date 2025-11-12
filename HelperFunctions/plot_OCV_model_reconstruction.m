function plot_OCV_model_reconstruction(meas_fullCell_SOC, meas_fullCell_voltage, calc_fullCell_SOC, full_cell_voltage, Q_DVA_meas, DVA_smooth_meas, Q_DVA_calc, DVA_smooth_calc, params, ...
    cathSOC, cathU_recon, anodeSOC, anodeU_recon, finalRMSE, i, cellName)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-10-10
%
% plot_OCV_model_reconstruction
% Tiled OCV plus DVA plotter with dynamic limits and parameterized styling

    TUM_colors;

    fontSize = 12;
    lineSize = 2;
    titleSize = 16;
    figSize1 = [3 2 25 15];
    fontName    = 'Times New Roman';

    % Hard-coded switch to show RMSE in title
    showRMSE = true; % Set to true to display RMSE in the title

    rmse_mV = finalRMSE * 1000;
    set(groot, 'defaultAxesTickLabelInterpreter','latex');

    % Full Cell
    [~, OCV_ICA_fullcell_meas, ICA_fullcell_meas] = calculate_ICA(meas_fullCell_SOC, meas_fullCell_voltage);
    [~, OCV_ICA_fullcell_calc, ICA_fullcell_calc] = calculate_ICA(calc_fullCell_SOC, full_cell_voltage);
    ICA_smooth_meas = smooth(ICA_fullcell_meas,  30, 'lowess');
    ICA_smooth_calc = smooth(ICA_fullcell_calc,  30, 'lowess');

    % Cathode
    [Q_DVA_cath, ~, DVA_cath] = calculate_DVA(cathSOC, cathU_recon);
    DVA_cath_smooth = smooth(DVA_cath, 30, 'lowess');

    % Anode
    [Q_DVA_anode, ~, DVA_anode] = calculate_DVA(anodeSOC, anodeU_recon);
    DVA_anode_smooth = smooth(DVA_anode, 30, 'lowess');

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
    
    % add extra plots for legend
    h1 = plot(NaN, NaN, 'k-', 'DisplayName', 'model','linewidth', lineSize);
    hold on;
    h2 = plot(NaN, NaN, 'k:', 'DisplayName', 'measurement','linewidth', lineSize);
    hold on;

    plot(calc_fullCell_SOC, full_cell_voltage, '-', 'Color', tumBlue, 'linewidth', lineSize);
    hold on;
    plot(meas_fullCell_SOC, meas_fullCell_voltage, ':', 'Color', darkGray, 'linewidth', lineSize);
    plot(cathSOC, cathU_recon, '-', 'Color', mediumGray, 'linewidth', lineSize);
    maxVoltage = max([max(meas_fullCell_voltage), max(full_cell_voltage), max(cathU_recon)]) + 0.05;
    minVoltage = min([min(meas_fullCell_voltage), min(full_cell_voltage)]) - 0.05;
    xlim([min(cathSOC)-0.1, max(cathSOC)+0.1]);
    ylim([minVoltage, maxVoltage]);
    set(gca, 'XTickLabel', [], 'Position', ax.Position + fig_offset);
    
    [ax, h] = addaxis(meas_fullCell_SOC, (meas_fullCell_voltage - full_cell_voltage)*1000, [-20, 20], ':', 'Color', green, 'linewidth', lineSize*0.7);
    addaxislabel(2, 'Full Cell Fitting Error / mV', 'Fontname', fontName, 'Interpreter', 'latex');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    set(ax, 'box', 'off', 'TickDir', 'out', 'linewidth', lineSize, 'fontSize', fontSize, 'Fontname', fontName, 'Ticklabelinterpreter', 'latex', 'Position', ax.Position + [+0.06 0 0 0]);
    
    yyaxis right;
    plot(anodeSOC, anodeU_recon, '-', 'Color', orange, 'linewidth', lineSize);
    hold on;
    ylabel('Anode Potential vs Li/Li$^+$ / V', 'fontSize', fontSize, 'FontName', fontName, 'Interpreter', 'latex');
    set(gca, 'box', 'off', 'TickDir', 'out', 'YColor', orange, 'linewidth', lineSize, 'fontSize', fontSize, 'FontName', fontName, 'TickLabelInterpreter', 'latex');
    
    ax = subplot(4,4, [1 5 9]);
    plot(ICA_smooth_calc, OCV_ICA_fullcell_calc, '-', 'Color', tumBlue, 'linewidth', lineSize);
    hold on;
    maxICA = max([max(ICA_smooth_meas), max(OCV_ICA_fullcell_meas)])+5;
    plot(ICA_smooth_meas, OCV_ICA_fullcell_meas, ':', 'Color', darkGray, 'linewidth', lineSize);
    plot(1000, 1000, '-', 'Color', mediumGray, 'linewidth', lineSize);
    plot(1000, 1000, '-', 'Color', orange, 'linewidth', lineSize);
    ylim([minVoltage, maxVoltage]);
    xlim([0, maxICA]);
    ylabel('Full Cell Voltage / V', 'fontSize', fontSize, 'FontName', fontName, 'Interpreter', 'latex');
    xlabel("$\frac{\Delta Q}{\Delta U}$ / $\frac{\mathrm{Ah}}{\mathrm{V}}$", 'fontSize', fontSize, 'FontName', fontName, 'Interpreter', 'latex');
    legend([h1, h2], 'Position', [0.1 0.15 0.08 0.05], 'Units', 'normalized', 'FontSize', fontSize, 'Interpreter', 'latex');
    
    % JE to MRe: yticks ok? Do we need ylabel againg for full-cell voltage?
    % -> ICA scaling is different to improve graphics
    set(gca, 'box', 'off', 'TickDir', 'out', 'YColor', tumBlue, 'linewidth', lineSize, 'fontSize', fontSize, 'FontName', fontName, 'TickLabelInterpreter', 'latex', 'Position', ax.Position + fig_offset + [0.02 0 -0.02 0]);
    
    [ax, ~] = addaxis([], []);
    set(ax, 'visible', false);
    
    [ax, h] = addaxis(NaN, NaN, [2.3, 4.4], '--', 'Color', mediumGray, 'linewidth', lineSize*0.7);
    addaxislabel(3, 'Cathode Potential vs Li/Li^+ / V');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    set(ax, 'box', 'off', 'TickDir', 'out', 'linewidth', lineSize, 'fontSize', fontSize, 'Fontname', fontName, 'Ticklabelinterpreter', 'latex', 'Position', ax.Position + [-0.05 0 0 0]); %'Fontname', fontName,
     
    ax = subplot(4,4, [14 15 16]);
    plot(Q_DVA_calc, DVA_smooth_calc, '-', 'Color', tumBlue, 'linewidth', lineSize);
    hold on;
    plot(Q_DVA_meas, DVA_smooth_meas, ':', 'Color', darkGray, 'linewidth', lineSize);
    plot(Q_DVA_cath, DVA_cath_smooth, '-', 'Color', mediumGray, 'linewidth', lineSize);
    plot(Q_DVA_anode, abs(DVA_anode_smooth), '-', 'Color', orange, 'linewidth', lineSize);
    % Determine overall limits from both datasets, then pad by 0.1
    % Robust one‑liner using nested mins/maxs
    xlim([min(min(Q_DVA_cath), min(Q_DVA_anode)) - 0.1,  max(max(Q_DVA_cath), max(Q_DVA_anode)) + 0.1]);
    ylim([-0.10, 3]);
    ylabel("$\frac{\Delta U}{\Delta Q}$ / $\frac{\mathrm{V}}{\mathrm{Ah}}$", 'fontSize', fontSize, 'FontName', fontName, 'Interpreter', 'latex');
    xlabel("Full cell SOC / -", 'fontSize', fontSize, 'FontName', fontName, 'Interpreter', 'latex');
    set(gca, 'box', 'off', 'TickDir', 'out', 'linewidth', lineSize, 'fontSize', fontSize, 'FontName', fontName, 'TickLabelInterpreter', 'latex', 'Position', ax.Position + fig_offset + [0 -0.02 0 0]);


end
