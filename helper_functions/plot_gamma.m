function plot_gamma(data, EFC, calendarOrCyclic, varargin)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%
% Input Arguments:
%   - varargin: Default material is Si. Can be replaced using a name-value pair.
%               Key: 'material'
%               Example: plot_gamma(data, EFC, calendarOrCyclic, 'material', 'Graphite')
%
% INPUTS:
%   data              - DMA result struct with fitted parameter vectors
%   EFC               - x-axis labels for each CU
%   calendarOrCyclic  - label mode selector (RPT number or EFC)
%   varargin          - optional name-value pairs, currently `material`
%
% OUTPUTS:
%   none              - creates a figure showing the blend-2 fraction trend

    fontSize = 14;
    lineSize = 1;
    markerSize = 8;
    figureSize = [10 5 20 15];
    fontName = 'Times New Roman';
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

    colors = tum_colors();

    % --- Input parser for optional arguments ---
    p = inputParser;
    addParameter(p, 'material', 'Si', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    material = p.Results.material;

    % Extract gamma values
    gamma = zeros(1, numel(EFC));
    for i = 1:numel(EFC)
        fieldName = sprintf('CU%d', i);
        gamma(i) = data.(fieldName).params(5) * 100;
    end

    % Plot
    fig = figure('Units', 'centimeters', 'Position', figureSize);
    ax = axes(fig);
    set(ax, 'ColorOrder', colors.colorOrder);
    plot(ax, EFC, gamma, '^', 'MarkerFaceColor', colors.tumBlue, 'Color', colors.tumBlue, ...
        'LineWidth', lineSize, 'LineStyle', '--', 'MarkerSize', markerSize);
    grid(ax, 'on')
    if calendarOrCyclic == 0    % Calendar aging
        xlabel(ax, 'RPT Number / -', 'Interpreter', 'latex', 'FontSize', fontSize, 'FontName', fontName);
    else                        % Cyclic aging
        xlabel(ax, 'EFC', 'Interpreter', 'latex', 'FontSize', fontSize, 'FontName', fontName);
    end
    ylabel(ax, ['$\gamma_{', material, '}\, / \,$ \%'], ...
        'FontSize', fontSize * 1.2, 'Interpreter', 'latex', 'FontName', fontName);

end
