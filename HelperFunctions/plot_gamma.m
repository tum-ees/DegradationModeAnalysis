function [] = plot_gamma(Data, EFC, CalendarOrCyclic, varargin)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%   
% Input Arguments:
%   - varargin: Default material is Si. Can be replaced using a Key-Value pair.
%               Key: 'Material'
%               Example: plot_gamma(Data, EFC, CalendarOrCyclic, 'Material', 'Graphite')

    fontSize = 14;
    lineSize    = 1;
    markerSize  = 8;
    figureSize = [10 5 20 15];
    fontName    = 'Times New Roman';
    set(groot, 'defaultAxesTickLabelInterpreter','latex');

    TUM_colors;

    % --- Input Parser für optionale Argumente ---
    p = inputParser;
    addParameter(p, 'Material', 'Si', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    material = p.Results.Material;

    % Extract Gamma values
    for i = 1:numel(EFC)
        fieldName = sprintf('CU%d', i);
        gamma(i) = Data.(fieldName).params(5) * 100;
    end

    % Plot
    figure('Units', 'centimeters', 'Position', figureSize);
    plot(EFC, gamma, '^', 'MarkerFaceColor', tumBlue, 'Color', tumBlue, 'LineWidth', lineSize, 'LineStyle', '--', 'MarkerSize', markerSize);
    grid on
    if CalendarOrCyclic == 0    % Calendar Aging
        xlabel('RPT Number / -', 'Interpreter','latex', 'FontSize', fontSize, 'FontName', fontName);
   else                         % Cyclic Aging
        xlabel('EFC', 'Interpreter','latex', 'FontSize', fontSize, 'FontName', fontName);
   end
    ylabel(['$\gamma_{', material, '}\, / \,$ \%'],  'FontSize', fontSize * 1.2, 'Interpreter', 'latex', 'FontName', fontName);

end
