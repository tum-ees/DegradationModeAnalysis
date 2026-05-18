%> regenerate_param_figures
%>
%> Regenerate per-CU parameter figures (under input_data/test_data/results/
%> parameter/) with fixed x-axis bounds, 20% wider figure, and higher-
%> resolution PNG output. Reads saved per-CU data from
%> results/P45B_parameters.mat - no need to rerun the fit.
%>
%> Run from the project root.

%% Setup
addpath(fullfile('.', 'helper_functions'));

%% Load saved data
matFile = fullfile('input_data', 'test_data', 'results', 'P45B_parameters.mat');
load(matFile, 'data');

cuIdx         = 1:9;                  % CU indices in the data struct
cellName      = data.s.cellName;
labelCfg      = data.s;
savePath      = fullfile('input_data', 'test_data', 'results', 'parameter');
resolutionDPI = 600;                  % higher than default 300 in save_figure
targetXLim    = [-0.32, 1.19];        % fixed bounds
widthScale    = 0.88;                 % previous 1.1 reduced by 20%
heightScale   = 0.832;                % previous 0.64 increased by 30%
legendCols    = 2;                    % wrap legend onto two columns

%% Helper: build a solution-shaped struct from data.CU<i>
buildSolution = @(cu) struct( ...
    'solverInput',   struct('qCell',  cu.measured.soc, ...
                            'ocvCell', cu.measured.voltage), ...
    'reconSOC',      cu.calculated.soc, ...
    'fullCellUModel',cu.calculated.voltage, ...
    'qDVAMeas',      cu.measured.qDVA, ...
    'dvaSmoothMeas', cu.measured.dva, ...
    'qDVACalc',      cu.calculated.qDVA, ...
    'dvaSmoothCalc', cu.calculated.dva, ...
    'params',        cu.params, ...
    'cathodeSOC',    cu.calculated.cathodeSOC, ...
    'cathodeURecon', cu.calculated.cathodeURecon, ...
    'anodeSOC',      cu.calculated.anodeSOC, ...
    'anodeURecon',   cu.calculated.anodeURecon ...
);

%% Render every CU, apply fixed xlim, widen, and save high-res
for i = cuIdx
    plot_ocv_model_param_show(buildSolution(data.(sprintf('CU%d', i))), ...
                              i, cellName, labelCfg);

    fig = gcf;

    % Scale figure width and height (centimeters), keep origin
    set(fig, 'Units', 'centimeters');
    pos    = get(fig, 'Position');
    pos(3) = pos(3) * widthScale;
    pos(4) = pos(4) * heightScale;
    set(fig, 'Position', pos);

    % Apply fixed xlim to both tile axes
    axs = findobj(fig, 'Type', 'axes');
    for k = 1:numel(axs)
        set(axs(k), 'XLim', targetXLim);
    end

    % Re-wrap legend onto multiple columns
    lgd = findobj(fig, 'Type', 'Legend');
    if ~isempty(lgd)
        set(lgd, 'NumColumns', legendCols, 'Orientation', 'vertical');
    end

    drawnow;                          % let listeners reposition arrow annotations

    label       = sprintf('%s_CU%d', cellName, i);
    fileNamePDF = fullfile(savePath, sprintf('%s.pdf', label));
    fileNameFIG = fullfile(savePath, sprintf('%s.fig', label));
    fileNamePNG = fullfile(savePath, sprintf('%s.png', label));
    exportgraphics(fig, fileNamePDF, 'ContentType', 'vector');
    savefig(fig, fileNameFIG);
    exportgraphics(fig, fileNamePNG, 'Resolution', resolutionDPI);
    close(fig);
end

fprintf('Regenerated %d figures: xlim = [%.3f, %.3f], width x%.2f, height x%.2f, legend cols %d, %d dpi.\n', ...
        numel(cuIdx), targetXLim(1), targetXLim(2), widthScale, heightScale, legendCols, resolutionDPI);
