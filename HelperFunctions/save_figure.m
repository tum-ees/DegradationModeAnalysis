function save_figure(outputFolder, fieldName, varargin)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%
% A helper function that helps to save the created figures in main_DMA
%
    % saveFigure(outputFolder, fieldName, 'savePNG', value, ...)

    % default value
    savePNG = false;

    % Evaluate options
    for v = 1:2:numel(varargin)
        key = lower(string(varargin{v}));
        val = varargin{v+1};
        switch key
            case {'savepng', 'save_png'}
                savePNG = logical(val);
            otherwise
                warning('saveFigure: unknown option "%s" – ignored.', varargin{v});
        end
    end

    % create directory if necessary
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    fileNamePDF = fullfile(outputFolder, sprintf('%s.pdf', fieldName));
    fileNameFIG = fullfile(outputFolder, sprintf('%s.fig', fieldName));
    exportgraphics(gcf, fileNamePDF, 'ContentType', 'vector');
    savefig(gcf, fileNameFIG);

    if savePNG
        fileNamePNG = fullfile(outputFolder, sprintf('%s.png', fieldName));
        exportgraphics(gcf, fileNamePNG, 'Resolution', 300);
    end
end