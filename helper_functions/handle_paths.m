function [savePathModel, savePathDMA] = handle_paths(settings)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%
% Build save paths for DMA results.
%
% INPUTS:
%   settings              - DMA settings struct with `pathSaveResults`
%
% OUTPUTS:
%   savePathModel         - directory for per-CU model reconstruction plots
%   savePathDMA           - directory for overall DMA summary figures
%
    pathSave      = settings.pathSaveResults;
    savePathModel = fullfile(pathSave, 'parameter');
    savePathDMA   = pathSave;
end
