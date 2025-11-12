%> Author: Mathias Rehm
%> Email: mathias.rehm@tum.de
%> Date: 2025-09-12

%% parsing_main_DMA

% script to handle any input of main_DMA; main function is to clean up
% main_DMA

% %% logic to handle optional varargin of this function %%%
p = inputParser;

% in settings all user settings are stored
addParameter(p, 'Settings', []); % struct with predefined settings
parse(p, varargin{:});

% give back variables of p
tableFilterNew    = p.Results.TableFilter;
fileNameFilterNew = p.Results.FileNameFilter;
cellName =          p.Results.CellName;
settings =          p.Results.Settings;

% provide robust defaults for the table filter
% Default key: 'battery_serial'
% only overwrite the tableFilter if user has given an input here!
if ~isempty(tableFilterNew)
    s.tableFilter = tableFilterNew;
end

% if user did not provide a fileNameFilter use the new one
if ~isempty(fileNameFilterNew)
    fileNameFilter = fileNameFilterNew;
    % if no tableFilter was given i) as input of main_DMA and ii) no
    % tableFilter was defined in main_DMA
end