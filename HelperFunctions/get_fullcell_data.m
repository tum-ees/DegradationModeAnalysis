function [fullCellData, toggleMissingData, fullTableOut] = get_fullcell_data(s, currentPath, ...
    cellIdentifier_keys, cellIdentifier_values, nameTableColumnOCV, i, varargin)
%> Authors : Mathias Rehm
%> E-mail : mathias.rehm@tum.de
%> Date    : 9/12/2025
%
% get_fullcell_data
% Fetch full-cell data for one CU/RPT using filters from main_DMA.
%
%> -------------------------------------------------------------------------
%>  * Original code by: Sebastian Karl (s.karl@tum.de)
%>    Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> -------------------------------------------------------------------------

%% 1) Setup and basic checks
basePath = regexprep(char(s.pathAgingStudy), '[\\/]', filesep);
toggleMissingData = false;
cuLabel = sprintf('CU%d', i);

% extract variables
inputIsAging_data_table = 1; % always 1; main_DMA only supports table mode
% -> future release may support folder mode as well

% Build Identifier argument (accept {}, 1xN or Nx1 cells as inputs paired to Nx2)
identifierArg = {};
if ~isempty(cellIdentifier_keys) || ~isempty(cellIdentifier_values)
    if isempty(cellIdentifier_keys) ~= isempty(cellIdentifier_values)
        error('TableFilter must include both Keys and Values, or neither.');
    end
    if ~iscell(cellIdentifier_keys) || ~iscell(cellIdentifier_values)
        error('CellIdentifier_Keys/Values must be cell arrays.');
    end
    if numel(cellIdentifier_keys) ~= numel(cellIdentifier_values)
        error('Keys and Values must have the same number of elements.');
    end
    K = cellfun(@char, cellIdentifier_keys(:),   'UniformOutput', false);
    V = cellfun(@char, cellIdentifier_values(:), 'UniformOutput', false);
    idPairs = [K V];  % Nx2 cell
    identifierArg = {'Identifier', idPairs};
end

%% --- Combine identifierArg with varargin ---
if numel(varargin)==1 && iscell(varargin{1}) && isempty(varargin{1})
    extraArgs = [identifierArg];
else
    extraArgs = [identifierArg, varargin];
end


% If CU folder missing, fall back to basePath
if ~exist(currentPath, 'dir')
    % fprintf('Folder "%s" not found. Using basePath instead.\n', cuLabel);
    currentPath = basePath;
else
    % If CU folders exist, basePath must NOT contain stray MAT files
    basePathFiles = dir(fullfile(basePath, '*.mat'));
    if ~isempty(basePathFiles)
        error(['Inconsistent data structure: pathAgingStudy contains both ' ...
            'CU folders and MAT-file(s). Please keep either CU subfolders ' ...
            'or a single MAT in the root, not both.']);
    end
end

%% 2) Discover files in currentPath
if isfolder(currentPath)
    files = dir(fullfile(currentPath, '*.mat'));
% case user gives aging_data_table as input
elseif isfile(currentPath)
    files{1} = currentPath;
    inputIsSingleFile = 1;
end

if isempty(files)
    if strcmp(currentPath, basePath)
        error('No MAT-file found in pathAgingStudy. Cannot continue.');
    else
        fprintf('No MAT-file found for %s. Skipping.\n', cuLabel);
        toggleMissingData = true;
        fullCellData = [];
        return;
    end
end

%% 3) BasePath mode (no CU folders): select file and import
% basePath and currentPath equal -> only one file as input (=
% agingDataTable)
if strcmp(currentPath, basePath)
    if numel(files) == 1
        if inputIsSingleFile
            matFilePath = currentPath;
        else
            matFilePath = fullfile(currentPath, files(1).name);
        end
    else
        error(['Multiple MAT-files found in pathAgingStudy. Currently, only one MAT-file (table) is supported. ' ...
        'Wait for future releases or add the missing functionality yourself.']);
        % fileNameFilter = s.fileNameFilter;
        % if isempty(fileNameFilter)
        %     error('Multiple MAT-files found in basePath. Please specify a FileNameFilter.');
        % end
        % filteredFiles = files(contains({files.name}, fileNameFilter));
        % if isempty(filteredFiles)
        %     error('No MAT-file matched the FileNameFilter in basePath.');
        % elseif numel(filteredFiles) > 1
        %     error('Multiple MAT-Files matched the FileNameFilter in basePath. Please refine the filter.');
        % end
        % matFilePath = fullfile(currentPath, filteredFiles(1).name);
    end

    % load data
    if isempty(identifierArg)
        [fullCellData, fullTableOut.agingDataTable, ~] = import_OCV(matFilePath, nameTableColumnOCV, inputIsAging_data_table, 'CheckUp', i, varargin);
    else
        [fullCellData, fullTableOut.agingDataTable, ~] = import_OCV(matFilePath, nameTableColumnOCV, inputIsAging_data_table, 'CheckUp', i, extraArgs{:});
    end

%% 4) CU-folder mode: import from the CU directory
else
    if numel(files) == 1
        if isempty(identifierArg)
            [fullCellData, fullTableOut.agingDataTable, ~] = import_OCV(currentPath, 'CheckUp', i, varargin);
        else
            [fullCellData, fullTableOut.agingDataTable, ~] = import_OCV(currentPath, 'CheckUp', i, extraArgs{:});
        end
    else
        if isempty(identifierArg)
            [fullCellData, fullTableOut.agingDataTable, ~] = import_OCV(currentPath, 'CheckUp', i, 'FileNameFilter', fileNameFilter, varargin);
        else
            [fullCellData, fullTableOut.agingDataTable, ~] = import_OCV(currentPath, 'CheckUp', i, ...
                identifierArg{:}, 'FileNameFilter', fileNameFilter, extraArgs{:});
        end
    end
end

%% 5) Final sanity
if isempty(fullCellData)
    fprintf('No full-cell data returned for %s. Skipping.\n', cuLabel);
    toggleMissingData = true;
end
end
