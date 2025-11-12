function [dataOut, data, varargout] = import_OCV(Path, nameTableColumnOCV, inputIsAging_data_table, varargin)
%> Authors : Mathias Rehm
%> E-mail : mathias.rehm@tum.de
%> Date    : 2025-09-15
%
% import_OCV
% Import OCV/OCP curve data for DMA_main. Supports:
%   • Parent dir with CU* subfolders  (pass 'CheckUp', N; optional 'FileNameFilter')
%   • Direct CUi folder               (same as above)
%   • Direct .mat file                (table "aging_data_table" or variables TestData[/TestInfo])
%
% Row selection for tables (aging_data_table):
%   • Optional 'Identifier'={ 'Col',Val; ... } AND-filter
%   • Then pick requested CheckUp via CU/CheckUp column; else row index fallback
%
% Payload extraction from selected row:
%   1) If nameTableColumnOCV is given and points to a table/timetable → wrap as .TestData
%   2) Else first table/timetable column → wrap as .TestData
%   3) Else first struct with .TestData (table/timetable) → use that
%   4) Else treat as flat table with SOC/U/MaxAhStep → normalize to struct
%
% varargout{1} (second output):
%   'ReturnColumn' (if provided & exists) → else EFC → RPT → CU → NaN/index fallback
%
% Notes:
%   • When wrapping tabular data, we also try to attach a TestInfo if a matching column
%     (starts with 'Testinfo') exists in the table, or if the MAT file provides TestInfo.
%   • SOC is computed and appended if missing.

%% 1) Parse NV-args
p = inputParser;
addParameter(p, 'Identifier', {});       % 1x2 or Nx2 {Col,Val}
addParameter(p, 'CheckUp', []);          % numeric
addParameter(p, 'FileNameFilter', '');   % substring to match a MAT in CUi
addParameter(p, 'ReturnColumn', '');     % optional override for second output
addParameter(p, 'agingDataTable', []);   % <-- map your injected NV-pair here
parse(p, varargin{:});

identifier     = normalizeIdentifier(p.Results.Identifier);
checkUpNum     = p.Results.CheckUp;
fileNameFilter = p.Results.FileNameFilter;
ReturnColName  = p.Results.ReturnColumn;
agingDataTable = p.Results.agingDataTable;

%% 2) Resolve Path kind (parent with CU*, direct CUi, or file)
% We accept all three to stay compatible with existing callers (incl. your current get_fullcell_data).
if isfolder(Path)
    [parentDir, leaf] = fileparts(Path);
    hasCUHere   = ~isempty(dir(fullfile(Path, 'CU*')));                       % parent with CU*
    isCUfolder  = ~isempty(regexp(leaf, '^CU\d+$', 'once'));                  % direct CUi

    if hasCUHere
        % 2a) Parent with CU* → dive into requested CUi and pick file
        requireCheckUp(checkUpNum, 'Path points to CU* structure but no "CheckUp" given.');
        targetFolder = fullfile(Path, sprintf('CU%d', checkUpNum));
        Path = pickMatInCU(targetFolder, fileNameFilter);

    elseif isCUfolder
        % 2b) Direct CUi folder → pick file
        Path = pickMatInCU(Path, fileNameFilter);

    else
        % 2c) Plain directory (no CU*). Choose a single MAT (optionally by filter).
        Path = pickMatInFolder(Path, fileNameFilter);
    end
end

%% 3) Load file
if ~isfile(Path), error('File "%s" not found.', Path); end
if isempty(agingDataTable)
    loaded   = load(Path);
    flagLoaded = false;
else
    loaded = agingDataTable;
    flagLoaded = true;
end
varNames = fieldnames(loaded);

% 3a) aging_data_table mode: take the first variable as the table
if inputIsAging_data_table
    if ~flagLoaded
        if numel(varNames) > 1
            warning('inputIsAging_data_table=1 but multiple vars in MAT. Using first: %s', varNames{1});
        end
        data = loaded.(varNames{1});
    else
        data = loaded;
    end

% 3b) Direct TestData[/TestInfo] files (common in your workflow)
elseif numel(varNames) == 2 && all(ismember({'TestData','TestInfo'}, varNames))
    dataOut = struct('TestData', loaded.TestData, 'TestInfo', loaded.TestInfo);
    dataOut = ensureSOC(dataOut, Path);
    varargout{1} = NaN;  % No table row => no EFC/RPT/CU to read here
    return;

elseif numel(varNames) == 1
    data = loaded.(varNames{1});

else
    error('MAT must contain either aging_data_table or [TestData (+TestInfo)]. Found: %s', strjoin(varNames, ', '));
end

%% 4) Interpret content
% 4a) If the loaded variable is already a struct with .TestData (rare but allowed)
if isstruct(data) && isfield(data, 'TestData') && istabular_local(data.TestData)
    dataOut = ensureSOC(data, Path);
    % Second output not well-defined here; try to read EFC from TestInfo if present
    if isfield(dataOut, 'TestInfo') && isfield(dataOut.TestInfo, 'EFC')
        varargout{1} = dataOut.TestInfo.EFC;
    else
        varargout{1} = NaN;
    end
    return;
end

% 4b) Main path: aging_data_table (table)
if ~istable(data)
    error('Unsupported content: expected a table (aging_data_table) or struct with TestData.');
end

selected = selectRowFromTable(data, identifier, checkUpNum);   % resolves CheckUp row
if isempty(selected)
    dataOut = {};
    varargout{1} = NaN;
    return;
end
[dataOut, vv] = extractPayloadFromRow(selected, nameTableColumnOCV, Path);  % your highlighted block, simplified & fixed

% 4c) Attach TestInfo if not set by extract step and present in table columns
if ~isfield(dataOut, 'TestInfo')
    try
        vnames = selected.Properties.VariableNames;
        idxTI  = find(startsWith(vnames, 'Testinfo', 'IgnoreCase', true), 1, 'first');
        if ~isempty(idxTI)
            dataOut.TestInfo = selected{1, idxTI};
        end
    catch
        % optional, ignore if incompatible
    end
end

% 4d) Second output selection (ReturnColumn → EFC → RPT → CU → fallback)
if ~isempty(ReturnColName) && ismember(ReturnColName, selected.Properties.VariableNames)
    varargout{1} = selected.(ReturnColName)(1);
else
    varargout{1} = chooseSecondOutput(selected, checkUpNum, ~isempty(ReturnColName));
end

% In case extractPayloadFromRow already found a preferred value (vv), prefer that only if ReturnColumn requested it
if ~isempty(vv) && ~isempty(ReturnColName)
    varargout{1} = vv;
end

end
%% ========================= Helper functions =========================

function identifier = normalizeIdentifier(identifier)
% Normalize Identifier to {} or Nx2 cell of {Key,Value} (char)
if isstruct(identifier)
    fn = fieldnames(identifier);
    if isempty(fn)
        identifier = {};
    else
        vals = struct2cell(identifier);
        identifier = [cellfun(@char, fn(:),  'UniformOutput', false), ...
                      cellfun(@char, vals(:),'UniformOutput', false)];
    end
elseif iscell(identifier)
    if isempty(identifier)
        % ok
    elseif isvector(identifier) && numel(identifier) == 2
        identifier = {char(identifier{1}), char(identifier{2})};
    elseif size(identifier,2) == 2
        identifier(:,1) = cellfun(@char, identifier(:,1), 'UniformOutput', false);
        identifier(:,2) = cellfun(@char, identifier(:,2), 'UniformOutput', false);
    else
        error('Identifier must be 1x2 or Nx2 cell array of {Key,Value} pairs, or a struct.');
    end
else
    if ~isempty(identifier)
        error('Identifier must be a cell array or struct.');
    end
end
end

function requireCheckUp(n, msg)
if isempty(n) || ~isscalar(n) || ~isfinite(n)
    error(msg);
end
end

function fpath = pickMatInCU(cuFolder, fileNameFilter)
% Pick a .mat inside a specific CUi folder (optionally filtered).
if ~isfolder(cuFolder), error('CheckUp folder "%s" not found.', cuFolder); end
matFiles = dir(fullfile(cuFolder, '*.mat'));
if isempty(matFiles), error('No MAT in "%s".', cuFolder); end
if ~isempty(fileNameFilter)
    matFiles = matFiles(contains({matFiles.name}, fileNameFilter));
    if isempty(matFiles)
        error('No MAT in "%s" matched filter "%s".', cuFolder, fileNameFilter);
    end
end
fpath = fullfile(cuFolder, matFiles(1).name);
end

function fpath = pickMatInFolder(folder, fileNameFilter)
% Choose a single MAT in a plain folder (no CU*).
files = dir(fullfile(folder, '*.mat'));
if isempty(files)
    error('No MAT-file found in folder "%s".', folder);
end
if numel(files) == 1 && isempty(fileNameFilter)
    fpath = fullfile(folder, files(1).name);
    return
end
if isempty(fileNameFilter)
    error('Multiple MAT-files in "%s". Please provide FileNameFilter.', folder);
end
filtered = files(contains({files.name}, fileNameFilter));
if isempty(filtered)
    error('No MAT-file in "%s" matched filter "%s".', folder, fileNameFilter);
elseif numel(filtered) > 1
    error('Multiple MAT-files in "%s" matched filter "%s". Refine the filter.', folder, fileNameFilter);
end
fpath = fullfile(folder, filtered(1).name);
end

function selected = selectRowFromTable(T, identifier, checkUpNum)
% Apply Identifier AND-filter (if given), then select the row for the requested CheckUp.
if ~isempty(identifier)
    mask = true(height(T),1);
    for r = 1:size(identifier,1)
        col = identifier{r,1}; val = identifier{r,2};
        if ~ismember(col, T.Properties.VariableNames)
            error('Identifier column "%s" not found.', col);
        end
        colData = T.(col);
        if iscellstr(colData) || isstring(colData) || ischar(colData) || iscategorical(colData)
            mask = mask & (string(colData) == string(val));
        elseif isnumeric(colData) || islogical(colData)
            if ischar(val) || isstring(val)
                v = str2double(val);
                if isnan(v), error('Identifier value for numeric column "%s" must be numeric(-string).', col); end
                val = v;
            end
            mask = mask & (colData == val);
        else
            error('Column "%s" has unsupported type %s.', col, class(colData));
        end
    end
    filtered = T(mask,:);
    if isempty(filtered)
        error('No rows matched the provided Identifier pairs (AND-combined).');
    end
else
    filtered = T;
end

if height(filtered) == 1 && isempty(checkUpNum)
    selected = filtered(1,:);
    return
end

% Try CU/CheckUp columns first
possibleCUCols = {'CU','CUs','CheckUp','CheckUps','CheckUpNumber','CUNumber','Check_Up'};
cuCol = '';
for c = possibleCUCols
    if ismember(c{1}, filtered.Properties.VariableNames), cuCol = c{1}; break; end
end

if ~isempty(cuCol)
    requireCheckUp(checkUpNum, 'Table has CU/CheckUp column but no "CheckUp" given.');
    rowIdx = find(filtered.(cuCol) == checkUpNum, 1, 'first');
    if isempty(rowIdx)
        warning('Requested CheckUp %d not found in column "%s".', checkUpNum, cuCol);
        selected = {};
    else
        selected = filtered(rowIdx,:);
    end
    
else
    requireCheckUp(checkUpNum, 'No CU column present; using row index requires "CheckUp".');
    if checkUpNum > height(filtered)
        warning('Requested CheckUp %d exceeds number of rows (%d).', checkUpNum, height(filtered));
        selected = {};
    else
        selected = filtered(checkUpNum,:);
    end
end
end

function [outStruct, retOverride] = extractPayloadFromRow(selected, nameTableColumnOCV, Path)
% Implement the functionality you highlighted, but simpler and correct.
% Strategy:
%   1) If nameTableColumnOCV exists and is tabular → use it.
%   2) Else first tabular column → use it.
%   3) Else first struct with .TestData(tabular) → use it.
%   4) Else treat selected as a "flat" row → normalize to struct with SOC/U/MaxAhStep.
%
% Returns:
%   outStruct: struct with at least .TestData (for 1-3) or normalized fields (4)
%   retOverride: [] always here (kept for compatibility with an optional future override)

retOverride = [];

rowData   = table2cell(selected(1,:));
varNamesS = selected.Properties.VariableNames;

% --- Prefer explicit column, if provided ---
candidates = 1:numel(rowData);
if ~isempty(nameTableColumnOCV)
    idx = find(strcmp(varNamesS, nameTableColumnOCV), 1);
    if ~isempty(idx)
        candidates = idx;
    else
        warning('nameTableColumnOCV "%s" not found. Falling back to auto-detect.', nameTableColumnOCV);
    end
end

% 1) tabular column?
idxTab = [];
for k = candidates
    if istabular_local(rowData{k}), idxTab = k; break; end
end
if isempty(idxTab)
    % 2) any tabular column in the row?
    for k = 1:numel(rowData)
        if istabular_local(rowData{k}), idxTab = k; break; end
    end
end
if ~isempty(idxTab)
    outStruct = struct('TestData', rowData{idxTab});
    outStruct = ensureSOC(outStruct, Path);
    return
end

% 3) struct with .TestData (table/timetable)?
idxStruct = [];
for k = candidates
    v = rowData{k};
    if isstruct(v) && isfield(v,'TestData') && istabular_local(v.TestData)
        idxStruct = k; break;
    end
end
if isempty(idxStruct)
    for k = 1:numel(rowData)
        v = rowData{k};
        if isstruct(v) && isfield(v,'TestData') && istabular_local(v.TestData)
            idxStruct = k; break;
        end
    end
end
if ~isempty(idxStruct)
    outStruct = ensureSOC(rowData{idxStruct}, Path);
    return
end

% 4) No tabular/struct payload → interpret as flat row with SOC/U/MaxAhStep
outStruct = normalizeTableOrStruct(selected);
end

function v = chooseSecondOutput(selected, checkUpNum, warnedForReturnCol)
% ReturnColumn → EFC → RPT → CU → fallback
E = {'EFC','EffCap','EffectiveCapacity'};
R = {'RPT','ReProTest','ReferencePerformanceTest'};
C = {'CU','CUs','CheckUp','CheckUps','CheckUpNumber','CUNumber','Check_Up'};

if istable(selected)
    x = firstHit(selected, E);
    if ~isempty(x), v = x; return; end

    x = firstHit(selected, R);
    if ~isempty(x)
        v = x;
        if ~warnedForReturnCol, warning('EFC column not found. Returned RPT instead.'); end
        return
    end

    x = firstHit(selected, C);
    if ~isempty(x)
        v = x;
        if ~warnedForReturnCol, warning('EFC or RPT column not found. Returned CU instead.'); end
        return
    end

    if height(selected) == 1
        v = NaN;
        if ~warnedForReturnCol, warning('No EFC/RPT/CU column found; returned NaN.'); end
    else
        v = checkUpNum;
        if ~warnedForReturnCol, warning('No EFC/RPT column found; returned CU (by index).'); end
    end
else
    v = NaN;
end
end

function x = firstHit(T, candidates)
x = [];
hit = intersect(candidates, T.Properties.VariableNames, 'stable');
if ~isempty(hit), x = T.(hit{1})(1); end
end

function S = ensureSOC(S, Path)
% Append SOC to S.TestData if missing. Handles table/timetable.
if ~isfield(S,'TestData') || ~istabular_local(S.TestData), return; end
if ismember('SOC', S.TestData.Properties.VariableNames), return; end
try
    AhStep = S.TestData.Ah_Step;
    U      = S.TestData.U;
    S.TestData.SOC = getSOC(AhStep, U);
catch ME
    warning('Could not calculate SOC for file "%s": %s', Path, ME.message);
end
end

function outStruct = normalizeTableOrStruct(data)
% Create a minimal uniform struct from a flat table (or struct) with SOC/U/MaxAhStep.
valid.SOC       = {'SOC','StateOfCharge','State_of_Charge'};
valid.U         = {'U','V','U_','V_','E_','Potential','OCP','OCV','Voltage'};
valid.MaxAhStep = {'MaxAhStep','AhStep','max_ahstepsize','CapacityStep','Capacity_Step'};

if isstruct(data), names = fieldnames(data); else, names = data.Properties.VariableNames; end
outStruct = struct();
for f = fieldnames(valid)'
    key = f{1}; cands = valid.(key); match = [];
    for c = cands
        cand = c{1};
        if numel(cand) == 1
            hit = find(strcmpi(names, cand), 1, 'first');
        else
            hit = find(startsWith(names, cand, 'IgnoreCase', true), 1, 'first');
        end
        if ~isempty(hit), match = hit; break; end
    end
    if isempty(match)
        if istable(data), tdesc = 'table column'; else, tdesc = 'struct field'; end
        error('No %s for "%s". Allowed: %s', tdesc, key, strjoin(cands, ', '));
    end
    outStruct.(key) = data.(names{match});
    if strcmp(key,'MaxAhStep')
        vals = outStruct.MaxAhStep; vals = vals(~isnan(vals));
        outStruct.MaxAhStep = tern(~isempty(vals), max(vals), NaN);
    end
end
end

function y = tern(cond, a, b)
if cond, y = a; else, y = b; end
end

function tf = istabular_local(x)
tf = istable(x) || istimetable(x);
end

function SOC = getSOC(AhStep, U, varargin)
% Works only for curves in one direction.
[minAh, maxAh] = bounds(abs(AhStep));
if mean(AhStep) < 0, AhStep = AhStep + maxAh; end
SOC = nan(length(AhStep),1);
for i = 1:length(AhStep)
    SOC(i) = (AhStep(i) - minAh) / maxAh;
end
plaus = mean(diff(SOC)) / mean(diff(U));
if plaus < 0 && (max(U)-min(U) > 0.5)
    warning('SOC is wrong. With increasing voltage SOC is decreasing or vice versa!');
end
end
