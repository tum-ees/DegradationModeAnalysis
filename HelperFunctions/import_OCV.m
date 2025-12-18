function [dataOut, data, varargout] = import_OCV(Path, nameTableColumnOCV, inputIsAging_data_table, varargin)
%> Authors : Mathias Rehm
%> E-mail : mathias.rehm@tum.de
%> Date    : 2025-12-03
%
% import_OCV
%
% Import OCV or OCP curve data for DMA_main
%
% Supported input types
%   • Parent directory containing CU* subfolders
%   • Direct CUi folder
%   • Direct MAT file with
%       - aging_data_table like table
%       - or variables TestData and TestInfo
%
% Row selection in aging_data_table
%   • Optional Identifier: AND filter on multiple columns
%   • Then select by CheckUp number using CU like columns
%   • If no CU like column exists, use CheckUp as row index
%
% Payload extraction from selected row
%   1) If nameTableColumnOCV points to a table or timetable
%   2) Else first table or timetable column in row
%   3) Else first struct with TestData that is table or timetable
%   4) Else treat row as flat SOC, U, MaxAhStep holder
%
% Outputs
%   dataOut      struct with TestData (and TestInfo if available)
%   data         raw table or struct used as input
%   varargout{1} tag or index
%                ReturnColumn (if given and found)
%                otherwise EFC, or RPT, or CU, or fallback

%% 1) Parse NV arguments
p = inputParser;
addParameter(p, 'Identifier', {});       % identifier pairs to pre filter rows
addParameter(p, 'CheckUp', []);          % numeric check up number
addParameter(p, 'FileNameFilter', '');   % optional substring filter for MAT name
addParameter(p, 'ReturnColumn', '');     % explicit column for third output
addParameter(p, 'agingDataTable', []);   % external aging_data_table
parse(p, varargin{:});

identifier     = normalizeIdentifier(p.Results.Identifier);
checkUpNum     = p.Results.CheckUp;
fileNameFilter = p.Results.FileNameFilter;
ReturnColName  = p.Results.ReturnColumn;
agingDataTable = p.Results.agingDataTable;

%% 2) Resolve Path type and decide which MAT file to use
% Path can be parent of CUi folders, a CUi folder itself, a plain folder or a MAT file

if isfolder(Path)
    [parentDir, leaf] = fileparts(Path); %#ok<ASGLU>

    % Does this directory contain CU* subfolders
    dirCU     = dir(fullfile(Path, 'CU*'));
    hasCUHere = ~isempty(dirCU);

    % Is this directory itself named like a CUi folder
    isCUfolder = ~isempty(regexp(leaf, '^CU\d+$', 'once'));

    if hasCUHere
        % Parent with CU* subfolders
        requireCheckUp(checkUpNum, 'Path points to CU* structure but no "CheckUp" given.');
        targetFolder = fullfile(Path, sprintf('CU%d', checkUpNum));
        Path = pickMatInCU(targetFolder, fileNameFilter);

    elseif isCUfolder
        % Direct CUi folder
        Path = pickMatInCU(Path, fileNameFilter);

    else
        % Plain folder
        Path = pickMatInFolder(Path, fileNameFilter);
    end
end

%% 3) Load MAT file or use injected table or struct

if ~isfile(Path) && isempty(agingDataTable)
    error('File "%s" not found.', Path);
end

if isempty(agingDataTable)
    loaded     = load(Path);
    flagLoaded = false;
else
    loaded     = agingDataTable;
    flagLoaded = true;
end

if isstruct(loaded)
    varNamesLoaded = fieldnames(loaded);
else
    varNamesLoaded = {};
end

% aging_data_table mode
if inputIsAging_data_table
    if ~flagLoaded
        if numel(varNamesLoaded) > 1
            warning('import_OCV:MultiVarMat', ...
                'inputIsAging_data_table=1 but MAT contains multiple variables. Using first: %s', varNamesLoaded{1});
        end
        if isempty(varNamesLoaded)
            error('inputIsAging_data_table=1 but loaded content is not a struct with variables.');
        end
        data = loaded.(varNamesLoaded{1});
    else
        % agingDataTable was passed directly as table
        data = loaded;
    end

% Direct TestData plus TestInfo
elseif ~isempty(varNamesLoaded) && numel(varNamesLoaded) == 2 && all(ismember({'TestData','TestInfo'}, varNamesLoaded))
    dataOut = struct('TestData', loaded.TestData, 'TestInfo', loaded.TestInfo);
    dataOut = ensureSOC(dataOut, Path);
    data    = loaded;
    varargout{1} = NaN;  % no EFC tag is defined here
    return;

% Single variable in MAT, assume it is aging_data_table like content
elseif ~isempty(varNamesLoaded) && numel(varNamesLoaded) == 1
    data = loaded.(varNamesLoaded{1});

else
    error('MAT must contain either aging_data_table or [TestData (+TestInfo)]. Found variables: %s', strjoin(varNamesLoaded, ', '));
end

%% 4) Interpret loaded content

% Case struct with TestData already present
if isstruct(data) && isfield(data, 'TestData') && istabular_local(data.TestData)
    dataOut       = ensureSOC(data, Path);
    data          = dataOut;  % second output
    if isfield(dataOut, 'TestInfo') && isfield(dataOut.TestInfo, 'EFC')
        varargout{1} = dataOut.TestInfo.EFC;
    else
        varargout{1} = NaN;
    end
    return;
end

% Main path, aging_data_table like table
if ~istable(data)
    error('Unsupported content: expected a table (aging_data_table) or struct with TestData.');
end

% First output will be derived from this selected row
selected = selectRowFromTable(data, identifier, checkUpNum);
if isempty(selected)
    dataOut      = {};
    data         = data;
    varargout{1} = NaN;
    return;
end

[dataOut, vv] = extractPayloadFromRow(selected, nameTableColumnOCV, Path);

% If we did not get TestInfo yet, try to attach from a column
if ~isfield(dataOut, 'TestInfo')
    try
        vnames = selected.Properties.VariableNames;
        idxTI  = find(startsWith(vnames, 'Testinfo', 'IgnoreCase', true), 1, 'first');
        if ~isempty(idxTI)
            dataOut.TestInfo = selected{1, idxTI};
        end
    catch
        % ignore if types are incompatible
    end
end

% Third output: priority ReturnColumn then EFC then RPT then CU then fallback
if ~isempty(ReturnColName) && ismember(ReturnColName, selected.Properties.VariableNames)
    varargout{1} = selected.(ReturnColName)(1);
else
    varargout{1} = chooseThirdOutput(selected, checkUpNum, ~isempty(ReturnColName));
end

% If extractPayloadFromRow ever returns its own preferred tag, allow override
if ~isempty(vv) && ~isempty(ReturnColName)
    varargout{1} = vv;
end

% Second output is the full aging table
data = data;

end

%% ========================= Helper functions =========================

function identifier = normalizeIdentifier(identifier)
% Normalize Identifier to {} or N by 2 cell array of {Key, Value}
% Keys are always chars, values are passed through

if isstruct(identifier)
    fn = fieldnames(identifier);
    if isempty(fn)
        identifier = {};
    else
        vals = struct2cell(identifier);
        n    = numel(fn);
        tmp  = cell(n,2);
        for k = 1:n
            tmp{k,1} = char(fn{k});
            tmp{k,2} = vals{k};
        end
        identifier = tmp;
    end

elseif iscell(identifier)
    if isempty(identifier)
        % ok, nothing to do
    elseif isvector(identifier) && numel(identifier) == 2
        identifier = {char(identifier{1}), identifier{2}};
    elseif size(identifier,2) == 2
        for k = 1:size(identifier,1)
            identifier{k,1} = char(identifier{k,1});
        end
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
% Check that a CheckUp number is a finite scalar

if isempty(n) || ~isscalar(n) || ~isfinite(n)
    error(msg);
end
end

function fpath = pickMatInCU(cuFolder, fileNameFilter)
% Pick a MAT file inside a specific CUi folder

if ~isfolder(cuFolder)
    error('CheckUp folder "%s" not found.', cuFolder);
end

matFiles = dir(fullfile(cuFolder, '*.mat'));
if isempty(matFiles)
    error('No MAT file in "%s".', cuFolder);
end

% sort by name for deterministic choice
[~, idxSort] = sort({matFiles.name});
matFiles = matFiles(idxSort);

if ~isempty(fileNameFilter)
    keep = false(size(matFiles));
    for k = 1:numel(matFiles)
        if contains(matFiles(k).name, fileNameFilter)
            keep(k) = true;
        end
    end
    matFiles = matFiles(keep);
    if isempty(matFiles)
        error('No MAT file in "%s" matched filter "%s".', cuFolder, fileNameFilter);
    end
end

fpath = fullfile(cuFolder, matFiles(1).name);
end

function fpath = pickMatInFolder(folder, fileNameFilter)
% Pick a MAT file in a plain folder without CU* structure

files = dir(fullfile(folder, '*.mat'));
if isempty(files)
    error('No MAT file found in folder "%s".', folder);
end

% deterministic order
[~, idxSort] = sort({files.name});
files = files(idxSort);

if numel(files) == 1 && isempty(fileNameFilter)
    fpath = fullfile(folder, files(1).name);
    return;
end

if isempty(fileNameFilter)
    error('Multiple MAT files in "%s". Please provide FileNameFilter.', folder);
end

keep = false(size(files));
for k = 1:numel(files)
    if contains(files(k).name, fileNameFilter)
        keep(k) = true;
    end
end
filtered = files(keep);

if isempty(filtered)
    error('No MAT file in "%s" matched filter "%s".', folder, fileNameFilter);
elseif numel(filtered) > 1
    error('Multiple MAT files in "%s" matched filter "%s". Refine the filter.', folder, fileNameFilter);
end

fpath = fullfile(folder, filtered(1).name);
end

function selected = selectRowFromTable(T, identifier, checkUpNum)
% Select one row from an aging_data_table
%
% Steps
%   1) apply Identifier pairs as AND filter
%   2) if multiple rows remain, use CU or CheckUp like column
%   3) if no CU like column exists, use CheckUp as row index
%
% CU like column can be numeric or string labels like "CU1" or "CU01"

% 1) apply Identifier filter if given
if ~isempty(identifier)
    mask = true(height(T),1);
    for r = 1:size(identifier,1)
        col = identifier{r,1};
        val = identifier{r,2};

        if ~ismember(col, T.Properties.VariableNames)
            error('Identifier column "%s" not found.', col);
        end

        colData = T.(col);

        % text like columns, including cell arrays of string or char
        if iscell(colData)
            % convert both sides to string for comparison
            colStr = string(colData);
            valStr = string(val);
            mask   = mask & (colStr == valStr);

        elseif isstring(colData) || ischar(colData) || iscategorical(colData)
            mask = mask & (string(colData) == string(val));

        % numeric or logical columns
        elseif isnumeric(colData) || islogical(colData)
            if ischar(val) || isstring(val) || iscategorical(val)
                v = str2double(string(val));
                if isnan(v)
                    error('Identifier value for numeric column "%s" must be numeric or numeric string.', col);
                end
                val = v;
            end
            mask = mask & (colData == val);

        else
            error('Identifier column "%s" has unsupported type %s.', col, class(colData));
        end
    end

    filtered = T(mask,:);
    if isempty(filtered)
        error('No rows matched the provided Identifier pairs.');
    end
else
    filtered = T;
end


% If there is exactly one row and no CheckUp was requested, return it
if height(filtered) == 1 && isempty(checkUpNum)
    selected = filtered(1,:);
    return;
end

% 2) look for CU like columns
possibleCUCols = {'CU','CUs','CheckUp','CheckUps','CheckUpNumber','CUNumber','Check_Up'};
cuCol = '';

for c = 1:numel(possibleCUCols)
    if ismember(possibleCUCols{c}, filtered.Properties.VariableNames)
        cuCol = possibleCUCols{c};
        break;
    end
end

if ~isempty(cuCol)
    % We need a CheckUp number in this case
    requireCheckUp(checkUpNum, 'Table has CU or CheckUp column but no "CheckUp" number given.');

    colData  = filtered.(cuCol);
    rowIdx   = [];

    if isnumeric(colData) || islogical(colData)
        % direct numeric compare
        rowIdx = find(colData == checkUpNum, 1, 'first');

    else
        % robust selection for string labels like "CU1", "CU01"
        labels = string(colData);
        % first try exact numeric string
        pattern1 = sprintf('^%d$', checkUpNum);
        % then patterns like CU1 or CU01 or CU_1
        pattern2 = sprintf('^CU[_ ]*0*%d$', checkUpNum);

        for i = 1:numel(labels)
            if ~isempty(regexp(labels(i), pattern1, 'once'))
                rowIdx = i;
                break;
            end
        end
        if isempty(rowIdx)
            for i = 1:numel(labels)
                if ~isempty(regexp(labels(i), pattern2, 'once'))
                    rowIdx = i;
                    break;
                end
            end
        end

        % final fallback: extract first integer from label and compare
        if isempty(rowIdx)
            for i = 1:numel(labels)
                digitsMatch = regexp(labels(i), '\d+', 'match', 'once');
                if ~isempty(digitsMatch)
                    valNum = str2double(digitsMatch);
                    if valNum == checkUpNum
                        rowIdx = i;
                        break;
                    end
                end
            end
        end
    end

    if isempty(rowIdx)
        warning('Requested CheckUp %d not found in column "%s".', checkUpNum, cuCol);
        selected = {};
    else
        selected = filtered(rowIdx,:);
    end

else
    % 3) no CU or CheckUp column, use CheckUp as row index
    requireCheckUp(checkUpNum, 'No CU or CheckUp column present; using row index requires "CheckUp".');
    if checkUpNum > height(filtered)
        warning('Requested CheckUp %d exceeds number of rows (%d).', checkUpNum, height(filtered));
        selected = {};
    else
        selected = filtered(checkUpNum,:);
    end
end
end

function [outStruct, retOverride] = extractPayloadFromRow(selected, nameTableColumnOCV, Path)
% Extract OCV payload from a single selected row
%
% Order
%   1) explicit OCV column if it exists and is table or timetable
%   2) first table or timetable column in row
%   3) first struct with TestData as table or timetable
%   4) interpret row as flat SOC or U or MaxAhStep container

retOverride = [];

rowData   = table2cell(selected(1,:));
varNamesS = selected.Properties.VariableNames;

% candidate columns to inspect first
candidates = 1:numel(rowData);
if ~isempty(nameTableColumnOCV)
    idx = find(strcmp(varNamesS, nameTableColumnOCV), 1);
    if ~isempty(idx)
        candidates = idx;
    else
        warning('import_OCV:MissingOCVColumn', ...
            'nameTableColumnOCV "%s" not found. Falling back to auto detect.', nameTableColumnOCV);
    end
end

% 1) explicit or first tabular column
idxTab = [];
for k = candidates
    if istabular_local(rowData{k})
        idxTab = k;
        break;
    end
end
if isempty(idxTab)
    for k = 1:numel(rowData)
        if istabular_local(rowData{k})
            idxTab = k;
            break;
        end
    end
end
if ~isempty(idxTab)
    outStruct = struct('TestData', rowData{idxTab});
    outStruct = ensureSOC(outStruct, Path);
    return;
end

% 2) struct with TestData field
idxStruct = [];
for k = candidates
    v = rowData{k};
    if isstruct(v) && isfield(v,'TestData') && istabular_local(v.TestData)
        idxStruct = k;
        break;
    end
end
if isempty(idxStruct)
    for k = 1:numel(rowData)
        v = rowData{k};
        if isstruct(v) && isfield(v,'TestData') && istabular_local(v.TestData)
            idxStruct = k;
            break;
        end
    end
end
if ~isempty(idxStruct)
    outStruct = ensureSOC(rowData{idxStruct}, Path);
    return;
end

% 3) flat row, treat as container of SOC and U and MaxAhStep type data
outStruct = normalizeFlatRow(selected);
end

function tag = chooseThirdOutput(selected, checkUpNum, warnedForReturnCol)
% Decide tag for third output
%
% Priority
%   EFC like
%   RPT like
%   CU like
%   Fallback NaN or CheckUp index

E = {'EFC','EffCap','EffectiveCapacity'};
R = {'RPT','ReProTest','ReferencePerformanceTest'};
C = {'CU','CUs','CheckUp','CheckUps','CheckUpNumber','CUNumber','Check_Up'};

if istable(selected)
    value = firstHit(selected, E);
    if ~isempty(value)
        tag = value;
        return;
    end

    value = firstHit(selected, R);
    if ~isempty(value)
        tag = value;
        if ~warnedForReturnCol
            warning('import_OCV:UseRPT', 'EFC column not found. Returned RPT instead.');
        end
        return;
    end

    value = firstHit(selected, C);
    if ~isempty(value)
        tag = value;
        if ~warnedForReturnCol
            warning('import_OCV:UseCU', 'EFC or RPT column not found. Returned CU instead.');
        end
        return;
    end

    if height(selected) == 1
        tag = NaN;
        if ~warnedForReturnCol
            warning('import_OCV:NoTag', 'No EFC, RPT or CU column found; returned NaN.');
        end
    else
        tag = checkUpNum;
        if ~warnedForReturnCol
            warning('import_OCV:UseIndex', 'No EFC or RPT column found; returned CheckUp number as tag.');
        end
    end
else
    tag = NaN;
end
end

function value = firstHit(T, candidates)
% Return first existing column from candidates list, first row

value = [];
for i = 1:numel(candidates)
    cand = candidates{i};
    if ismember(cand, T.Properties.VariableNames)
        value = T.(cand)(1);
        return;
    end
end
end

function S = ensureSOC(S, Path)
% Ensure that TestData has SOC column
%
% Steps
%   1) normalize OCV column names (aliases for U, Ah_Step, SOC)
%   2) if SOC is missing and U and Ah_Step exist, compute SOC

if ~isfield(S,'TestData') || ~istabular_local(S.TestData)
    return;
end

% rename aliases and warn if non canonical names were used
S.TestData = unifyOCVColumnNames(S.TestData, Path);

if ismember('SOC', S.TestData.Properties.VariableNames)
    return;
end

if ~ismember('Ah_Step', S.TestData.Properties.VariableNames) || ...
   ~ismember('U', S.TestData.Properties.VariableNames)
    warning('import_OCV:MissingColumns', ...
        'File "%s": could not find both "Ah_Step" and "U" after name normalization. SOC will not be added.', Path);
    return;
end

try
    AhStep = S.TestData.Ah_Step;
    U      = S.TestData.U;
    S.TestData.SOC = getSOC(AhStep, U);
catch ME
    warning('import_OCV:SocFailure', 'Could not calculate SOC for file "%s": %s', Path, ME.message);
end
end

function T = unifyOCVColumnNames(T, Path)
% Map OCV related column names to canonical U, Ah_Step, SOC
%
% Aliases for U
%   V, U_, V_, E_, Potential, OCP, OCV, Voltage, U_V
%
% Aliases for Ah_Step
%   AhStep, Ahstep, AH_Step, Ah_Step_, Capacity_Step, CapacityStep,
%   ah_step, ah_Step, ahstep
%
% Aliases for SOC
%   StateOfCharge, State_of_Charge

if ~istabular_local(T)
    return;
end

vnames = T.Properties.VariableNames;

alias.U       = {'V','U_','V_','E_','Potential','OCP','OCV','Voltage','U_V'};
alias.Ah_Step = {'AhStep','Ahstep','AH_Step','Ah_Step_','Capacity_Step','CapacityStep', ...
                 'ah_step','ah_Step','ahstep'};
alias.SOC     = {'StateOfCharge','State_of_Charge'};

canonList = fieldnames(alias);

for idx = 1:numel(canonList)
    canon = canonList{idx};
    if any(strcmp(canon, vnames))
        continue;
    end
    candList = alias.(canon);
    hitIdx   = [];
    hitName  = '';
    for j = 1:numel(candList)
        cand = candList{j};
        idxHit = find(strcmpi(vnames, cand), 1, 'first');
        if ~isempty(idxHit)
            hitIdx  = idxHit;
            hitName = vnames{idxHit};
            break;
        end
    end
    if ~isempty(hitIdx)
        T.Properties.VariableNames{hitIdx} = canon;
        vnames{hitIdx} = canon;
        warning('import_OCV:AliasColumnName', ...
            'In file "%s", column "%s" is treated as "%s".', Path, hitName, canon);
    end
end
end

function outStruct = normalizeFlatRow(data)
% Interpret a one row table or struct without nested tabular payload
% and normalize names for SOC, U and MaxAhStep

valid.SOC       = {'SOC','soc','StateOfCharge','State_of_Charge'};
valid.U         = {'U','V','U_','V_','E_','Potential','potential','OCP','ocp', ...
                   'OCV','ocv','Voltage','voltage','U_V'};
valid.MaxAhStep = {'MaxAhStep','maxahstep','AhStep','Ah_Step','ah_step','ah_Step', ...
                   'max_ahstepsize','CapacityStep','capacitystep', 'Capacity_Step', ...
                   'capacity_step', 'capStep', 'capstep', 'CapStep'};

if isstruct(data)
    names = fieldnames(data);
    isTab = false;
else
    names = data.Properties.VariableNames;
    isTab = true;
end

outStruct = struct();
fieldsValid = fieldnames(valid);

for idxKey = 1:numel(fieldsValid)
    key   = fieldsValid{idxKey};
    cands = valid.(key);
    matchIdx  = [];
    usedName  = '';

    for idxCand = 1:numel(cands)
        cand = cands{idxCand};
        hit  = find(strcmpi(names, cand), 1, 'first');
        if ~isempty(hit)
            matchIdx = hit;
            usedName = names{hit};
            break;
        end
    end

    if isempty(matchIdx)
        if isTab
            tdesc = 'table column';
        else
            tdesc = 'struct field';
        end
        error('No %s for "%s". Allowed names are %s', tdesc, key, strjoin(cands, ', '));
    end

    if ~strcmp(usedName, key)
        warning('import_OCV:AliasVarName', ...
            'Variable "%s" is treated as "%s".', usedName, key);
    end

    if isTab
        value = data.(usedName);
    else
        value = data.(usedName);
    end

    outStruct.(key) = value;

    if strcmp(key,'MaxAhStep')
        vals = outStruct.MaxAhStep;
        if iscell(vals)
            % Look for first numeric entry in cell array
            numericVals = [];
            for k = 1:numel(vals)
                if isnumeric(vals{k})
                    numericVals = vals{k};
                    break;
                end
            end
            vals = numericVals;
        end
        if isnumeric(vals)
            valsLocal = vals(~isnan(vals));
            if ~isempty(valsLocal)
                outStruct.MaxAhStep = max(valsLocal);
            else
                outStruct.MaxAhStep = NaN;
            end
        else
            warning('import_OCV:MaxAhStepNonNumeric', ...
                'MaxAhStep values are not numeric; storing NaN.');
            outStruct.MaxAhStep = NaN;
        end
    end
end
end

function tf = istabular_local(x)
% True for tables and timetables

tf = istable(x) || istimetable(x);
end

function SOC = getSOC(AhStep, U, varargin)
% Simple SOC from Ah_Step
%
% This assumes a single direction curve
% SOC is mapped linearly from min(AhStep) to max(AhStep)

[minAh, maxAh] = bounds(AhStep);
if maxAh == minAh
    SOC = zeros(size(AhStep));
else
    SOC = (AhStep - minAh) ./ (maxAh - minAh);
end

% Sanity check: if voltage span is large and SOC and U move in opposite directions, warn
du = diff(U);
ds = diff(SOC);
if ~isempty(du) && ~isempty(ds)
    meanDu = mean(du);
    meanDs = mean(ds);
    if meanDu ~= 0
        plaus = meanDs / meanDu;
        if plaus < 0 && (max(U) - min(U) > 0.5)
            warning('import_OCV:SocDirection', ...
                'SOC is possibly wrong. With increasing voltage SOC is decreasing or vice versa.');
        end
    end
end
end
