%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  load_blend.m                                                   v3.1
%  ------------------------------------------------------------------------
%  Purpose
%  -------
%  Robustly extract a blend curve, i.e. a pair of vectors
%      * voltage              [V]  -> S.voltage
%      * normalized capacity  [-]  -> S.normalizedCapacity
%
%  from almost any MATLAB variable that lives inside a *.mat file.
%  The routine
%      1) loads the first variable found in the file,
%      2) recursively walks through structs, tables, and cell arrays,
%      3) tries to match sensible field or column names such as
%         "U / Voltage / v", "SOC / q / Capacity", ...
%      4) asks once via `inputdlg` if it cannot find a hit.
%
%  Input
%  -----
%  matPath   char | string
%            Absolute or relative path to a *.mat file.
%
%  Output
%  ------
%  S         struct with two equally long vectors
%              * S.voltage            - raw voltage [V]
%              * S.normalizedCapacity - 0...1, 0 = fully discharged,
%                                        1 = fully charged
%
%  Normalization rules
%  -------------------
%  * Vectors whose label contains "soc" are inverted (1 - SOC) so that
%    the curve grows with discharged capacity.
%  * Vectors spanning a range > 1.1 (~ raw Ah) are linearly rescaled.
%  * Everything else is assumed to be 0...1 already, but is clamped.
%
%  Example
%  -------
%      S = load_blend('calendarAging_001.mat');
%      plot(S.normalizedCapacity, S.voltage), grid on
%
%  Version history
%  ---------------
%  v3.0  2025-07-10  "robust-as-hell" original release - M. Rehm
%  v3.1  2025-07-10  * Added full header and detailed docs
%                    * Fixed undefined-variable bug in table branch
%                    * Guarded optional args in dig_for_curve
%                    * Minor name-list cleanup
%
%  Author
%  ------
%  Mathias Rehm  <mathias.rehm@tum.de>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = load_blend(matPath)
%LOAD_BLEND Load a voltage-capacity blend curve from an arbitrary *.mat file.

% ------------------------------------------------------------------ 1) load
raw = load(matPath);
vars = fieldnames(raw);
top = raw.(vars{1});                        % first variable in the file

% ------------------------------------------------------------------ 2) automatic hunt
try
    S = dig_for_curve(top);                 % may recurse many levels
catch
    % fallback -> ask user once
    [vName, qName, ok] = ask_user_for_names(top);
    if ~ok
        error('load_blend:Abort', 'User canceled field-name dialog.');
    end
    S = dig_for_curve(top, vName, qName);
end
end

% ===================================================================
%  dig_for_curve  - recursive extractor
% ===================================================================
function out = dig_for_curve(data, vNameUser, qNameUser)
% Optional-argument sanity (prevents "undefined variable" errors)
if nargin < 2 || isempty(vNameUser), vNameUser = ''; end
if nargin < 3 || isempty(qNameUser), qNameUser = ''; end

% -- unwrap "TestData" layer if present (case-insensitive) -------------
if isstruct(data)
    tdField = first_hit(fieldnames(data), {'TestData'}, 'exact');
    if ~isempty(tdField)
        out = dig_for_curve(data.(tdField), vNameUser, qNameUser);
        return
    end
end

% -- 1) struct --------------------------------------------------------
if isstruct(data)
    fn = fieldnames(data);

    % voltage keys (user string first -> highest priority)
    voltKeys = [cellstr(vNameUser), {'voltage', 'u', 'v'}, ...
                fn(contains(lower(fn), 'volt'))];
    % capacity keys
    capKeys = [cellstr(qNameUser), ...
                {'normalizedcapacity', 'soc', 'capacity', 'q'}, ...
                fn(contains(lower(fn), {'ah', 'cap'}))];

    vField = first_hit(fn, voltKeys);
    qField = first_hit(fn, capKeys);

    if ~isempty(vField) && ~isempty(qField)
        vecV = data.(vField)(:);
        vecQ = data.(qField)(:);
        out = pack_result(vecV, vecQ, qField);
        return
    end
end

% -- 2) table ---------------------------------------------------------
if istabular(data)
    vn = data.Properties.VariableNames;

    vCand = [{'voltage', 'u', 'v'}, cellstr(vNameUser)];
    qCand = [{'normalizedcapacity', 'soc', 'capacity', 'q', 'ah'}, ...
             cellstr(qNameUser)];

    vCol = first_hit(vn, vCand, 'contains');
    qCol = first_hit(vn, qCand, 'contains');

    if ~isempty(vCol) && ~isempty(qCol)
        vecV = data.(vCol)(:);
        vecQ = data.(qCol)(:);
        out = pack_result(vecV, vecQ, qCol);
        return
    end
end

% -- 3) cell array (scan every element) -------------------------------
if iscell(data)
    for idx = 1:numel(data)
        try
            out = dig_for_curve(data{idx}, vNameUser, qNameUser);
            return                                  % first success wins
        catch
        end
    end
end

error('load_blend:NotFound', 'Voltage and capacity vectors not found.');
end

% ===================================================================
%  pack_result  - build output struct, normalize capacity
% ===================================================================
function S = pack_result(voltageVec, capVec, capLabel)
% Normalize capacity depending on its label and numeric range.
lowlab = lower(capLabel);

if contains(lowlab, 'soc')                   % invert SOC
    cap = 1 - capVec;

elseif max(capVec) - min(capVec) > 1.1       % raw Ah -> rescale 0-1
    cap = rescale(capVec, 0, 1);

else                                         % assume already 0-1
    cap = min(max(capVec, 0), 1);            % clamp just in case
end

S.voltage = voltageVec;
S.normalizedCapacity = cap;
end

% ===================================================================
%  first_hit  - return first matching name from a pool
% ===================================================================
function hit = first_hit(pool, keys, mode)
if nargin < 3, mode = 'exact'; end
pool = cellstr(pool);            % normalize to char cell
keys = cellstr(keys);

hit = '';
for k = 1:numel(keys)
    key = strtrim(keys{k});
    if isempty(key), continue, end

    switch mode
        case 'exact'
            idx = find(strcmpi(pool, key), 1);

        case 'contains'
            idx = find(contains(lower(pool), lower(key)), 1);
    end

    if ~isempty(idx)
        hit = pool{idx};
        break
    end
end
end

% ===================================================================
%  ask_user_for_names  - interactive fallback
% ===================================================================
function [vName, qName, ok] = ask_user_for_names(obj)

if isstruct(obj)
    inventory = strjoin(fieldnames(obj), ', ');
    infoTxt = sprintf('Struct fields:\n%s', inventory);
elseif istable(obj)
    inventory = strjoin(obj.Properties.VariableNames, ', ');
    infoTxt = sprintf('Table columns:\n%s', inventory);
else
    infoTxt = 'Enter names manually';
end

prompt = {[infoTxt newline newline 'Voltage field / column:'], ...
          'Capacity field / column (normCap, SOC, Ah):'};

answer = inputdlg(prompt, 'Specify names', [1 70; 1 70]);
ok = ~isempty(answer);

if ok
    vName = strtrim(answer{1});
    qName = strtrim(answer{2});
else
    vName = '';
    qName = '';
end
end
