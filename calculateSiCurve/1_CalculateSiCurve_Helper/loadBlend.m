%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  loadBlend.m                                                   v3.1
%  ------------------------------------------------------------------------
%  Purpose
%  -------
%  Robustly extract a *blend curve* – i.e. a pair of vectors
%      · voltage              [V]  → S.voltage
%      · normalised capacity  [–]  → S.normalizedCapacity
%
%  from (almost) any MATLAB variable that lives inside a *.mat file.
%  The routine
%      1) loads the **first** variable found in the file,
%      2) recursively walks through structs, tables, and cell-arrays,
%      3) tries to match sensible field/column names such as
%         “U / Voltage / v”, “SOC / q / Capacity”, …
%      4) asks once via `inputdlg` if it cannot find a hit.
%
%  Input
%  -----
%  matPath   char | string  
%            Absolute or relative path to a *.mat file.
%
%  Output
%  ------
%  S         struct with two equally-long vectors  
%              • S.voltage            – *raw* voltage [V]  
%              • S.normalizedCapacity – 0…1, **0 = fully discharged**,
%                                        1 = fully charged
%
%  Normalisation rules
%  -------------------
%  • Vectors whose label contains “soc” are **inverted** (1 – SOC) so that
%    the curve grows with discharged capacity.  
%  • Vectors spanning a range > 1.1 (≈ raw Ah) are linearly rescaled.  
%  • Everything else is assumed to be 0…1 already (but is clamped).
%
%  Example
%  -------
%      S = loadBlend('calendarAging_001.mat');
%      plot(S.normalizedCapacity, S.voltage), grid on
%
%  Version history
%  ---------------
%  v3.0  2025-07-10  “robust-as-hell” original release – M. Rehm  
%  v3.1  2025-07-10  • Added full header & detailed docs  
%                    • Fixed *undefined-variable* bug in table branch  
%                    • Guarded optional args in digForCurve  
%                    • Minor name-list clean-up
%
%  Author
%  ------
%  Mathias Rehm  <mathias.rehm@tum.de>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = loadBlend(matPath)
%LOADBLEND Load a voltage-capacity blend curve from an arbitrary *.mat file.

% ------------------------------------------------------------------ 1) load
raw  = load(matPath);
vars = fieldnames(raw);
top  = raw.(vars{1});                        % first variable in the file

% ------------------------------------------------------------------ 2) automatic hunt
try
    S = digForCurve(top);                   % may recurse many levels
catch
    % fallback → ask user once
    [vName,qName,ok] = askUserForNames(top);
    if ~ok
        error('loadBlend:Abort','User cancelled field-name dialog.');
    end
    S = digForCurve(top, vName, qName);
end
end

% ===================================================================
%  digForCurve  – recursive extractor
% ===================================================================
function out = digForCurve(data, vNameUser, qNameUser)
% Optional-argument sanity (prevents “undefined variable” errors)
if nargin < 2 || isempty(vNameUser), vNameUser = ''; end
if nargin < 3 || isempty(qNameUser), qNameUser = ''; end

% -- unwrap "TestData" layer if present (case-insensitive) -------------
if isstruct(data)
    tdField = firstHit(fieldnames(data),{'TestData'},'exact');
    if ~isempty(tdField)
        out = digForCurve(data.(tdField), vNameUser, qNameUser);
        return
    end
end

% -- 1) struct --------------------------------------------------------
if isstruct(data)
    fn = fieldnames(data);

    % voltage keys (user string first ⇒ highest priority)
    voltKeys = [cellstr(vNameUser), {'voltage','u','v'}, ...
                fn(contains(lower(fn),'volt'))];
    % capacity keys
    capKeys  = [cellstr(qNameUser), ...
                {'normalizedcapacity','soc','capacity','q'}, ...
                fn(contains(lower(fn),{'ah','cap'}))];

    vField = firstHit(fn, voltKeys);
    qField = firstHit(fn, capKeys );

    if ~isempty(vField) && ~isempty(qField)
        vecV  = data.(vField)(:);
        vecQ  = data.(qField)(:);
        out   = packResult(vecV, vecQ, qField);
        return
    end
end

% -- 2) table ---------------------------------------------------------
if istabular(data)
    vn = data.Properties.VariableNames;

    vCand = [{'voltage','u','v'}, cellstr(vNameUser)];
    qCand = [{'normalizedcapacity','soc','capacity','q','ah'}, ...
             cellstr(qNameUser)];

    vCol = firstHit(vn, vCand,'contains');
    qCol = firstHit(vn, qCand,'contains');

    if ~isempty(vCol) && ~isempty(qCol)
        vecV  = data.(vCol)(:);
        vecQ  = data.(qCol)(:);
        out   = packResult(vecV, vecQ, qCol);
        return
    end
end

% -- 3) cell array (scan every element) -------------------------------
if iscell(data)
    for idx = 1:numel(data)
        try
            out = digForCurve(data{idx}, vNameUser, qNameUser);
            return                                  % first success wins
        catch
        end
    end
end

error('loadBlend:NotFound','Voltage & capacity vectors not found.');
end

% ===================================================================
%  packResult  – build output struct, normalise capacity
% ===================================================================
function S = packResult(voltageVec, capVec, capLabel)
% normalise capacity depending on its label & numeric range
lowlab = lower(capLabel);

if contains(lowlab,'soc')                    % invert SOC
    cap = 1 - capVec;

elseif max(capVec)-min(capVec) > 1.1         % raw Ah → rescale 0–1
    cap = rescale(capVec, 0, 1);

else                                         % assume already 0–1
    cap = min(max(capVec,0),1);              % clamp just in case
end

S.voltage            = voltageVec;
S.normalizedCapacity = cap;
end

% ===================================================================
%  firstHit  – return first matching name from a pool
% ===================================================================
function hit = firstHit(pool, keys, mode)
if nargin<3, mode='exact'; end
pool = cellstr(pool);            % normalise to char cell
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
%  askUserForNames  – interactive fallback
% ===================================================================
function [vName,qName,ok] = askUserForNames(obj)

if isstruct(obj)
    inventory = strjoin(fieldnames(obj),', ');
    infoTxt   = sprintf('Struct fields:\n%s', inventory);
elseif istable(obj)
    inventory = strjoin(obj.Properties.VariableNames,', ');
    infoTxt   = sprintf('Table columns:\n%s', inventory);
else
    infoTxt   = 'Enter names manually';
end

prompt = { [infoTxt newline newline 'Voltage field / column:'], ...
           'Capacity field / column (normCap, SOC, Ah):' };

answer = inputdlg(prompt,'Specify names',[1 70; 1 70]);
ok = ~isempty(answer);

if ok
    vName = strtrim(answer{1});
    qName = strtrim(answer{2});
else
    vName = ''; qName = '';
end
end
