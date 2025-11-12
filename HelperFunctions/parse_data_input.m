function [rawSOC, rawU, Capa] = parse_data_input(inData, dataType)
%> Authors : Mathias Rehm
%> E-mail  : mathias.rehm@tum.de
%> Date    : 2025-09-12
%> - Additional code by Moritz Guenthner (moritz.guenthner@tum.de)
%
% parse_data_input
% Parse input data (struct/cell) and extract SOC, voltage, and capacity;
% warn if duplicate SOC values are detected (exact equality).
%
% This function:
%   * Provides a flexible way to extract SOC (rawSOC), voltage (rawU), and 
%     capacity (MaxAhStep) from various possible data input formats (struct 
%     vs. cell array vs. loaded .mat struct).
%   * Always assumes that if the input is a "dataList" (cell array), the 
%     relevant entry is at {1,1}. 
%   * For the voltage field: 
%       - If dataType == 'cathode', we first look for 'UCa'. If not found, 
%         we fallback to 'U'. 
%       - If dataType == 'anode', we first look for 'UAn'. If not found, 
%         we fallback to 'U'. 
%       - If dataType == 'anodeBlend2', we first look for 'UAn'. If not found, 
%         we fallback to 'U'.  -> in this case: we do not check whether SOC
%         values are unique (as a calculated blend2 curve might have non-unique
%         SOC values, see Rehm et al. (2026))
%       - If dataType == 'fullcell', we directly look for 'U' only.
%       - Additionally, if the above are not found, 'voltage' is allowed as a 
%         second fallback.
%   * For the SOC field:
%       - We first look for 'SOC'. If not found, we fallback to 
%         'normalizedCapacity'.
%   * If we cannot find capacity, it defaults to NaN.
%   * If essential fields (SOC, voltage) are not found, an error is thrown.
%
% Usage:
%   [rawSOC, rawU, Capa] = parse_data_input(inData, dataType)
%
% Inputs:
%   inData  - Data in one of the recognized formats:
%             (1) cell array, with inData{1,1}.TestData.(SOC & U) or 
%                 .TestData.(SOC & UCa) or .TestData.(SOC & UAn),
%             (2) struct with a .TestData sub-struct, containing 
%                 .TestData.SOC and .TestData.(U / UCa / UAn),
%             (3) struct with top-level fields .SOC and .(U / UCa/UAn).
%   dataType - A string specifying the nature of this data:
%              'cathode' => tries 'UCa' first, fallback to 'U'
%              'anode'   => tries 'UAn' first, fallback to 'U'
%              'fullcell'=> uses 'U' only
%
% Outputs:
%   rawSOC - The vector of SOC data extracted from inData
%   rawU   - The vector of Voltage data extracted from inData
%   Capa   - The capacity (MaxAhStep) if available, else NaN

%% 1) Initialize outputs and decide voltage fields
    rawSOC = [];
    rawU   = [];
    Capa   = NaN;  % If we cannot find capacity, leave as NaN
    curveIsBlend2 = false; % flag for blend2 curves (no duplicate SOC check)

    switch lower(dataType)
        case 'cathode'
            primaryUField   = 'UCa';
            fallbackUField  = 'U';
        case 'anodeblend2'
            primaryUField   = 'UAn';
            fallbackUField  = 'U';  
        case 'anode'
            primaryUField   = 'UAn';
            fallbackUField  = 'U';          
        case 'fullcell'
            primaryUField   = 'U';   % no fallback for full-cell
            fallbackUField  = '';
        otherwise
            error('parseDataInput: unknown dataType "%s". Use "cathode", "anode", or "fullcell".', dataType);
    end

%% 2) Extract data from input container
    if iscell(inData)
        if size(inData,1) < 1 || size(inData,2) < 1
            error('parseDataInput: inData cell array is empty.');
        end
        thisData = inData{1,1};
        
        if isfield2(thisData, 'TestData')
            rawSOC = chooseSOCField(thisData.TestData);
            rawU   = chooseVoltageField(thisData.TestData, primaryUField, fallbackUField);
            Capa   = getMaxAhStep(thisData.TestData);
        else
            rawSOC = chooseSOCField(thisData);
            if ~isempty(rawSOC)
                rawU = chooseVoltageField(thisData, primaryUField, fallbackUField);
            end
            Capa     = getMaxAhStep(thisData);
        end

    elseif isstruct(inData)
        if isfield2(inData, 'TestData')
            rawSOC = chooseSOCField(inData.TestData);
            rawU   = chooseVoltageField(inData.TestData, primaryUField, fallbackUField);
            Capa   = getMaxAhStep(inData.TestData);
        else
            rawSOC = chooseSOCField(inData);
            rawU   = chooseVoltageField(inData, primaryUField, fallbackUField);
            Capa   = getMaxAhStep(inData);
        end

    else
        error('parseDataInput: unrecognized input type. Must be cell array or struct.');
    end

%% 3) Duplicate SOC check (efficient: sort + diff on finite values)
%    NOTE: Exact equality is used here. If you prefer a tolerance-based
%    check for quantized data, replace "diff(vv) == 0" with, e.g.,
%    "abs(diff(vv)) <= 1e-9".
    if isnumeric(rawSOC) && ~isempty(rawSOC) && ~curveIsBlend2
        v = rawSOC(:);
        v = v(isfinite(v));              % ignore NaN/Inf
        if numel(v) > 1
            vv = sort(v);                % O(n log n)
            dupMask = diff(vv) == 0;     % exact duplicates
            if any(dupMask)
                dupVals = unique(vv(dupMask));
                nDup = numel(dupVals);
                % Show up to 5 example duplicate values
                k = min(5, nDup);
                ex = dupVals(1:k).';
                exStr = sprintf('%.8g, ', ex);
                exStr = exStr(1:end-2);  % strip trailing ", "
                if nDup > k
                    exStr = [exStr, ', ...'];
                end
                % Strong warning: duplicates may cause huge downstream problems
                warning('parseDataInput:DuplicateSOC', ...
                    ['Duplicate SOC values detected (%d unique duplicate value%s). ', ...
                     'This may cause HUGE problems in interpolation, fitting, and lookup-table operations ', ...
                     '(e.g., non-monotonic axes, NaNs from interp1, ill-conditioned regressions). ', ...
                     'Please ensure SOC is strictly unique. Examples: %s'], ...
                     nDup, char('s'*(nDup>1)), exStr);
            end
        end
    end

%% 4) Final checks
    if isempty(rawSOC) || isempty(rawU)
        error('parseDataInput: Could not find valid SOC or Voltage data for dataType="%s".', dataType);
    end

end

% =====================================================================
% Local helper function: chooseVoltageField
%   Tries primaryUField in the struct, then fallback if needed.
%   Additionally, allows 'voltage' as a second fallback.
% =====================================================================
function val = chooseVoltageField(s, primaryField, fallbackField)
    if isfield2(s, primaryField)
        tmp = s.(primaryField);
        if ~isempty(tmp)
            val = tmp;
            return;
        end
    end
    if ~isempty(fallbackField) && isfield2(s, fallbackField)
        tmp = s.(fallbackField);
        if ~isempty(tmp)
            val = tmp;
            return;
        end
    end
    if isfield2(s, 'voltage')
        tmp = s.voltage;
        if ~isempty(tmp)
            val = tmp;
            return;
        end
    end
    val = [];
end

% =====================================================================
% Local helper function: chooseSOCField
%   Tries 'SOC' in the struct, then falls back to 'normalizedCapacity'.
% =====================================================================
function val = chooseSOCField(s)
    if isfield2(s, 'SOC')
        tmp = s.SOC;
        if ~isempty(tmp)
            val = tmp;
            return;
        end
    end
    if isfield2(s, 'normalizedCapacity')
        tmp = s.normalizedCapacity;
        if ~isempty(tmp)
            val = tmp;
            return;
        end
    end
    val = [];
end

% =====================================================================
% Local helper: getMaxAhStep
%   Returns |last-first| from the best available Ah-like column:
%   'Ah_Step' > 'Ah' > any name containing 'Ah' (case-insensitive; if
%   multiple, choose the one with larger |last-first|). Returns NaN if none.
% =====================================================================
function Capa = getMaxAhStep(s)
    if istable(s) || istimetable(s), names = s.Properties.VariableNames; get = @(n) s.(n);
    elseif isstruct(s),              names = fieldnames(s);             get = @(n) s.(n);
    else, Capa = NaN; return; 
    end

    % 1) Exact 'Ah_Step'
    if any(strcmp(names,'Ah_Step')), v = get('Ah_Step'); if isnumeric(v)&&~isempty(v), Capa = abs(v(end)-v(1)); return, end, end
    % 2) Exact 'Ah'
    if any(strcmp(names,'Ah')),      v = get('Ah');      if isnumeric(v)&&~isempty(v), Capa = abs(v(end)-v(1)); return, end, end
    % 3) Any name containing 'Ah' → pick larger amplitude
    idx = find(contains(upper(string(names)),'AH')); Capa = NaN; best = -inf;
    for k = idx(:).'
        v = get(names{k}); if ~isnumeric(v) || isempty(v), continue, end
        amp = abs(v(end)-v(1)); if amp > best, best = amp; Capa = amp; end
    end
end

% =====================================================================
% Local helper function: isfield2 
%   A more robust variant of isfield, returning true only if the 
%   field genuinely exists among fieldnames(s). Helps avoid 
%   certain corner cases with overshadowing or 0x0 struct arrays.
% =====================================================================
function tf = isfield2(s, fName)
    if ~isstruct(s) && ~istabular(s)
        tf = false;
        return;
    end
    fn = fieldnames(s);
    tf = any(strcmp(fn, fName));
end
