function [rawSOC, rawU, Capa] = parse_data_input(inData, dataType)
%> Authors : Mathias Rehm
%> E mail  : mathias.rehm@tum.de
%> Date    : 2025-09-12
%> Additional code by Moritz Guenthner (moritz.guenthner@tum.de)
%
% parse_data_input
% Parse input data struct or cell and extract SOC, voltage, and capacity
% warn if duplicate SOC values are detected exact equality
%
% This function
%   Provides a flexible way to extract SOC rawSOC, voltage rawU, and 
%   capacity MaxAhStep from various possible data input formats struct 
%   versus cell array versus loaded mat struct
%   Always assumes that if the input is a dataList cell array the 
%   relevant entry is at 1,1
%   For the voltage field 
%       If dataType is cathode or cathodeBlend2 we first look for UCa If not found 
%         we fallback to U
%       If dataType is anode we first look for UAn; if not found 
%         we fallback to U
%       If dataType is anodeBlend2 we first look for UAn; if not found 
%         we fallback to U. In this case we do not check whether SOC
%         values are unique as a calculated blend2 curve might have non unique
%         SOC values
%       If dataType is fullcell we directly look for U only
%       Additionally if the above are not found voltage is allowed as a 
%         second fallback
%   For the SOC field
%       We first look for SOC If not found we fallback to 
%       normalizedCapacity
%   If we cannot find capacity it defaults to NaN
%   If essential fields SOC voltage are not found an error is thrown
%
% SOC ordering and monotonicity
%   For non-blend2 curves SOC is cleaned so it increases from start to end
%   by removing SOC/U pairs that would break a non-decreasing sequence
%   For blend2 curves only a global flip is applied if needed and
%   internally non-monotonic SOC is left untouched but a warning is given
%
% Usage
%   rawSOC, rawU, Capa  parse_data_input(inData, dataType)
%
% Inputs
%   inData  Data in one of the recognized formats
%           1 cell array with inData1,1 TestData SOC and U or 
%             TestData SOC and UCa or TestData SOC and UAn
%           2 struct with a TestData sub struct containing 
%             TestData SOC and TestData U or UCa or UAn
%           3 struct with top level fields SOC and U or UCa or UAn
%   dataType  A string specifying the nature of this data
%             cathode        tries UCa first fallback to U
%             cathodeBlend2  same as cathode
%             anode          tries UAn first fallback to U
%             anodeBlend2    same as anode but duplicate SOC check is skipped
%             fullcell       uses U only
%
% Outputs
%   rawSOC  Vector of SOC data
%   rawU    Vector of voltage data
%   Capa    Capacity MaxAhStep if available else NaN

%% 1) Initialize outputs and decide voltage fields
rawSOC = [];
rawU   = [];
Capa   = NaN;  % If we cannot find capacity, leave as NaN
curveIsBlend2 = false; % flag for blend2 curves no duplicate SOC check and no point removal

switch lower(dataType)
    case 'cathode'
        primaryUField   = 'UCa';
        fallbackUField  = 'U';
    case 'cathodeblend2'
        primaryUField   = 'UCa';
        fallbackUField  = 'U';
        curveIsBlend2   = true;
    case 'anodeblend2'
        primaryUField   = 'UAn';
        fallbackUField  = 'U';
        curveIsBlend2   = true;
    case 'anode'
        primaryUField   = 'UAn';
        fallbackUField  = 'U';
    case 'fullcell'
        primaryUField   = 'U';   % no fallback for full cell
        fallbackUField  = '';
    otherwise
        error('parseDataInput: unknown dataType "%s". Use "cathode", "cathodeBlend2", "anode", "anodeBlend2" or "fullcell".', dataType);
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
        Capa   = getMaxAhStep(thisData);
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

%% 3) SOC ordering and monotonicity handling
if ~isempty(rawSOC) && ~isempty(rawU)
    % First handle global direction so SOC starts low and ends high
    if rawSOC(1) > rawSOC(end)
        rawSOC = flip(rawSOC);
        rawU   = flip(rawU);
    end

    % Mixed sign differences indicate truly non-monotonic SOC
    socVec   = rawSOC(:);
    d        = diff(socVec);
    dFinite  = d(isfinite(d));
    dNonZero = dFinite(dFinite ~= 0);

    if ~isempty(dNonZero)
        sgn = sign(dNonZero);
        if ~(all(sgn > 0) || all(sgn < 0))
            if curveIsBlend2
                % Blend2 curves keep all points but must still warn
                warning('parseDataInput:NonMonotonicSOC', ...
                    ['SOC for "%s" data is non-monotonic; order kept. ' ...
                    '-> This is not a problem in case of blend electrodes!'], dataType);
            else
                % For all other curves remove offending points
                [rawSOC, rawU, nRemoved] = removeNonMonotonicSOCPoints(rawSOC, rawU, dataType);
                if nRemoved > 0
                    warning('parseDataInput:NonMonotonicSOCCleaned', ...
                        ['SOC for "%s" data is non-monotonic. ', ...
                         '%d SOC/voltage point%s removed to enforce non-decreasing SOC.'], ...
                         dataType, nRemoved, char('s'*(nRemoved > 1)));
                end
            end
        end
    end
end

%% 4) Duplicate SOC check
%    Exact equality is used here for numeric SOC
if isnumeric(rawSOC) && ~isempty(rawSOC) && ~curveIsBlend2
    v = rawSOC(:);
    v = v(isfinite(v));              % ignore NaN and Inf
    if numel(v) > 1
        vv = sort(v);                % n log n
        dupMask = diff(vv) == 0;     % exact duplicates
        if any(dupMask)
            dupVals = unique(vv(dupMask));
            nDup = numel(dupVals);
            k = min(5, nDup);
            ex = dupVals(1:k).';
            exStr = sprintf('%.8g, ', ex);
            exStr = exStr(1:end-2);  % strip trailing ", "
            if nDup > k
                exStr = [exStr, ', ...'];
            end
            warning('parseDataInput:DuplicateSOC', ...
                ['Duplicate SOC values detected (%d unique duplicate value%s). ', ...
                 'This may cause large problems in interpolation, fitting, and lookup table operations ', ...
                 'for example non monotonic axes or NaNs from interp1. ', ...
                 'Please ensure SOC is strictly unique. Examples: %s'], ...
                 nDup, char('s'*(nDup>1)), exStr);
        end
    end
end

%% 5) Final checks
if isempty(rawSOC) || isempty(rawU)
    error('parseDataInput: Could not find valid SOC or Voltage data for dataType="%s".', dataType);
end

end

% =====================================================================
% Local helper function chooseVoltageField
%   Tries primaryUField in the struct then fallback if needed
%   Additionally allows voltage as a second fallback
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
% Local helper function chooseSOCField
%   Tries SOC in the struct then falls back to normalizedCapacity
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
% Local helper getMaxAhStep
%   Returns absolute last minus first from the best available Ah like column
%   Ah_Step  Ah  any name containing Ah case insensitive
%   Returns NaN if none
% =====================================================================
function Capa = getMaxAhStep(s)
if istable(s) || istimetable(s)
    names = s.Properties.VariableNames;
    get = @(n) s.(n);
elseif isstruct(s)
    names = fieldnames(s);
    get = @(n) s.(n);
else
    Capa = NaN;
    return;
end

% 1 exact Ah_Step
if any(strcmp(names,'Ah_Step'))
    v = get('Ah_Step');
    if isnumeric(v) && ~isempty(v)
        Capa = abs(v(end) - v(1));
        return;
    end
end
% 2 exact Ah
if any(strcmp(names,'Ah'))
    v = get('Ah');
    if isnumeric(v) && ~isempty(v)
        Capa = abs(v(end) - v(1));
        return;
    end
end
% 3 any name containing Ah choose larger amplitude
idx = find(contains(upper(string(names)),'AH'));
Capa = NaN;
best =  -inf;
for k = idx(:).'
    v = get(names{k});
    if ~isnumeric(v) || isempty(v)
        continue;
    end
    amp = abs(v(end) - v(1));
    if amp > best
        best = amp;
        Capa = amp;
    end
end
end

% =====================================================================
% Local helper function isfield2 
%   Robust isfield variant for struct and tabular types
% =====================================================================
function tf = isfield2(s, fName)
if ~(isstruct(s) || istabular(s))
    tf = false;
    return;
end
fn = fieldnames(s);
tf = any(strcmp(fn, fName));
end

% =====================================================================
% Local helper removeNonMonotonicSOCPoints
%   For non-blend2 data types remove SOC/U pairs that would break
%   a non-decreasing SOC sequence after the global flip step
% =====================================================================
function [socOut, uOut, nRemoved] = removeNonMonotonicSOCPoints(socIn, uIn, dataType)
if numel(socIn) ~= numel(uIn)
    error('parseDataInput:LengthMismatch', ...
        'SOC and voltage length mismatch for "%s" data (SOC=%d, U=%d).', ...
        dataType, numel(socIn), numel(uIn));
end

keepRow = isrow(socIn);
socVec = socIn(:);
uVec   = uIn(:);

n = numel(socVec);
if n <= 2
    socOut = socIn;
    uOut   = uIn;
    nRemoved = 0;
    return;
end

keep = false(n,1);

% start from the first finite SOC value
idx0 = find(isfinite(socVec), 1, 'first');
if isempty(idx0)
    socOut = [];
    uOut   = [];
    nRemoved = n;
    return;
end

keep(idx0) = true;
lastSOC = socVec(idx0);

for i = idx0+1:n
    s = socVec(i);
    if ~isfinite(s)
        % drop NaN or Inf entries
        keep(i) = false;
    else
        if s >= lastSOC
            keep(i) = true;
            lastSOC = s;
        else
            % this point would make SOC decrease so we drop it
            keep(i) = false;
        end
    end
end

socClean = socVec(keep);
uClean   = uVec(keep);
nRemoved = n - sum(keep);

if keepRow
    socOut = socClean.';
    uOut   = uClean.';
else
    socOut = socClean;
    uOut   = uClean;
end
end
