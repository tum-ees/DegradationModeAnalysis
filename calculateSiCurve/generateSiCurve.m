%% 01 HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author  : Mathias Rehm
%  E-mail  : mathias.rehm@tum.de
%  Date    : 2025-07-28
%
%  Function : generateSiCurve
%  ------------------------------------------------------------------------
%  Generates an artificial Si half-cell OCV curve from
%     * a measured graphite-Si blend (blendPath)
%     * a pure graphite reference   (graphitePath | graphiteSource)
%
%  Q_blend = gamma*Q_Si + (1-gamma)*Q_Gr  ->  Q_Si = ...
%
%  GUI (single dialog):
%     - Browse for blend file           (*.mat)
%     - Browse for save target          (*.mat)
%     - Direction (lithiation | delithiation)
%     - Graphite reference dropdown
%     - gamma-Si field
%     - Filter checkbox                 (LOWESS + deduplication on inputs)
%     - Optional PAV (Pool Adjacent Violators) isotonic regression on output
%       silicon curve
%
%  Optional name-value pairs (case-insensitive):
%     'blendPath'        path\to\blend.mat
%     'savePath'         where\to\save\Si.mat
%     'graphitePath'     explicit graphite file (overrides source)
%     'graphiteSource'   Rehm2026 | Rehm2025 | Hossain | Wetjen | Schmitt
%     'lithDirection'    lithiation | delithiation
%     'gammaSi'          silicon fraction gamma  (0 < gamma < 1)
%     'filterInputData'  true | false          (default false)
%                        LOWESS smoothing + deduplication on input curves.
%     'pavOutput'        true | false          (default true)
%                        Apply Pool Adjacent Violators (PAV) isotonic
%                        regression to the output silicon curve to enforce
%                        monotonicity.
%     'plotFlag'         true | false          (default true)
%
%  Returns
%     siliconStruct.voltage            [V]
%     siliconStruct.normalizedCapacity [0...1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function siliconStruct = generateSiCurve(varargin)

%% 02 SETUP
set(groot,'defaultAxesTickLabelInterpreter','latex')   % LaTeX ticks axis-wide

% ---- Defaults -----------------------------------------------------------
blendPath         = "";
savePath          = "";
graphitePath      = "";
graphiteSource    = "Rehm2026";
lithDirection     = "lithiation";
gammaSi           = NaN;
filterInputData   = false;
pavOutput         = true;
plotFlag          = true;

% ---- Parse varargin -----------------------------------------------------
if mod(numel(varargin), 2) ~= 0
    error('Optional inputs must be provided as name-value pairs.');
end

for k = 1:2:numel(varargin)
    key = varargin{k};
    val = varargin{k+1};

    keyStr = lower(char(string(key)));

    switch keyStr
        case 'blendpath'
            blendPath = string(val);

        case 'savepath'
            savePath = string(val);

        case 'graphitepath'
            graphitePath = string(val);

        case 'graphitesource'
            graphiteSource = capitalize(val);

        case 'lithdirection'
            lithDirection = lower(string(val));

        case 'gammasi'
            gammaSi = double(val);

        case 'filterinputdata'
            filterInputData = parseLogicalFlag(val, keyStr);

        case 'pavoutput'
            pavOutput = parseLogicalFlag(val, keyStr);

        case 'plotflag'
            plotFlag = parseLogicalFlag(val, keyStr);

        otherwise
            warning('Unknown parameter "%s" ignored.', keyStr);
    end
end


%% 03 PATHS
thisDir  = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
graphDir = fullfile(projRoot,'InputData','Graphite');

if ~contains(builtin('path'), projRoot)
    addpath(genpath(projRoot));
end

%% 04 GUI
if strlength(blendPath)==0 || isnan(gammaSi) || strlength(savePath)==0
    d = dialog('Name','Generate Silicon Curve','Units','normalized', ...
               'Position',[0.30 0.38 0.50 0.46]);

    y0=0.84; dy=0.12; w=0.70;

    % -- Blend file -------------------------------------------------------
    uicontrol(d,'Style','text','Units','normalized', ...
        'Position',[0.02 y0 0.25 0.08], 'String','Path to blend curve *.mat');
    hBlend = uicontrol(d,'Style','edit','Units','normalized', ...
        'Position',[0.28 y0 w 0.08], 'HorizontalAlignment','left', ...
        'String',char(blendPath));
    uicontrol(d,'Style','push','String','Browse ...','Units','normalized', ...
        'Position',[0.02 y0-0.02 0.20 0.06], ...
        'Callback',@(~,~) setFullPath(hBlend,@uigetfile, ...
            {'*.mat','MAT-files (*.mat)'}, 'Select blend *.mat'));

    % -- Save target ------------------------------------------------------
    y0 = y0 - dy;
    uicontrol(d,'Style','text','Units','normalized', ...
        'Position',[0.02 y0 0.25 0.08], 'String','Save extracted curve');
    hSave = uicontrol(d,'Style','edit','Units','normalized', ...
        'Position',[0.28 y0 w 0.08], 'HorizontalAlignment','left', ...
        'String',char(savePath));
    uicontrol(d,'Style','push','String','Browse ...','Units','normalized', ...
        'Position',[0.02 y0-0.02 0.20 0.06], ...
        'Callback',@(~,~) setFullPath(hSave,@uiputfile, ...
            {'*.mat','MAT-files (*.mat)'}, 'Save extracted Si curve as'));

    % -- Direction --------------------------------------------------------
    y0 = y0 - dy;
    uicontrol(d,'Style','text','Units','normalized', ...
        'Position',[0.02 y0 0.25 0.08], 'String','Direction');
    hDir = uicontrol(d,'Style','popupmenu','Units','normalized', ...
        'Position',[0.28 y0 0.30 0.08], ...
        'String',{'lithiation','delithiation'}, ...
        'Value',strcmpi(lithDirection,'delithiation')+1);

    % -- Graphite reference ----------------------------------------------
    y0 = y0 - dy;
    refs = {'Rehm2026','Schmitt','Rehm2025','Hossain','Wetjen'};
    defaultIdx = find(strcmpi(refs,graphiteSource)); if isempty(defaultIdx); defaultIdx=1; end
    uicontrol(d,'Style','text','Units','normalized', ...
        'Position',[0.02 y0 0.25 0.08], 'String','Graphite ref.');
    hRef = uicontrol(d,'Style','popupmenu','Units','normalized', ...
        'Position',[0.28 y0 0.30 0.08], ...
        'String',refs,'Value',defaultIdx);

    % Schmitt only for lithiation
    hDir.Callback = @(src,~) set(hRef,'String', ...
        iif(src.Value==2, setdiff(refs, {'Schmitt'}, 'stable'), refs));

    % -- gamma-Si ---------------------------------------------------------
    y0 = y0 - dy;
    uicontrol(d,'Style','text','Units','normalized', ...
        'Position',[0.02 y0 0.25 0.08], 'String','gamma-Si');
    hGamma = uicontrol(d,'Style','edit','Units','normalized', ...
        'Position',[0.28 y0 0.20 0.08], 'String',num2str(gammaSi));
    uicontrol(d,'Style','text','Units','normalized', ...
        'Position',[0.50 y0 0.48 0.10],'HorizontalAlignment','left', ...
        'String',['(use Si2-peak value if you have' newline ...
                  'the blend electrode in delithiation direction)']);

    % -- Filter & monotonic output checkboxes -----------------------------
    y0 = y0 - dy;
    hChkFilter = uicontrol(d,'Style','checkbox','Units','normalized', ...
        'Position',[0.02 y0 0.50 0.08], ...
        'String','Filter input data (LOWESS)', 'Value',filterInputData);
    hChkPAV = uicontrol(d,'Style','checkbox','Units','normalized', ...
        'Position',[0.52 y0 0.45 0.08], ...
        'String','Monotone filter output (PAV)', 'Value',pavOutput);

    % -- Info line --------------------------------------------------------
    y0 = y0 - dy;
    uicontrol(d,'Style','text','Units','normalized', ...
        'Position',[0.02 y0 0.96 0.08], ...
        'String','Rehm2026, Schmitt or Rehm2025 recommended - P45B approx. 0.245 gamma-Si', ...
        'FontAngle','italic');

    % -- OK / Cancel ------------------------------------------------------
    uicontrol(d,'Style','push','String','OK','Units','normalized', ...
        'Position',[0.62 0.04 0.15 0.08], 'Callback',@(~,~) uiresume(d));
    uicontrol(d,'Style','push','String','Cancel','Units','normalized', ...
        'Position',[0.80 0.04 0.15 0.08], 'Callback',@(~,~) delete(d));

    uiwait(d);
    if ~isvalid(d); error('Canceled by user.'); end

    % -- Read values (no chained {} right after get!) ---------------------
    blendPath       = string(get(hBlend,'String'));
    savePath        = string(get(hSave ,'String'));

    dirList         = get(hDir,'String');
    lithDirection   = string(dirList{get(hDir,'Value')});

    refList         = get(hRef,'String');
    graphiteSource  = string(refList{get(hRef,'Value')});

    gammaSi         = str2double(get(hGamma,'String'));
    filterInputData = logical(get(hChkFilter,'Value'));
    pavOutput       = logical(get(hChkPAV,'Value'));

    delete(d)
end

%% 05 RESOLVE PATHS
if strlength(graphitePath)==0
    if strcmpi(graphiteSource,'Schmitt') && strcmpi(lithDirection,'delithiation')
        warning('Schmitt reference only exists for lithiation - switching to Rehm2026.');
        graphiteSource = "Rehm2026";
    end
    graphitePath = fullfile(graphDir, ...
        sprintf('Gr_%s_%s.mat', capitalize(lithDirection), capitalize(graphiteSource)));
end

assert(isfile(graphitePath),'Graphite file not found.');
assert(isfile(blendPath)   ,'Blend file not found.');
if gammaSi<=0 || gammaSi>=1 || isnan(gammaSi)
    error('\gamma_{Si} must be a number between 0 and 1.');
end

%% 06 LOAD & PRE-CLEAN
G  = loadOCV(graphitePath);
BL = loadBlend(blendPath);

G  = smoothUnique(G ,filterInputData);
BL = smoothUnique(BL,filterInputData);

%% 07 ALIGN TO COMMON VOLTAGE WINDOW
Vmin = max([min(G.voltage),  min(BL.voltage)]);
Vmax = min([max(G.voltage),  max(BL.voltage)]);

G  = trimAndRenorm(G ,Vmin,Vmax);
BL = trimAndRenorm(BL,Vmin,Vmax);

N  = max(numel(G.voltage), numel(BL.voltage));
V  = linspace(Vmin,Vmax,N).';

QG = interp1(G.voltage , G.normalizedCapacity , V ,'linear',0);
QB = interp1(BL.voltage, BL.normalizedCapacity, V ,'linear',0);

maskFirst = false(size(V));  maskFirst(1) = true;
maskFlat  = (QG==min(QG)) & (QG==max(QG));
maskKeep  = ~(maskFirst | maskFlat);
V  = V(maskKeep);  QG = QG(maskKeep);  QB = QB(maskKeep);

%% 08 CALCULATE SILICON CURVE
QS = (QB - (1-gammaSi).*QG) ./ gammaSi;
QS(QS<0) = 0;   QS(QS>1) = 1;

% --- Optionally enforce monotonicity on output via PAV --------------------
if pavOutput
    if QS(end) < QS(1)
        QS = pavIsotonic(QS, 'nonincreasing');
    else
        QS = pavIsotonic(QS, 'nondecreasing');
    end
end

siliconStruct.voltage            = V;
siliconStruct.normalizedCapacity = QS;

%% 09 SAVE RESULT
if strlength(savePath)==0
    [file,path] = uiputfile({'*.mat','MAT-files (*.mat)'}, ...
                            'Save extracted Si curve as','SiCurve.mat');
    if isequal(file,0); savePath=""; else; savePath=fullfile(path,file); end
end
if strlength(savePath)>0
    if ~endsWith(savePath,'.mat','IgnoreCase',true)
        savePath = savePath + ".mat";
    end
    save(savePath,'siliconStruct','-mat');
end

%% 10 PLOT
if plotFlag
    tumBlue   = [  0 101 189]/255;
    tumOrange = [227 114  34]/255;

    figure('Name','Graphite / Silicon / Blend'); hold on; grid on; box on;
    plot(QG,V,'-','Color',tumBlue  ,'LineWidth',1.6, ...
        'DisplayName',['Graphite (' char(graphiteSource) ')']);
    plot(QS,V,'-','Color',tumOrange,'LineWidth',1.6, ...
        'DisplayName','Silicon');
    plot(QB,V,'--','Color',[0 0 0],'LineWidth',1.6, ...
        'DisplayName','SiGr-Blend');

    xlabel('normalized capacity $Q/Q_{\max}$','Interpreter','latex');
    ylabel('$U$ / V','Interpreter','latex');
    title(sprintf('Artificial silicon curve  (\\gamma_{Si}=%.2f)',gammaSi));
    legend('Interpreter','latex','Location','best');
    ylim([-0.05 0.85]);  xlim([-0.05 1.05]);
end

end   % ------------------- end main function ------------------------------


%% H1 loadOCV --------------------------------------------------------------
function S = loadOCV(matPath)
% Load first struct in a .mat file; verify required fields exist.
    raw = load(matPath);
    fn  = fieldnames(raw);
    S   = raw.(fn{1});
    assert(all(isfield(S,{'voltage','normalizedCapacity'})), ...
        'OCV file must contain: voltage & normalizedCapacity.');
end

%% H2 loadBlend ------------------------------------------------------------
function OUT = loadBlend(matPath)
% Load blend; supports struct or cell array with TestData.
    raw = load(matPath);  fn = fieldnames(raw);  tmp = raw.(fn{1});
    if isstruct(tmp) && all(isfield(tmp,{'voltage','normalizedCapacity'}))
        OUT = tmp;
    elseif iscell(tmp) && isfield(tmp{1},'TestData')
        T = tmp{1}.TestData;
        OUT.voltage            = T.voltage(:);
        OUT.normalizedCapacity = T.normalizedCapacity(:);
    else
        error('Unknown blend file format.');
    end
end

%% H3 smoothUnique ---------------------------------------------------------
function out = smoothUnique(inStruct,doSmooth)
% LOWESS (window 30) on voltage + remove duplicate voltages.
    if doSmooth
        inStruct.voltage = smooth(inStruct.voltage,30,'lowess');
    end
    [uVolt,idx] = unique(inStruct.voltage);
    out.voltage            = uVolt;
    out.normalizedCapacity = inStruct.normalizedCapacity(idx);
end

%% H3b pavIsotonic ---------------------------------------------------------
function q_mono = pavIsotonic(q_raw, direction)
% Pool Adjacent Violators (PAV) isotonic regression.
% Merges adjacent monotonicity-violating blocks by replacing them with
% their weighted mean.
% Finds q_mono minimizing sum((q_mono - q_raw)^2) subject to monotonicity.
%   direction : 'nondecreasing' (default) | 'nonincreasing'
    if nargin < 2, direction = 'nondecreasing'; end

    q = q_raw(:);
    n = numel(q);

    if strcmpi(direction, 'nonincreasing')
        q = -q;
    end

    % Block arrays: value = weighted mean, sz = number of points
    val = zeros(n,1);
    sz  = zeros(n,1);
    nb  = 0;                              % number of active blocks

    for i = 1:n
        nb      = nb + 1;
        val(nb) = q(i);
        sz(nb)  = 1;

        % Merge while the last two blocks violate non-decreasing order
        while nb > 1 && val(nb-1) > val(nb)
            total     = sz(nb-1) + sz(nb);
            val(nb-1) = (val(nb-1)*sz(nb-1) + val(nb)*sz(nb)) / total;
            sz(nb-1)  = total;
            nb        = nb - 1;
        end
    end

    % Write block values back to output vector
    q_mono = zeros(n,1);
    idx = 1;
    for k = 1:nb
        q_mono(idx : idx+sz(k)-1) = val(k);
        idx = idx + sz(k);
    end

    if strcmpi(direction, 'nonincreasing')
        q_mono = -q_mono;
    end
end

%% H3c parseLogicalFlag ----------------------------------------------------
function tf = parseLogicalFlag(val, paramName)
% Parse logical option from logical, numeric, string, or char scalar input.
    if islogical(val) && isscalar(val)
        tf = val;
        return
    end

    if isnumeric(val) && isscalar(val)
        if ~ismember(val, [0 1])
            error('Parameter "%s" must be true/false or 0/1.', paramName);
        end
        tf = logical(val);
        return
    end

    if (isstring(val) && isscalar(val)) || ischar(val)
        valStr = lower(strtrim(char(string(val))));
        if any(strcmp(valStr, {'true','1','on','yes'}))
            tf = true;
            return
        end
        if any(strcmp(valStr, {'false','0','off','no'}))
            tf = false;
            return
        end
    end

    error('Parameter "%s" must be a logical scalar.', paramName);
end

%% H4 trimAndRenorm --------------------------------------------------------
function out = trimAndRenorm(inStruct,Vmin,Vmax)
% Trim to [Vmin,Vmax] and rescale capacity to [0...1].
    mask = inStruct.voltage>=Vmin & inStruct.voltage<=Vmax;
    inStruct.voltage            = inStruct.voltage(mask);
    inStruct.normalizedCapacity = inStruct.normalizedCapacity(mask);
    inStruct.normalizedCapacity = rescale(inStruct.normalizedCapacity,0,1);
    out = inStruct;
end

%% H5 smoothLocalDescent ---------------------------------------------------
% -> functionality removed

%% H6 capitalize -----------------------------------------------------------
function str = capitalize(s)
% Capitalize first letter; rest lower-case.
    s = char(s);
    str = [upper(s(1)) lower(s(2:end))];
end

%% H7 iif ------------------------------------------------------------------
function y = iif(cond,a,b)
% Inline if: returns a if cond else b.
    if cond, y=a; else, y=b; end
end

%% H8 setFullPath ----------------------------------------------------------
function setFullPath(hEdit,dlgFcn,filterSpec,dlgTitle)
% Helper for Browse buttons -> write full path into edit field.
    [file,path] = dlgFcn(filterSpec,dlgTitle);
    if isequal(file,0), return; end
    full = fullfile(path,file);
    if isequal(dlgFcn,@uiputfile) && ~endsWith(full,'.mat','IgnoreCase',true)
        full = full + ".mat";
    end
    set(hEdit,'String', full);
end
