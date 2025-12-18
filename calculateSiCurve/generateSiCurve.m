%% 01 HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author  : Mathias Rehm
%  E-mail  : mathias.rehm@tum.de
%  Date    : 2025-07-28
%
%  Function : generateSiCurve   ·   Version 3.7
%  ------------------------------------------------------------------------
%  Generates an artificial Si half-cell OCV curve from
%     · a measured graphite-Si blend (BlendPath)
%     · a pure graphite reference   (GraphitePath | GraphiteSource)
%
%  Q_blend = γ·Q_Si + (1-γ)·Q_Gr   →   Q_Si = (Q_blend − (1-γ)Q_Gr)/γ
%
%  GUI (single dialog):
%     – Browse for blend file           (*.mat)
%     – Browse for save target          (*.mat)
%     – Direction (lithiation | delithiation)
%     – Graphite reference dropdown
%     – γ-Si field
%     – Filter checkbox                 (LOWESS + deduplication)
%     – SmoothBadQS checkbox            (Savitzky–Golay on jump regions)
%
%  Optional name–value pairs (case-insensitive):
%     'BlendPath'        path\to\blend.mat
%     'SavePath'         where\to\save\Si.mat
%     'GraphitePath'     explicit graphite file (overrides source)
%     'GraphiteSource'   Kuecher | Hossain | Wetjen | Schmitt | Rehm
%     'LithDirection'    lithiation | delithiation
%     'GammaSi'          silicon fraction γ  (0 < γ < 1)
%     'FilterInputData'  true | false          (default true)
%     'SmoothBadQS'      true | false          (default true)
%     'SmoothWindow'     odd frame length for sgolay      (default 11)
%     'PlotFlag'         true | false          (default true)
%
%  Returns
%     siliconStruct.voltage            [V]
%     siliconStruct.normalizedCapacity [0…1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function siliconStruct = generateSiCurve(varargin)

%% 02 SETUP
set(groot,'defaultAxesTickLabelInterpreter','latex')   % LaTeX ticks axis-wide

% ---- Defaults -----------------------------------------------------------
BlendPath         = "";              
SavePath          = "";              
GraphitePath      = "";              
GraphiteSource    = "Kuecher";       
LithDirection     = "lithiation";    
GammaSi           = NaN;             
filterInputData   = true;            
smoothBadQS       = true;            
smoothWindow      = 71;              
PlotFlag          = true;            

% ---- Parse varargin -----------------------------------------------------
for k = 1:2:numel(varargin)
    key = varargin{k};
    val = varargin{k+1};

    keyStr = lower(char(key));  % work with lower case char

    if contains(keyStr, 'blend')
        BlendPath = string(val);

    elseif contains(keyStr, 'save') && contains(keyStr, 'path')
        % requires both "save" and "path" in the key, e.g. "SavePath"
        SavePath = string(val);

    elseif contains(keyStr, 'graphite') && contains(keyStr, 'path')
        % explicit graphite file
        GraphitePath = string(val);

    elseif contains(keyStr, 'graphite')
        % graphite source selector
        GraphiteSource = capitalize(val);

    elseif contains(keyStr, 'lith') && contains(keyStr, 'dir')
        % e.g. "LithDirection"
        LithDirection = lower(string(val));

    elseif contains(keyStr, 'gamma')
        GammaSi = double(val);

    elseif contains(keyStr, 'filter')
        filterInputData = logical(val);

    elseif contains(keyStr, 'smooth') && contains(keyStr, 'window')
        % smooth window length for sgolay
        smoothWindow = double(val);

    elseif contains(keyStr, 'smooth') && contains(keyStr, 'bad') && contains(keyStr, 'qs')
        % SmoothBadQS flag
        smoothBadQS = logical(val);

    elseif contains(keyStr, 'plot')
        PlotFlag = logical(val);

    else
        warning('Unknown parameter "%s" ignored.', keyStr);
    end
end


%% 03 PATHS
thisDir  = fileparts(mfilename('fullpath'));          
projRoot = fileparts(thisDir);                        
graphDir = fullfile(projRoot,'InputData','Graphite'); 

if ~contains(feval('path'), projRoot)
    addpath(genpath(projRoot));
end

%% 04 GUI
if strlength(BlendPath)==0 || isnan(GammaSi) || strlength(SavePath)==0
    d = dialog('Name','Generate Silicon Curve','Units','normalized', ...
               'Position',[0.30 0.38 0.50 0.46]);

    y0=0.84; dy=0.12; w=0.70;

    % -- Blend file -------------------------------------------------------
    uicontrol(d,'Style','text','Units','normalized', ...
        'Position',[0.02 y0 0.25 0.08], 'String','Path to blend curve *.mat');
    hBlend = uicontrol(d,'Style','edit','Units','normalized', ...
        'Position',[0.28 y0 w 0.08], 'HorizontalAlignment','left', ...
        'String',char(BlendPath));
    uicontrol(d,'Style','push','String','Browse …','Units','normalized', ...
        'Position',[0.02 y0-0.02 0.20 0.06], ...
        'Callback',@(~,~) setFullPath(hBlend,@uigetfile, ...
            {'*.mat','MAT-files (*.mat)'}, 'Select blend *.mat'));

    % -- Save target ------------------------------------------------------
    y0 = y0 - dy;
    uicontrol(d,'Style','text','Units','normalized', ...
        'Position',[0.02 y0 0.25 0.08], 'String','Save extracted curve');
    hSave = uicontrol(d,'Style','edit','Units','normalized', ...
        'Position',[0.28 y0 w 0.08], 'HorizontalAlignment','left', ...
        'String',char(SavePath));
    uicontrol(d,'Style','push','String','Browse …','Units','normalized', ...
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
        'Value',strcmpi(LithDirection,'delithiation')+1);

    % -- Graphite reference ----------------------------------------------
    y0 = y0 - dy;
    refs = {'Kuecher','Schmitt','Hossain','Wetjen','Rehm'};
    defaultIdx = find(strcmpi(refs,GraphiteSource)); if isempty(defaultIdx); defaultIdx=1; end
    uicontrol(d,'Style','text','Units','normalized', ...
        'Position',[0.02 y0 0.25 0.08], 'String','Graphite ref.');
    hRef = uicontrol(d,'Style','popupmenu','Units','normalized', ...
        'Position',[0.28 y0 0.30 0.08], ...
        'String',refs,'Value',defaultIdx);

    % Schmitt only for lithiation
    hDir.Callback = @(src,~) set(hRef,'String', ...
        iif(src.Value==2, setdiff(refs, {'Schmitt'}, 'stable'), refs));

    % -- γ-Si -------------------------------------------------------------
    y0 = y0 - dy;
    gammaChar = char(hex2dec('03B3'));   % γ
    uicontrol(d,'Style','text','Units','normalized', ...
        'Position',[0.02 y0 0.25 0.08], 'String',[gammaChar '-Si']);
    hGamma = uicontrol(d,'Style','edit','Units','normalized', ...
        'Position',[0.28 y0 0.20 0.08], 'String',num2str(GammaSi));
    uicontrol(d,'Style','text','Units','normalized', ...
        'Position',[0.50 y0 0.48 0.10],'HorizontalAlignment','left', ...
        'String',['(use Si2-peak value if you have' newline ...
                  'the blend electrode in delithiation direction)']);

    % -- Filter & Smooth checkboxes --------------------------------------
    y0 = y0 - dy;
    hChkFilter = uicontrol(d,'Style','checkbox','Units','normalized', ...
        'Position',[0.02 y0 0.50 0.08], ...
        'String','Filter input data (LOWESS)', 'Value',filterInputData);
    hChkSmooth = uicontrol(d,'Style','checkbox','Units','normalized', ...
        'Position',[0.52 y0 0.45 0.08], ...
        'String','Smooth bad QS (sgolay)', 'Value',smoothBadQS);

    % -- Info line --------------------------------------------------------
    y0 = y0 - dy;
    uicontrol(d,'Style','text','Units','normalized', ...
        'Position',[0.02 y0 0.96 0.08], ...
        'String','Kuecher or Schmitt recommended – P45B ≈ 0.245 γ-Si', ...
        'FontAngle','italic');

    % -- OK / Cancel ------------------------------------------------------
    uicontrol(d,'Style','push','String','OK','Units','normalized', ...
        'Position',[0.62 0.04 0.15 0.08], 'Callback',@(~,~) uiresume(d));
    uicontrol(d,'Style','push','String','Cancel','Units','normalized', ...
        'Position',[0.80 0.04 0.15 0.08], 'Callback',@(~,~) delete(d));

    uiwait(d);
    if ~isvalid(d); error('Cancelled by user.'); end

    % -- Read values (no chained {} right after get!) ---------------------
    BlendPath      = string(get(hBlend,'String'));
    SavePath       = string(get(hSave ,'String'));

    dirList        = get(hDir,'String');
    LithDirection  = string(dirList{get(hDir,'Value')});

    refList        = get(hRef,'String');
    GraphiteSource = string(refList{get(hRef,'Value')});

    GammaSi        = str2double(get(hGamma,'String'));
    filterInputData= logical(get(hChkFilter,'Value'));
    smoothBadQS    = logical(get(hChkSmooth,'Value'));

    delete(d)
end

%% 05 RESOLVE PATHS
if strlength(GraphitePath)==0
    if strcmpi(GraphiteSource,'Schmitt') && strcmpi(LithDirection,'delithiation')
        warning('Schmitt reference only exists for lithiation – switching to Kuecher.');
        GraphiteSource = "Kuecher";
    end
    GraphitePath = fullfile(graphDir, ...
        sprintf('Gr_%s_%s.mat', capitalize(LithDirection), capitalize(GraphiteSource)));
end

assert(isfile(GraphitePath),'Graphite file not found.');
assert(isfile(BlendPath)   ,'Blend file not found.');
if GammaSi<=0 || GammaSi>=1 || isnan(GammaSi)
    error('\gamma_{Si} must be a number between 0 and 1.');
end
if mod(smoothWindow,2)==0
    warning('SmoothWindow must be odd – incremented by 1.');
    smoothWindow = smoothWindow + 1;
end

%% 06 LOAD & PRE-CLEAN
G  = loadOCV(GraphitePath);          
BL = loadBlend(BlendPath);           

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
QS = (QB - (1-GammaSi).*QG) ./ GammaSi;
QS(QS<0) = 0;   QS(QS>1) = 1;

% --- Smooth only where needed (no monotonic forcing) ---------------------
if smoothBadQS
    for i = 1 : 1
        QS = smoothdata(QS, 'movmedian', smoothWindow);
    end
end

siliconStruct.voltage            = V;

siliconStruct.normalizedCapacity = QS;

%% 09 SAVE RESULT
if strlength(SavePath)==0
    [file,path] = uiputfile({'*.mat','MAT-files (*.mat)'}, ...
                            'Save extracted Si curve as','SiCurve.mat');
    if isequal(file,0); SavePath=""; else; SavePath=fullfile(path,file); end
end
if strlength(SavePath)>0
    if ~endsWith(SavePath,'.mat','IgnoreCase',true)
        SavePath = SavePath + ".mat";
    end
    save(SavePath,'siliconStruct','-mat');
end

%% 10 PLOT
if PlotFlag
    tumBlue   = [  0 101 189]/255;
    tumOrange = [227 114  34]/255;

    figure('Name','Graphite · Silicon · Blend'); hold on; grid on; box on;
    plot(QG,V,'-','Color',tumBlue  ,'LineWidth',1.6, ...
        'DisplayName',['Graphite (' char(GraphiteSource) ')']);
    plot(QS,V,'-','Color',tumOrange,'LineWidth',1.6, ...
        'DisplayName','Silicon');
    plot(QB,V,'--','Color',[0 0 0],'LineWidth',1.6, ...
        'DisplayName','SiC-Blend');

    xlabel('normalised capacity $Q/Q_{\max}$','Interpreter','latex');
    ylabel('$U$ / V','Interpreter','latex');
    title(sprintf('Artificial silicon curve  (\\gamma_{Si}=%.2f)',GammaSi));
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

%% H4 trimAndRenorm --------------------------------------------------------
function out = trimAndRenorm(inStruct,Vmin,Vmax)
% Trim to [Vmin,Vmax] and rescale capacity to [0…1].
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
% Capitalise first letter; rest lower-case.
    s = char(s); 
    str = [upper(s(1)) lower(s(2:end))];
end

%% H7 iif ------------------------------------------------------------------
function y = iif(cond,a,b)
% Inline “if”: returns a if cond else b.
    if cond, y=a; else, y=b; end
end

%% H8 setFullPath ----------------------------------------------------------
function setFullPath(hEdit,dlgFcn,filterSpec,dlgTitle)
% Helper for Browse buttons → write full path into edit field.
    [file,path] = dlgFcn(filterSpec,dlgTitle);
    if isequal(file,0), return; end
    full = fullfile(path,file);
    if isequal(dlgFcn,@uiputfile) && ~endsWith(full,'.mat','IgnoreCase',true)
        full = full + ".mat";
    end
    set(hEdit,'String', full);
end
