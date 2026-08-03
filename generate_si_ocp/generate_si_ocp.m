%% 01 HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author  : Mathias Rehm
%  E-mail  : mathias.rehm@tum.de
%  Date    : 2025-07-28
%
%  Function : generate_si_ocp
%  ------------------------------------------------------------------------
%  Generates an artificial Si half-cell OCV curve from
%     * a measured graphite-Si blend (blendPath)
%     * a pure graphite reference   (graphitePath | graphiteSource)
%
%  qBlend = gamma*qSi + (1-gamma)*qGr  ->  qSi = ...
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
%     'filterBlend'      true | false           (default false)
%                        LOWESS smoothing + deduplication on the measured
%                        blend curve.
%     'filterGraphite'   true | false           (default false)
%                        LOWESS smoothing + deduplication on the graphite
%                        reference.
%     'filterInputData'  deprecated alias. If set, applied to both curves.
%     'pavOutput'        true | false          (default true)
%                        Apply Pool Adjacent Violators (PAV) isotonic
%                        regression to the output silicon curve to enforce
%                        monotonicity.
%     'interpolationReadyOutput' true | false  (default false)
%                        Collapse the PAV plateaus (keep both endpoints,
%                        drop the interior) and resample the result onto a
%                        uniform, strictly increasing capacity grid using
%                        PCHIP. The table then has unique, strictly
%                        monotone capacity coordinates, preserves the
%                        overall voltage support and represents plateau
%                        regions as steep continuous ramps about one grid
%                        cell wide. Meant for consumers that read the
%                        table directly as a look-up table over capacity,
%                        such as PyBaMM. It does not smooth the DVA of a
%                        q(U) consumer.
%     'plotFlag'         true | false          (default true)
%
%  Returns
%     siliconStruct.voltage            [V]
%     siliconStruct.normalizedCapacity [0...1]
%
%  Optional second output, the graphite reference as it enters the
%  reconstruction, interpolated onto the common voltage grid:
%     graphiteStruct.voltage            [V]  (the common grid)
%     graphiteStruct.normalizedCapacity [0...1]
%     graphiteStruct.meta               source, direction, path, Vwindow
%  meta.source names the built-in reference the graphite curve was resolved
%  from, or 'custom' when 'graphitePath' was passed explicitly, because no
%  built-in reference is involved then.
%  With 'interpolationReadyOutput' true the silicon table is resampled onto
%  the uniform 1001-point capacity grid while graphiteStruct stays on the
%  common voltage grid, so the two outputs no longer share one grid.
%
%  Saving
%     A non-empty 'savePath' writes both structs to that file. An empty
%     'savePath' skips saving.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [siliconStruct, graphiteStruct] = generate_si_ocp(varargin)

%% 02 SETUP
set(groot,'defaultAxesTickLabelInterpreter','latex')   % LaTeX ticks axis-wide

% ---- Defaults -----------------------------------------------------------
blendPath         = "";
savePath          = "";
graphitePath      = "";
graphiteSource    = "Rehm2026";
lithDirection     = "lithiation";
gammaSi           = NaN;
filterBlend       = false;
filterGraphite    = false;
pavOutput         = true;
interpolationReadyOutput = false;
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

        case 'filterblend'
            filterBlend = parseLogicalFlag(val, keyStr);

        case 'filtergraphite'
            filterGraphite = parseLogicalFlag(val, keyStr);

        case 'filterinputdata'
            warning('generate_si_ocp:DeprecatedOption', ...
                ['''filterInputData'' is deprecated; use ''filterBlend'' ', ...
                 'and ''filterGraphite'' instead.']);
            filterInputData = parseLogicalFlag(val, keyStr);
            filterBlend    = filterInputData;
            filterGraphite = filterInputData;

        case 'pavoutput'
            pavOutput = parseLogicalFlag(val, keyStr);

        case 'interpolationreadyoutput'
            interpolationReadyOutput = parseLogicalFlag(val, keyStr);

        case 'plotflag'
            plotFlag = parseLogicalFlag(val, keyStr);

        otherwise
            warning('Unknown parameter "%s" ignored.', keyStr);
    end
end

if interpolationReadyOutput && ~pavOutput
    error('interpolationReadyOutput requires pavOutput=true.');
end


%% 03 PATHS
thisDir  = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
graphDir = fullfile(projRoot,'input_data','graphite');

% Add project tree to the MATLAB path if needed.
if ~contains(path, projRoot)
    addpath(genpath(projRoot));
end

%% 04 GUI
if strlength(blendPath)==0 || isnan(gammaSi)
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
        'Callback',@(~,~) set_full_path(hBlend,@uigetfile, ...
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
        'Callback',@(~,~) set_full_path(hSave,@uiputfile, ...
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
    hChkBlend = uicontrol(d,'Style','checkbox','Units','normalized', ...
        'Position',[0.02 y0 0.32 0.08], ...
        'String','Filter blend (LOWESS)', 'Value',filterBlend);
    hChkGraphite = uicontrol(d,'Style','checkbox','Units','normalized', ...
        'Position',[0.34 y0 0.32 0.08], ...
        'String','Filter graphite (LOWESS)', 'Value',filterGraphite);
    hChkPAV = uicontrol(d,'Style','checkbox','Units','normalized', ...
        'Position',[0.66 y0 0.32 0.08], ...
        'String','Monotone output (PAV)', 'Value',pavOutput);

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
    filterBlend     = logical(get(hChkBlend   ,'Value'));
    filterGraphite  = logical(get(hChkGraphite,'Value'));
    pavOutput       = logical(get(hChkPAV     ,'Value'));

    delete(d)
end

%% 05 RESOLVE PATHS
% Remember whether the graphite file came from one of the built-in
% references or straight from the caller, so the metadata does not claim a
% source that was never used.
graphiteFromSource = strlength(graphitePath)==0;
if graphiteFromSource
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
graphiteData = load_ocv(graphitePath);
blendData = load_blend(blendPath);

graphiteData = smoothUnique(graphiteData, filterGraphite);
blendData = smoothUnique(blendData, filterBlend);

%% 07 ALIGN TO COMMON VOLTAGE WINDOW
voltageMin = max([min(graphiteData.voltage), min(blendData.voltage)]);
voltageMax = min([max(graphiteData.voltage), max(blendData.voltage)]);

graphiteData = trimAndRenorm(graphiteData, voltageMin, voltageMax);
blendData = trimAndRenorm(blendData, voltageMin, voltageMax);

numPoints = max(numel(graphiteData.voltage), numel(blendData.voltage));
voltageCommon = linspace(voltageMin, voltageMax, numPoints).';

qGraphite = interp1(graphiteData.voltage, graphiteData.normalizedCapacity, voltageCommon, 'linear');
qBlend = interp1(blendData.voltage, blendData.normalizedCapacity, voltageCommon, 'linear');

assert(all(isfinite(qGraphite)) && all(isfinite(qBlend)), ...
    'Common voltage grid is outside the aligned input support.');

%% 08 CALCULATE SILICON CURVE
qSilicon = (qBlend - (1 - gammaSi) .* qGraphite) ./ gammaSi;
qSilicon(qSilicon < 0) = 0;
qSilicon(qSilicon > 1) = 1;

% --- Optionally enforce monotonicity on output via PAV --------------------
if pavOutput
    if qSilicon(end) < qSilicon(1)
        qSilicon = pavIsotonic(qSilicon, 'nonincreasing');
    else
        qSilicon = pavIsotonic(qSilicon, 'nondecreasing');
    end
end

outputVoltage = voltageCommon;
outputCapacity = qSilicon;
if interpolationReadyOutput
    [outputCapacity, outputVoltage] = makeInterpolationReady( ...
        voltageCommon, qSilicon);
end

siliconStruct.voltage            = outputVoltage;
siliconStruct.normalizedCapacity = outputCapacity;

% The graphite reference on the common grid, so callers can reuse the
% aligned reconstruction basis instead of repeating the trimming and the
% interpolation.
graphiteStruct.voltage            = voltageCommon;
graphiteStruct.normalizedCapacity = qGraphite;
if graphiteFromSource
    graphiteStruct.meta.source    = char(graphiteSource);
else
    graphiteStruct.meta.source    = 'custom';
end
graphiteStruct.meta.direction     = char(lithDirection);
graphiteStruct.meta.path          = char(graphitePath);
graphiteStruct.meta.Vwindow       = [voltageMin voltageMax];

%% 09 SAVE RESULT
% An empty savePath means: do not save.
if strlength(savePath)>0
    if ~endsWith(savePath,'.mat','IgnoreCase',true)
        savePath = savePath + ".mat";
    end
    save(savePath,'siliconStruct','graphiteStruct','-mat');
end

%% 10 PLOT
if plotFlag
    tumBlue   = [  0 101 189]/255;
    tumOrange = [227 114  34]/255;

    figure('Name','Graphite / Silicon / Blend'); hold on; grid on; box on;
    plot(qGraphite, voltageCommon, '-','Color',tumBlue  ,'LineWidth',1.6, ...
        'DisplayName',['Graphite (' char(graphiteSource) ')']);
    plot(outputCapacity, outputVoltage, '-','Color',tumOrange,'LineWidth',1.6, ...
        'DisplayName','Silicon');
    plot(qBlend, voltageCommon, '--','Color',[0 0 0],'LineWidth',1.6, ...
        'DisplayName','SiGr-Blend');

    xlabel('normalized capacity $Q/Q_{\max}$','Interpreter','latex');
    ylabel('$U$ / V','Interpreter','latex');
    title(sprintf('Artificial silicon curve  (\\gamma_{Si}=%.2f)',gammaSi));
    legend('Interpreter','latex','Location','best');
    ylim([-0.05 0.85]);  xlim([-0.05 1.05]);
end

end   % ------------------- end main function ------------------------------


%% H1 load_ocv --------------------------------------------------------------
function S = load_ocv(matPath)
% Load first struct in a .mat file; verify required fields exist.
    raw = load(matPath);
    fn  = fieldnames(raw);
    S   = raw.(fn{1});
    assert(all(isfield(S,{'voltage','normalizedCapacity'})), ...
        'OCV file must contain: voltage & normalizedCapacity.');
end

%% H2 load_blend ------------------------------------------------------------
function OUT = load_blend(matPath)
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

%% H3 smoothUnique ----------------------------------------------------------
function out = smoothUnique(inStruct, doSmooth)
% LOWESS (window 30) on voltage + remove duplicate voltages.
    if doSmooth
        inStruct.voltage = smooth(inStruct.voltage,30,'lowess');
    end
    [uVolt,idx] = unique(inStruct.voltage);
    out.voltage            = uVolt;
    out.normalizedCapacity = inStruct.normalizedCapacity(idx);
end

%% H3b pavIsotonic ----------------------------------------------------------
function qMono = pavIsotonic(qRaw, direction)
% Pool Adjacent Violators (PAV) isotonic regression.
% Merges adjacent monotonicity-violating blocks by replacing them with
% their weighted mean.
% Finds qMono minimizing sum((qMono - qRaw)^2) subject to monotonicity.
%   direction : 'nondecreasing' (default) | 'nonincreasing'
    if nargin < 2, direction = 'nondecreasing'; end

    q = qRaw(:);
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
    qMono = zeros(n,1);
    idx = 1;
    for k = 1:nb
        qMono(idx : idx+sz(k)-1) = val(k);
        idx = idx + sz(k);
    end

    if strcmpi(direction, 'nonincreasing')
        qMono = -qMono;
    end
end

%% H3c makeInterpolationReady -----------------------------------------------
function [qGrid, voltageGrid] = makeInterpolationReady(voltage, capacity)
% Resample U(q) onto a uniform, strictly increasing capacity grid.
%
% PAV pools monotonicity violators onto their common mean, so the isotonic
% output carries runs of exactly equal capacity (plateaus). A plateau is a
% vertical segment of U(q): several samples share one capacity coordinate
% while their voltages differ. Discarding plateaus with unique(...,'stable')
% keeps only the first sample of every run and therefore throws the
% plateau's voltage extent away.
%
% collapsePlateaus keeps BOTH plateau endpoints and separates them by a
% tiny shift instead, so the capacity coordinates become unique and
% strictly monotone while the overall voltage support of the curve is
% preserved exactly. It serves the same purpose as the strict_sto export in
% PyDMA. In the resampled table a plateau region shows up as a steep but
% continuous ramp roughly one output grid cell wide, and two plateaus that
% fall into the same grid cell merge into a single ramp. The mode is meant
% for consumers that read the table directly as a look-up table over
% capacity, such as PyBaMM, where it avoids the voltage jump that several
% samples stacked on one capacity coordinate produce. It does not smooth
% the curve: a consumer that differentiates q(U) still sees the plateau
% edges of the PAV result.
    voltage  = voltage(:);
    capacity = capacity(:);

    if max(capacity) <= min(capacity)
        error(['Interpolation-ready output requires a non-degenerate ', ...
               'capacity range, but the PAV curve is constant.']);
    end

    % Collapse on the native (voltage-ordered) samples. collapsePlateaus
    % detects the monotonicity direction itself and needs the plateau
    % endpoints in their original order; sorting by capacity beforehand
    % would reverse the voltage order inside each plateau and fold U(q).
    [voltageCollapsed, capacityCollapsed] = collapsePlateaus(voltage, capacity);

    % The collapse only drops plateau interiors, so the first and the last
    % sample of the input must survive it unchanged.
    assert(voltageCollapsed(1) == voltage(1) && ...
           voltageCollapsed(end) == voltage(end), ...
        'Collapsing the plateaus changed the voltage support of the curve.');

    [qSorted, order] = sort(capacityCollapsed, 'ascend');
    voltageSorted = voltageCollapsed(order);

    if numel(qSorted) < 2
        error('Interpolation-ready output requires at least two capacity values.');
    end

    % No deduplication here. A residual tie means collapsePlateaus failed to
    % separate two samples, and silently dropping one of them would discard
    % its voltage together with a part of the voltage support, which is
    % exactly the failure mode this output mode has to avoid.
    tieIdx = find(diff(qSorted) <= 0, 1, 'first');
    if ~isempty(tieIdx)
        error(['Collapsed capacity is not strictly increasing: the sorted ', ...
               'samples %d and %d both sit at q = %.17g (voltages %.17g ', ...
               'and %.17g).'], tieIdx, tieIdx+1, qSorted(tieIdx), ...
               voltageSorted(tieIdx), voltageSorted(tieIdx+1));
    end

    qGrid = linspace(qSorted(1), qSorted(end), 1001).';
    voltageGrid = interp1(qSorted, voltageSorted, qGrid, 'pchip');

    assert(all(isfinite(voltageGrid)), 'PCHIP output contains non-finite values.');
    assert(all(diff(qGrid) > 0), 'Output capacity grid is not strictly increasing.');
    assert(min(diff(qGrid)) >= 1e-5, ...
        'Output capacity spacing is below the required minimum of 1e-5.');
    assert(all(diff(voltageGrid) > 0) || all(diff(voltageGrid) < 0), ...
        'PCHIP voltage output is not strictly monotone.');
end

%% H3c2 collapsePlateaus ----------------------------------------------------
function [voltageOut, capacityOut] = collapsePlateaus(voltage, capacity, epsShift)
% Replace plateaus in CAPACITY (runs of equal values) by their two voltage
% endpoints, separated by a tiny shift so the result is strictly monotone.
% (The shift argument is named epsShift rather than eps so it does not
%  shadow the MATLAB builtin eps inside this function.)
%
% For every run of length L >= 2 at level q:
%   - keep only the first and last index of the run (drop the interior),
%   - move the two kept samples to q-s and q+s (reversed for a
%     non-increasing curve).
% Samples that are not part of a plateau keep their exact capacity, so
% linear interpolation across the original plateau range is preserved to
% within s and voltage -> capacity consumers see the same curve.
%
% The shift s is bounded twice, and both bounds are what keeps the output
% strictly monotone:
%   - s <= epsShift, the global tie-breaking scale, and
%   - s <= a quarter of the distance to the neighbouring plateau levels, so
%     the shifted endpoints of two adjacent runs can never meet, however
%     close the two levels are.
% Levels that are closer than a few ulps are pooled into one run before the
% shift is computed, because at that distance a quarter of the gap is no
% longer representable (see below).
% A run that sits exactly on the boundary of the capacity range is shifted
% INWARD only: the boundary sample keeps its exact value and its partner
% moves into the range. The output therefore stays inside the input range
% by construction and never has to be clamped back. Clamping was the origin
% of a silent data loss: it pushed shifted samples back onto the boundary
% value, re-created exact ties there and let the caller drop the tied
% samples together with their voltage support.
%
%   voltage  : voltage samples, ordered together with capacity
%   capacity : monotone (post-PAV) capacity, possibly containing plateaus
%   epsShift : optional tie-breaking shift; default 1e-5 of the capacity
%              range (tiny, yet robust under downstream cubic interpolation)
    voltage  = voltage(:);
    capacity = capacity(:);
    n = numel(capacity);
    if n < 2
        voltageOut = voltage;
        capacityOut = capacity;
        return
    end

    if capacity(end) >= capacity(1)
        direction = 1;
    else
        direction = -1;
    end

    qRange = max(capacity) - min(capacity);
    if nargin < 3 || isempty(epsShift)
        if qRange > 0
            epsShift = qRange * 1e-5;
        else
            epsShift = 1e-9;
        end
    end
    epsShift = max(double(epsShift), eps('double'));

    % PAV-pooled means can drift in floating point so adjacent buckets that
    % should have merged end up slightly violating monotonicity. The
    % cumulative min/max snap forces (non-)monotonicity exactly, so the
    % plateau detection below works on a clean signal. Flipping the sign for
    % a falling curve lets both directions share one code path: qWork is
    % non-decreasing, and the mirrored result is flipped back at the end.
    if direction > 0
        qWork = cummax(capacity);
    else
        qWork = -cummin(capacity);
    end
    qWorkMin = qWork(1);
    qWorkMax = qWork(end);

    % Runs of equal capacity. Their levels are strictly increasing.
    runStart = [1; find(diff(qWork) > 0) + 1];
    runEnd   = [runStart(2:end) - 1; n];
    levels   = qWork(runStart);

    % Pool levels that are numerically indistinguishable. The shift applied
    % below is a quarter of the distance to the neighbouring level, so with
    % a gap under 4 ulps of the level magnitude the shifted endpoint rounds
    % straight back onto its own level and the pair ties again. Requiring
    % 8 ulps keeps every shifted endpoint at least 2 ulps away from its
    % level under round-to-nearest, on both sides of the run. Levels that
    % close are the same level at machine precision, so their runs are
    % merged and the merged run keeps the first and the last sample of the
    % whole group. This is the same kind of floating point pooling as the
    % cumulative min/max snap above.
    keepRun = true(numel(levels), 1);
    reference = levels(1);
    for r = 2:numel(levels)
        if levels(r) - reference <= 8*max(eps(levels(r)), eps(reference))
            keepRun(r) = false;
        else
            reference = levels(r);
        end
    end
    if ~all(keepRun)
        kept = find(keepRun);
        mergedEnd = zeros(numel(kept), 1);
        mergedEnd(1:end-1) = runEnd(kept(2:end) - 1);
        mergedEnd(end) = runEnd(end);
        mergedLevels = levels(kept);
        % The first group starts at the lower end of the range anyway. The
        % last group is represented by the upper end instead of by its own
        % first level, so both boundary values stay exact and the
        % inward-only shift below still recognises them.
        mergedLevels(end) = levels(end);
        runStart = runStart(kept);
        runEnd   = mergedEnd;
        levels   = mergedLevels;
    end
    numRuns  = numel(levels);

    % Shift budget per run: the global scale, capped at a quarter of the
    % distance to either neighbouring level.
    shift = repmat(epsShift, numRuns, 1);
    if numRuns > 1
        gaps = diff(levels);
        shift(1:end-1) = min(shift(1:end-1), 0.25*gaps);
        shift(2:end)   = min(shift(2:end)  , 0.25*gaps);
    end

    % Keep the first and the last sample of every run and separate them.
    voltageOut  = zeros(2*numRuns, 1);
    capacityOut = zeros(2*numRuns, 1);
    p = 0;
    for r = 1:numRuns
        q = levels(r);
        s = shift(r);

        p = p + 1;
        voltageOut(p) = voltage(runStart(r));
        if runEnd(r) == runStart(r)
            capacityOut(p) = q;      % isolated sample: keep it untouched
            continue
        end

        if q == qWorkMin
            capacityOut(p) = q;      % boundary run: shift inward only
        else
            capacityOut(p) = q - s;
        end

        p = p + 1;
        voltageOut(p) = voltage(runEnd(r));
        if q == qWorkMax
            capacityOut(p) = q;      % boundary run: shift inward only
        else
            capacityOut(p) = q + s;
        end
    end
    voltageOut  = voltageOut(1:p);
    capacityOut = capacityOut(1:p);

    % Safety net. After the pooling above every shift is representable and
    % neighbouring output values differ by at least two ulps, so this sweep
    % is not expected to change anything. It steps by a single ulp and the
    % result is not clamped back into the range: a clamp would re-create
    % exactly the ties that this function exists to prevent.
    for k = 2:p
        if capacityOut(k) <= capacityOut(k-1)
            capacityOut(k) = capacityOut(k-1) + eps(capacityOut(k-1));
        end
    end

    % Contract of this function, guaranteed by construction rather than
    % repaired: strictly monotone capacity that stays inside the range of
    % the snapped input curve.
    assert(all(diff(capacityOut) > 0), ...
        'Collapsed capacity is not strictly monotone.');
    assert(capacityOut(1) >= qWorkMin && capacityOut(end) <= qWorkMax, ...
        'Collapsed capacity left the range of the input curve.');

    capacityOut = direction * capacityOut;
end

%% H3d parseLogicalFlag -----------------------------------------------------
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

%% H4 trimAndRenorm ---------------------------------------------------------
function out = trimAndRenorm(inStruct, Vmin, Vmax)
% Interpolate exact overlap boundaries, retain interior samples, and
% rescale capacity to [0...1]. This keeps the subsequent common grid fully
% inside the aligned support instead of filling dropped edge samples with 0.
    voltage = inStruct.voltage(:);
    capacity = inStruct.normalizedCapacity(:);
    [voltage, order] = sort(voltage, 'ascend');
    capacity = capacity(order);
    [voltage, uniqueIdx] = unique(voltage, 'stable');
    capacity = capacity(uniqueIdx);

    boundaryCapacity = interp1(voltage, capacity, [Vmin; Vmax], 'linear');
    assert(all(isfinite(boundaryCapacity)), ...
        'Overlap boundaries are outside the input support.');

    interior = voltage > Vmin & voltage < Vmax;
    out.voltage = [Vmin; voltage(interior); Vmax];
    out.normalizedCapacity = [boundaryCapacity(1); ...
        capacity(interior); boundaryCapacity(2)];
    out.normalizedCapacity = rescale(out.normalizedCapacity,0,1);
end

%% H5 capitalize -----------------------------------------------------------
function str = capitalize(s)
% Capitalize first letter; rest lower-case.
    s = char(s);
    str = [upper(s(1)) lower(s(2:end))];
end

%% H6 iif ------------------------------------------------------------------
function y = iif(cond,a,b)
% Inline if: returns a if cond else b.
    if cond, y=a; else, y=b; end
end

%% H7 set_full_path ----------------------------------------------------------
function set_full_path(hEdit,dlgFcn,filterSpec,dlgTitle)
% Helper for Browse buttons -> write full path into edit field.
    [file,path] = dlgFcn(filterSpec,dlgTitle);
    if isequal(file,0), return; end
    full = fullfile(path,file);
    if isequal(dlgFcn,@uiputfile) && ~endsWith(full,'.mat','IgnoreCase',true)
        full = full + ".mat";
    end
    set(hEdit,'String', full);
end
