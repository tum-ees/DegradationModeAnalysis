function [halfAndFullCellData] = calculate_half_cell_data(halfAndFullCellData, settings)
%> Author: Josef Eizenhammer, (josef.eizenhammer@tum.de), Moritz Guenthner (moritz.guenthner@tum.de) 
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-09-22
%
% Preprocess half cell data anode and cathode, interpolate to the given data
% length, and provide blend ready curves for anode and cathode.
%
% INPUTS:
%   halfAndFullCellData  - struct used as the shared half/full-cell container
%   settings             - DMA settings struct with file paths, blend flags,
%                          smoothing parameters, and interpolation length
%
% OUTPUTS:
%   halfAndFullCellData  - input struct enriched with normalized half-cell
%                          curves and blend interpolation grids
%

% -----------------------------------------------------------------------
% 1) Load raw data containers
% -----------------------------------------------------------------------
useCathodeBlend = settings.useCathodeBlend;
useAnodeBlend   = settings.useAnodeBlend;

anodeBlend1Raw = convertIfMat(settings.blendAn1DataPath);
anodeBlend2Raw = convertIfMat(settings.blendAn2DataPath);

if useCathodeBlend
    cathodeBlend1Raw = convertIfMat(settings.blendCa1DataPath);
    cathodeBlend2Raw = convertIfMat(settings.blendCa2DataPath);
else
    cathodeDataRaw = convertIfMat(settings.blendCa1DataPath);
end

% -----------------------------------------------------------------------
% 2) Cathode side handling (Blend1 as reference, optional Blend2)
% -----------------------------------------------------------------------
if useCathodeBlend
    % Blend1 and Blend2 cathode curves
    [blendCa1SOC, blendCa1U, ~] = parse_data_input(cathodeBlend1Raw, 'cathode');
    blendCa1U = smooth(blendCa1U, settings.smoothingPoints, 'lowess');

    [blendCa2SOC, blendCa2U, ~] = parse_data_input(cathodeBlend2Raw, 'cathodeBlend2');
    blendCa2U = smooth(blendCa2U, settings.smoothingPoints, 'lowess');

    % Common voltage window for cathode blends
    commonVoltageCa = linspace( ...
        max([min(blendCa2U), min(blendCa1U)]), ...
        min([max(blendCa2U), max(blendCa1U)]), ...
        settings.dataLength);
    try
        qCathodeBlend1Interp = interp1(blendCa1U, blendCa1SOC, commonVoltageCa, 'linear', 0);
        qCathodeBlend2Interp = interp1(blendCa2U, blendCa2SOC, commonVoltageCa, 'linear', 0);
    catch
        error("If you are using LFP, the error is very likely due to flat cathode potential. For blend electrodes neighboring voltage values must not be equal!")
    end
    
    % Normalized cathode OCV for the Blend1 reference
    cathodeSOCSingle = linspace(0, 1, settings.dataLength);
    cathodeUSingle   = interp1(blendCa1SOC, blendCa1U, cathodeSOCSingle, 'linear', 'extrap');

    normCathodeSOC    = cathodeSOCSingle;
    normCathodeU      = cathodeUSingle;
else
    % Single cathode curve (no blend2)
    [rawCatSOC, rawCatU, ~] = parse_data_input(cathodeDataRaw, 'cathode');
    catUSmooth             = smooth(rawCatU, settings.smoothingPoints, 'lowess');

    % SOC based single curve representation
    cathodeSOCSingle = linspace(0, 1, settings.dataLength);
    cathodeUSingle   = interp1(rawCatSOC, catUSmooth, cathodeSOCSingle, 'linear', 'extrap');

    normCathodeSOC    = cathodeSOCSingle;
    normCathodeU      = cathodeUSingle;

    % No blend in use
    qCathodeBlend1Interp = [];
    qCathodeBlend2Interp = [];
    commonVoltageCa         = [];
end

% -----------------------------------------------------------------------
% 3) Anode side handling  Blend1 for example graphite  Blend2 for example Si
%     Same pattern as cathode:
%     blend mode  common voltage grid  plus single curve from Blend1
%     non blend  single curve only
% -----------------------------------------------------------------------
if useAnodeBlend
    % Blend path: use common voltage grid to mix Blend1 and Blend2
    [blend1SOC, blend1U, ~] = parse_data_input(anodeBlend1Raw, 'anode');
    blend1U                 = smooth(blend1U, settings.smoothingPoints, 'lowess');

    [blend2SOC, blend2U, ~] = parse_data_input(anodeBlend2Raw, 'anodeBlend2');
    blend2U                 = smooth(blend2U, settings.smoothingPoints, 'lowess');

    % Common voltage window for anode blends
    commonVoltageAn = linspace( ...
        max([min(blend2U), min(blend1U)]), ...
        min([max(blend2U), max(blend1U)]), ...
        settings.dataLength);

    try
        qAnodeBlend2Interp = interp1(blend2U, blend2SOC, commonVoltageAn, 'linear', 0);
    catch
        % Fallback: drop duplicate U values so interp1 gets unique x
        [blend2UUnique, idxBlend2] = unique(blend2U, 'stable');
        qAnodeBlend2Interp = interp1(blend2UUnique, blend2SOC(idxBlend2), ...
            commonVoltageAn, 'linear', 0);
    end
    try
        qAnodeBlend1Interp = interp1(blend1U, blend1SOC, commonVoltageAn, 'linear', 0);
    catch
        % Fallback: drop duplicate U values so interp1 gets unique x
        [blend1UUnique, idxBlend1] = unique(blend1U, 'stable');
        qAnodeBlend1Interp = interp1(blend1UUnique, blend1SOC(idxBlend1), ...
            commonVoltageAn, 'linear', 0);
    end

    % Normalized single curve representation from Blend1
    anodeSOCSingle = linspace(0, 1, settings.dataLength);
    anodeUSingle   = interp1(blend1SOC, blend1U, anodeSOCSingle, 'linear', 'extrap');
else
    % Non blend path: direct SOC to U curve
    [rawAnSOC, rawAnU, ~] = parse_data_input(anodeBlend1Raw, 'anode');
    anUSmooth             = smooth(rawAnU, settings.smoothingPoints, 'lowess');

    % SOC based single curve representation
    anodeSOCSingle = linspace(0, 1, settings.dataLength);
    anodeUSingle   = interp1(rawAnSOC, anUSmooth, anodeSOCSingle, 'linear', 'extrap');

    % No blend in use
    commonVoltageAn       = [];
    qAnodeBlend1Interp = [];
    qAnodeBlend2Interp = [];
end

% For symmetry with cathode also provide normAnode fields as aliases
normAnodeSOC = anodeSOCSingle;
normAnodeU   = anodeUSingle;

% -----------------------------------------------------------------------
% 4) Store results in halfAndFullCellData struct
% -----------------------------------------------------------------------
% Cathode
halfAndFullCellData.normCathodeSOC         = normCathodeSOC;
halfAndFullCellData.normCathodeU           = normCathodeU;
halfAndFullCellData.cathodeSOCSingle      = cathodeSOCSingle;
halfAndFullCellData.cathodeUSingle        = cathodeUSingle;

halfAndFullCellData.commonVoltageCathode   = commonVoltageCa;
halfAndFullCellData.qCathodeBlend2Interp = qCathodeBlend2Interp;
halfAndFullCellData.qCathodeBlend1Interp = qCathodeBlend1Interp;

% Anode
halfAndFullCellData.normAnodeSOC           = normAnodeSOC;
halfAndFullCellData.normAnodeU             = normAnodeU;
halfAndFullCellData.anodeSOCSingle        = anodeSOCSingle;
halfAndFullCellData.anodeUSingle          = anodeUSingle;

halfAndFullCellData.commonVoltageAnode     = commonVoltageAn;
halfAndFullCellData.qAnodeBlend2Interp   = qAnodeBlend2Interp;
halfAndFullCellData.qAnodeBlend1Interp   = qAnodeBlend1Interp;

% -----------------------------------------------------------------------
% Nested helper convertIfMat
% -----------------------------------------------------------------------
    function dataOut = convertIfMat(dataIn)
        if ischar(dataIn) || isstring(dataIn)
            tmp = load(dataIn);
            fn  = fieldnames(tmp);
            dataOut = tmp.(fn{1});
        else
            dataOut = dataIn;
        end
    end

end
