function [half_and_full_cell_data] = calculate_half_cell_data(half_and_full_cell_data, settings)
%> Author: Josef Eizenhammer, (josef.eizenhammer@tum.de), Moritz Guenthner (moritz.guenthner@tum.de) 
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-09-22
%
% Preprocess half cell data anode and cathode, interpolate to the given data
% length, and provide blend ready curves for anode and cathode.
%

% -----------------------------------------------------------------------
% 1) Load raw data containers
% -----------------------------------------------------------------------
useCathodeBlend = settings.useCathodeBlend;
useAnodeBlend   = settings.useAnodeBlend;

anodeBlend1_raw = convertIfMat(settings.blendAn1Data_path);
anodeBlend2_raw = convertIfMat(settings.blendAn2Data_path);

if useCathodeBlend
    cathBlend1_raw = convertIfMat(settings.blendCa1Data_path);
    cathBlend2_raw = convertIfMat(settings.blendCa2Data_path);
else
    cathodeData_raw = convertIfMat(settings.blendCa1Data_path);
end

% -----------------------------------------------------------------------
% 2) Cathode side handling (Blend1 as reference, optional Blend2)
% -----------------------------------------------------------------------
if useCathodeBlend
    % Blend1 and Blend2 cathode curves
    [BlendCa1SOC, BlendCa1U, ~] = parse_data_input(cathBlend1_raw, 'cathode');
    BlendCa1U = smooth(BlendCa1U, settings.smoothingPoints, 'lowess');

    [BlendCa2SOC, BlendCa2U, ~] = parse_data_input(cathBlend2_raw, 'cathodeBlend2');
    BlendCa2U = smooth(BlendCa2U, settings.smoothingPoints, 'lowess');

    % Common voltage window for cathode blends
    commonVoltageCa = linspace( ...
        max([min(BlendCa2U), min(BlendCa1U)]), ...
        min([max(BlendCa2U), max(BlendCa1U)]), ...
        settings.dataLength);
    try
        Q_cathode_blend1_interp = interp1(BlendCa1U, BlendCa1SOC, commonVoltageCa, 'linear', 0);
        Q_cathode_blend2_interp = interp1(BlendCa2U, BlendCa2SOC, commonVoltageCa, 'linear', 0);
    catch
        error("If you are using LFP, the error is very likely due to flat cathode potential. For blend electrodes neighboring voltage values must not be equal!")
    end
    
    % Normalized cathode OCV for the Blend1 reference
    cathode_SOC_single = linspace(0, 1, settings.dataLength);
    cathode_U_single   = interp1(BlendCa1SOC, BlendCa1U, cathode_SOC_single, 'linear', 'extrap');

    normCathode_SOC    = cathode_SOC_single;
    normCathode_U      = cathode_U_single;
else
    % Single cathode curve (no blend2)
    [rawCatSOC, rawCatU, ~] = parse_data_input(cathodeData_raw, 'cathode');
    catU_smooth            = smooth(rawCatU, settings.smoothingPoints, 'lowess');

    % SOC based single curve representation
    cathode_SOC_single = linspace(0, 1, settings.dataLength);
    cathode_U_single   = interp1(rawCatSOC, catU_smooth, cathode_SOC_single, 'linear', 'extrap');

    normCathode_SOC    = cathode_SOC_single;
    normCathode_U      = cathode_U_single;

    % No blend in use
    Q_cathode_blend1_interp = [];
    Q_cathode_blend2_interp = [];
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
    [Blend1SOC, Blend1U, ~] = parse_data_input(anodeBlend1_raw, 'anode');
    Blend1U                 = smooth(Blend1U, settings.smoothingPoints, 'lowess');

    [Blend2SOC, Blend2U, ~] = parse_data_input(anodeBlend2_raw, 'anodeBlend2');
    Blend2U                 = smooth(Blend2U, settings.smoothingPoints, 'lowess');

    % Common voltage window for anode blends
    commonVoltageAn = linspace( ...
        max([min(Blend2U), min(Blend1U)]), ...
        min([max(Blend2U), max(Blend1U)]), ...
        settings.dataLength);

    Q_anode_blend2_interp = interp1(Blend2U, Blend2SOC, commonVoltageAn, 'linear', 0);
    Q_anode_blend1_interp = interp1(Blend1U, Blend1SOC, commonVoltageAn, 'linear', 0);

    % Normalized single curve representation from Blend1
    anode_SOC_single = linspace(0, 1, settings.dataLength);
    anode_U_single   = interp1(Blend1SOC, Blend1U, anode_SOC_single, 'linear', 'extrap');
else
    % Non blend path: direct SOC to U curve
    [rawAnSOC, rawAnU, ~] = parse_data_input(anodeBlend1_raw, 'anode');
    anU_smooth            = smooth(rawAnU, settings.smoothingPoints, 'lowess');

    % SOC based single curve representation
    anode_SOC_single = linspace(0, 1, settings.dataLength);
    anode_U_single   = interp1(rawAnSOC, anU_smooth, anode_SOC_single, 'linear', 'extrap');

    % No blend in use
    commonVoltageAn       = [];
    Q_anode_blend1_interp = [];
    Q_anode_blend2_interp = [];
end

% For symmetry with cathode also provide normAnode fields as aliases
normAnode_SOC = anode_SOC_single;
normAnode_U   = anode_U_single;

% -----------------------------------------------------------------------
% 4) Store results in half_and_full_cell_data struct
% -----------------------------------------------------------------------
% Cathode
half_and_full_cell_data.normCathode_SOC         = normCathode_SOC;
half_and_full_cell_data.normCathode_U           = normCathode_U;
half_and_full_cell_data.cathode_SOC_single      = cathode_SOC_single;
half_and_full_cell_data.cathode_U_single        = cathode_U_single;

half_and_full_cell_data.commonVoltage_cathode   = commonVoltageCa;
half_and_full_cell_data.Q_cathode_blend2_interp = Q_cathode_blend2_interp;
half_and_full_cell_data.Q_cathode_blend1_interp = Q_cathode_blend1_interp;

% Anode
half_and_full_cell_data.normAnode_SOC           = normAnode_SOC;
half_and_full_cell_data.normAnode_U             = normAnode_U;
half_and_full_cell_data.anode_SOC_single        = anode_SOC_single;
half_and_full_cell_data.anode_U_single          = anode_U_single;

half_and_full_cell_data.commonVoltage_anode     = commonVoltageAn;
half_and_full_cell_data.Q_anode_blend2_interp   = Q_anode_blend2_interp;
half_and_full_cell_data.Q_anode_blend1_interp   = Q_anode_blend1_interp;

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
