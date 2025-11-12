function [half_and_full_cell_data] = calculate_half_cell_data(half_and_full_cell_data, settings)
%> Author: Josef Eizenhammer, (josef.eizenhammer@tum.de), Moritz Guenthner (moritz.guenthner@tum.de) 
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-09-22
%
% This function preprocesses input data for half cell and interpolates the
% the data to the given data length
% 

% -----------------------------------------------------------------------
% 1) Prepare data for processing
% -----------------------------------------------------------------------

anodeBlend1 = settings.Blend1Data_path;
anodeBlend2 = settings.Blend2Data_path;
cathodeData = settings.cathodeData_path;

cathodeData     = convertIfMat(cathodeData);
anodeBlend2     = convertIfMat(anodeBlend2);
anodeBlend1     = convertIfMat(anodeBlend1);

Umin = max([min(anodeBlend2.voltage) , min(anodeBlend1.voltage)]);
Umax = min([max(anodeBlend2.voltage) , max(anodeBlend1.voltage)]);

anodeBlend2  = norm_blend_voltage(anodeBlend2 ,Umin,Umax);
anodeBlend1 = norm_blend_voltage(anodeBlend1,Umin,Umax);

% -----------------------------------------------------------------------
% 2) Parse the cathode data (rawCatSOC, rawCatU, etc.) -> dataType='cathode'
% -----------------------------------------------------------------------
[rawCatSOC, rawCatU, ~] = parse_data_input(cathodeData, 'cathode');
catU_smooth       = smooth(rawCatU, settings.smoothingPoints, 'lowess');
normCathode_SOC   = linspace(0, 1, settings.dataLength);
normCathode_U     = interp1(rawCatSOC, catU_smooth, normCathode_SOC, 'linear','extrap');

% -----------------------------------------------------------------------
% 3) Parse the anode: Blend1 (e.g graphite) & Blend2 data (always via parse_data_input)
% -----------------------------------------------------------------------
% ----- Blend2 (e.g. silicon) -----
[Blend2SOC, Blend2U, ~]  = parse_data_input(anodeBlend2, 'anodeBlend2');
Blend2U                  = smooth(Blend2U, settings.smoothingPoints, 'lowess');
[uniqueBlend2U, iBlend2] = unique(Blend2U);
rawBlend2Voltage         = uniqueBlend2U;
rawBlend2Capacity        = Blend2SOC(iBlend2);

% ----- Blend1 (e.g. graphite) -----
[Blend1SOC, Blend1U, ~] = parse_data_input(anodeBlend1, 'anode');
Blend1U                 = smooth(Blend1U, settings.smoothingPoints, 'lowess');
[uniqueBlend1U, iBlend1]= unique(Blend1U);
rawBlend1Voltage        = uniqueBlend1U;
rawBlend1Capacity       = Blend1SOC(iBlend1);

% define a common voltage axis for blending (normed to dataLength points)
commonVoltage = linspace( ...
    max([min(rawBlend2Voltage), min(rawBlend1Voltage)]), ...
    min([max(rawBlend2Voltage), max(rawBlend1Voltage)]), ...
    settings.dataLength);

Q_Blend2_interp = interp1(rawBlend2Voltage, rawBlend2Capacity, commonVoltage, 'linear', 0);
Q_Blend1_interp = interp1(rawBlend1Voltage, rawBlend1Capacity, commonVoltage, 'linear', 0);

half_and_full_cell_data.normCathode_SOC = normCathode_SOC;
half_and_full_cell_data.normCathode_U   = normCathode_U;
half_and_full_cell_data.commonVoltage   = commonVoltage;
half_and_full_cell_data.Q_Blend2_interp = Q_Blend2_interp;
half_and_full_cell_data.Q_Blend1_interp = Q_Blend1_interp;

% -----------------------------------------------------------------------
% Nested helper: convertIfMat
%   Checks if the input is a string/char path; if so, loads it. Otherwise
%   returns as-is. This ensures everything is eventually a struct/cell
%   and can be parsed by parse_data_input.
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

