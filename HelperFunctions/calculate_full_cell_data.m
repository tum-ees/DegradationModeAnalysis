function [half_and_full_cell_data] = calculate_full_cell_data(fullCellData, half_and_full_cell_data, settings)
%> Author: Moritz Guenthner (moritz.guenthner@tum.de), Josef Eizenhammer, (josef.eizenhammer@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-09-22
%
% This function preprocesses input data for full cell and interpolates the
% the data to the given data length
% 


% -----------------------------------------------------------------------
% parse the full-cell data (fcSOC_raw, fcU_raw, capa_act) -> dataType='fullcell'
% -----------------------------------------------------------------------
fullCellData    = convertIfMat(fullCellData);

[fcSOC_raw, fcU_raw, capa_act] = parse_data_input(fullCellData, 'fullcell');
fcU_smooth   = smooth(fcU_raw, settings.smoothingPoints, 'lowess');
fullCell_SOC = linspace(0, 1, settings.dataLength);
try
    fullCell_U   = interp1(fcSOC_raw, fcU_smooth, fullCell_SOC, 'linear','extrap');
catch    
    warning('Rare interp1 failure; enforcing unique, sorted SOC and retrying.');
    [xu,ia] = unique(fcSOC_raw(:)); 
    fullCell_U = interp1(xu, fcU_smooth(ia), fullCell_SOC, 'linear','extrap');
end

% For Q0, we measure the raw SOC range
Q0 = max(fcSOC_raw) - min(fcSOC_raw);

half_and_full_cell_data.fullCell_SOC      = fullCell_SOC;
half_and_full_cell_data.fullCell_U      = fullCell_U;
half_and_full_cell_data.capa_act        = capa_act;
half_and_full_cell_data.Q0              = Q0;

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
