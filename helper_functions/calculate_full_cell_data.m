function [halfAndFullCellData] = calculate_full_cell_data(fullCellData, halfAndFullCellData, settings)
%> Author: Moritz Guenthner (moritz.guenthner@tum.de), Josef Eizenhammer, (josef.eizenhammer@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-09-22
%
% This function preprocesses input data for full cell and interpolates the
% the data to the given data length
%
% INPUTS:
%   fullCellData         - raw full-cell curve source as struct, table-like
%                          payload, or path already resolved upstream
%   halfAndFullCellData  - struct that receives the processed full-cell data
%   settings             - DMA settings struct with smoothing and grid size
%
% OUTPUTS:
%   halfAndFullCellData  - input struct enriched with full-cell SOC, voltage,
%                          capacity, and precomputed DVA/ICA curves


% -----------------------------------------------------------------------
% parse the full-cell data (fcSOCRaw, fcURaw, capaAct) -> dataType='fullCell'
% -----------------------------------------------------------------------
fullCellData    = convertIfMat(fullCellData);

[fcSOCRaw, fcURaw, capaAct] = parse_data_input(fullCellData, 'fullCell');
fcUSmooth = smooth(fcURaw, settings.smoothingPoints, 'lowess');
fullCellSOC = linspace(0, 1, settings.dataLength);
try
    fullCellU = interp1(fcSOCRaw, fcUSmooth, fullCellSOC, 'linear','extrap');
catch    
    warning('Rare interp1 failure; enforcing unique, sorted SOC and retrying.');
    [xu, ia] = unique(fcSOCRaw(:));
    fullCellU = interp1(xu, fcUSmooth(ia), fullCellSOC, 'linear','extrap');
end

% q0 is the raw SOC span of the full-cell curve
q0 = max(fcSOCRaw) - min(fcSOCRaw);

halfAndFullCellData.fullCellSOC = fullCellSOC;
halfAndFullCellData.fullCellU   = fullCellU;
halfAndFullCellData.capaAct     = capaAct;
halfAndFullCellData.q0          = q0;

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
