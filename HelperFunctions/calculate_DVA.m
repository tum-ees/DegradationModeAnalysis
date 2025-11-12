function [Q_DVA, OCV_DVA, DVA] = calculate_DVA(Q, OCV, steps)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%
%   * This function calculates the DVA with the corresponding OCV and charge
%     vector
%   * size of the output vectors can be adjusted by setting steps
%
    if nargin < 3
        steps = 1001; % Default number of steps
    end
    
    % Remove NaN values
    validIdx = ~isnan(Q);
    Q = Q(validIdx);
    OCV = OCV(validIdx);
    
    % Ensure unique values for interpolation
    [Q, index] = unique(Q);
    
    % Interpolate OCV over a fixed number of points
    OCV_DVA = interp1(Q, OCV(index), linspace(min(Q), max(Q), steps));
    Q_DVA = linspace(min(Q), max(Q), steps);
    
    % Initialize DVA and calculate differential voltage analysis
    DVA = zeros(length(OCV_DVA), 1);
    for i = 2:length(DVA)
        dU = OCV_DVA(i) - OCV_DVA(i - 1);
        dQ = Q_DVA(i) - Q_DVA(i - 1);
        DVA(i - 1) = dU / dQ;
    end
    % last value will be set to 0 otherwise -> makes curves more smooth
    DVA(end) = DVA(end-1);

end
