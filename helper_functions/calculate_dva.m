function [qDVA, ocvDVA, dva] = calculate_dva(Q, OCV, steps)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%
%   * This function calculates the DVA with the corresponding OCV and charge
%     vector
%   * size of the output vectors can be adjusted by setting steps
%
% INPUTS:
%   Q      - charge or SOC vector used as interpolation axis
%   OCV    - voltage vector corresponding to Q
%   steps  - optional number of interpolation points (default: 1001)
%
% OUTPUTS:
%   qDVA   - uniformly sampled Q axis used for the derivative
%   ocvDVA - interpolated voltage on qDVA
%   dva    - differential voltage analysis curve dU/dQ
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
    ocvDVA = interp1(Q, OCV(index), linspace(min(Q), max(Q), steps));
    qDVA = linspace(min(Q), max(Q), steps);
    
    % Initialize DVA and calculate differential voltage analysis
    dva = zeros(length(ocvDVA), 1);
    for i = 2:length(dva)
        dU = ocvDVA(i) - ocvDVA(i - 1);
        dQ = qDVA(i) - qDVA(i - 1);
        dva(i - 1) = dU / dQ;
    end
    % last value will be set to 0 otherwise -> makes curves more smooth
    dva(end) = dva(end-1);

end
