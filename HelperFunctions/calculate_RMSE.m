function rmse = calculate_RMSE(measuredCurve, calculatedCurve, socRange)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-01-01
%> Last changes: 2025-03-12
%
% Calculates the Root Mean Square Error (RMSE) between a measured curve and a 
% calculated curve over a specified State-of-Charge (SOC) range.
%
% Requirements:
%   - measuredCurve and calculatedCurve must be vectors of the same length.
%   - The data points are assumed to be uniformly distributed in SOC from 0 to 1.
%   - socRange is an optional two-element vector [soc_min, soc_max] that specifies 
%     the SOC range over which to compute the RMSE. Default is [0, 1].
%
% Example:
%   rmse = calculateRMSE(measuredCurve, calculatedCurve, [0.2, 0.8]);
%

% Set default SOC range if not provided
if nargin < 3
    socRange = [0, 1];
end

% Ensure the curves are of the same length
if length(measuredCurve) ~= length(calculatedCurve)
    error('Curves must be of the same length to calculate RMSE.');
end

% Create a vector representing SOC values uniformly distributed between 0 and 1
socVec = linspace(0, 1, length(measuredCurve));

% Find indices corresponding to the specified SOC range
indices = (socVec >= socRange(1)) & (socVec <= socRange(2));

if sum(indices) == 0
    error('No data points found within the specified SOC range.');
end

% Calculate the differences only for the selected SOC range
differences = measuredCurve(indices) - calculatedCurve(indices);

% Square the differences
squaredDifferences = differences .^ 2;

% Calculate the mean of the squared differences
meanSquaredDifferences = mean(squaredDifferences);

% Take the square root of the mean squared differences to obtain the RMSE
rmse = sqrt(meanSquaredDifferences);

end
