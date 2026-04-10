function uMean = calculate_inhomogeneity(soc, u, inhomMax, inhomOffsetFraction)
%> Author: Moritz Guenthner (moritz.guenthner@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-11-24
%
% This function calculates the inhomogeneity effect on the OCV curve
% based on a Gaussian distribution of local SOCs around the mean SOC.
% Inhomogeneity is zero at 0 percent full-cell SOC and maximum at 100
% percent full-cell SOC (default behavior, inhomOffsetFraction = 0).
%
% Optional 4th argument inhomOffsetFraction (default 0):
%   If > 0, the inhomogeneity spread at SOC = 0 is
%   inhomOffsetFraction * maxSpread instead of zero.
%   At SOC = 1 the spread is always 100 percent.
%
% INPUTS:
%   soc                  - SOC grid on which the electrode curve is defined
%   u                    - electrode potential corresponding to soc
%   inhomMax             - maximum SOC spread used for the Gaussian average
%   inhomOffsetFraction  - optional fraction of the maximum spread already
%                          present at soc = 0
%
% OUTPUTS:
%   uMean                - inhomogeneity-averaged electrode potential

if nargin < 4 || isempty(inhomOffsetFraction)
    inhomOffsetFraction = 0;
end

if inhomMax <= 0
    uMean = u;
    return
end

% Precompute the normalized Gaussian support once.
persistent x mu lastSigma lastWeights
if isempty(x)
    x = linspace(0.5, 1.5, 61);
    mu = 1;
    lastSigma = NaN;
    lastWeights = [];
end

sigma = inhomMax;

% Recompute weights only if sigma changed.
if isempty(lastWeights) || abs(sigma - lastSigma) > 1e-8
    z = (x - mu) ./ sigma;
    weights = exp(-0.5 * z.^2);
    weights = weights / sum(weights);
    lastSigma = sigma;
    lastWeights = weights;
else
    weights = lastWeights;
end

soc = soc(:);
u = u(:);

try
    uMean = computeUMean(soc, u, x, weights, inhomOffsetFraction);
catch ME
    diagMsg = sprintf(['griddedInterpolant failed in calculate_inhomogeneity.\n', ...
        'SOC size: %s | U size: %s\n', ...
        'NaN/Inf SOC: %d | NaN/Inf U: %d\n', ...
        'Non-decreasing SOC: %d | Non-increasing SOC: %d\n', ...
        'Unique SOC count: %d | minSOC: %.6g | maxSOC: %.6g\n', ...
        'Example SOC head: %s\n', ...
        'Original error: %s'], ...
        mat2str(size(soc)), mat2str(size(u)), ...
        sum(~isfinite(soc)), sum(~isfinite(u)), ...
        all(diff(soc) >= 0), all(diff(soc) <= 0), ...
        numel(unique(soc(isfinite(soc)))), ...
        safeStat(@min, soc), safeStat(@max, soc), ...
        mat2str(soc(1:min(end, 5)).'), ...
        ME.message);

    newME = MException('calculate_inhomogeneity:InterpolantFailure', diagMsg);
    newME = addCause(newME, ME);
    throw(newME);
end

end

% -----------------------------------------------------------------------
% Local helper: computeUMean
% -----------------------------------------------------------------------
function uMean = computeUMean(soc, u, x, weights, inhomOffsetFraction)
% Build query grid.
% With inhomOffsetFraction = 0 (default), spread is proportional to SOC:
% zero spread at SOC = 0 and full spread at SOC = 1.
% With inhomOffsetFraction > 0, the spread at SOC = 0 starts at a finite
% fraction and grows linearly to full spread at SOC = 1.
xDev = x - 1;
alphaEff = inhomOffsetFraction + (1 - inhomOffsetFraction) .* soc;
xQuery = soc + alphaEff .* xDev;

% 'nearest' extrapolation clamps out-of-range queries to the nearest
% boundary value, which is appropriate for electrode OCV curves.
interpolant = griddedInterpolant(soc, u, 'linear', 'nearest');
electrodeOCVDist = interpolant(xQuery);

% Weighted average across columns.
uMean = electrodeOCVDist * weights(:);
end

% -----------------------------------------------------------------------
% Local helper: safeStat
% -----------------------------------------------------------------------
function val = safeStat(funHandle, vec)
if isempty(vec) || all(~isfinite(vec))
    val = NaN;
else
    val = funHandle(vec(isfinite(vec)));
end
end
