function U_mean = calculate_inhomogeneity(SOC, U, inhom_max, inhom_offset_fraction)
%> Author: Moritz Guenthner (moritz.guenthner@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-11-24
%
% This function calculates the inhomogeneity effect on the OCV curve
% based on a Gaussian distribution of local SOCs around the mean SOC.
% Inhomogeneity is zero at 0 percent full cell SOC and maximum at 100 percent full
% cell SOC (default behavior, inhom_offset_fraction = 0).
%
% Optional 4th argument inhom_offset_fraction (default 0):
%   If > 0, the inhomogeneity spread at SOC=0 is inhom_offset_fraction * max_spread
%   instead of zero. E.g. 0.5 means 50% of max inhomogeneity already at SOC=0.
%   At SOC=1 the spread is always 100% (unchanged).

    if nargin < 4 || isempty(inhom_offset_fraction)
        inhom_offset_fraction = 0;
    end

    debugInhom = true;

    if inhom_max <= 0
        U_mean = U;
        return
    end

    % Precompute x once
    persistent x mu last_sigma last_w
    if isempty(x)
        x = linspace(0.5, 1.5, 61);
        mu = 1;
        last_sigma = NaN;
        last_w = [];
    end

    sigma = inhom_max;

    % Recompute weights only if sigma changed
    if isempty(last_w) || abs(sigma - last_sigma) > 1e-8
        z = (x - mu) ./ sigma;
        w = exp(-0.5 * z.^2);
        w = w / sum(w);
        last_sigma = sigma;
        last_w = w;
    else
        w = last_w;
    end

    SOC = SOC(:);
    U   = U(:);

    if ~debugInhom
        U_mean = computeUmean(SOC, U, x, w, inhom_offset_fraction);
        return;
    end

    try
        U_mean = computeUmean(SOC, U, x, w, inhom_offset_fraction);
    catch ME
        diagMsg = sprintf(['griddedInterpolant failed in calculate_inhomogeneity (debug mode).\n', ...
            'SOC size: %s | U size: %s\n', ...
            'NaN/Inf SOC: %d | NaN/Inf U: %d\n', ...
            'Non-decreasing SOC: %d | Non-increasing SOC: %d\n', ...
            'Unique SOC count: %d | minSOC: %.6g | maxSOC: %.6g\n', ...
            'Example SOC head: %s\n', ...
            'Original error: %s'], ...
            mat2str(size(SOC)), mat2str(size(U)), ...
            sum(~isfinite(SOC)), sum(~isfinite(U)), ...
            all(diff(SOC) >= 0), all(diff(SOC) <= 0), ...
            numel(unique(SOC(isfinite(SOC)))), ...
            safeStat(@min, SOC), safeStat(@max, SOC), ...
            mat2str(SOC(1:min(end,5)).'), ...
            ME.message);

        newME = MException('calculate_inhomogeneity:InterpolantFailure', diagMsg);
        newME = addCause(newME, ME);
        throw(newME);
    end

end

% -----------------------------------------------------------------------
% Local helper: computeUmean
% -----------------------------------------------------------------------
function U_mean = computeUmean(SOC, U, x, w, inhom_offset_fraction)
    % Build query grid.
    % With inhom_offset_fraction = 0 (default): spread is proportional to SOC
    %   (zero spread at SOC=0, full spread at SOC=1).
    % With inhom_offset_fraction > 0: spread at SOC=0 is offset_fraction * max_spread,
    %   growing linearly to full spread at SOC=1.
    x_dev = x - 1;  % deviations from centre: linspace(-0.5, 0.5, 61)
    alpha_eff = inhom_offset_fraction + (1 - inhom_offset_fraction) .* SOC;  % N×1
    Xq = SOC + alpha_eff .* x_dev;  % N×61 via implicit broadcasting

    % Faster interpolation using griddedInterpolant.
    % 'nearest' extrapolation clamps out-of-range queries to the nearest
    % boundary value (U at socMin for Xq<0, U at socMax for Xq>1), which is
    % physically correct for electrode OCV curves and avoids the artifact
    % from the old U(end) clamping when the offset produces negative Xq values.
    F = griddedInterpolant(SOC, U, 'linear', 'nearest');
    E_OC_dist = F(Xq);

    % Weighted average across columns
    U_mean = E_OC_dist * w(:);
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
