function U_mean = calculate_inhomogeneity(SOC, U, inhom_max)
%> Author: Moritz Guenthner (moritz.guenthner@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-11-24
%
% This function calculates the inhomogeneity effect on the OCV curve
% based on a Gaussian distribution of local SOCs around the mean SOC.
% Inhomogeneity is zero at 0 percent full cell SOC and maximum at 100 percent full
% cell SOC.

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
        U_mean = computeUmean(SOC, U, x, w);
        return;
    end

    try
        U_mean = computeUmean(SOC, U, x, w);
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
function U_mean = computeUmean(SOC, U, x, w)
    % Build query grid
    Xq = SOC * x;

    % Faster interpolation using griddedInterpolant
    F = griddedInterpolant(SOC, U, 'linear', 'linear');
    E_OC_dist = F(Xq);

    % Match old behavior: outside range gives U(end) (handle non-monotonic SOC order)
    socMin = min(SOC);
    socMax = max(SOC);
    outMask = (Xq < socMin) | (Xq > socMax);
    E_OC_dist(outMask) = U(end);

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

