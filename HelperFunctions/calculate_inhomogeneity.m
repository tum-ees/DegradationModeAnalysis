function U_mean = calculate_inhomogeneity(SOC, U, inhomo_max)
%> Author: Moritz Guenthner (moritz.guenthner@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-10-05
%
% This function calculates the inhomogeneity effect on the OCV curve
% based on a Gaussian distribution of local SOCs around the mean SOC.
% Inhomogeneity is zero at 0% full-cell SOC and maximum at 100% full-
% cell SOC.

    % Skip the calculation if no inhomogeneity should be applied.
    if inhomo_max <= 0
        U_mean = U;
        return;
    end

    x     = linspace(0.5, 1.5, 61);
    mu    = 1;
    sigma = inhomo_max;

    % Normalized Gaussian weights
    w = normpdf(x, mu, sigma);
    w = w / sum(w);

    % Build query grid: each column = SOC scaled by x(i)
    Xq = SOC(:) * x;   % size [length(SOC) - length(x)]

    % Interpolate all at once; fill missing with U(end)
    E_OC_dist = interp1(SOC, U, Xq, 'linear', U(end));

    % Weighted average across columns
    U_mean = E_OC_dist * w(:);   % vector of length(SOC)

end
