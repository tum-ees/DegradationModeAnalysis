function Diff_ICA = fit_ICA(X, myData, Q0, ROI_ICA_min, ROI_ICA_max, precompICA)
%> Author: Josef Eizenhammer (josef.eizenhammer@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> additional code by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-05-19
%
% x is optimized parameters with fixed length
% [alpha_an, beta_an, alpha_cat, beta_cat, gamma_an_blend2, gamma_ca_blend2, inhom_an, inhom_ca]
%
% This function
%   1) Builds the anode curve based on gamma_an_blend2
%   2) Interpolates anode, cathode, and measured OCV on the same Q grid
%   3) Computes discrete dQ dU ICA for measured and modeled OCV
%   4) Compares ICA curves only inside ROI_ICA_min and ROI_ICA_max
%      by summing squared differences and normalizing by the ROI length

% 1) Ensure X is a single row vector
if size(X,1) > 1
    error('fit_ICA is non vectorized and expects a single parameter vector.')
end
if iscolumn(X)
    X = X(:).';
end
X = expand_params_full(X);

% 2) Parameter unpacking from fixed layout
alpha_an       = X(1);
beta_an        = X(2);
alpha_cat      = X(3);
beta_cat       = X(4);
gamma_an_blend2   = X(5);
gamma_ca_blend2   = X(6);
inhom_mag_an  = X(7);
inhom_mag_ca  = X(8);

Q  = myData.Q_cell(:);   

% 3) Reuse precomputed mask and measured ICA if provided
if nargin >= 6 && ~isempty(precompICA)
    mask    = precompICA.mask;
    ICA_OCV = precompICA.measuredICA;
else
    mask    = build_ROI_mask(Q, ROI_ICA_min, ROI_ICA_max);
    OCV_ICA = interp1(Q, myData.OCV_cell, Q, 'linear', 0);
end

% 5) Prepare anode source curve; split between blend and non-blend paths
if myData.useAnodeBlend && ~isempty(myData.Q_anode_blend1_interp)
    [anodeSOC_src, anodeU_src] = calculate_blend_curve(gamma_an_blend2, myData, 'anode');
else
    anodeSOC_src = myData.anode_SOC_single;
    anodeU_src   = myData.anode_U_single;
end

% Apply anode inhomogeneity if needed
if inhom_mag_an ~= 0
    anodeU_src = calculate_inhomogeneity(anodeSOC_src, anodeU_src, inhom_mag_an);
end

% Interpolate anode potential on Q grid
anodePot = interp1(alpha_an * anodeSOC_src + beta_an, ...
                   anodeU_src, Q, 'linear', 0);

% 6) Build cathode (blend optional) and apply inhomogeneity
if myData.useCathodeBlend && ~isempty(myData.Q_cathode_blend1_interp)
    [cathSOC_src, cathU_src] = calculate_blend_curve(gamma_ca_blend2, myData, 'cathode');
else
    cathSOC_src = myData.cathode_SOC_single;
    cathU_src   = myData.cathode_U_single;
end
if inhom_mag_ca ~= 0
    cathU_src = calculate_inhomogeneity( ...
        cathSOC_src, cathU_src, inhom_mag_ca);
end

% 7) Interpolate cathode potential on Q grid
cathodePot = interp1(alpha_cat * cathSOC_src + beta_cat, ...
                     cathU_src, Q, 'linear', 0);

% 8) Compute ICA for the measurement if not precomputed
if nargin < 6 || isempty(precompICA)
    [~, ~, ICA_OCV] = calculate_ICA(Q, OCV_ICA);
    ICA_OCV = ICA_OCV / Q0;
    ICA_OCV = apply_filter(ICA_OCV, 'filtermethod', 'sgolay');
end

% 9) Compute ICA for the modeled curve
OCV_sum = cathodePot - anodePot;
[~, ~, ICA_calc] = calculate_ICA(Q, OCV_sum);
ICA_calc = ICA_calc / Q0;
ICA_calc = apply_filter(ICA_calc, 'filtermethod', 'sgolay');

% 10) Compare only within ROI
diffArray = (ICA_calc - ICA_OCV).^2;
diffArray(~mask) = 0;

Diff_ICA = sum(diffArray) / sum(mask);

end
