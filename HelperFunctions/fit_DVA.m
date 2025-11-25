function Diff_DVA = fit_DVA(X, myData, Q0, ROI_DVA_min, ROI_DVA_max, precompDVA)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> additional code by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-03-11
%
% x is optimized parameters with fixed length
% [alpha_an, beta_an, alpha_cat, beta_cat, gamma_an_blend2, gamma_ca_blend2, inhom_an, inhom_ca]
%
% This function
%   1) Builds the anode curve based on gamma_an_blend2
%   2) Interpolates anode, cathode, and measured OCV on the same Q grid
%   3) Computes discrete dU dQ DVA for anode, cathode, and OCV
%   4) Forms DVA_Sum = DVA_cathode minus DVA_anode and compares it to DVA_OCV
%      only within ROI_DVA_min and ROI_DVA_max by summing squared differences
%      Diff_DVA is that sum normalized by the ROI length

% 1) Ensure X is a single row vector
if size(X,1) > 1
    error('fit_DVA is non vectorized and expects a single parameter vector.')
end
if iscolumn(X)
    X = X(:).';
end
X = expand_params_full(X);

% 2) Read parameters from fixed layout
alpha_an      = X(1);
beta_an       = X(2);
alpha_cat     = X(3);
beta_cat      = X(4);
gamma_an_blend2  = X(5);
gamma_ca_blend2  = X(6);
inhom_mag_an  = X(7);
inhom_mag_ca  = X(8);

Q  = myData.Q_cell(:);
nQ = length(Q);

% 3) Reuse precomputed mask and measured DVA if provided
if nargin >= 6 && ~isempty(precompDVA)
    mask    = precompDVA.mask;
    dva_ocv = precompDVA.measuredDVA;
else
    mask    = build_ROI_mask(Q, ROI_DVA_min, ROI_DVA_max);
    dva_ocv = precompute_measured_DVA(Q, myData.OCV_cell, Q0);
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

% 7) Build cathode potential
cathodePot = interp1(alpha_cat * cathSOC_src + beta_cat, ...
                     cathU_src, Q, 'linear', 0);

% 8) Compute discrete DVA for measured OCV if not precomputed
if nargin < 6 || isempty(precompDVA)
    OCV_DVA = interp1(Q, myData.OCV_cell, Q, 'linear', 0);
    dva_ocv = zeros(nQ, 1);
    for idx = 2:nQ
        dU = OCV_DVA(idx) - OCV_DVA(idx-1);
        dQ = Q(idx) - Q(idx-1);
        dva_ocv(idx-1) = (dU / dQ) * Q0;
    end
    dva_ocv = apply_filter(dva_ocv, 'filtermethod', 'sgolay');
end

% 9) Compute discrete DVA for anode and cathode
dva_anode   = zeros(nQ, 1);
dva_cathode = zeros(nQ, 1);
for idx = 2:nQ
    dU_an  = anodePot(idx)   - anodePot(idx-1);
    dU_cat = cathodePot(idx) - cathodePot(idx-1);
    dQ     = Q(idx) - Q(idx-1);

    dva_anode(idx-1)   = (dU_an  / dQ) * Q0;
    dva_cathode(idx-1) = (dU_cat / dQ) * Q0;
end

% 10) DVA sum and smoothing
dva_sum = dva_cathode - dva_anode;
dva_sum = apply_filter(dva_sum, 'filtermethod', 'sgolay');

% 11) Compare only within ROI
diffArray = (dva_sum - dva_ocv).^2;
diffArray(~mask) = 0;

Diff_DVA = sum(diffArray) / sum(mask);

end
