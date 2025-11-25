function Diff_OCV = fit_OCV(X, myData, ROI_OCV_min, ROI_OCV_max)
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
%   1) Calculates the blended cell OCV based on the given parameters x
%   2) Compares the blended OCV OCV_Calc with the measured OCV myData.OCV_cell
%      only within the intervals specified by ROI_OCV_min and ROI_OCV_max
%   3) The comparison is done by summing the squared differences
%      OCV_Calc minus OCV_measured squared for all points in each ROI interval
%      ignoring all points outside these intervals

% 1) Ensure X is a single row vector
if size(X,1) > 1
    error('fit_OCV is non vectorized and expects a single parameter vector.')
end
if iscolumn(X)
    X = X(:).';
end
X = expand_params_full(X);

% 2) Ensure myData.OCV_cell is a column vector for consistent operations
myData.OCV_cell = myData.OCV_cell(:);

% 3) Build a single logical mask for the region(s) of interest
Q = myData.Q_cell(:);
mask = false(size(Q));

if isscalar(ROI_OCV_min) && isscalar(ROI_OCV_max)
    % Case 1: one value
    mask(:) = (Q >= ROI_OCV_min) & (Q <= ROI_OCV_max);
elseif numel(ROI_OCV_min) == 2 && numel(ROI_OCV_max) == 2
    % Case 2: two values
    mask(:) = (Q >= ROI_OCV_min(1) & Q <= ROI_OCV_min(2)) | ...
              (Q >= ROI_OCV_max(1) & Q <= ROI_OCV_max(2));
else
    error('ROI_OCV_min and ROI_OCV_max must either consist of one or two values.');
end

% 4) Parameter unpacking from fixed layout
alpha_an       = X(1);
beta_an        = X(2);
alpha_cat      = X(3);
beta_cat       = X(4);
gamma_an_blend2   = X(5);
gamma_ca_blend2   = X(6);
inhom_mag_an  = X(7);
inhom_mag_ca  = X(8);

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
    cathU_src = calculate_inhomogeneity(cathSOC_src, cathU_src, inhom_mag_ca);
end

% 7) Construct cathode potential
cathPot = interp1(alpha_cat * cathSOC_src + beta_cat, ...
                  cathU_src, Q, 'linear', 0);

% 8) Full cell OCV
OCV_Calc = cathPot - anodePot;

% 9) Compare to measured OCV within ROI
OCV_errors = (OCV_Calc(:) - myData.OCV_cell).^2;
OCV_errors(~mask) = 0;

% 10) Sum along Q dimension and normalize by ROI length
Diff_OCV = sum(OCV_errors) / sum(mask);

% 11) OPEN QUESTION: Do we want to take the square root and normalize by
% dividing by the region? If so comment out the following section
% if numel(ROI_OCV_min) == 1 && numel(ROI_OCV_max) == 1
%     Diff_OCV = sqrt(Diff_OCV/(ROI_OCV_max-ROI_OCV_min));
% elseif numel(ROI_OCV_min) == 2 && numel(ROI_OCV_max) == 2
%     Diff_OCV = sqrt(Diff_OCV/(ROI_OCV_max(2)-ROI_OCV_max(1)+ROI_OCV_min(2)-ROI_OCV_min(1)));
% else
%     error('ROI_OCV_min and ROI_OCV_max must either consist of one or two values');
% end

end
