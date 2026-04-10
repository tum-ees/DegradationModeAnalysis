function diffICA = fit_ica(X, solverInput, q0, roiICAMin, roiICAMax, precompICA)
%> Author: Josef Eizenhammer (josef.eizenhammer@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> additional code by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-05-19
%
% x is optimized parameters with fixed length
% [alphaAn, betaAn, alphaCat, betaCat, gammaAnBlend2, gammaCaBlend2, inhomAn, inhomCa]
%
% This function
%   1) Builds the anode curve based on gammaAnBlend2
%   2) Interpolates anode, cathode, and measured OCV on the same Q grid
%   3) Computes discrete dQ dU ICA for measured and modeled OCV
%   4) Compares ICA curves only inside roiICAMin and roiICAMax
%      by summing squared differences and normalizing by the ROI length
%
% INPUTS:
%   X            - DMA parameter vector in the expanded fixed ordering
%   solverInput  - struct with measured full-cell data and half-cell curves
%   q0           - active full-cell capacity used for ICA scaling
%   roiICAMin    - lower ICA ROI bound or bounds
%   roiICAMax    - upper ICA ROI bound or bounds
%   precompICA   - optional precomputed measured ICA curve on qCell
%
% OUTPUTS:
%   diffICA      - normalized summed squared ICA mismatch inside the ROI

% 1) Ensure X is a single row vector
if size(X,1) > 1
    error('fit_ica is non vectorized and expects a single parameter vector.')
end
if iscolumn(X)
    X = X(:).';
end
X = expand_params_full(X);

% 2) Parameter unpacking from fixed layout
alphaAn       = X(1);
betaAn        = X(2);
alphaCat      = X(3);
betaCat       = X(4);
gammaAnBlend2 = X(5);
gammaCaBlend2 = X(6);
inhomMagAn    = X(7);
inhomMagCa    = X(8);

Q  = solverInput.qCell(:);   

% 3) Reuse precomputed mask and measured ICA if provided
if nargin >= 6 && ~isempty(precompICA)
    mask    = precompICA.mask;
else
    mask    = build_roi_mask(Q, roiICAMin, roiICAMax);
    ocvICA = interp1(Q, solverInput.ocvCell, Q, 'linear', 0);
end

% 5) Prepare anode source curve; split between blend and non-blend paths
if solverInput.useAnodeBlend && ~isempty(solverInput.qAnodeBlend1Interp)
    [anodeSOCSrc, anodeUSrc] = calculate_blend_curve(gammaAnBlend2, solverInput, 'anode');
else
    anodeSOCSrc = solverInput.anodeSOCSingle;
    anodeUSrc   = solverInput.anodeUSingle;
end

% Apply anode inhomogeneity if needed
if inhomMagAn ~= 0
    anodeUSrc = calculate_inhomogeneity(anodeSOCSrc, anodeUSrc, inhomMagAn, solverInput.inhomAnodeOffset);
end

% Interpolate anode potential on Q grid
anodePot = interp1(alphaAn * anodeSOCSrc + betaAn, ...
                   anodeUSrc, Q, 'linear', 0);

% 6) Build cathode (blend optional) and apply inhomogeneity
if solverInput.useCathodeBlend && ~isempty(solverInput.qCathodeBlend1Interp)
    [cathodeSOCSrc, cathUSrc] = calculate_blend_curve(gammaCaBlend2, solverInput, 'cathode');
else
    cathodeSOCSrc = solverInput.cathodeSOCSingle;
    cathUSrc      = solverInput.cathodeUSingle;
end
if inhomMagCa ~= 0
    cathUSrc = calculate_inhomogeneity( ...
        cathodeSOCSrc, cathUSrc, inhomMagCa, solverInput.inhomCathodeOffset);
end

% 7) Interpolate cathode potential on Q grid
cathodePot = interp1(alphaCat * cathodeSOCSrc + betaCat, ...
                     cathUSrc, Q, 'linear', 0);

% 8) Compute ICA for the measurement if not precomputed
if nargin < 6 || isempty(precompICA)
    [~, ~, icaOCV] = calculate_ica(Q, ocvICA);
    icaOCV = icaOCV / q0;
    icaOCV = apply_filter(icaOCV, 'filtermethod', 'sgolay');
else
    icaOCV = precompICA.measuredICA;
end

% 9) Compute ICA for the modeled curve
ocvSum = cathodePot - anodePot;
[~, ~, icaCalc] = calculate_ica(Q, ocvSum);
icaCalc = icaCalc / q0;
icaCalc = apply_filter(icaCalc, 'filtermethod', 'sgolay');

% 10) Compare only within ROI
diffArray = (icaCalc - icaOCV).^2;
diffArray(~mask) = 0;

diffICA = sum(diffArray) / sum(mask);

end
