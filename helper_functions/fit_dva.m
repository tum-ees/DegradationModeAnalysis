function diffDVA = fit_dva(X, solverInput, q0, roiDVAMin, roiDVAMax, precompDVA)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> additional code by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-03-11
%
% x is optimized parameters with fixed length
% [alphaAn, betaAn, alphaCat, betaCat, gammaAnBlend2, gammaCaBlend2, inhomAn, inhomCa]
%
% This function
%   1) Builds the anode curve based on gammaAnBlend2
%   2) Interpolates anode, cathode, and measured OCV on the same Q grid
%   3) Computes discrete dU dQ DVA for anode, cathode, and OCV
%   4) Forms dvaSum = dvaCathode minus dvaAnode and compares it to measured DVA
%      only within roiDVAMin and roiDVAMax by summing squared differences
%      diffDVA is that sum normalized by the ROI length
%
% INPUTS:
%   X            - DMA parameter vector in the expanded fixed ordering
%   solverInput  - struct with measured full-cell data and half-cell curves
%   q0           - active full-cell capacity used for DVA scaling
%   roiDVAMin    - lower DVA ROI bound or bounds
%   roiDVAMax    - upper DVA ROI bound or bounds
%   precompDVA   - optional precomputed measured DVA curve on qCell
%
% OUTPUTS:
%   diffDVA      - normalized summed squared DVA mismatch inside the ROI

% 1) Ensure X is a single row vector
if size(X,1) > 1
    error('fit_dva is non vectorized and expects a single parameter vector.')
end
if iscolumn(X)
    X = X(:).';
end
X = expand_params_full(X);

% 2) Read parameters from fixed layout
alphaAn      = X(1);
betaAn       = X(2);
alphaCat     = X(3);
betaCat      = X(4);
gammaAnBlend2 = X(5);
gammaCaBlend2 = X(6);
inhomMagAn   = X(7);
inhomMagCa   = X(8);

Q  = solverInput.qCell(:);
nQ = length(Q);

% 3) Reuse precomputed mask and measured DVA if provided
if nargin >= 6 && ~isempty(precompDVA)
    mask    = precompDVA.mask;
    dva_ocv = precompDVA.measuredDVA;
else
    mask    = build_roi_mask(Q, roiDVAMin, roiDVAMax);
    dva_ocv = precompute_measured_dva(Q, solverInput.ocvCell, q0);
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

% 7) Build cathode potential
cathodePot = interp1(alphaCat * cathodeSOCSrc + betaCat, ...
                     cathUSrc, Q, 'linear', 0);

% 8) Compute discrete DVA for measured OCV if not precomputed
if nargin < 6 || isempty(precompDVA)
    ocvDVA = interp1(Q, solverInput.ocvCell, Q, 'linear', 0);
    dva_ocv = zeros(nQ, 1);
    for idx = 2:nQ
        dU = ocvDVA(idx) - ocvDVA(idx-1);
        dQ = Q(idx) - Q(idx-1);
        dva_ocv(idx-1) = (dU / dQ) * q0;
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

    dva_anode(idx-1)   = (dU_an  / dQ) * q0;
    dva_cathode(idx-1) = (dU_cat / dQ) * q0;
end

% 10) DVA sum and smoothing
dva_sum = dva_cathode - dva_anode;
dva_sum = apply_filter(dva_sum, 'filtermethod', 'sgolay');

% 11) Compare only within ROI
diffArray = (dva_sum - dva_ocv).^2;
diffArray(~mask) = 0;

diffDVA = sum(diffArray) / sum(mask);

end
