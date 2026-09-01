function diffOCV = fit_ocv(X, solverInput, roiOCVMin, roiOCVMax)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> additional code by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-03-11
%
% x is optimized parameters with fixed length
% [alphaAn, betaAn, alphaCat, betaCat, gammaAnBlend2, gammaCaBlend2, inhomAn,
%  inhomCa, rOffset]
%
% This function
%   1) Calculates the blended cell OCV based on the given parameters x
%   2) Compares the blended OCV ocvCalc with the measured OCV solverInput.ocvCell
%      only within the intervals specified by roiOCVMin and roiOCVMax
%   3) The comparison is done by summing the squared differences
%      ocvCalc minus the measured OCV squared for all points in each ROI interval
%      ignoring all points outside these intervals
%
% INPUTS:
%   X            - DMA parameter vector in the expanded fixed ordering
%   solverInput  - struct with measured full-cell data and half-cell curves
%   roiOCVMin    - lower OCV ROI bound or bounds
%   roiOCVMax    - upper OCV ROI bound or bounds
%
% OUTPUTS:
%   diffOCV      - normalized summed squared OCV mismatch inside the ROI

% 1) Ensure X is a single row vector
if size(X,1) > 1
    error('fit_ocv is non vectorized and expects a single parameter vector.')
end
if iscolumn(X)
    X = X(:).';
end
X = expand_params_full(X);

% 2) Ensure solverInput.ocvCell is a column vector for consistent operations
solverInput.ocvCell = solverInput.ocvCell(:);

% 3) Build a single logical mask for the region(s) of interest
Q = solverInput.qCell(:);
mask = false(size(Q));

if isscalar(roiOCVMin) && isscalar(roiOCVMax)
    % Case 1: one value
    mask(:) = (Q >= roiOCVMin) & (Q <= roiOCVMax);
elseif numel(roiOCVMin) == 2 && numel(roiOCVMax) == 2
    % Case 2: two values
    mask(:) = (Q >= roiOCVMin(1) & Q <= roiOCVMin(2)) | ...
              (Q >= roiOCVMax(1) & Q <= roiOCVMax(2));
else
    error('roiOCVMin and roiOCVMax must each contain one or two values.');
end

% 4) Parameter unpacking from fixed layout
alphaAn       = X(1);
betaAn        = X(2);
alphaCat      = X(3);
betaCat       = X(4);
gammaAnBlend2 = X(5);
gammaCaBlend2 = X(6);
inhomMagAn    = X(7);
inhomMagCa    = X(8);
rOffset       = X(9);

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
    cathUSrc = calculate_inhomogeneity(cathodeSOCSrc, cathUSrc, inhomMagCa, solverInput.inhomCathodeOffset);
end

% 7) Construct cathode potential
cathPot = interp1(alphaCat * cathodeSOCSrc + betaCat, ...
                  cathUSrc, Q, 'linear', 0);

% 8) Full cell OCV, lifted onto the measured pOCV by the ohmic drop.
% The term is a constant over the curve, so it moves the level without
% touching the shape. That is deliberate: fit_dva and fit_ica never read
% the offset slot, and even where a DVA is derived from the lifted curve a
% constant vanishes under the derivative, so the offset can only earn its
% place on the OCV term.
ocvCalc = cathPot - anodePot + resistance_offset_voltage(rOffset, solverInput);

% 9) Compare to measured OCV within ROI
ocvErrors = (ocvCalc(:) - solverInput.ocvCell).^2;
ocvErrors(~mask) = 0;

% 10) Sum along Q dimension and normalize by ROI length
diffOCV = sum(ocvErrors) / sum(mask);

end
