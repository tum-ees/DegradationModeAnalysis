function errVal = calculate_dva_error( ...
    gammaSi, ...
    qSiInterp, qGrInterp, qMeasuredInterp, ...
    commonVoltage, ...
    roiMin, roiMax ...
)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-10-10
%
% OBJECTIVE FUNCTION: Minimizes DVA error in [ROI_min, ROI_max].
%
% INPUTS:
%   gammaSi         - silicon fraction used to blend graphite and silicon
%   qSiInterp       - interpolated silicon capacity curve on commonVoltage
%   qGrInterp       - interpolated graphite capacity curve on commonVoltage
%   qMeasuredInterp - measured full-cell capacity curve on commonVoltage
%   commonVoltage   - common voltage grid used for all DVA calculations
%   roiMin          - lower ROI bound in Q/SOC space
%   roiMax          - upper ROI bound in Q/SOC space
%
% OUTPUTS:
%   errVal          - normalized summed squared DVA error inside the ROI

    % Blend the capacities
    qBlend = gammaSi * qSiInterp + (1 - gammaSi) * qGrInterp;

    % Calculate the DVAs
    [~, ~, dvaBlend]     = calculate_dva(qBlend, commonVoltage);
    [qDVAMeas, ~, dvaMeasured] = calculate_dva(qMeasuredInterp, commonVoltage);

    % Smooth the DVAs
    dvaBlendSmooth = smooth(dvaBlend, 5, 'lowess');
    dvaMeasuredSmooth = smooth(dvaMeasured, 5, 'lowess');

    % =========== Apply the Region of Interest ===========
    % We only compare within [roiMin, roiMax] for the midpoints qDVA
    roiIdx = (qDVAMeas >= roiMin) & (qDVAMeas <= roiMax);

    % Mean-squared error in the ROI
    diffROI = (dvaBlendSmooth(roiIdx) - dvaMeasuredSmooth(roiIdx)).^2;
    errVal = sum(diffROI);
end
