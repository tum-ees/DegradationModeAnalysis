function errVal = calculate_DVA_error( ...
    gamma_Si, ...
    Q_Si_interp, Q_Gr_interp, Q_measured_interp, ...
    commonVoltage, ...
    ROI_min, ROI_max ...
)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-10-10
%
% OBJECTIVE FUNCTION: Minimizes DVA error in [ROI_min, ROI_max].

    % Blend the capacities
    Q_blend = gamma_Si * Q_Si_interp + (1 - gamma_Si) * Q_Gr_interp;

    % Calculate the DVAs
    [Q_DVA_blend, ~, DVA_blend]     = calculate_DVA(Q_blend,           commonVoltage);
    [Q_DVA_meas,  ~, DVA_measured]  = calculate_DVA(Q_measured_interp, commonVoltage);

    % Smooth the DVAs
    DVA_blend_smooth    = smooth(DVA_blend,    5, 'lowess');
    DVA_measured_smooth = smooth(DVA_measured, 5, 'lowess');

    % =========== Apply the Region of Interest ===========
    % We only compare within [ROI_min, ROI_max] for the midpoints Q_DVA
    roiIdx = (Q_DVA_meas >= ROI_min) & (Q_DVA_meas <= ROI_max);

    % Mean-squared error in the ROI
    diffROI  = (DVA_blend_smooth(roiIdx) - DVA_measured_smooth(roiIdx)).^2;
    errVal   = sum(diffROI);
end
