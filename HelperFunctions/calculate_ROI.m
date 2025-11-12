function [lowestROI, highestROI] = calculate_ROI(s)
%> Author: Josef Eizenhammer (josef.eizenhammer@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-09-11
%
% This function performs:
%   * Aproach of region-based RMSE. Define a single bounding region.
%   * We'll take the min across all ROI arrays.
%   * Use the corresponding region according to selected fitting method.
    

%% ----------- LOAD DATA FROM SETTINGS ------------------------------------
ROI_OCV_min = s.ROI_OCV_min;
ROI_DVA_min = s.ROI_DVA_min;
ROI_ICA_min = s.ROI_ICA_min;
ROI_OCV_max = s.ROI_OCV_max;
ROI_DVA_max = s.ROI_DVA_max;
ROI_ICA_max = s.ROI_ICA_max;

weightDVA = s.weightDVA;
weightICA = s.weightICA;

%% ----------- CALCULATE ROI ----------------------------------------------
    if (weightDVA ~= 0) && (weightICA == 0)
        lowestROI  = min([ROI_OCV_min(:); ROI_DVA_min]);
        highestROI = max([ROI_OCV_max(:); ROI_DVA_max]);
    elseif (weightDVA == 0) && (weightICA ~= 0)
        lowestROI  = min([ROI_OCV_min(:); ROI_ICA_min]);
        highestROI = max([ROI_OCV_max(:); ROI_ICA_max]);
    elseif (weightDVA == 0) && (weightICA == 0)
        lowestROI  = min(ROI_OCV_min(:));
        highestROI = max(ROI_OCV_max(:));
    else
        lowestROI  = min([ROI_OCV_min(:); ROI_DVA_min; ROI_ICA_min]);
        highestROI = max([ROI_OCV_max(:); ROI_DVA_max; ROI_ICA_max]);
    end
end

