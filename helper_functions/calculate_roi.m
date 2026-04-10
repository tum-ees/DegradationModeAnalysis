function [lowestROI, highestROI] = calculate_roi(settings)
%> Author: Josef Eizenhammer (josef.eizenhammer@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-09-11
%
% This function performs:
%   * Aproach of region-based RMSE. Define a single bounding region.
%   * We'll take the min across all ROI arrays.
%   * Use the corresponding region according to selected fitting method.
%
% INPUTS:
%   settings    - DMA settings struct containing ROI bounds and objective
%                 weights for OCV, DVA, and ICA fitting
%
% OUTPUTS:
%   lowestROI   - global lower ROI bound used for plotting and comparisons
%   highestROI  - global upper ROI bound used for plotting and comparisons
    

%% ----------- LOAD DATA FROM SETTINGS ------------------------------------
roiOCVMin = settings.roiOCVMin;
roiDVAMin = settings.roiDVAMin;
roiICAMin = settings.roiICAMin;
roiOCVMax = settings.roiOCVMax;
roiDVAMax = settings.roiDVAMax;
roiICAMax = settings.roiICAMax;

weightDVA = settings.weightDVA;
weightICA = settings.weightICA;

%% ----------- CALCULATE ROI ----------------------------------------------
    if (weightDVA ~= 0) && (weightICA == 0)
        lowestROI  = min([roiOCVMin(:); roiDVAMin]);
        highestROI = max([roiOCVMax(:); roiDVAMax]);
    elseif (weightDVA == 0) && (weightICA ~= 0)
        lowestROI  = min([roiOCVMin(:); roiICAMin]);
        highestROI = max([roiOCVMax(:); roiICAMax]);
    elseif (weightDVA == 0) && (weightICA == 0)
        lowestROI  = min(roiOCVMin(:));
        highestROI = max(roiOCVMax(:));
    else
        lowestROI  = min([roiOCVMin(:); roiDVAMin; roiICAMin]);
        highestROI = max([roiOCVMax(:); roiDVAMax; roiICAMax]);
    end
end

