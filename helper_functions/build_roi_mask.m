function mask = build_roi_mask(Q, ROI_min, ROI_max)
%> Author: Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-11-24
%
%build_roi_mask Build logical mask for one or two ROI intervals.
%
% INPUTS:
%   Q       - column or row vector used as the ROI axis
%   ROI_min - scalar lower bound or first two-element ROI interval
%   ROI_max - scalar upper bound or second two-element ROI interval
%
% OUTPUTS:
%   mask    - logical vector marking samples inside the requested ROI

Q = Q(:);
if isscalar(ROI_min) && isscalar(ROI_max)
    mask = (Q >= ROI_min) & (Q <= ROI_max);
elseif numel(ROI_min) == 2 && numel(ROI_max) == 2
    mask = (Q >= ROI_min(1) & Q <= ROI_min(2)) | ...
           (Q >= ROI_max(1) & Q <= ROI_max(2));
else
    error('ROI_min and ROI_max must either consist of one or two values.');
end
end
