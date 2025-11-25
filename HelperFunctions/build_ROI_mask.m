function mask = build_ROI_mask(Q, ROI_min, ROI_max)
%> Author: Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-11-24
%
%BUILD_ROI_MASK Build logical mask for one or two ROI intervals.
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
