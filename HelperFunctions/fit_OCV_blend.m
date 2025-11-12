function Diff_OCV = fit_OCV_blend(X, myData, ROI_OCV_min, ROI_OCV_max)
%{
%% Original Comments
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-03-11
%
% x = [alpha_an, beta_an, alpha_cat, beta_cat] or 
%     [alpha_an, beta_an, alpha_cat, beta_cat, gamma_Blend2]
%
% If x has only 4 elements, then gamma_Blend2 is assumed to be 0 (pure Blend1 mode).
%
% This function:
%   1) Calculates the blended cell OCV based on the given parameters x.
%   2) Compares the blended OCV (OCV_Calc) with the measured OCV (myData.OCV_cell),
%      but **only** within the intervals specified by [ROI_OCV_min, ROI_OCV_max].
%   3) The comparison is done by summing the squared differences 
%      (OCV_Calc - OCV_measured)^2 for all points in each ROI interval, 
%      ignoring all points outside these intervals.
% 
% In this vectorized version:
%   - X is [N×4] or [N×5], one row per solution.
%   - We return an N×1 vector of errors.
%}

    % 1) Ensure X is at least 2D if it's just a single row
    if size(X,1) == 1 && size(X,2) <= 7
        X = reshape(X,1,[]); 
    end
    N = size(X,1);

    % 2) Ensure myData.OCV_cell is a column vector for consistent broadcasting
    myData.OCV_cell = myData.OCV_cell(:);

    % 3) Build a single logical mask for the region(s) of interest

    Q = myData.Q_cell(:); 
    mask = false(size(Q));

    if numel(ROI_OCV_min) == 1 && numel(ROI_OCV_max) == 1
        % Case 1: one value
         mask = (Q >= ROI_OCV_min) & (Q <= ROI_OCV_max);
    elseif numel(ROI_OCV_min) == 2 && numel(ROI_OCV_max) == 2
        % Case 2: Two values (lower and upper region)
        mask = (Q >= ROI_OCV_min(1) & Q <= ROI_OCV_min(2)) | ...
            (Q >= ROI_OCV_max(1) & Q <= ROI_OCV_max(2));
    else
        error('ROI_OCV_min and ROI_OCV_max must either consist of one or two values.');
    end

    % 4) Prepare output and intermediate storage
    Diff_OCV = zeros(N,1);
    OCV_Calc = zeros(length(Q), N); % each column => OCV for one row of X

    % 5) Compute modeled OCV for each row of X
    for i = 1:N
        x_i = X(i,:);
        OCV_Calc(:,i) = calculate_OCV_blend(x_i, myData);
    end

    % 6) Compare to measured OCV within ROI
    % elementwise difference => squared => mask out
    OCV_errors = (OCV_Calc - myData.OCV_cell).^2;  % [length(Q), N] - [length(Q),1] => broadcasts to [length(Q),N]
    OCV_errors(~mask, :) = 0;                      % zero outside the ROI

    % 7) Sum along Q dimension => one scalar per solution
    Diff_OCV = sum(OCV_errors,1).';                % shape => [N x 1]

    Diff_OCV = Diff_OCV / sum(mask);

    % 8) OPEN QUESTION: Do we want to take the square root and normalize by
    % dividing by the region? If so comment out the following section
    % if numel(ROI_OCV_min) == 1 && numel(ROI_OCV_max) == 1
    %     Diff_OCV = sqrt(Diff_OCV/(ROI_OCV_max-ROI_OCV_min));
    % elseif numel(ROI_OCV_min) == 2 && numel(ROI_OCV_max) == 2
    %     Diff_OCV = sqrt(Diff_OCV/(ROI_OCV_max(2)-ROI_OCV_max(1)+ROI_OCV_min(2)-ROI_OCV_min(1)));
    % else
    %     error('ROI_OCV_min and ROI_OCV_max must either consist of one or two values');
    % end

end
