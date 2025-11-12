function [Q_ICA, OCV_ICA, ICA] = calculate_ICA(Q, OCV, steps)
%% ========================================================================
%> Author: Josef Eizenhammer (josef.eizenhammer@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-06-18
%
%  FUNCTION: calculate_ICA
%  --------------------------------------------------------------------
%  Calculates an Incremental Capacity Analysis (ICA) curve from input
%  charge (Q) and open-circuit voltage (OCV) vectors.
%
%  INPUTS:
%     Q      – charge vector  [Ah]
%     OCV    – voltage vector [V]   (must be same length as Q)
%     steps  – (optional) number of interpolation points (default = 1001)
%
%  OUTPUTS:
%     Q_ICA  – uniformly sampled charge vector        [Ah]
%     OCV_ICA– interpolated voltage vector            [V]
%     ICA    – incremental capacity dQ/dU             [Ah V⁻¹]
%
%  NOTES / ASSUMPTIONS:
%     • NaNs are stripped before processing.
%     • Duplicate Q entries are removed via unique().
%     • ICA is smoothed with LOWESS (span = 30 samples).
%     • No explicit guard against dU = 0; caller should validate results.
% ========================================================================

    % ---------------------- Default argument handling --------------------
    if nargin < 3
        steps = 1001;                      % use 1001 points if not provided
    end

    % ---------------------------- Sanitization ---------------------------
    validIdx = ~isnan(Q);                  % ignore NaNs in Q/OCV pairs
    Q   = Q(validIdx);
    OCV = OCV(validIdx);

    [Q, idxUnique] = unique(Q);            % ensure monotonic Q for interp
    OCV = OCV(idxUnique);

    % --------------------- Uniform re-sampling of curve ------------------
    Q_ICA  = linspace(min(Q), max(Q), steps);
    OCV_ICA = interp1(Q, OCV, Q_ICA);      % linear interpolation

    OCV_ICA = assure_non_zero_dV(OCV_ICA);  % fix flat segments

    % ------------------ Incremental Capacity Calculation -----------------
    ICA = zeros(size(OCV_ICA));            % preallocate for speed
    for i = 2:numel(ICA)
        dU = OCV_ICA(i) - OCV_ICA(i-1);
        dQ = Q_ICA(i)  - Q_ICA(i-1);

        % NOTE: dU can be ~0 if voltage plateau is perfectly flat
        %       User may want to impose a minimum dU threshold.
        ICA(i-1) = dQ / dU;
    end
end
