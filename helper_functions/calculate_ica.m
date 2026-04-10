function [qICA, ocvICA, ica] = calculate_ica(Q, OCV, steps)
%% ========================================================================
%> Author: Josef Eizenhammer (josef.eizenhammer@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-06-18
%
%  FUNCTION: calculate_ica
%  --------------------------------------------------------------------
%  Calculates an Incremental Capacity Analysis (ICA) curve from input
%  charge (Q) and open-circuit voltage (OCV) vectors.
%
%  INPUTS:
%     Q      - charge vector  [Ah]
%     OCV    - voltage vector [V]   (must be same length as Q)
%     steps  - optional number of interpolation points (default = 1001)
%
%  OUTPUTS:
%     qICA   - uniformly sampled charge vector        [Ah]
%     ocvICA - interpolated voltage vector            [V]
%     ica    - incremental capacity dQ/dU             [Ah V^-1]
%
%  NOTES / ASSUMPTIONS:
%     * NaNs are stripped before processing.
%     * Duplicate Q entries are removed via unique().
%     * ICA is smoothed with LOWESS (span = 30 samples).
%     * No explicit guard against dU = 0; caller should validate results.
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
    qICA   = linspace(min(Q), max(Q), steps);
    ocvICA = interp1(Q, OCV, qICA);        % linear interpolation

    ocvICA = assure_non_zero_dv(ocvICA);   % fix flat segments

    % ------------------ Incremental Capacity Calculation -----------------
    ica = zeros(size(ocvICA));             % preallocate for speed
    for i = 2:numel(ica)
        dU = ocvICA(i) - ocvICA(i-1);
        dQ = qICA(i)  - qICA(i-1);

        % NOTE: dU can be ~0 if voltage plateau is perfectly flat
        %       User may want to impose a minimum dU threshold.
        ica(i-1) = dQ / dU;
    end
end
