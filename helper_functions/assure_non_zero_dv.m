function vecOut = assure_non_zero_dv(vecIn)
%% ========================================================================
%> Author: Josef Eizenhammer (josef.eizenhammer@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-06-18
%
%  FUNCTION: assure_non_zero_dv
%  --------------------------------------------------------------------
%  Replaces runs of identical values in a vector with a linear
%  interpolation between neighboring distinct points.
%
%  INPUTS:
%     vecIn  - 1-D numeric vector
%
%  OUTPUTS:
%     vecOut - vector with plateaus smoothed into ramps
%
%  NOTES:
%     * Intended to avoid division-by-zero in ICA where dU = 0.
%     * Works in-place except for plateau regions.
% ========================================================================

    vecOut = vecIn;                 % initialize output
    n = length(vecIn);
    i = 2;                          % start from second element

    while i < n
        if vecIn(i) == vecIn(i-1)   % plateau detected
            startIdx = i;
            plateauValue = vecIn(i);

            % ------------------- locate end of plateau -------------------
            while i < n && vecIn(i) == plateauValue
                i = i + 1;
            end
            endIdx = i - 1;

            prevValue = vecIn(startIdx - 1);        % value before plateau
            count     = endIdx - startIdx + 1;      % plateau length

            if i <= n && vecIn(i) ~= plateauValue   % there is a next peak
                nextValue = vecIn(i);
                delta = (nextValue - prevValue) / (count + 1);

                % linearly fill plateau
                for k = 1:count
                    vecOut(startIdx + k - 1) = prevValue + delta * k;
                end
            else
                % plateau extends to or near end of vector ----------------
                nextValue = vecIn(end);

                % move backwards to last differing value
                j = startIdx - 1;
                while j > 0 && vecIn(j) == plateauValue
                    j = j - 1;
                end

                if j > 0
                    prevValue = vecIn(j);
                    totalSteps = n - j;
                    delta = (nextValue - prevValue) / totalSteps;

                    for k = 1:(n - j - 1)
                        vecOut(j + k) = prevValue + delta * k;
                    end
                end
                break;  % reached end of input
            end
        else
            i = i + 1; % advance if no plateau
        end
    end
end
