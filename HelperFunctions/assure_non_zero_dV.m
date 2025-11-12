function vec_out = assure_non_zero_dV(vec_in)
%% ========================================================================
%> Author: Josef Eizenhammer (josef.eizenhammer@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-06-18
%
%  FUNCTION: assure_non_zero_dV
%  --------------------------------------------------------------------
%  Replaces runs of identical values in a vector with a linear
%  interpolation between neighboring distinct points.
%
%  INPUT:
%     vec_in  – 1-D numeric vector
%
%  OUTPUT:
%     vec_out – vector with plateaus smoothed into ramps
%
%  NOTES:
%     • Intended to avoid division-by-zero in ICA where dU = 0.
%     • Works in-place except for plateau regions.
% ========================================================================

    vec_out = vec_in;                 % initialize output
    n = length(vec_in);
    i = 2;                            % start from second element

    while i < n
        if vec_in(i) == vec_in(i-1)   % plateau detected
            start_idx = i;
            plateau_val = vec_in(i);

            % ------------------- locate end of plateau -------------------
            while i < n && vec_in(i) == plateau_val
                i = i + 1;
            end
            end_idx = i - 1;

            prev_val = vec_in(start_idx - 1);        % value before plateau
            count    = end_idx - start_idx + 1;      % plateau length

            if i <= n && vec_in(i) ~= plateau_val    % there is a next peak
                next_val = vec_in(i);
                delta = (next_val - prev_val) / (count + 1);

                % linearly fill plateau
                for k = 1:count
                    vec_out(start_idx + k - 1) = prev_val + delta * k;
                end
            else
                % plateau extends to or near end of vector ----------------
                next_val = vec_in(end);

                % move backwards to last differing value
                j = start_idx - 1;
                while j > 0 && vec_in(j) == plateau_val
                    j = j - 1;
                end

                if j > 0
                    prev_val = vec_in(j);
                    total_steps = n - j;
                    delta = (next_val - prev_val) / total_steps;

                    for k = 1:(n - j - 1)
                        vec_out(j + k) = prev_val + delta * k;
                    end
                end
                break;  % reached end of input
            end
        else
            i = i + 1; % advance if no plateau
        end
    end
end
