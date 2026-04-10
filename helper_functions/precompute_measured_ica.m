function icaMeasured = precompute_measured_ica(Q, OCV, q0)
% precompute_measured_ica
% Precompute the measured ICA curve on the supplied Q grid.
%
% INPUTS:
%   Q            - measured full-cell Q or SOC vector
%   OCV          - measured voltage vector corresponding to Q
%   q0           - active full-cell capacity used for normalization
%
% OUTPUTS:
%   icaMeasured  - filtered measured ICA curve on the input grid

Q  = Q(:);
OCV = OCV(:);

% Interpolate OCV on itself to ensure column format
ocvInterp = interp1(Q, OCV, Q, 'linear', 0);
[~, ~, icaMeasured] = calculate_ica(Q, ocvInterp);
icaMeasured = icaMeasured / q0;
icaMeasured = apply_filter(icaMeasured, 'filtermethod', 'sgolay');
end
