function dvaMeasured = precompute_measured_dva(Q, OCV, q0)
% precompute_measured_dva
% Precompute the measured DVA curve on the supplied Q grid.
%
% INPUTS:
%   Q            - measured full-cell Q or SOC vector
%   OCV          - measured voltage vector corresponding to Q
%   q0           - active full-cell capacity used for normalization
%
% OUTPUTS:
%   dvaMeasured  - filtered measured DVA curve on the input grid

Q  = Q(:);
OCV = OCV(:);

nQ = length(Q);
dvaMeasured = zeros(nQ, 1);
for idx = 2:nQ
    dU = OCV(idx) - OCV(idx-1);
    dQ = Q(idx) - Q(idx-1);
    dvaMeasured(idx-1) = (dU / dQ) * q0;
end
dvaMeasured = apply_filter(dvaMeasured, 'filtermethod', 'sgolay');
end
