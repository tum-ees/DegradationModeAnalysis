function ICA_measured = precompute_measured_ICA(Q, OCV, Q0)
% Precompute measured ICA on the supplied Q grid.
Q  = Q(:);
OCV = OCV(:);

% Interpolate OCV on itself to ensure column format
OCV_interp = interp1(Q, OCV, Q, 'linear', 0);
[~, ~, ICA_measured] = calculate_ICA(Q, OCV_interp);
ICA_measured = ICA_measured / Q0;
ICA_measured = apply_filter(ICA_measured, 'filtermethod', 'sgolay');
end
