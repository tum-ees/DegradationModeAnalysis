function DVA_measured = precompute_measured_DVA(Q, OCV, Q0)
% Precompute measured DVA on the supplied Q grid.
Q  = Q(:);
OCV = OCV(:);

nQ = length(Q);
DVA_measured = zeros(nQ, 1);
for idx = 2:nQ
    dU = OCV(idx) - OCV(idx-1);
    dQ = Q(idx) - Q(idx-1);
    DVA_measured(idx-1) = (dU / dQ) * Q0;
end
DVA_measured = apply_filter(DVA_measured, 'filtermethod', 'sgolay');
end
