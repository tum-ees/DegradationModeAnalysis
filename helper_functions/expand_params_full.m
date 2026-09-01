function Xfull = expand_params_full(X)
%> Author: Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-11-24
%
%EXPAND_PARAMS_FULL Pad legacy parameter vectors to fixed length 9.
% Fixed ordering (future-proof):
% [alphaAn, betaAn, alphaCat, betaCat, gammaAnBlend2, gammaCaBlend2, inhomAn,
%  inhomCa, rOffset]
%
% Accepted lengths (single vector only):
%   9 -> already full order.
%   8 -> no resistance offset; interprets X = [1:8].
%   7 -> missing cathode blend placeholder; interprets X = [1:5, inhom_an, inhom_ca].
%   6 -> non-blend with inhomogeneity: X = [1:4, inhom_an, inhom_ca].
%   5 -> blend without inhomogeneity:  [1:5].
%   4 -> non-blend, no inhomogeneity:  [1:4].
%
% A padded slot is zero, and a zero resistance offset contributes no voltage,
% so every vector stored before slot 9 existed still expands to the same
% reconstruction it did then.
%
% INPUTS:
%   X      - legacy or current DMA parameter vector with length 4 to 9
%
% OUTPUTS:
%   Xfull  - 1x9 parameter vector in the current fixed ordering

if size(X,1) > 1
    error('expand_params_full expects a single parameter vector.');
end
if iscolumn(X)
    X = X(:).';
end

switch numel(X)
    case 9
        Xfull = X;
    case 8
        Xfull = [X(1:8), 0];
    case 7
        Xfull = [X(1:5), 0, X(6:7), 0];
    case 6
        Xfull = [X(1:4), 0, 0, X(5:6), 0];
    case 5
        Xfull = [X(1:5), 0, 0, 0, 0];
    case 4
        Xfull = [X(1:4), 0, 0, 0, 0, 0];
    otherwise
        error('Unsupported parameter vector length %d. Expected 4-9.', numel(X));
end

end
