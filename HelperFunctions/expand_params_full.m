function Xfull = expand_params_full(X)
%> Author: Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-11-24
%
%EXPAND_PARAMS_FULL Pad legacy parameter vectors to fixed length 8.
% Fixed ordering (future-proof):
% [alpha_an, beta_an, alpha_cat, beta_cat, gamma_an_blend2, gamma_ca_blend2, inhom_an, inhom_ca]
%
% Accepted lengths (single vector only):
%   8 -> already full order.
%   7 -> missing cathode blend placeholder; interprets X = [1:5, inhom_an, inhom_ca].
%   6 -> non-blend with inhomogeneity: X = [1:4, inhom_an, inhom_ca].
%   5 -> blend without inhomogeneity:  [1:5].
%   4 -> non-blend, no inhomogeneity:  [1:4].

if size(X,1) > 1
    error('expand_params_full expects a single parameter vector.');
end
if iscolumn(X)
    X = X(:).';
end

switch numel(X)
    case 8
        Xfull = X;
    case 7
        Xfull = [X(1:5), 0, X(6:7)];
    case 6
        Xfull = [X(1:4), 0, 0, X(5:6)];
    case 5
        Xfull = [X(1:5), 0, 0, 0];
    case 4
        Xfull = [X(1:4), 0, 0, 0, 0];
    otherwise
        error('Unsupported parameter vector length %d. Expected 4-8.', numel(X));
end

end
