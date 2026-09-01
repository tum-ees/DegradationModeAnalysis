function uOffset = resistance_offset_voltage(rOffset, solverInput)
%> Author: Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2026-08-28
%
%RESISTANCE_OFFSET_VOLTAGE Ohmic drop that separates a pOCV from the OCV.
%
% A pOCV is measured at a small but finite current, so it does not sit on the
% equilibrium OCV. In charge the terminal voltage is above it, in discharge
% below it, by the ohmic drop of the cell. The reconstruction builds an OCV,
% so without this term the model is compared against a curve that is shifted
% away from it by a constant.
%
% The current is constant over a CC pOCV, so the drop is one number for the
% whole curve. It is kept as a resistance rather than a voltage because a
% resistance can be checked against a measured DCIR, and because expressing
% it that way lets the sign follow the direction instead of being fitted.
%
% One caveat on how to read the fitted value: everything else that displaces
% the two curves from one another by a constant lands in it as well, the
% voltage hysteresis included. It is not a substitute for a DCIR measurement.
%
% INPUTS:
%   rOffset      - fitted series resistance in Ohm
%   solverInput  - solver bundle carrying rOffsetCurrent in A (magnitude of
%                  the pOCV current) and rOffsetSign (+1 charge, -1 discharge)
%
% OUTPUTS:
%   uOffset      - voltage to add to the reconstructed full-cell OCV, in V

if rOffset == 0
    % Also the path taken whenever the offset is switched off, so a solver
    % bundle without the two fields stays valid as long as nothing asks for
    % a non-zero offset.
    uOffset = 0;
    return;
end

if ~isfield(solverInput, 'rOffsetCurrent') || ~isfield(solverInput, 'rOffsetSign')
    error('resistanceOffsetVoltage:MissingCurrent', ...
        ['A non-zero resistance offset was requested but the solver input ' ...
        'carries no pOCV current. Set s.pOCVCurrent in main_dma.']);
end

uOffset = solverInput.rOffsetSign * rOffset * solverInput.rOffsetCurrent;
end
