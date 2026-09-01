function [rBounds, rInit] = resistance_offset_bounds(settings, capaAct)
%> Author: Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2026-08-28
%
%RESISTANCE_OFFSET_BOUNDS Turn the capacity-normalised limit into a resistance.
%
% The limit is given in Ohm*Ah rather than in Ohm, because a resistance on
% its own is not comparable between cells of different size. The voltage the
% offset can reach is
%
%     R * I = (limit / capaAct) * I = limit * C-rate
%
% so one limit gives every cell in a study the same room in volts as long as
% the check-up runs at the same rate, and gives a cell run at twice the rate
% twice the room, which is what an ohmic drop does. The rate here is the
% momentary one, current over the actual capacity of the check-up: a study
% that fixes its current from the nominal capacity sees the room in volts
% grow as the cell fades (a quarter more at 80 percent state of health) and
% differ a few percent between cells. That is a property of holding the
% current fixed, not of the bound.
%
% The sign is a separate decision. A pOCV lies above the OCV in charge and
% below it in discharge, so a positive resistance is the physical case and is
% all that is allowed by default. Half-cell references are pseudo-OCPs
% measured at a finite rate themselves and carry their own drop into the
% reconstruction, which can exceed the drop in the full-cell measurement; a
% negative net resistance is then the correct answer, and
% allowNegativeResistanceOffset opens that half of the range.
%
% INPUTS:
%   settings  - DMA settings struct, read for resistanceOffsetLimit and
%               allowNegativeResistanceOffset
%   capaAct   - active full-cell capacity of this check-up, in Ah
%
% OUTPUTS:
%   rBounds   - [lower, upper] for the fitted resistance, in Ohm
%   rInit     - starting value inside those bounds, in Ohm

if isfield(settings, 'resistanceOffsetLimit')
    limitOhmAh = settings.resistanceOffsetLimit;
else
    limitOhmAh = 0.25;
end
if ~isscalar(limitOhmAh) || ~isfinite(limitOhmAh) || limitOhmAh < 0
    error('resistanceOffsetBounds:BadLimit', ...
        's.resistanceOffsetLimit must be a non-negative scalar in Ohm*Ah.');
end
if ~isscalar(capaAct) || ~isfinite(capaAct) || capaAct <= 0
    error('resistanceOffsetBounds:BadCapacity', ...
        ['The active capacity is %g Ah, so the limit cannot be turned into ' ...
        'a resistance.'], capaAct);
end

rMax = limitOhmAh / capaAct;

if isfield(settings, 'allowNegativeResistanceOffset') && ...
        settings.allowNegativeResistanceOffset
    rBounds = [-rMax, rMax];
else
    rBounds = [0, rMax];
end

rInit = mean(rBounds);
end
