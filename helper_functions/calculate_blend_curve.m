function [blendSOC, blendVoltage] = calculate_blend_curve(gammaBlend2, solverInput, electrode)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%
% This function calculates the SOC and voltage of a blended electrode curve (anode or cathode).
%
% INPUTS:
%   gammaBlend2  - blend share of component 2
%   solverInput  - struct containing commonVoltage / q* fields
%   electrode    - 'anode' (default) or 'cathode'
%
% OUTPUTS:
%   blendSOC      - blended electrode SOC curve on the common voltage grid
%   blendVoltage  - voltage grid corresponding to the blended SOC curve
%
    if nargin < 3 || isempty(electrode)
        electrode = 'anode';
    end

    switch lower(electrode)
        case 'anode'
            commonVoltage = solverInput.commonVoltageAnode;
            qBlend2Interp = solverInput.qAnodeBlend2Interp;
            qBlend1Interp = solverInput.qAnodeBlend1Interp;
        case 'cathode'
            commonVoltage = solverInput.commonVoltageCathode;
            qBlend2Interp = solverInput.qCathodeBlend2Interp;
            qBlend1Interp = solverInput.qCathodeBlend1Interp;
        otherwise
            error('Unsupported electrode type "%s". Use ''anode'' or ''cathode''.', electrode);
    end

    if isempty(commonVoltage) || isempty(qBlend1Interp) || isempty(qBlend2Interp)
        error('calculate_blend_curve: missing blend data for electrode "%s".', electrode);
    end

    % 1) Weighted sum
    qBlend = gammaBlend2 * qBlend2Interp + (1 - gammaBlend2) * qBlend1Interp;

    % 2) Normalize to 0..1 so it behaves like an SOC
    minQ = min(qBlend);
    maxQ = max(qBlend);
    qNorm = (qBlend - minQ) / (maxQ - minQ);

    % 3) Sort qNorm so we can invert cleanly
    [qSorted, idx] = sort(qNorm);
    vSorted = commonVoltage(idx);

    % 4) Create a uniform 0..1 SOC, and map to the corresponding voltage
    blendSOC    = linspace(0, 1, length(qSorted));
    blendVoltage= interp1(qSorted, vSorted, blendSOC, 'linear','extrap');

end
