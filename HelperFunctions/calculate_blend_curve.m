function [blendSOC, blendVoltage] = calculate_blend_curve(gamma_blend2, myData, electrode)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%
% This function calculates the SOC and voltage of a blended electrode curve (anode or cathode).
%
% Inputs:
%   gamma_blend2 : blend share of component 2
%   myData       : struct containing commonVoltage / Q_* fields
%   electrode    : 'anode' (default) or 'cathode'
%
    if nargin < 3 || isempty(electrode)
        electrode = 'anode';
    end

    switch lower(electrode)
        case 'anode'
            commonVoltage   = myData.commonVoltage_anode;
            Q_blend2_interp = myData.Q_anode_blend2_interp;
            Q_blend1_interp = myData.Q_anode_blend1_interp;
        case 'cathode'
            commonVoltage   = myData.commonVoltage_cathode;
            Q_blend2_interp = myData.Q_cathode_blend2_interp;
            Q_blend1_interp = myData.Q_cathode_blend1_interp;
        otherwise
            error('Unsupported electrode type "%s". Use ''anode'' or ''cathode''.', electrode);
    end

    if isempty(commonVoltage) || isempty(Q_blend1_interp) || isempty(Q_blend2_interp)
        error('calculate_blend_curve: missing blend data for electrode "%s".', electrode);
    end

    % 1) Weighted sum
    Q_blend = gamma_blend2 * Q_blend2_interp + (1 - gamma_blend2)*Q_blend1_interp;

    % 2) Normalize to 0..1 so it behaves like an SOC
    minQ = min(Q_blend);
    maxQ = max(Q_blend);
    Q_norm = (Q_blend - minQ) / (maxQ - minQ);

    % 3) Sort Q_norm so we can invert cleanly
    [Q_sorted, idx] = sort(Q_norm);
    V_sorted = commonVoltage(idx);

    % 4) Create a uniform 0..1 SOC, and map to the corresponding voltage
    blendSOC    = linspace(0, 1, length(Q_sorted));
    blendVoltage= interp1(Q_sorted, V_sorted, blendSOC, 'linear','extrap');

end
