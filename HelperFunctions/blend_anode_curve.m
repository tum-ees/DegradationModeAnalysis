function [blendSOC, blendVoltage] = blend_anode_curve(gamma_Blend2, myData, inhomo_mag)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%
% This function calculates the SOC and voltage of the blend model
% 
    % Unpack
    commonVoltage = myData.commonVoltage;
    Q_Blend2_interp   = myData.Q_Blend2_interp;
    Q_Blend1_interp   = myData.Q_Blend1_interp;

    % 1) Weighted sum
    Q_blend = gamma_Blend2 * Q_Blend2_interp + (1 - gamma_Blend2)*Q_Blend1_interp;

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

    % 5) Optional: compute inhomgenity
    if ~(inhomo_mag == 0)
        blendVoltage = calculate_inhomogeneity(blendSOC, blendVoltage, inhomo_mag);
    end 

end