function out = norm_blend_voltage(inStruct,Vmin,Vmax)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%
% This function trims Blend1 & Blend2 to equal voltage max/min
%
    % Trim to overlapping voltage region and rescale capacity to [0…1]
    mask = inStruct.voltage >= Vmin & inStruct.voltage <= Vmax;

    inStruct.voltage            = inStruct.voltage(mask);
    inStruct.normalizedCapacity = inStruct.normalizedCapacity(mask);

    inStruct.normalizedCapacity = rescale(inStruct.normalizedCapacity,0,1);
    out = inStruct;
end


