function [LAM_anode, LAM_cathode, LI, LAM_anode_blend2, LAM_anode_blend1, LAM_cathode_blend2, LAM_cathode_blend1] = calculate_degradation_modes(params, capa_act, capa_anode_init, capa_cathode_init, capa_inventory_init, gamma_an_blend2_init, gamma_ca_blend2_init, varargin)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%
% This function calculates the degradation modes for the anode, cathode and
% LI and additionally for the two blend materials based on the obatined
% parameters.
%
    % Parse optional flag
    fitReverse = false;
    if ~isempty(varargin)
        fitReverse = logical(varargin{1});
    end
    
    alpha_anode   = params(1);
    beta_anode    = params(2);
    alpha_cathode = params(3);
    beta_cathode  = params(4);
    
    % if ~blend mode -> params(5/6) == 0
    gamma_an_blend2 = params(5);
    gamma_ca_blend2 = params(6); % cathode blend share

    % Compute current anode capacity
    capa_anode = alpha_anode * capa_act;
    capa_cathode = alpha_cathode * capa_act;

    % Initial sub-capacities based on reference gamma_an_blend2_init
    capa_anode_blend2_init = capa_anode_init * gamma_an_blend2_init;
    capa_anode_blend1_init = capa_anode_init * (1 - gamma_an_blend2_init);
    capa_cathode_blend2_init = capa_cathode_init * gamma_ca_blend2_init;
    capa_cathode_blend1_init = capa_cathode_init * (1 - gamma_ca_blend2_init);
    
    % Current sub-capacities using the optimized gamma_an_blend2 (or zero in pure Blend1 mode)
    capa_anode_blend2 = capa_anode * gamma_an_blend2;
    capa_anode_blend1 = capa_anode * (1 - gamma_an_blend2);
    capa_cathode_blend2 = capa_cathode * gamma_ca_blend2;
    capa_cathode_blend1 = capa_cathode * (1 - gamma_ca_blend2);
    
    % Calculate loss for Blend2 and Blend1 parts separately
    LAM_anode_blend2 = safe_loss(capa_anode_blend2_init, capa_anode_blend2);
    LAM_anode_blend1 = safe_loss(capa_anode_blend1_init, capa_anode_blend1);
    LAM_cathode_blend2 = safe_loss(capa_cathode_blend2_init, capa_cathode_blend2);
    LAM_cathode_blend1 = safe_loss(capa_cathode_blend1_init, capa_cathode_blend1);
    
    % Overall anode and cathode degradation
    LAM_anode = (capa_anode_init - capa_anode) / capa_anode_init;
    LAM_cathode = (capa_cathode_init - capa_cathode) / capa_cathode_init;
    
    % Inventory loss calculation
    capa_inventory = (alpha_cathode + beta_cathode - beta_anode) * capa_act;
    LI = (capa_inventory_init - capa_inventory) / capa_inventory_init;

    if fitReverse
        LAM_cathode = -LAM_cathode*capa_cathode_init/capa_cathode;
        LAM_anode = -LAM_anode*capa_anode_init/capa_anode;
        LAM_anode_blend1 = -LAM_anode_blend1*capa_anode_blend1_init/capa_anode_blend1;
        LAM_anode_blend2 = -LAM_anode_blend2*capa_anode_blend2_init/capa_anode_blend2;
        if capa_cathode_blend1 ~= 0
            LAM_cathode_blend1 = -LAM_cathode_blend1*capa_cathode_blend1_init/capa_cathode_blend1;
        end
        if capa_cathode_blend2 ~= 0
            LAM_cathode_blend2 = -LAM_cathode_blend2*capa_cathode_blend2_init/capa_cathode_blend2;
        end
    end

    % nested helper to avoid divide-by-zero
    function lossVal = safe_loss(initVal, currentVal)
        if initVal == 0
            lossVal = 0;
        else
            lossVal = (initVal - currentVal) / initVal;
        end
    end
    
end
