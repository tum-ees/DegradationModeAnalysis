function [LAM_Anode, LAM_Cathode, LI, LAM_Anode_Blend2, LAM_Anode_Blend1] = calculate_degradation_modes(params, Capa_act, Capa_Anode_init, Capa_Cathode_init, Capa_Inventory_init, gamma_Blend2_init, varargin)
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
    
    alpha_Anode   = params(1);
    beta_Anode    = params(2);
    alpha_Cathode = params(3);
    beta_Cathode  = params(4);
    
    if numel(params) < 5
        gamma_Blend2 = 0;
    else
        gamma_Blend2 = params(5);
    end

    % Compute current anode capacity
    Capa_Anode = alpha_Anode * Capa_act;

    % Initial sub-capacities based on reference gamma_Blend2_init
    Capa_Anode_Blend2_init = Capa_Anode_init * gamma_Blend2_init;
    Capa_Anode_Blend1_init = Capa_Anode_init * (1 - gamma_Blend2_init);
    
    % Current sub-capacities using the optimized gamma_Blend2 (or zero in pure Blend1 mode)
    Capa_Anode_Blend2 = Capa_Anode * gamma_Blend2;
    Capa_Anode_Blend1 = Capa_Anode * (1 - gamma_Blend2);
    
    % Calculate loss for Blend2 and Blend1 parts separately
    LAM_Anode_Blend2 = (Capa_Anode_Blend2_init - Capa_Anode_Blend2) / Capa_Anode_Blend2_init;
    LAM_Anode_Blend1 = (Capa_Anode_Blend1_init - Capa_Anode_Blend1) / Capa_Anode_Blend1_init;
    
    % Overall anode and cathode degradation
    LAM_Anode = (Capa_Anode_init - Capa_Anode) / Capa_Anode_init;
    Capa_Cathode = alpha_Cathode * Capa_act;
    LAM_Cathode = (Capa_Cathode_init - Capa_Cathode) / Capa_Cathode_init;
    
    % Inventory loss calculation
    Capa_Inventory = (alpha_Cathode + beta_Cathode - beta_Anode) * Capa_act;
    LI = (Capa_Inventory_init - Capa_Inventory) / Capa_Inventory_init;

    if fitReverse
        LAM_Cathode = -LAM_Cathode*Capa_Cathode_init/Capa_Cathode;
        LAM_Anode = -LAM_Anode*Capa_Anode_init/Capa_Anode;
        LAM_Anode_Blend1 = -LAM_Anode_Blend1*Capa_Anode_Blend1_init/Capa_Anode_Blend1;
        LAM_Anode_Blend2 = -LAM_Anode_Blend2*Capa_Anode_Blend2_init/Capa_Anode_Blend2;
    end

    
end