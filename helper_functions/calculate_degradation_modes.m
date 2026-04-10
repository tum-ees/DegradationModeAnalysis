function [lamAnode, lamCathode, li, lamAnodeBlend2, lamAnodeBlend1, lamCathodeBlend2, lamCathodeBlend1] = calculate_degradation_modes(params, capaAct, capaAnodeInit, capaCathodeInit, capaInventoryInit, gammaAnBlend2Init, gammaCaBlend2Init, varargin)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%
% This function calculates the degradation modes for the anode, cathode and
% li and additionally for the two blend materials based on the obtained
% parameters.
%
% INPUTS:
%   params             - fitted DMA parameter vector
%   capaAct            - fitted active full-cell capacity
%   capaAnodeInit      - reference beginning-of-life anode capacity
%   capaCathodeInit    - reference beginning-of-life cathode capacity
%   capaInventoryInit  - reference beginning-of-life lithium inventory
%   gammaAnBlend2Init  - reference anode blend-2 fraction
%   gammaCaBlend2Init  - reference cathode blend-2 fraction
%   varargin           - optional logical flag for reverse-fit handling
%
% OUTPUTS:
%   lamAnode           - aggregate anode loss of active material
%   lamCathode         - aggregate cathode loss of active material
%   li                 - loss of lithium inventory
%   lamAnodeBlend2     - blend-2 contribution to anode LAM
%   lamAnodeBlend1     - blend-1 contribution to anode LAM
%   lamCathodeBlend2   - blend-2 contribution to cathode LAM
%   lamCathodeBlend1   - blend-1 contribution to cathode LAM
%
    % Parse optional flag
    fitReverse = false;
    if ~isempty(varargin)
        fitReverse = logical(varargin{1});
    end
    
    alphaAnode   = params(1);
    betaAnode    = params(2);
    alphaCathode = params(3);
    betaCathode  = params(4);
    
    % if ~blend mode -> params(5/6) == 0
    gammaAnBlend2 = params(5);
    gammaCaBlend2 = params(6); % cathode blend share

    % Compute current anode capacity
    capaAnode = alphaAnode * capaAct;
    capaCathode = alphaCathode * capaAct;

    % Initial sub-capacities based on reference gammaAnBlend2Init
    capaAnodeBlend2Init = capaAnodeInit * gammaAnBlend2Init;
    capaAnodeBlend1Init = capaAnodeInit * (1 - gammaAnBlend2Init);
    capaCathodeBlend2Init = capaCathodeInit * gammaCaBlend2Init;
    capaCathodeBlend1Init = capaCathodeInit * (1 - gammaCaBlend2Init);
    
    % Current sub-capacities using the optimized gammaAnBlend2 (or zero in pure Blend1 mode)
    capaAnodeBlend2 = capaAnode * gammaAnBlend2;
    capaAnodeBlend1 = capaAnode * (1 - gammaAnBlend2);
    capaCathodeBlend2 = capaCathode * gammaCaBlend2;
    capaCathodeBlend1 = capaCathode * (1 - gammaCaBlend2);
    
    % Calculate loss for Blend2 and Blend1 parts separately
    lamAnodeBlend2 = safeLoss(capaAnodeBlend2Init, capaAnodeBlend2);
    lamAnodeBlend1 = safeLoss(capaAnodeBlend1Init, capaAnodeBlend1);
    lamCathodeBlend2 = safeLoss(capaCathodeBlend2Init, capaCathodeBlend2);
    lamCathodeBlend1 = safeLoss(capaCathodeBlend1Init, capaCathodeBlend1);
    
    % Overall anode and cathode degradation
    lamAnode = (capaAnodeInit - capaAnode) / capaAnodeInit;
    lamCathode = (capaCathodeInit - capaCathode) / capaCathodeInit;
    
    % Inventory loss calculation
    capaInventory = (alphaCathode + betaCathode - betaAnode) * capaAct;
    li = (capaInventoryInit - capaInventory) / capaInventoryInit;

    if fitReverse
        lamCathode = -lamCathode * capaCathodeInit / capaCathode;
        lamAnode = -lamAnode * capaAnodeInit / capaAnode;
        lamAnodeBlend1 = -lamAnodeBlend1 * capaAnodeBlend1Init / capaAnodeBlend1;
        lamAnodeBlend2 = -lamAnodeBlend2 * capaAnodeBlend2Init / capaAnodeBlend2;
        if capaCathodeBlend1 ~= 0
            lamCathodeBlend1 = -lamCathodeBlend1 * capaCathodeBlend1Init / capaCathodeBlend1;
        end
        if capaCathodeBlend2 ~= 0
            lamCathodeBlend2 = -lamCathodeBlend2 * capaCathodeBlend2Init / capaCathodeBlend2;
        end
    end

    % nested helper to avoid divide-by-zero
    function lossVal = safeLoss(initVal, currentVal)
        if initVal == 0
            lossVal = 0;
        else
            lossVal = (initVal - currentVal) / initVal;
        end
    end
    
end
