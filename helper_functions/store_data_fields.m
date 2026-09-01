function data = store_data_fields(data, fName, solution, lamAnode, lamCathode, li, lamAnodeBlend2, lamAnodeBlend1, lamCathodeBlend2, lamCathodeBlend1)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%
% A helper function that packages results into one data struct.  
% That way, we can easily compare solutions or store the best.
%
% INPUTS:
%   data               - aggregate DMA output struct
%   fName              - target CU field name, for example `CU1`
%   solution           - best-solution struct for the current CU
%   lamAnode           - aggregate anode LAM
%   lamCathode         - aggregate cathode LAM
%   li                 - loss of lithium inventory
%   lamAnodeBlend2     - blend-2 contribution to anode LAM
%   lamAnodeBlend1     - blend-1 contribution to anode LAM
%   lamCathodeBlend2   - blend-2 contribution to cathode LAM
%   lamCathodeBlend1   - blend-1 contribution to cathode LAM
%
% OUTPUTS:
%   data               - updated aggregate DMA output struct

    data.(fName).params                   = solution.params;
    % The active capacity of this check-up. Every degradation mode is
    % relative to it, and the resistance offset is bounded per capacity, so
    % without it a stored result cannot be read back and checked.
    data.(fName).capaAct                  = solution.capaAct;
    data.(fName).measured.soc             = solution.solverInput.qCell;
    data.(fName).measured.voltage         = solution.solverInput.ocvCell;
    data.(fName).measured.qDVA            = solution.qDVAMeas;
    data.(fName).measured.dva             = solution.dvaSmoothMeas;
    data.(fName).measured.qICA            = solution.qICAMeas;
    data.(fName).measured.ica             = solution.icaSmoothMeas;
    data.(fName).calculated.soc           = solution.reconSOC;
    data.(fName).calculated.voltage       = solution.fullCellUModel;
    data.(fName).calculated.qDVA          = solution.qDVACalc;
    data.(fName).calculated.dva           = solution.dvaSmoothCalc;
    data.(fName).calculated.qICA          = solution.qICACalc;
    data.(fName).calculated.ica           = solution.icaSmoothCalc;
    data.(fName).calculated.cathodeSOC    = solution.cathodeSOC;
    data.(fName).calculated.cathodeURecon = solution.cathodeURecon;
    data.(fName).calculated.anodeSOC      = solution.anodeSOC;
    data.(fName).calculated.anodeURecon   = solution.anodeURecon;
    data.(fName).lamAnode                 = lamAnode;
    data.(fName).lamCathode               = lamCathode;
    data.(fName).li                       = li;
    data.(fName).lamAnodeBlend2           = lamAnodeBlend2;
    data.(fName).lamAnodeBlend1           = lamAnodeBlend1;
    data.(fName).lamCathodeBlend2         = lamCathodeBlend2;
    data.(fName).lamCathodeBlend1         = lamCathodeBlend1;
    data.(fName).rmseFitRegion            = solution.rmseFitRegion;
    data.(fName).rmseFullRange            = solution.rmseFullRange;
    data.(fName).rmseDVAFull              = solution.rmseDVAFull;
end

