function solution = store_solution_struct(solverInput, reconSOC, fullCellUModel, ...
    qDVAMeas, dvaSmoothMeas, qDVACalc, dvaSmoothCalc, ...
    qICAMeas, icaSmoothMeas, qICACalc, icaSmoothCalc, ...
    capaAct, params, algorithm, weightOCV, weightDVA, weightICA, ...
    cathodeSOC, cathodeURecon, anodeSOC, anodeURecon,...
    rmseFitRegion, rmseFullRange, rmseDVAFull)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%
% A helper function that packages results into one struct.  
% That way, we can easily compare solutions or store the best.
%
% INPUTS:
%   solverInput      - preprocessed solver input bundle for the current CU
%   reconSOC         - reconstructed full-cell SOC vector
%   fullCellUModel   - reconstructed full-cell voltage vector
%   qDVAMeas         - measured DVA x-axis
%   dvaSmoothMeas    - measured smoothed DVA curve
%   qDVACalc         - calculated DVA x-axis
%   dvaSmoothCalc    - calculated smoothed DVA curve
%   qICAMeas         - measured ICA x-axis
%   icaSmoothMeas    - measured smoothed ICA curve
%   qICACalc         - calculated ICA x-axis
%   icaSmoothCalc    - calculated smoothed ICA curve
%   capaAct          - fitted active full-cell capacity
%   params           - fitted DMA parameter vector
%   algorithm        - solver label used to obtain the solution
%   weightOCV        - OCV objective weight
%   weightDVA        - DVA objective weight
%   weightICA        - ICA objective weight
%   cathodeSOC       - cathode SOC axis for electrode reconstruction
%   cathodeURecon    - cathode potential curve
%   anodeSOC         - anode SOC axis for electrode reconstruction
%   anodeURecon      - anode potential curve
%   rmseFitRegion    - RMSE over the configured fit region
%   rmseFullRange    - RMSE over the full SOC range
%   rmseDVAFull      - full-range DVA RMSE
%
% OUTPUTS:
%   solution         - packaged solution struct used throughout the framework

    solution.solverInput   = solverInput;
    solution.reconSOC      = reconSOC;
    solution.fullCellUModel = fullCellUModel;
    solution.qDVAMeas      = qDVAMeas;
    solution.dvaSmoothMeas = dvaSmoothMeas;
    solution.qDVACalc      = qDVACalc;
    solution.dvaSmoothCalc = dvaSmoothCalc;
    solution.qICAMeas      = qICAMeas;
    solution.icaSmoothMeas = icaSmoothMeas;
    solution.qICACalc      = qICACalc;
    solution.icaSmoothCalc = icaSmoothCalc;
    solution.capaAct       = capaAct;
    solution.params        = params;
    solution.algorithm     = algorithm;
    solution.weightOCV     = weightOCV;
    solution.weightDVA     = weightDVA;
    solution.weightICA     = weightICA;
    solution.cathodeSOC    = cathodeSOC;
    solution.cathodeURecon = cathodeURecon;
    solution.anodeSOC      = anodeSOC;
    solution.anodeURecon   = anodeURecon;
    solution.rmseFitRegion = rmseFitRegion;
    solution.rmseFullRange = rmseFullRange;
    solution.rmseDVAFull   = rmseDVAFull;
end
