function solStruct = store_solution_struct(myData, reconSOC, fcU_model, ...
    Q_DVA_meas, DVA_smooth_meas, Q_DVA_calc, DVA_smooth_calc, ...
    Q_ICA_meas, ICA_smooth_meas, Q_ICA_calc, ICA_smooth_calc, ...
    Capa_act, params, Algorithm, weightOCV, weightDVA, weightICA, ...
    cathSOC, cathU_recon, anodeSOC, anodeU_recon,...
    rmse_fitRegion, rmse_fullRange, rmse_DVA_full)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%
% A helper function that packages results into one struct.  
% That way, we can easily compare solutions or store the best.
%
    solStruct.myData           = myData;
    solStruct.reconSOC        = reconSOC;
    solStruct.fcU_model       = fcU_model;
    solStruct.Q_DVA_meas      = Q_DVA_meas;
    solStruct.DVA_smooth_meas = DVA_smooth_meas;
    solStruct.Q_DVA_calc      = Q_DVA_calc;
    solStruct.DVA_smooth_calc = DVA_smooth_calc;
    solStruct.Q_ICA_meas      = Q_ICA_meas;
    solStruct.ICA_smooth_meas = ICA_smooth_meas;
    solStruct.Q_ICA_calc      = Q_ICA_calc;
    solStruct.ICA_smooth_calc = ICA_smooth_calc;
    solStruct.Capa_act        = Capa_act;
    solStruct.params          = params;
    solStruct.Algorithm       = Algorithm;
    % solStruct.window          = window;
    solStruct.weightOCV       = weightOCV;
    solStruct.weightDVA       = weightDVA;
    solStruct.weightICA       = weightICA;
    solStruct.cathSOC         = cathSOC;
    solStruct.cathU_recon     = cathU_recon;
    solStruct.anodeSOC        = anodeSOC;
    solStruct.anodeU_recon    = anodeU_recon;
    solStruct.rmse_fitRegion  = rmse_fitRegion;
    solStruct.rmse_fullRange  = rmse_fullRange;
    solStruct.rmse_DVA_full   = rmse_DVA_full;
end