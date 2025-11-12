function [Data] = store_data_fields(Data, fName, sol, LAM_Anode, LAM_Cathode, LI, LAM_Anode_Blend2, LAM_Anode_Blend1)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Date: 2025-09-11
%
% A helper function that packages results into one data struct.  
% That way, we can easily compare solutions or store the best.
%
    Data.(fName).params                  = sol.params;
    Data.(fName).measured.SOC            = sol.myData.Q_cell;
    Data.(fName).measured.voltage        = sol.myData.OCV_cell;
    Data.(fName).measured.Q_DVA          = sol.Q_DVA_meas;
    Data.(fName).measured.DVA            = sol.DVA_smooth_meas;
    Data.(fName).measured.Q_ICA          = sol.Q_ICA_meas;
    Data.(fName).measured.ICA            = sol.ICA_smooth_meas;
    Data.(fName).calculated.SOC          = sol.reconSOC;
    Data.(fName).calculated.voltage      = sol.fcU_model;
    Data.(fName).calculated.Q_DVA        = sol.Q_DVA_calc;
    Data.(fName).calculated.DVA          = sol.DVA_smooth_calc;
    Data.(fName).calculated.Q_ICA        = sol.Q_ICA_calc;
    Data.(fName).calculated.ICA          = sol.ICA_smooth_calc;
    Data.(fName).calculated.cathSOC      = sol.cathSOC;
    Data.(fName).calculated.cathU_recon  = sol.cathU_recon;
    Data.(fName).calculated.anodeSOC     = sol.anodeSOC;
    Data.(fName).calculated.anodeU_recon = sol.anodeU_recon;
    Data.(fName).LAM_Anode               = LAM_Anode;
    Data.(fName).LAM_Cathode             = LAM_Cathode;
    Data.(fName).LI                      = LI;
    Data.(fName).LAM_Anode_Blend2        = LAM_Anode_Blend2;
    Data.(fName).LAM_Anode_Blend1        = LAM_Anode_Blend1;
    Data.(fName).RMSE_fitRegion          = sol.finalRMSE_fit;
    Data.(fName).RMSE_0_9                = sol.finalRMSE_full;
    Data.(fName).RMSE_DVA                = sol.rmse_DVA_full;
end

