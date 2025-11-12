function Diff_DVA = fit_DVA_mae_blend(X, myData, Q0, ROI_DVA_min, ROI_DVA_max)
%{
%% Original Comments
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-03-11
%
% x = [alpha_an, beta_an, alpha_cat, beta_cat] or 
%     [alpha_an, beta_an, alpha_cat, beta_cat, gamma_Blend2]
%
% If x has only 4 elements, then gamma_Blend2 is assumed to be 0 (pure Blend1 mode).
%
% This function:
%   1) Builds the anode curve (possibly a blend) based on x(5) if present.
%   2) Interpolates the potentials of anode, cathode, and measured OCV on the 
%      same Q-grid (myData.Q_cell).
%   3) Computes the discrete dU/dQ (DVA) for anode, cathode, and OCV.
%   4) Forms DVA_Sum = DVA_Cathode - DVA_Anode and compares it to DVA_OCV
%      only within the specified region of interest [ROI_DVA_min, ROI_DVA_max]
%      by summing the squared differences. The output Diff_DVA is that sum.
% 
% In this vectorized version:
%   - X is [N×4] or [N×5], one row per solution
%   - We return an N×1 vector of errors
%}
    % Handle dimension
    if size(X,1)==1 && size(X,2) <= 7
        X = reshape(X,1,[]);
    end
    N = size(X,1);

    Q = myData.Q_cell;   
    nQ = length(Q);

    % We want to produce N errors
    Diff_DVA = zeros(N,1);

    % Build a single mask for the ROI
    mask = (Q >= ROI_DVA_min & Q <= ROI_DVA_max);

    % Precompute measured OCV at Q-grid
    OCV_DVA = interp1(Q, myData.OCV_cell, Q, 'linear', 0);

    % We'll store anode potentials in a [nQ x N] array, likewise cathode
    anodePotAll   = zeros(nQ, N);
    cathodePotAll = zeros(nQ, N);

    %----- 1) For each solution, build the anode curve and do interpolation -----
    for i = 1:N
        x_i = X(i,:);
        alpha_an  = x_i(1);
        beta_an   = x_i(2);
        alpha_cat = x_i(3);
        beta_cat  = x_i(4);

        gamma_Blend2  = 0;
        inhomo_mag_an = 0;
        inhomo_mag_ca = 0;
        if length(x_i) == 5
            gamma_Blend2 = x_i(5);
        elseif length (x_i) > 5
            gamma_Blend2 = x_i(5);
            inhomo_mag_an = x_i(6);
            inhomo_mag_ca = x_i(7);
        end

        % Build or blend the anode curve
        if gamma_Blend2 ~= 0
            [blendSOC, blendU] = blend_anode_curve(gamma_Blend2, myData, inhomo_mag_an);
            anodePotAll(:,i) = interp1(alpha_an*blendSOC + beta_an, blendU, Q, 'linear', 0);
        else
            % Pure Blend1 approach
            anodePotAll(:,i) = interp1(alpha_an*myData.Q_Blend2_interp + beta_an, ...
                                       myData.commonVoltage, Q, 'linear', 0);
        end
        
        % Calculate inhomogeneities
        if ~(inhomo_mag_ca == 0)
            myData.normCathode_U = calculate_inhomogeneity(myData.normCathode_SOC, myData.normCathode_U, inhomo_mag_ca);
        end

        % Build the cathode potential
        cathodePotAll(:,i) = interp1(alpha_cat*myData.normCathode_SOC + beta_cat, ...
                                     myData.normCathode_U, Q, 'linear', 0);
    end

    %----- 2) Compute discrete DVA for each solution -----
    % We'll store DVA in a [nQ x N] array as well
    dva_anode   = zeros(nQ, N);
    dva_cathode = zeros(nQ, N);
    dva_ocv     = zeros(nQ, 1);  % same for all solutions

    % The measured OCV is only 1 curve, so we do this once
    for idx = 2:nQ
        dU = OCV_DVA(idx) - OCV_DVA(idx-1);
        dQ = Q(idx) - Q(idx-1);
        dva_ocv(idx-1) = (dU / dQ) * Q0;
    end

    dva_ocv = filterIt(dva_ocv, 'filtermethod', 'sgolay');

    % Now do anode & cathode in a vectorized for-loop style
    for i = 2:nQ
        dU_an = anodePotAll(i,:) - anodePotAll(i-1,:);
        dU_cat= cathodePotAll(i,:) - cathodePotAll(i-1,:);
        dQ    = Q(i) - Q(i-1);

        dva_anode(i-1,:)   = (dU_an  / dQ) * Q0;  % row vector
        dva_cathode(i-1,:) = (dU_cat / dQ) * Q0;
    end

    % Now DVA_Sum = Cathode - Anode, shape [nQ x N]
    dva_sum = dva_cathode - dva_anode;

    dva_sum = filterIt(dva_sum, 'filtermethod', 'sgolay');

    %----- 3) Compare only within [ROI_DVA_min, ROI_DVA_max] -----
    % We form squared differences for each solution
    %   - replicate dva_ocv across columns to subtract from dva_sum
    %   - or do a simple loop
    ocvRep = repmat(dva_ocv, 1, N);  % [nQ x N]
    diffArray = (dva_sum - ocvRep).^2;  % [nQ x N]

    % Mask out points outside ROI
    diffArray(~mask,:) = 0;

    % Sum over Q dimension
    DiffVals = sum(diffArray,1);  % 1xN
    Diff_DVA = DiffVals(:);       % Nx1
    Diff_DVA = Diff_DVA / sum(mask);
end
