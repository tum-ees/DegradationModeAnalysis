function Diff_ICA = fit_ICA_mae_blend(X, myData, Q0, ROI_ICA_min, ROI_ICA_max)
%> Author: Josef Eizenhammer (josef.eizenhammer@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-05-19
%
% x = [alpha_an, beta_an, alpha_cat, beta_cat] or 
%     [alpha_an, beta_an, alpha_cat, beta_cat, gamma_Blend2] or
%      [alpha_an, beta_an, alpha_cat, beta_cat, gamma_Blend2, anodeInhomogenity]
%
% If x has only 4 elements, then gamma_Blend2 is assumed to be 0 (pure Blend1 mode).
%
% This function:
%   1) Builds the anode curve (possibly a blend) based on x(5) if present.
%   2) Interpolates the potentials of anode, cathode, and measured OCV on the 
%      same Q-grid (myData.Q_cell).
%   3) Computes the discrete dQ/dU (ICA) for anode, cathode, and OCV.
%   4) Forms ICA_Sum = ICA_Cathode - ICA_Anode and compares it to ICA_OCV
%      only within the specified region of interest [ROI_ICA_min, ROI_ICA_max]
%      by summing the squared differences. The output Diff_ICA is that sum.
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


    % Build a single mask for the ROI
    mask = (Q >= ROI_ICA_min & Q <= ROI_ICA_max);

    % Precompute measured OCV at Q-grid
    OCV_ICA = interp1(Q, myData.OCV_cell, Q, 'linear', 0);

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

    %----- 2) Compute discrete ICA for the measurement -----

    [~, ~, ICA_OCV] = calculate_ICA(Q, OCV_ICA);
    ICA_OCV = ICA_OCV/Q0;
    ICA_OCV = filterIt(ICA_OCV, 'filtermethod', 'sgolay');

    %----- 3) Compute discrete ICA for the calculated curves -----
    
    % initilaize and calculate OCV of the calcuted anode and cathode
    % potential curves
    OCV_sum = cathodePotAll(:,:)-anodePotAll(:,:);
    [~, ~, ICA_calc] = calculate_ICA(Q, OCV_sum);
    ICA_calc = ICA_calc/Q0;
    ICA_calc = filterIt(ICA_calc, 'filtermethod', 'sgolay');

    %----- 4) Compare only within [ROI_ICA_min, ROI_ICA_max] -----
    % We form squared differences for each solution
    %   - replicate ICA_ocv across columns to subtract from ICA_sum
    %   - or do a simple loop
    ICA_OCV = repmat(ICA_OCV, 1, N);  % [nQ x N]
    ICA_calc = repmat(ICA_calc, 1, N);
    diffArray = (ICA_calc - ICA_OCV).^2;  % [nQ x N]

    % Mask out points outside ROI
    diffArray(~mask,:) = 0;

    % Sum over Q dimension
    DiffVals = sum(diffArray,1);  % 1xN
    Diff_ICA = DiffVals(:);
    Diff_ICA = Diff_ICA / sum(mask);
    
end
