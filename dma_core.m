function [myData, reconSOC, fcU_model, Q_DVA_meas, DVA_smooth_meas, ...
    Q_DVA_calc, DVA_smooth_calc, Q_ICA_meas, ICA_smooth_meas, ...
    Q_ICA_calc, ICA_smooth_calc, capa_act, params, AlgorithmOut, ...
    weightOCV_Out, cathSOC, normCathode_U, anodeSOC, blendU] = ...
    dma_core(half_and_full_cell_data, settings, refData, ...
    LAM_Anode_prev, LAM_Cathode_prev, LAM_Anode_Blend1_prev, ...
    LAM_Anode_Blend2_prev, Inhmg_An_prev, Inhmg_Ca_prev, ...
    allowAnodeInhomogeneity, allowCathodeInhomogeneity, fitReverse)
%> -------------------------------------------------------------------------
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Sebastian Karl (s.karl@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Additional code by Moritz Guenthner (moritz.guenthner@tum.de)
%> Additional code by Mathias Rehm (mathias.rehm@tum.de)

%> Date: 2025-03-11
%
% OVERVIEW:
%   * Combines (blends) anode and cathode data to reconstruct a full-cell OCV.
%   * Uses user-defined or pre-computed half-cell and full-cell data (either
%     in a struct or cell-array format, or a .mat path).
%   * Always parses data via parse_data_input, even for Blend1 and Blend2,
%     so that user inputs are fully flexible (any data could be in any format).
%   * Defines objective functions (OCV- and DVA-based) for optimization, then
%     runs a chosen optimization algorithm (e.g., 'ga', 'fmincon', etc.).
%   * Returns the fitted parameters, reconstructed curves (OCV, DVA), and
%     relevant data used for further analysis.
%
% Inputs:
%   half_and_full_cell_data
%     - struct containing the preprocessed half- and full-cell data.
%   settings
%     - settings struct defined in main_DMA.
%   refData
%     - optional struct with reference capacities for penalties in
%       objectiveWithPenalty; must include .Capa_Anode_init,
%       .Capa_Cathode_init, .Capa_Inventory_init, and .gamma_Blend2_init.
%   LAM_Anode_prev
%     - optional previous anode loss for penalties; empty disables the
%       negative-anode penalty.
%   LAM_Cathode_prev
%     - optional previous cathode loss; empty disables the negative-cathode
%       penalty.
%   allowAnodeInhomogeneity / allowCathodeInhomogeneity
%     - logical flags that permit inhomogeneity parameters.
%   fitReverse
%     - if true, process CUs in reverse order (last CU first).
%
% Outputs:
%   myData
%     - struct with required OCV/DVA data (normCathode_SOC, Q_cell, etc.).
%   reconSOC
%     - uniform SOC array (0..1) used for reconstruction.
%   fcU_model
%     - reconstructed full-cell OCV from the model.
%   Q_DVA_meas
%     - charge grid for measured DVA.
%   DVA_smooth_meas
%     - smoothed measured DVA array.
%   Q_DVA_calc
%     - charge grid for calculated (modeled) DVA.
%   DVA_smooth_calc
%     - smoothed calculated DVA array.
%   Capa_act
%     - actual capacity from the input data (e.g., fullCellData).
%   params
%     - optimised parameters [alpha_an, beta_an, alpha_cat, beta_cat,
%       (gamma_Blend2)].
%   AlgorithmOut
%     - echo of the chosen algorithm.
%   weightOCV_Out
%     - echo of the chosen OCV weight.
%   cathSOC           - Cathode SOC axis used in final reconstruction.
%   normCathode_U     - Interpolated cathode potential (smooth).
%   anodeSOC          - Anode SOC axis used in final reconstruction.
%   blendU            - The blended anode voltage array.
%

% ======= SETTINGS / CONSTANTS / PREPROCESSING =======
warning('off', 'all');

% Get Settings from settings struct
smoothingPoints         = settings.smoothingPoints;
dataLength              = settings.dataLength;
Algorithm               = settings.Algorithm;
weightOCV               = settings.weightOCV;
weightDVA               = settings.weightDVA;
weightICA               = settings.weightICA;
ROI_OCV_min             = settings.ROI_OCV_min;
ROI_OCV_max             = settings.ROI_OCV_max;
ROI_DVA_min             = settings.ROI_DVA_min;
ROI_DVA_max             = settings.ROI_DVA_max;
ROI_ICA_min             = settings.ROI_ICA_min;
ROI_ICA_max             = settings.ROI_ICA_max;
useBlendElectrodeModel  = settings.useBlend;
maxAnodeGain            = settings.maxAnodeGain;
maxCathodeGain          = settings.maxCathodeGain;
maxBlend1Gain           = settings.maxBlend1Gain;
maxBlend2Gain           = settings.maxBlend2Gain;
maxAnodeLoss            = settings.maxAnodeLoss;
maxCathodeLoss          = settings.maxCathodeLoss;
maxBlend1Loss           = settings.maxBlend1Loss;
maxBlend2Loss           = settings.maxBlend2Loss;
lowerBoundaries         = settings.LowerBoundaries;
upperBoundaries         = settings.UpperBoundaries;
gammaBlend2_upperBound  = settings.gammaBlend2_upperBound;
capa_act                = half_and_full_cell_data.Capa_act;

% Build per-electrode inhomogeneity limits upfront so the first CU can vary.
inhomoUpperBoundBase = settings.maxInhomogeneity;
if isscalar(inhomoUpperBoundBase)
    inhomoUpperBoundBase = repmat(inhomoUpperBoundBase, 1, 2);
else
    inhomoUpperBoundBase = inhomoUpperBoundBase(:).';
    if numel(inhomoUpperBoundBase) == 1
        inhomoUpperBoundBase = repmat(inhomoUpperBoundBase, 1, 2);
    else
        inhomoUpperBoundBase = inhomoUpperBoundBase(1:2);
    end
end

inhomoDeltaPerCU = settings.maxInhomogeneityDelta;
if isscalar(inhomoDeltaPerCU)
    inhomoDeltaPerCU = repmat(inhomoDeltaPerCU, 1, 2);
else
    inhomoDeltaPerCU = inhomoDeltaPerCU(:).';
    if numel(inhomoDeltaPerCU) == 1
        inhomoDeltaPerCU = repmat(inhomoDeltaPerCU, 1, 2);
    else
        inhomoDeltaPerCU = inhomoDeltaPerCU(1:2);
    end
end

% If user did not provide maxCathodeGain, default to 0
if ~exist('maxCathodeGain','var') || isempty(maxCathodeGain)
    maxCathodeGain = 0;
end

% If user did not provide LAM_Cathode_prev, default to []
if ~exist('LAM_Cathode_prev','var')
    LAM_Cathode_prev = [];
end

% -----------------------------------------------------------------------
% 4) Pack everything into myData struct, to pass to objective functions
% -----------------------------------------------------------------------
myData.normCathode_SOC = half_and_full_cell_data.normCathode_SOC;
myData.normCathode_U   = half_and_full_cell_data.normCathode_U;
myData.Q_cell          = half_and_full_cell_data.fullCell_SOC;
myData.OCV_cell        = half_and_full_cell_data.fullCell_U;
myData.commonVoltage   = half_and_full_cell_data.commonVoltage;
myData.Q_Blend2_interp = half_and_full_cell_data.Q_Blend2_interp;
myData.Q_Blend1_interp = half_and_full_cell_data.Q_Blend1_interp;
myData.Q0              = half_and_full_cell_data.Q0;

% -----------------------------------------------------------------------
% 5) Define vectorized objective functions and run optimization
% -----------------------------------------------------------------------
% Sub-objectives (each must accept [N×4 or N×5], return [N×1])
funOCV = @(X) fit_OCV_blend(X, myData, ROI_OCV_min, ROI_OCV_max);
funDVA = @(X) fit_DVA_mae_blend(X, myData, myData.Q0, ROI_DVA_min, ROI_DVA_max);
funICA = @(X) fit_ICA_mae_blend(X, myData, myData.Q0, ROI_ICA_min, ROI_ICA_max);

% If reference data AND either previous anode or cathode loss are provided,
% add a penalty. (You can adjust logic if you require both to be non-empty.)
% After penalty decision make decision which fititing method
% should be used via evaluating the weighting factors
if ~isempty(refData) && ...
        (~isempty(LAM_Anode_prev) || ~isempty(LAM_Cathode_prev))
    % Vectorized objective with penalty
    funMulti = @(X) objectiveWithPenalty( ...
        X, myData, myData.Q0, ROI_OCV_min, ROI_OCV_max, ...
        ROI_DVA_min, ROI_DVA_max, ROI_ICA_min, ROI_ICA_max, ...
        weightOCV, weightDVA, weightICA, refData, ...
        LAM_Anode_prev, LAM_Cathode_prev, LAM_Anode_Blend1_prev, ...
        LAM_Anode_Blend2_prev, capa_act, useBlendElectrodeModel, ...
        maxAnodeGain, maxCathodeGain, maxBlend1Gain, maxBlend2Gain, ...
        maxAnodeLoss, maxCathodeLoss, maxBlend1Loss, maxBlend2Loss, ...
        fitReverse); % 0 disallows negative cathode loss
else
    % No penalty: sum OCV + DVA + ICA with corresponding weights.
    % Must return Nx1, so skip mean(...).
    % Dynamic handling of the function handle
    funList = {};
    if weightDVA ~= 0
        funList{end+1} = @(X) weightDVA * funDVA(X);
    end
    if weightOCV ~= 0
        funList{end+1} = @(X) weightOCV * funOCV(X);
    end
    if weightICA ~= 0
        funList{end+1} = @(X) weightICA * funICA(X);
    end
    funMulti = @(X) sum(cellfun(@(f) f(X), funList));
end

% -----------------------------------------------------------------------
% 6) Base parameters
% -----------------------------------------------------------------------
if useBlendElectrodeModel
    % params = [alpha_an, beta_an, alpha_cat, beta_cat, gamma_Blend2]
    init_params = [1.05, -0.005, 1.1, -0.01, 0.2];
    lb          = [lowerBoundaries, 0.02];
    ub          = [upperBoundaries, gammaBlend2_upperBound];
else
    % params = [alpha_an, beta_an, alpha_cat, beta_cat]
    init_params = [1.2, 0.0, 1.1, -0.1];
    lb          =  lowerBoundaries;
    ub          =  upperBoundaries;
end

% -----------------------------------------------------------------------
% 7) Inhomogeneity parameters (optional)
% -----------------------------------------------------------------------
% Boolean mask (1 = active, 0 = inactive)
inhomoMask = [allowAnodeInhomogeneity, allowCathodeInhomogeneity];

% Start with the global per-electrode limits; tighten only if a previous CU exists.
maxInhomUB = inhomoUpperBoundBase;
if isempty(Inhmg_An_prev)
    maxInhomUB(1) = inhomoUpperBoundBase(1);    % first CU: keep full range
else
    maxInhomUB(1) = min(inhomoUpperBoundBase(1), ...
        inhomoDeltaPerCU(1) + Inhmg_An_prev);
end

if isempty(Inhmg_Ca_prev)
    maxInhomUB(2) = inhomoUpperBoundBase(2);    % first CU: keep full range
else
    maxInhomUB(2) = min(inhomoUpperBoundBase(2), ...
        inhomoDeltaPerCU(2) + Inhmg_Ca_prev);
end

% Initial values, lower/upper bounds
inhomoInit = 0.03       * inhomoMask;        % [0.02 0.02] or [0.02 0] ...
inhomoLB   = 0          * inhomoMask;        % [0 0] or [0 0]
inhomoUB   = maxInhomUB .* inhomoMask;       % [max max] or [max 0]
inhomoInit = min(inhomoInit, inhomoUB);      % clamp initial guess inside bounds

% -----------------------------------------------------------------------
% 8) Concatenate vectors 
% -----------------------------------------------------------------------
if any(inhomoMask)
    init_params = [init_params, inhomoInit];
    lb          = [lb,          inhomoLB];
    ub          = [ub,          inhomoUB];
end

% Run chosen optimization
switch Algorithm
    case 'patternsearch'
        params = patternsearch(funMulti, init_params, [], [], [], [], lb, ub);

    case 'particleswarm'
        options = optimoptions('particleswarm', ...
            'SwarmSize', 1000, ...
            'UseParallel', true);
        params = particleswarm(funMulti, numel(init_params), lb, ub, options);

    case 'ga'
        % GA in vectorized mode: calls funMulti(X) with X = [popSize x nParams]
        options = optimoptions('ga', ...
            'UseParallel', true, ...
            'PopulationSize', 500, ...
            'MaxGenerations', 100, ...
            'MaxStallGenerations', 50, ...
            'EliteCount', round(0.05 * 300), ...
            'CrossoverFraction', 0.8, ...
            'MutationFcn', {@mutationadaptfeasible, 0.2}, ...
            'SelectionFcn', @selectiontournament, ...
            'CrossoverFcn', @crossoverscattered, ...
            'UseVectorized', false, ...
            'Display', 'off');

        params = ga(funMulti, numel(init_params), [], [], [], [], ...
            lb, ub, [], options);

    case 'lsqnonlin'
        params = lsqnonlin(funMulti, init_params, lb, ub);

    case 'fmincon'
        opts = optimoptions(@fmincon, 'Algorithm','sqp', ...
            'MaxFunEvals',1e5, 'TolFun',1e-8);
        params = fmincon(funMulti, init_params, [], [], [], [], ...
            lb, ub, [], opts);

    case 'GlobalSearch'
        opts = optimoptions(@fmincon, 'Algorithm','sqp', ...
            'MaxFunEvals',1e5, 'TolFun',1e-10, 'TolCon', 1e-10);
        problem = createOptimProblem('fmincon', ...
            'objective', funMulti, 'x0', init_params, ...
            'lb', lb, 'ub', ub, 'options', opts);
        gs = GlobalSearch;
        params = run(gs, problem);

    otherwise
        error('Invalid solver chosen');
end

% -----------------------------------------------------------------------
% 9) Build the final reconstruction from the optimized parameters
% -----------------------------------------------------------------------
if useBlendElectrodeModel
    alpha_an        = params(1);
    beta_an         = params(2);
    alpha_cat       = params(3);
    beta_cat        = params(4);
    gamma_Blend2    = params(5);
    if allowAnodeInhomogeneity || allowCathodeInhomogeneity
        inhomo_val_an = params(6);
        inhomo_val_ca = params(7);
    else
        inhomo_val_an = 0;
        inhomo_val_ca = 0;
    end
    
else
    alpha_an  = params(1);
    beta_an   = params(2);
    alpha_cat = params(3);
    beta_cat  = params(4);
    gamma_Blend2  = 0;
end

% Compute anode curve (blend or pure Blend1)
[blendSOC, blendU] = blend_anode_curve(gamma_Blend2, myData, inhomo_val_an);
% Compute cathode curve with inhomogeneities
if allowCathodeInhomogeneity
    half_and_full_cell_data.normCathode_U = ...
        calculate_inhomogeneity( ...
            half_and_full_cell_data.normCathode_SOC, ...
            half_and_full_cell_data.normCathode_U, inhomo_val_ca);
end

% Shift/scale SOC for anode & cathode
anodeSOC        = alpha_an  * blendSOC      + beta_an;
cathSOC         = alpha_cat * ...
    half_and_full_cell_data.normCathode_SOC + beta_cat;
normCathode_U   = half_and_full_cell_data.normCathode_U;

% Evaluate on a uniform 0..1 axis
reconSOC     = linspace(0, 1, dataLength);
anodeU_recon = interp1(anodeSOC, blendU,         reconSOC, 'linear', 0);
cathU_recon  = interp1(cathSOC, ...
    half_and_full_cell_data.normCathode_U, reconSOC, 'linear', 0);
fcU_model    = cathU_recon - anodeU_recon;

% Compute measured vs. calculated DVA
[Q_DVA_meas, ~, DVA_meas]  = ...
    calculate_DVA(myData.Q_cell, myData.OCV_cell, dataLength + 1);
[Q_DVA_calc, ~, DVA_calc]  = ...
    calculate_DVA(reconSOC, fcU_model, dataLength + 1);
DVA_smooth_meas = smooth(DVA_meas,  smoothingPoints, 'lowess');
DVA_smooth_calc = smooth(DVA_calc,  smoothingPoints, 'lowess');

% Compute measured vs. calculated ICA
[Q_ICA_meas, ~, ICA_meas]  = ...
    calculate_ICA(myData.Q_cell, myData.OCV_cell, dataLength + 1);
[Q_ICA_calc, ~, ICA_calc]  = ...
    calculate_ICA(reconSOC, fcU_model, dataLength + 1);
ICA_smooth_meas = smooth(ICA_meas,  smoothingPoints, 'lowess');
ICA_smooth_calc = smooth(ICA_calc,  smoothingPoints, 'lowess');

% Return the optimization settings
AlgorithmOut      = Algorithm;
weightOCV_Out = weightOCV;

% -----------------------------------------------------------------------
% Nested function: objectiveWithPenalty (vectorized)
%   Must return [Nx1] if X is [Nx4 or Nx5].
%   Applies penalties when anode/cathode LAM exceed allowed limits.
% -----------------------------------------------------------------------
    function f = objectiveWithPenalty( ...
            X, myData, Q0, ROI_OCV_min, ROI_OCV_max, ROI_DVA_min, ...
            ROI_DVA_max, ROI_ICA_min, ROI_ICA_max, weightOCV, weightDVA, ...
            weightICA, refDataLoc, LAM_prevAnLocal, LAM_prevCathLocal, ...
            LAM_prevAnBlend1Local, LAM_prevAnBlend2Local, Capa_actLoc, ...
            useBlend, aAnodeLossLocal, aCathodeLossLocal, ...
            aAnodeBlend1LossLocal, aAnodeBlend2LossLocal, ...
            limitPositiveAnodeLossLocal, limitPositiveCathodeLossLocal, ...
            limitPositiveBlend1LossLocal, limitPositiveBlend2LossLocal, ...
            fitReverseLocal)

        baseVal = zeros(size(X,1),1);
        if weightOCV ~= 0
            baseVal = baseVal + weightOCV * ...
                fit_OCV_blend(X, myData, ROI_OCV_min, ROI_OCV_max);
        end
        if weightDVA ~= 0
            baseVal = baseVal + weightDVA * ...
                fit_DVA_mae_blend(X, myData, Q0, ROI_DVA_min, ROI_DVA_max);
        end
        if weightICA ~= 0
            baseVal = baseVal + weightICA * ...
                fit_ICA_mae_blend(X, myData, Q0, ROI_ICA_min, ROI_ICA_max);
        end

        Npop = size(X,1);
        penalty = zeros(Npop,1);
        scale = 1e8;

        for i = 1:Npop
            x_i = X(i,:);

            if useBlend
                paramsLoc = x_i;
            else
                paramsLoc = [x_i, 0];
            end
            
            [LAM_An, LAM_Cath, ~, LAM_An_Blend2, LAM_An_Blend1] = ...
                calculate_degradation_modes( ...
                paramsLoc, Capa_actLoc, refDataLoc.Capa_Anode_init, ...
                refDataLoc.Capa_Cathode_init, refDataLoc.Capa_Inventory_init, ...
                refDataLoc.gamma_Blend2_init, fitReverseLocal);
            

            LAM_currentAn = LAM_An;
            LAM_currentCath = LAM_Cath;
            LAM_currentAnBlend1 = LAM_An_Blend1;
            LAM_currentAnBlend2 = LAM_An_Blend2;

            tmpPenalty = 0;

            if ~isempty(LAM_prevAnLocal)
                neg = (LAM_prevAnLocal - aAnodeLossLocal) - LAM_currentAn;
                pos = LAM_currentAn - ...
                    (LAM_prevAnLocal + limitPositiveAnodeLossLocal);
                tmpPenalty = tmpPenalty + ...
                    scale * max(neg,0)^2 + scale * max(pos,0)^2;
            end

            if ~isempty(LAM_prevCathLocal)
                neg = (LAM_prevCathLocal - aCathodeLossLocal) - LAM_currentCath;
                pos = LAM_currentCath - ...
                (LAM_prevCathLocal + limitPositiveCathodeLossLocal);
                tmpPenalty = tmpPenalty + ...
                scale * max(neg,0)^2 + scale * max(pos,0)^2;
            end

            if ~isempty(LAM_prevAnBlend1Local)
                neg = (LAM_prevAnBlend1Local - aAnodeBlend1LossLocal) - ...
                    LAM_currentAnBlend1;
                pos = LAM_currentAnBlend1 - ...
                    (LAM_prevAnBlend1Local + limitPositiveBlend1LossLocal);
                tmpPenalty = tmpPenalty + ...
                    scale * max(neg,0)^2 + scale * max(pos,0)^2;
            end

            if ~isempty(LAM_prevAnBlend2Local)
                neg = (LAM_prevAnBlend2Local - aAnodeBlend2LossLocal) - ...
                    LAM_currentAnBlend2;
                pos = LAM_currentAnBlend2 - ...
                    (LAM_prevAnBlend2Local + limitPositiveBlend2LossLocal);
                tmpPenalty = tmpPenalty + ...
                    scale * max(neg,0)^2 + scale * max(pos,0)^2;
            end

            penalty(i) = tmpPenalty;
        end

        f = baseVal + penalty;
    end
end
