function [myData, reconSOC, fcU_model, Q_DVA_meas, DVA_smooth_meas, ...
    Q_DVA_calc, DVA_smooth_calc, Q_ICA_meas, ICA_smooth_meas, ...
    Q_ICA_calc, ICA_smooth_calc, capa_act, params, algorithmOut, ...
    weightOCV_Out, cathSOC, normCathode_U, anodeSOC, blendU] = ...
    dma_core(half_and_full_cell_data, settings, refData, ...
        LAM_anode_prev, LAM_cathode_prev, LAM_anode_blend1_prev, ...
        LAM_anode_blend2_prev, inhom_An_prev, inhom_Ca_prev, ...
        allowAnodeInhomogeneity, allowCathodeInhomogeneity, fitReverse)
%> -------------------------------------------------------------------------
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Additional code by Sebastian Karl (s.karl@tum.de)
%> Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%> Additional code by Moritz Guenthner (moritz.guenthner@tum.de)
%> Additional code by Mathias Rehm (mathias.rehm@tum.de)
%
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
%       objectiveWithPenalty; must include .capa_anode_init,
%       .capa_cathode_init, .capa_inventory_init, and .gamma_an_blend2_init.
%   LAM_anode_prev
%     - optional previous anode loss for penalties; empty disables the
%       negative-anode penalty.
%   LAM_cathode_prev
%     - optional previous cathode loss for penalties; empty disables the
%       negative-cathode penalty.
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
%   capa_act
%     - actual capacity from the input data (e.g., fullCellData).
%   params
%     - optimised parameters [alpha_an, beta_an, alpha_cat, beta_cat,
%       (gamma_an_blend2), (gamma_ca_blend2 future), (inhom_an), (inhom_cat)].
%   algorithmOut
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
algorithm               = settings.algorithm;
weightOCV               = settings.weightOCV;
weightDVA               = settings.weightDVA;
weightICA               = settings.weightICA;
ROI_OCV_min             = settings.ROI_OCV_min;
ROI_OCV_max             = settings.ROI_OCV_max;
ROI_DVA_min             = settings.ROI_DVA_min;
ROI_DVA_max             = settings.ROI_DVA_max;
ROI_ICA_min             = settings.ROI_ICA_min;
ROI_ICA_max             = settings.ROI_ICA_max;
useAnodeBlendModel    = settings.useAnodeBlend;
useCathodeBlendModel  = settings.useCathodeBlend;
maxAnodeGain            = settings.maxAnodeGain;
maxCathodeGain          = settings.maxCathodeGain;
maxAnBlend1Gain         = settings.maxAnBlend1Gain;
maxAnBlend2Gain         = settings.maxAnBlend2Gain;
maxAnodeLoss            = settings.maxAnodeLoss;
maxCathodeLoss          = settings.maxCathodeLoss;
maxAnBlend1Loss         = settings.maxAnBlend1Loss;
maxAnBlend2Loss         = settings.maxAnBlend2Loss;
lowerBoundaries         = settings.lowerBoundaries;
upperBoundaries         = settings.upperBoundaries;
gammaAnBlend2_upperBound  = settings.gammaAnBlend2_upperBound;
capa_act                = half_and_full_cell_data.capa_act;

% Build per-electrode inhomogeneity limits upfront so the first CU can vary.
inhomUpperBoundBase = settings.maxInhomogeneity;
if isscalar(inhomUpperBoundBase)
    inhomUpperBoundBase = repmat(inhomUpperBoundBase, 1, 2);
else
    inhomUpperBoundBase = inhomUpperBoundBase(:).';
    if isscalar(inhomUpperBoundBase)
        inhomUpperBoundBase = repmat(inhomUpperBoundBase, 1, 2);
    else
        inhomUpperBoundBase = inhomUpperBoundBase(1:2);
    end
end

inhomDeltaPerCU = settings.maxInhomogeneityDelta;
if isscalar(inhomDeltaPerCU)
    inhomDeltaPerCU = repmat(inhomDeltaPerCU, 1, 2);
else
    inhomDeltaPerCU = inhomDeltaPerCU(:).';
    if isscalar(inhomDeltaPerCU)
        inhomDeltaPerCU = repmat(inhomDeltaPerCU, 1, 2);
    else
        inhomDeltaPerCU = inhomDeltaPerCU(1:2);
    end
end

% If user did not provide maxCathodeGain, default to 0
if ~exist('maxCathodeGain','var') || isempty(maxCathodeGain)
    maxCathodeGain = 0;
end

% If user did not provide LAM_cathode_prev, default to []
if ~exist('LAM_cathode_prev','var')
    LAM_cathode_prev = [];
end

% -----------------------------------------------------------------------
% 4) Pack everything into myData struct, to pass to objective functions
% -----------------------------------------------------------------------
myData.normCathode_SOC = half_and_full_cell_data.normCathode_SOC;
myData.normCathode_U   = half_and_full_cell_data.normCathode_U;
myData.cathode_SOC_single = half_and_full_cell_data.cathode_SOC_single;
myData.cathode_U_single   = half_and_full_cell_data.cathode_U_single;
myData.Q_cell          = half_and_full_cell_data.fullCell_SOC;
myData.OCV_cell        = half_and_full_cell_data.fullCell_U;
myData.commonVoltage_anode   = half_and_full_cell_data.commonVoltage_anode;
myData.Q_anode_blend2_interp = half_and_full_cell_data.Q_anode_blend2_interp;
myData.Q_anode_blend1_interp = half_and_full_cell_data.Q_anode_blend1_interp;
myData.Q_cathode_blend2_interp = half_and_full_cell_data.Q_cathode_blend2_interp;
myData.Q_cathode_blend1_interp = half_and_full_cell_data.Q_cathode_blend1_interp;
myData.commonVoltage_cathode   = half_and_full_cell_data.commonVoltage_cathode;
myData.anode_SOC_single        = half_and_full_cell_data.anode_SOC_single;
myData.anode_U_single          = half_and_full_cell_data.anode_U_single;
myData.useCathodeBlend         = useCathodeBlendModel;
myData.useAnodeBlend           = useAnodeBlendModel;
myData.Q0              = half_and_full_cell_data.Q0;

% Precompute static masks and measured derivatives to avoid per-iteration recomputation
dvaPrecomp.mask        = build_ROI_mask(myData.Q_cell, ROI_DVA_min, ROI_DVA_max);
dvaPrecomp.measuredDVA = precompute_measured_DVA(myData.Q_cell, myData.OCV_cell, myData.Q0);
icaPrecomp.mask        = build_ROI_mask(myData.Q_cell, ROI_ICA_min, ROI_ICA_max);
icaPrecomp.measuredICA = precompute_measured_ICA(myData.Q_cell, myData.OCV_cell, myData.Q0);

% -----------------------------------------------------------------------
% 5) Define vectorized objective functions and run optimization
% -----------------------------------------------------------------------
% Sub-objectives (each must accept fixed-length 8 params, return [N×1])
funOCV = @(X) fit_OCV(X, myData, ROI_OCV_min, ROI_OCV_max);
funDVA = @(X) fit_DVA(X, myData, myData.Q0, ROI_DVA_min, ROI_DVA_max, dvaPrecomp);
funICA = @(X) fit_ICA(X, myData, myData.Q0, ROI_ICA_min, ROI_ICA_max, icaPrecomp);

% If reference data AND either previous anode or cathode loss are provided,
% add a penalty. (You can adjust logic if you require both to be non-empty.)
% After penalty decision make decision which fititing method
% should be used via evaluating the weighting factors
if ~isempty(refData) && ...
        (~isempty(LAM_anode_prev) || ~isempty(LAM_cathode_prev))
    % Vectorized objective with penalty
    funMulti = @(X) objectiveWithPenalty( ...
        X, myData, myData.Q0, ROI_OCV_min, ROI_OCV_max, ...
        ROI_DVA_min, ROI_DVA_max, ROI_ICA_min, ROI_ICA_max, ...
        weightOCV, weightDVA, weightICA, refData, ...
        LAM_anode_prev, LAM_cathode_prev, LAM_anode_blend1_prev, ...
        LAM_anode_blend2_prev, capa_act, useAnodeBlendModel, useCathodeBlendModel, ...
                maxAnodeGain, maxCathodeGain, maxAnBlend1Gain, maxAnBlend2Gain, ...
                maxAnodeLoss, maxCathodeLoss, maxAnBlend1Loss, maxAnBlend2Loss, ...
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

funMultiFull = funMulti;

% -----------------------------------------------------------------------
% 6) Inhomogeneity parameters (optional)
% -----------------------------------------------------------------------
% Boolean mask (1 = active, 0 = inactive)
inhomMask = [allowAnodeInhomogeneity, allowCathodeInhomogeneity];

% Start with the global per-electrode limits; tighten only if a previous CU exists.
maxInhomUB = inhomUpperBoundBase;
if isempty(inhom_An_prev)
    maxInhomUB(1) = inhomUpperBoundBase(1);    % first CU: keep full range
else
    maxInhomUB(1) = min(inhomUpperBoundBase(1), ...
        inhomDeltaPerCU(1) + inhom_An_prev);
end

if isempty(inhom_Ca_prev)
    maxInhomUB(2) = inhomUpperBoundBase(2);    % first CU: keep full range
else
    maxInhomUB(2) = min(inhomUpperBoundBase(2), ...
        inhomDeltaPerCU(2) + inhom_Ca_prev);
end

% Initial values, lower/upper bounds
inhomInit = 0.03       * inhomMask;        % [0.02 0.02] or [0.02 0] ...
inhomLB   = 0          * inhomMask;        % [0 0] or [0 0]
inhomUB   = maxInhomUB .* inhomMask;       % [max max] or [max 0]
inhomInit = min(inhomInit, inhomUB);      % clamp initial guess inside bounds

% -----------------------------------------------------------------------
% 7) Build full 8 parameter vectors, then reduce to active subset
% -----------------------------------------------------------------------
% Full fixed order:
% [alpha_an, beta_an, alpha_cat, beta_cat, gamma_an_blend2, gamma_ca_blend2, inhom_an, inhom_ca]
% gamma_ca_blend2 (slot 6) is reserved for a future cathode blend release.

fullInit = zeros(1, 8);
fullLB   = zeros(1, 8);
fullUB   = zeros(1, 8);

% Base 4 parameters always exist
if useAnodeBlendModel
    baseInit = [1.05, -0.005, 1.1, -0.01];
else
    baseInit = [1.2, 0.0, 1.1, -0.1];
end

fullInit(1:4) = baseInit;
fullLB(1:4)   = lowerBoundaries;
fullUB(1:4)   = upperBoundaries;

% Anode gamma slot is only active if blend model is used
if useAnodeBlendModel
    fullInit(5) = 0.2;
    fullLB(5)   = 0.02;
    fullUB(5)   = gammaAnBlend2_upperBound;
else
    fullInit(5) = 0;
    fullLB(5)   = 0;
    fullUB(5)   = 0;
end

% Cathode gamma (Blend2) if enabled
if useCathodeBlendModel
    fullInit(6) = 0.2;
    fullLB(6)   = 0.02;
    fullUB(6)   = settings.gammaCaBlend2_upperBound;
else
    fullInit(6) = 0;
    fullLB(6)   = 0;
    fullUB(6)   = 0;
end

% Inhomogeneity slots (already masked)
fullInit(7:8) = inhomInit;
fullLB(7:8)   = inhomLB;
fullUB(7:8)   = inhomUB;

% Active mask and reduced vectors for the solver
activeMask = [true true true true useAnodeBlendModel useCathodeBlendModel ...
    allowAnodeInhomogeneity allowCathodeInhomogeneity];
freeIdx = find(activeMask);

init_params = fullInit(freeIdx);
lb          = fullLB(freeIdx);
ub          = fullUB(freeIdx);

% Expander for free to full (fixed) layout
expandParamsFixed  = @(Xfree) local_expand_params_fixed(Xfree, freeIdx);

% Wrap objective so solvers only see free variables,
% while keeping fixed ordering for existing fit functions
funMulti = @(Xfree) funMultiFull(expandParamsFixed(Xfree));

% Run chosen optimization
switch algorithm
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

% Expand to fixed full length 8 for output and reconstruction
params = expandParamsFixed(params);

% -----------------------------------------------------------------------
% 9) Build the final reconstruction from the optimized parameters
% -----------------------------------------------------------------------
alpha_an          = params(1);
beta_an           = params(2);
alpha_cat         = params(3);
beta_cat          = params(4);
gamma_an_blend2   = params(5);
gamma_ca_blend2   = params(6);
inhom_val_an      = params(7);
inhom_val_ca      = params(8);

% Compute anode curve; split between blend and non-blend paths
if myData.useAnodeBlend && ~isempty(myData.Q_anode_blend1_interp)
    [blendSOC, blendU] = calculate_blend_curve(gamma_an_blend2, myData, 'anode');
else
    blendSOC = half_and_full_cell_data.anode_SOC_single;
    blendU   = half_and_full_cell_data.anode_U_single;
end

% Apply anode inhomogeneity to the source curve if enabled
if allowAnodeInhomogeneity
    blendU = calculate_inhomogeneity(blendSOC, blendU, inhom_val_an);
end

% Compute cathode curve (blend optional) with inhomogeneities
if myData.useCathodeBlend && ~isempty(myData.Q_cathode_blend1_interp)
    [cathSOC_src, cathU_src] = calculate_blend_curve(gamma_ca_blend2, myData, 'cathode');
else
    cathSOC_src = half_and_full_cell_data.normCathode_SOC;
    cathU_src   = half_and_full_cell_data.normCathode_U;
end
if allowCathodeInhomogeneity
    cathU_src = calculate_inhomogeneity( ...
        cathSOC_src, cathU_src, inhom_val_ca);
end

% Shift/scale SOC for anode & cathode
anodeSOC        = alpha_an  * blendSOC      + beta_an;
cathSOC         = alpha_cat * cathSOC_src + beta_cat;
normCathode_U   = cathU_src;

% Evaluate on a uniform 0..1 axis
reconSOC     = linspace(0, 1, dataLength);
anodeU_recon = interp1(anodeSOC, blendU,         reconSOC, 'linear', 0);
cathU_recon  = interp1(cathSOC, cathU_src, reconSOC, 'linear', 0);
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
algorithmOut      = algorithm;
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
            LAM_prevAnBlend1Local, LAM_prevAnBlend2Local, capa_actLoc, ...
            useAnodeBlendLocal, useCathodeBlendLocal, aAnodeLossLocal, aCathodeLossLocal, ...
            aAnodeBlend1LossLocal, aAnodeBlend2LossLocal, ...
            limitPositiveAnodeLossLocal, limitPositiveCathodeLossLocal, ...
            limitPositiveBlend1LossLocal, limitPositiveBlend2LossLocal, ...
            fitReverseLocal)

        baseVal = zeros(size(X,1),1);
        if weightOCV ~= 0
            baseVal = baseVal + weightOCV * ...
                fit_OCV(X, myData, ROI_OCV_min, ROI_OCV_max);
        end
        if weightDVA ~= 0
            baseVal = baseVal + weightDVA * ...
                fit_DVA(X, myData, Q0, ROI_DVA_min, ROI_DVA_max, dvaPrecomp);
        end
        if weightICA ~= 0
            baseVal = baseVal + weightICA * ...
                fit_ICA(X, myData, Q0, ROI_ICA_min, ROI_ICA_max, icaPrecomp);
        end

        Npop = size(X,1);
        penalty = zeros(Npop,1);
        scale = 1e8;

        for i = 1:Npop
            x_i = X(i,:);

            paramsLoc = x_i;
            if ~useAnodeBlendLocal
                paramsLoc(5) = 0; % enforce zero anode blend fraction
            end
            if ~useCathodeBlendLocal
                paramsLoc(6) = 0; % enforce zero cathode blend fraction
            end
            
            
            [LAM_An, LAM_Cath, ~, LAM_An_blend2, LAM_An_blend1, ~, ~] = ...
                calculate_degradation_modes( ...
                paramsLoc, capa_actLoc, refDataLoc.capa_anode_init, ...
                refDataLoc.capa_cathode_init, refDataLoc.capa_inventory_init, ...
                refDataLoc.gamma_an_blend2_init, refDataLoc.gamma_ca_blend2_init, fitReverseLocal);
            

            LAM_currentAn = LAM_An;
            LAM_currentCath = LAM_Cath;
            LAM_currentAnBlend1 = LAM_An_blend1;
            LAM_currentAnBlend2 = LAM_An_blend2;

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

% -----------------------------------------------------------------------
% Nested function: local_expand_params_fixed
%   Expands free parameter vectors to full length 8 fixed order
% -----------------------------------------------------------------------
    function Xfull = local_expand_params_fixed(Xfree, freeIdxLoc)

        if isvector(Xfree)
            Xfree = Xfree(:).';
        end

        Xfull = zeros(size(Xfree, 1), 8);
        Xfull(:, freeIdxLoc) = Xfree;
    end

end
