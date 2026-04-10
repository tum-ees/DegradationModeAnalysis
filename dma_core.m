function [solverInput, reconSOC, fullCellUModel, qDVAMeas, dvaSmoothMeas, ...
    qDVACalc, dvaSmoothCalc, qICAMeas, icaSmoothMeas, ...
    qICACalc, icaSmoothCalc, capaAct, params, algorithmOut, ...
    weightOCVOut, cathodeSOC, cathodeURecon, anodeSOC, anodeURecon] = ...
    dma_core(halfAndFullCellData, settings, refData, ...
        lamAnodePrev, lamCathodePrev, lamAnodeBlend1Prev, ...
        lamAnodeBlend2Prev, inhomAnPrev, inhomCaPrev, ...
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
%   * Core solver used by the DMA framework to fit electrode parameters per CU.
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
%   halfAndFullCellData
%     - struct containing the preprocessed half- and full-cell data.
%   settings
%     - settings struct defined in main_dma.
%   refData
%     - optional struct with reference capacities for penalties in
%       objectiveWithPenalty; must include .capaAnodeInit,
%       .capaCathodeInit, .capaInventoryInit, and .gammaAnBlend2Init.
%   lamAnodePrev
%     - optional previous anode loss for penalties; empty disables the
%       negative-anode penalty.
%   lamCathodePrev
%     - optional previous cathode loss for penalties; empty disables the
%       negative-cathode penalty.
%   allowAnodeInhomogeneity / allowCathodeInhomogeneity
%     - logical flags that permit inhomogeneity parameters.
%   fitReverse
%     - if true, process CUs in reverse order (last CU first).
%
% Outputs:
%   solverInput
%     - struct with required OCV/DVA data (normCathodeSOC, qCell, etc.).
%   reconSOC
%     - uniform SOC array (0..1) used for reconstruction.
%   fullCellUModel
%     - reconstructed full-cell OCV from the model.
%   qDVAMeas
%     - charge grid for measured DVA.
%   dvaSmoothMeas
%     - smoothed measured DVA array.
%   qDVACalc
%     - charge grid for calculated (modeled) DVA.
%   dvaSmoothCalc
%     - smoothed calculated DVA array.
%   capaAct
%     - actual capacity from the input data (e.g., fullCellData).
%   params
%     - optimised parameters [alphaAn, betaAn, alphaCat, betaCat,
%       (gammaAnBlend2), (gammaCaBlend2), (inhomAn), (inhomCat)].
%   algorithmOut
%     - echo of the chosen algorithm.
%   weightOCVOut
%     - echo of the chosen OCV weight.
%   cathodeSOC
%     - cathode SOC axis after blending and shifting.
%   cathodeURecon
%     - cathode potential curve paired with cathodeSOC.
%   anodeSOC
%     - anode SOC axis after blending and shifting.
%   anodeURecon
%     - anode potential curve paired with anodeSOC.
%

% ======= SETTINGS / CONSTANTS / PREPROCESSING =======
warnState = warning('off', 'all');
warnCleanup = onCleanup(@() warning(warnState));

% Get Settings from settings struct
smoothingPoints         = settings.smoothingPoints;
dataLength              = settings.dataLength;
algorithm               = settings.algorithm;
weightOCV               = settings.weightOCV;
weightDVA               = settings.weightDVA;
weightICA               = settings.weightICA;
roiOCVMin               = settings.roiOCVMin;
roiOCVMax               = settings.roiOCVMax;
roiDVAMin               = settings.roiDVAMin;
roiDVAMax               = settings.roiDVAMax;
roiICAMin               = settings.roiICAMin;
roiICAMax               = settings.roiICAMax;
useAnodeBlendModel      = settings.useAnodeBlend;
useCathodeBlendModel    = settings.useCathodeBlend;
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
gammaAnBlend2UpperBound = settings.gammaAnBlend2UpperBound;
capaAct                 = halfAndFullCellData.capaAct;

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

if isfield(settings, 'inhomAnodeOffset')
    inhomAnodeOffset = settings.inhomAnodeOffset;
else
    inhomAnodeOffset = 0;
end
if isfield(settings, 'inhomCathodeOffset')
    inhomCathodeOffset = settings.inhomCathodeOffset;
else
    inhomCathodeOffset = 0;
end
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

% If user did not provide lamCathodePrev, default to []
if ~exist('lamCathodePrev','var')
    lamCathodePrev = [];
end

% -----------------------------------------------------------------------
% 4) Pack everything into solverInput struct, to pass to objective functions
% -----------------------------------------------------------------------
solverInput.normCathodeSOC       = halfAndFullCellData.normCathodeSOC;
solverInput.normCathodeU         = halfAndFullCellData.normCathodeU;
solverInput.cathodeSOCSingle     = halfAndFullCellData.cathodeSOCSingle;
solverInput.cathodeUSingle       = halfAndFullCellData.cathodeUSingle;
solverInput.qCell                = halfAndFullCellData.fullCellSOC;
solverInput.ocvCell              = halfAndFullCellData.fullCellU;
solverInput.commonVoltageAnode   = halfAndFullCellData.commonVoltageAnode;
solverInput.qAnodeBlend2Interp   = halfAndFullCellData.qAnodeBlend2Interp;
solverInput.qAnodeBlend1Interp   = halfAndFullCellData.qAnodeBlend1Interp;
solverInput.qCathodeBlend2Interp = halfAndFullCellData.qCathodeBlend2Interp;
solverInput.qCathodeBlend1Interp = halfAndFullCellData.qCathodeBlend1Interp;
solverInput.commonVoltageCathode = halfAndFullCellData.commonVoltageCathode;
solverInput.anodeSOCSingle       = halfAndFullCellData.anodeSOCSingle;
solverInput.anodeUSingle         = halfAndFullCellData.anodeUSingle;
solverInput.useCathodeBlend      = useCathodeBlendModel;
solverInput.useAnodeBlend        = useAnodeBlendModel;
solverInput.inhomAnodeOffset     = inhomAnodeOffset;
solverInput.inhomCathodeOffset   = inhomCathodeOffset;
solverInput.q0                   = halfAndFullCellData.q0;

% Precompute static masks and measured derivatives to avoid per-iteration recomputation
dvaPrecomp.mask        = build_roi_mask(solverInput.qCell, roiDVAMin, roiDVAMax);
dvaPrecomp.measuredDVA = precompute_measured_dva(solverInput.qCell, solverInput.ocvCell, solverInput.q0);
icaPrecomp.mask        = build_roi_mask(solverInput.qCell, roiICAMin, roiICAMax);
icaPrecomp.measuredICA = precompute_measured_ica(solverInput.qCell, solverInput.ocvCell, solverInput.q0);

% -----------------------------------------------------------------------
% 5) Define vectorized objective functions and run optimization
% -----------------------------------------------------------------------
% Sub-objectives (each must accept fixed-length 8 params, return [Nx1])
funOCV = @(X) fit_ocv(X, solverInput, roiOCVMin, roiOCVMax);
funDVA = @(X) fit_dva(X, solverInput, solverInput.q0, roiDVAMin, roiDVAMax, dvaPrecomp);
funICA = @(X) fit_ica(X, solverInput, solverInput.q0, roiICAMin, roiICAMax, icaPrecomp);

% If reference data AND either previous anode or cathode loss are provided,
% add a penalty. (You can adjust logic if you require both to be non-empty.)
% After penalty decision make decision which fititing method
% should be used via evaluating the weighting factors
if ~isempty(refData) && ...
        (~isempty(lamAnodePrev) || ~isempty(lamCathodePrev))
    % Vectorized objective with penalty
    funMulti = @(X) objectiveWithPenalty( ...
        X, solverInput, solverInput.q0, roiOCVMin, roiOCVMax, ...
        roiDVAMin, roiDVAMax, roiICAMin, roiICAMax, ...
        weightOCV, weightDVA, weightICA, refData, ...
        lamAnodePrev, lamCathodePrev, lamAnodeBlend1Prev, ...
        lamAnodeBlend2Prev, capaAct, useAnodeBlendModel, useCathodeBlendModel, ...
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
if isempty(inhomAnPrev)
    maxInhomUB(1) = inhomUpperBoundBase(1);    % first CU: keep full range
else
    maxInhomUB(1) = min(inhomUpperBoundBase(1), ...
        inhomDeltaPerCU(1) + inhomAnPrev);
end

if isempty(inhomCaPrev)
    maxInhomUB(2) = inhomUpperBoundBase(2);    % first CU: keep full range
else
    maxInhomUB(2) = min(inhomUpperBoundBase(2), ...
        inhomDeltaPerCU(2) + inhomCaPrev);
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
% [alphaAn, betaAn, alphaCat, betaCat, gammaAnBlend2, gammaCaBlend2, ...
% inhomAn, inhomCa];

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
    fullUB(5)   = gammaAnBlend2UpperBound;
else
    fullInit(5) = 0;
    fullLB(5)   = 0;
    fullUB(5)   = 0;
end

% Cathode gamma (Blend2) if enabled
if useCathodeBlendModel
    fullInit(6) = 0.2;
    fullLB(6)   = 0.02;
    fullUB(6)   = settings.gammaCaBlend2UpperBound;
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

initParams = fullInit(freeIdx);
lb         = fullLB(freeIdx);
ub         = fullUB(freeIdx);

% Expander for free to full (fixed) layout
expandParamsFixed  = @(Xfree) local_expand_params_fixed(Xfree, freeIdx);

% Wrap objective so solvers only see free variables,
% while keeping fixed ordering for existing fit functions
funMulti = @(Xfree) funMultiFull(expandParamsFixed(Xfree));

% Run chosen optimization
switch algorithm
    case 'patternsearch'
        params = patternsearch(funMulti, initParams, [], [], [], [], lb, ub);

    case 'particleswarm'
        options = optimoptions('particleswarm', ...
            'SwarmSize', 1000, ...
            'UseParallel', true);
        params = particleswarm(funMulti, numel(initParams), lb, ub, options);

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

        params = ga(funMulti, numel(initParams), [], [], [], [], ...
            lb, ub, [], options);

    case 'lsqnonlin'
        params = lsqnonlin(funMulti, initParams, lb, ub);

    case 'fmincon'
        opts = optimoptions(@fmincon, 'Algorithm','sqp', ...
            'MaxFunEvals',1e5, 'TolFun',1e-8);
        params = fmincon(funMulti, initParams, [], [], [], [], ...
            lb, ub, [], opts);

    case 'GlobalSearch'
        opts = optimoptions(@fmincon, 'Algorithm','sqp', ...
            'MaxFunEvals',1e5, 'TolFun',1e-10, 'TolCon', 1e-10);
        problem = createOptimProblem('fmincon', ...
            'objective', funMulti, 'x0', initParams, ...
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
alphaAn       = params(1);
betaAn        = params(2);
alphaCat      = params(3);
betaCat       = params(4);
gammaAnBlend2 = params(5);
gammaCaBlend2 = params(6);
inhomValAn    = params(7);
inhomValCa    = params(8);

% Compute anode curve; split between blend and non-blend paths
if solverInput.useAnodeBlend && ~isempty(solverInput.qAnodeBlend1Interp)
    [blendSOC, blendU] = calculate_blend_curve(gammaAnBlend2, solverInput, 'anode');
else
    blendSOC = halfAndFullCellData.anodeSOCSingle;
    blendU   = halfAndFullCellData.anodeUSingle;
end

% Apply anode inhomogeneity to the source curve if enabled
if allowAnodeInhomogeneity
    blendU = calculate_inhomogeneity(blendSOC, blendU, inhomValAn, inhomAnodeOffset);
end
anodeURecon = blendU;

% Compute cathode curve (blend optional) with inhomogeneities
if solverInput.useCathodeBlend && ~isempty(solverInput.qCathodeBlend1Interp)
    [cathodeSOCSrc, cathUSrc] = calculate_blend_curve(gammaCaBlend2, solverInput, 'cathode');
else
    cathodeSOCSrc = halfAndFullCellData.normCathodeSOC;
    cathUSrc      = halfAndFullCellData.normCathodeU;
end
if allowCathodeInhomogeneity
    cathUSrc = calculate_inhomogeneity( ...
        cathodeSOCSrc, cathUSrc, inhomValCa, inhomCathodeOffset);
end

% Shift/scale SOC for anode & cathode
anodeSOC      = alphaAn * blendSOC + betaAn;
cathodeSOC    = alphaCat * cathodeSOCSrc + betaCat;
cathodeURecon = cathUSrc;

% Evaluate on a uniform 0..1 axis
reconSOC = linspace(0, 1, dataLength);
anodeUModel = interp1(anodeSOC, anodeURecon, reconSOC, 'linear', 0);
cathodeUModel = interp1(cathodeSOC, cathodeURecon, reconSOC, 'linear', 0);
fullCellUModel = cathodeUModel - anodeUModel;

% Compute measured vs. calculated DVA
[qDVAMeas, ~, dvaMeas]  = ...
    calculate_dva(solverInput.qCell, solverInput.ocvCell, dataLength + 1);
[qDVACalc, ~, dvaCalc]  = ...
    calculate_dva(reconSOC, fullCellUModel, dataLength + 1);
dvaSmoothMeas = smooth(dvaMeas, smoothingPoints, 'lowess');
dvaSmoothCalc = smooth(dvaCalc, smoothingPoints, 'lowess');

% Compute measured vs. calculated ICA
[qICAMeas, ~, icaMeas]  = ...
    calculate_ica(solverInput.qCell, solverInput.ocvCell, dataLength + 1);
[qICACalc, ~, icaCalc]  = ...
    calculate_ica(reconSOC, fullCellUModel, dataLength + 1);
icaSmoothMeas = smooth(icaMeas, smoothingPoints, 'lowess');
icaSmoothCalc = smooth(icaCalc, smoothingPoints, 'lowess');

% Return the optimization settings
algorithmOut      = algorithm;
weightOCVOut = weightOCV;

% -----------------------------------------------------------------------
% Nested function: objectiveWithPenalty (vectorized)
%   Must return [Nx1] if X is [Nx4 or Nx5].
%   Applies penalties when anode/cathode LAM exceed allowed limits.
% -----------------------------------------------------------------------
    function f = objectiveWithPenalty( ...
            X, solverInput, q0, roiOCVMinLocal, roiOCVMaxLocal, roiDVAMinLocal, ...
            roiDVAMaxLocal, roiICAMinLocal, roiICAMaxLocal, weightOCV, weightDVA, ...
            weightICA, refDataLoc, lamPrevAnLocal, lamPrevCathLocal, ...
            lamPrevAnBlend1Local, lamPrevAnBlend2Local, capaActLoc, ...
            useAnodeBlendLocal, useCathodeBlendLocal, aAnodeLossLocal, aCathodeLossLocal, ...
            aAnodeBlend1LossLocal, aAnodeBlend2LossLocal, ...
            limitPositiveAnodeLossLocal, limitPositiveCathodeLossLocal, ...
            limitPositiveBlend1LossLocal, limitPositiveBlend2LossLocal, ...
            fitReverseLocal)

        baseVal = zeros(size(X,1),1);
        if weightOCV ~= 0
            baseVal = baseVal + weightOCV * ...
                fit_ocv(X, solverInput, roiOCVMinLocal, roiOCVMaxLocal);
        end
        if weightDVA ~= 0
            baseVal = baseVal + weightDVA * ...
                fit_dva(X, solverInput, q0, roiDVAMinLocal, roiDVAMaxLocal, dvaPrecomp);
        end
        if weightICA ~= 0
            baseVal = baseVal + weightICA * ...
                fit_ica(X, solverInput, q0, roiICAMinLocal, roiICAMaxLocal, icaPrecomp);
        end

        Npop = size(X,1);
        penalty = zeros(Npop,1);
        scale = 1e8; % large penalty weight so constraint violations dominate the objective

        for i = 1:Npop
            xCandidate = X(i,:);

            paramsLoc = xCandidate;
            if ~useAnodeBlendLocal
                paramsLoc(5) = 0; % enforce zero anode blend fraction
            end
            if ~useCathodeBlendLocal
                paramsLoc(6) = 0; % enforce zero cathode blend fraction
            end
            
            
            [lamAn, lamCath, ~, lamAnBlend2, lamAnBlend1, ~, ~] = ...
                calculate_degradation_modes( ...
                paramsLoc, capaActLoc, refDataLoc.capaAnodeInit, ...
                refDataLoc.capaCathodeInit, refDataLoc.capaInventoryInit, ...
                refDataLoc.gammaAnBlend2Init, refDataLoc.gammaCaBlend2Init, fitReverseLocal);
            

            lamCurrentAn = lamAn;
            lamCurrentCath = lamCath;
            lamCurrentAnBlend1 = lamAnBlend1;
            lamCurrentAnBlend2 = lamAnBlend2;

            tmpPenalty = 0;

            if ~isempty(lamPrevAnLocal)
                neg = (lamPrevAnLocal - aAnodeLossLocal) - lamCurrentAn;
                pos = lamCurrentAn - ...
                    (lamPrevAnLocal + limitPositiveAnodeLossLocal);
                tmpPenalty = tmpPenalty + ...
                    scale * max(neg,0)^2 + scale * max(pos,0)^2;
            end

            if ~isempty(lamPrevCathLocal)
                neg = (lamPrevCathLocal - aCathodeLossLocal) - lamCurrentCath;
                pos = lamCurrentCath - ...
                (lamPrevCathLocal + limitPositiveCathodeLossLocal);
                tmpPenalty = tmpPenalty + ...
                scale * max(neg,0)^2 + scale * max(pos,0)^2;
            end

            if ~isempty(lamPrevAnBlend1Local)
                neg = (lamPrevAnBlend1Local - aAnodeBlend1LossLocal) - ...
                    lamCurrentAnBlend1;
                pos = lamCurrentAnBlend1 - ...
                    (lamPrevAnBlend1Local + limitPositiveBlend1LossLocal);
                tmpPenalty = tmpPenalty + ...
                    scale * max(neg,0)^2 + scale * max(pos,0)^2;
            end

            if ~isempty(lamPrevAnBlend2Local)
                neg = (lamPrevAnBlend2Local - aAnodeBlend2LossLocal) - ...
                    lamCurrentAnBlend2;
                pos = lamCurrentAnBlend2 - ...
                    (lamPrevAnBlend2Local + limitPositiveBlend2LossLocal);
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

