function [Data, s] = main_DMA(userSettingsOutside)
%> -------------------------------------------------------------------------
%>  * Original code by:
%>    - Can Korkmaz (can.korkmaz@tum.de), supervised by Mathias Rehm
%>    - Additional code by Mathias Rehm (mathias.rehm@tum.de)
%>    - Additional code by Josef Eizenhammer (josef.eizenhammer@tum.de)
%>    - Additional code by Moritz Guenthner (moritz.guenthner@tum.de)
%>    - Additional code by Sebastian Karl (s.karl@tum.de)

%> Date: 2025-11-12
%>
%> OVERVIEW:
%>   * Entry point of the DMA framework; orchestrates analysis across CUs.
%>   * Performs a Degradation Mode Analysis (DMA) across checkup points (CU). 
%>   * Loads/constructs half-cell (cathode) data + full-cell data,
%>     then calls `dma_core` for each CU to fit parameters.
%>   * Applies constraints (RMSE threshold, maximum tries, etc.)
%>     to find acceptable solutions and keep the best (lowest RMSE) per CU.
%>   * Includes logic for partial or missing CUs and region-based RMSE.
%>     Adds reference-data checks for LAM constraints.
%>   * Plots OCV/DVA reconstructions, relevant degradation modes, etc.
%>
%> USAGE:
%>   - All user-defined settings are stored in the struct "s" (below).
%>      -> just follow the steps from 1) to 16) below (and give input in
%       a data format fitting the requirements)
%    - You can overwrite all settings WITHIN main_DMA by giving them as
%       input of settings (e.g. myUserSetting.tableFilter = {'Battery_serial',
%       '26'} and main_DMA(myUserSetting) will overwrite only the 
%       battery_serial hard coded within main_DMA.
%>   - The function returns a struct `Data` containing results for each CU:
%>       Data.CU1.params, Data.CU1.RMSE, Data.CU1.LAM_anode, ...
%>   - At the end, optional overall plots (DMA, Si content) can be created
%>       if toggles in `s` are set to true.
%>   - pOCVs (fullCellData) must be stored in one table (.mat) or in folders
%>       per CU. Use directories named CU1, CU2, ... in the basePath.
%>   - Half-cell data (e.g. cathodeData) has to be stored in a single file 
%>       as a table; for half-cell data, exact naming and correct SOC/UOCV 
%>       columns are required (in contrast to full-cell data no support in 
%>       data handling; see folder .\inputData to copy the structure and SOC 
%>       convention).
%>
%> DEPENDENCIES:
%>   - `dma_core`, `calculate_RMSE`, `calculate_degradation_modes`,
%>     `plot_OCV_model_param_show`,
%>     `plot_DMA`, `plot_gamma`.
%>   - Make sure you have all relevant half-cell data files accessible.
%>   - Ensure correct path additions for your local environment.
%>
%> -------------------------------------------------------------------------


%% -------------- SETUP ALL USER SETTINGS IN STRUCT s --------------
% Settings are overwritten if varargin contains a settings struct (see below)
s = struct();

% 1) Define path to your Aging Study (in case you use .mat table to store
% your pOCVs: path to your agingDataTable; in case of singular files: CU
% folder struct);
s.pathAgingStudy = ".\InputData" + ...
    "\TestData\P45B_serial23_aging_data_table.mat";

s.pathSaveResults = ".\InputData" + ...
    "\TestData\Results";

% 2) Select calendarOrCyclic
% s.CalendarOrCyclic selects the label semantics used in figures/titles only:
%   - 0 = Calendar study -> labels show "Check-up number".
%   - 1 = Cyclic study  -> labels show "EFC" (effective full cycles).
s.CalendarOrCyclic = 1;

% 3) Define number of CUs / vector for EFCs 
% s.nCUs is a vector whose length defines the number of check-ups processed.
%   - The i-th element of s.nCUs is a label/value attached to CU i:
%       * Calendar (0): label is the check-up number shown in plots.
%       * Cyclic   (1): label is the EFC shown in plots.
%   - Important mapping:
%       * A table containing all pOCVs over aging is needed 
%         (see the minimum working example as reference).
%       * The i-th CU row is selected from the table by CU index i 
%         (via a CU/CheckUp/RPT column if present, else by row position). 
%         The numeric contents of s.nCUs(i) are not used for lookup;
%         they are used later only for labelling/plotting.
%
% Fitting order (forward vs. reverse) is inferred from the order of s.nCUs:
%   - If s.nCUs is ascending (e.g., [0 100 200 300]),
%     fitting proceeds CU1->CU2->... .
%   - If s.nCUs is descending (e.g., [1400 1200 1000 ... 0]), fitting proceeds
%     in reverse (last->first). This helps when early CUs are noisy and you
%     prefer to seed from the most aged state.

if s.CalendarOrCyclic == 1
    % specify the logic of the CUs here (might be e.g. 0:1:10 or 0:100:800 for
    % calendar or cyclic aging; use e.g. 800:-100:0 to use reverse fitting)
    s.nCUs = 0:100:800;
else
    s.nCUs = 1:1:8; % e.g. 1:1:8;
end

% 4) Choose whether you use pOCVs in charge or discharge direction
% set which direction for the DMA; can be 'charge' or 'discharge'
s.direction                 = 'charge';
if exist('userSettingsOutside','var') && ...
        isfield(userSettingsOutside,'direction')
    s.direction = userSettingsOutside.direction;
end

% 5) Define how the algorithm can find your cell
% tell main_DMA how to find all lines in the table of this battery with the 
% right settings (e.g. 'Battery_serial', '102'; and TBegin = 25)
s.tableFilter(1,:) = {'Battery_serial'; '23'};
% optional: several filters
% tableFilter(2,:) = {'TBegin', '25'};
% optional: in case you have several tables in one table -> give the name
% of the columnname you need (e.g. Testdata_pOCV_DCH; default: empty ->
% there is only one column having structure or tabular format)
if strcmp(s.direction, 'charge')
    s.nameTableColumnOCV = 'Testdata_pOCV_CH';
else
    s.nameTableColumnOCV = 'Testdata_pOCV_DCH';
end

% 6) Specify how many (successfull) optimizations should pe performed per
% cell and per CU -> reqAccepted is the number of optimizations with an
% RMSE between measured and reconstructed curve below the threshold
s.reqAccepted           = 3;     % required number of solutions per CU
s.maxTriesOverall       = 10;     % maximum overall tries per CU
s.rmseThreshold         = 0.01; % acceptable RMSE threshold, in Volt
s.measureRuntimePerCU   = 1;     % toggle to measure runtime of a single CU run
% 7) Decide, which algorithm you want to use to find the parameters for the
% reconstruction of the OCV. Possible algorithms are:
% -> 'ga' (Genetic Algorithm, required: Global Optimization Toolbox)
%     -> ga is recommended: High computation costs, but best results
% -> 'patternsearch' (required: Global Optimization Toolbox)
% -> 'particleswarm' (required: Global Optimization Toolbox)
% -> 'lsqnonlin' (Nonlinear Least Squares, required: Optimization Toolbox)
% -> 'fmincon' (required: Optimization Toolbox)
% -> 'GlobalSearch' (required: Global Optimization Toolbox)
s.algorithm = 'ga';

% 8) Set the data length for the input data (data gets resampled to 
% e.g. 1000 dataPoints)
% Shorter length results in shorter runtime
% Default case: s.dataLength = 1000
s.dataLength  = 1000;

% 9) Set the number of smoothing points, which are used as an input for the 
% smooth function.
% This function is applied to the input data before conducting the DMA
% Default case: s.smoothingPoints = 30; LOWESS filter is used here
s.smoothingPoints = 30;

% 10) Decide, which type of fitting (OCV, DVA, or ICA) you want to do and
% how they should be weighted against each other. In most cases OCV and DVA
% together lead to best results.
% Weighting options: Default case (DVA fitting with OCV fitting
% -> s.weightOCV = 100, s.weightDVA = 1, s.weightICA = 0; 0 means not used)
% You use LFP? -> we recommend weightOCV / weightDVA = 10 / 3;
s.weightOCV        = 100;   % weighting for the OCV term
s.weightDVA        = 1;     % weighting for the DVA region
s.weightICA        = 0;     % weighting for the ICA region

% 11) Decide, in which SOC regions the fitting of DVA, ICA, and OCV should
% be done. 
% -> Good, robust results for fitting between 0 and 100% SOC for OCV 
% and middle region (e.g. between 5 and 95%) for DVA.
% You use LFP? -> we recommend to include only the high and low SOC region
% for the OCV fitting (e.g. 0-15% and 85-100% SOC) and middle region (e.g.
% 10-90% SOC) for DVA.
% Region-of-interest (ROI) settings:
% lower boundary of DVA fitting region (e.g. 0.1 = 10% SOC)
s.ROI_DVA_min = 0.1;
% upper boundary of DVA fitting region (e.g. 0.9 = 90% SOC)
s.ROI_DVA_max = 0.9; 
% lower boundary of ICA fitting region (e.g. 0.10 = 10% SOC)
s.ROI_ICA_min = 0.13;
% upper boundary of ICA fitting region (e.g. 0.9 = 90% SOC)
s.ROI_ICA_max = 0.9; 
%  lower and upper boundary of OCV fitting region (e.g. 0.02 = 2% SOC);
% can be also given as regions (e.g. s.ROI_OCV_min = [0.02; 0.25] and 
% s.ROI_OCV_max = [0.80; 0.98]) -> in this case only the regions 2-25% SOC 
% and 80-98% SOC are used for the OCV fitting.
s.ROI_OCV_min = 0;
s.ROI_OCV_max = 1;

% 12) Change the boundaries for the parameters if you are very familiar with 
% the algorithm. If not, keep them! 
% they are saved as: params = [alpha_an, beta_an, alpha_cat, beta_cat]
% Use these to set lower and upper boundaries for the fitting process
% recommendation: to start with the fitting, just keep those very wide
% boundaries -> this assures the physically correct solution is within the
% allowed boundaries
% Default case:
%   s.lowerBoundaries          = [0.8, -1.0, 0.8, -1.0];
%   s.upperBoundaries          = [2.0,  0.2, 2.0,  0.2];
s.lowerBoundaries          = [1, -1.0, 1, -1.0];
s.upperBoundaries          = [2.0,  0, 2.1,  0];

% 13) Define whether your anode is blend material (e.g. silicon-graphite) -> 
% set to true in this case. Assure that gammaAnBlend2_upperBound is not lower 
% than actual capacity share!
s.useAnodeBlend            = true;  % Turn on anode blend-electrode model by default
s.gammaAnBlend2_init       = 0.25;  % initial guess for gamma(Blend2)
s.gammaAnBlend2_upperBound = 0.30;  % e.g. 0.5 if Blend2 can reach 50% share

% Cathode blend option (disabled by default; enable and set paths below if needed)
s.useCathodeBlend          = false;
s.gammaCaBlend2_init       = 0.5;    % initial guess for cathode blend2 share
s.gammaCaBlend2_upperBound = 1.0; % example upper bound if cathode blend2 exists

% 14) Decide, whether to model inhomogeneity or not. Default is false.
% Setting to include or exclude inhomogeneity in the modeling.
% Do not use cathodeInhomogeneity for LFP cells!
s.allowAnodeInhomogeneity   = true;
s.allowCathodeInhomogeneity = true;
% Toggle to enable inhomogeneity already for the very first CU when active
s.allowFirstCycleInhomogeneity = true;
% set maximum allowed Inhomogeneity for both the anode and cathode
s.maxInhomogeneity          = 0.3;           
% maximum relative increase (Δ per CU), e.g. 0.1 for 10 % per CU
s.maxInhomogeneityDelta     = 0.1;          

% 15) Decide whether negative LAM per electrode or blend is allowed.
% The user can allow some small negative anode/cathode loss (gain). 
s.maxCathodeGain   = 0.010; % e.g. 0.01 -> allows 1% cathode gain per CU
s.maxAnodeGain     = 0.010; % e.g. 0.01 -> allows 1% anode gain per CU
s.maxAnBlend1Gain  = 0.005; % e.g. 0.005 -> allows 0.5% anode Blend1 gain per CU
s.maxAnBlend2Gain  = 0.010; % e.g. 0.01 -> allows 1% anode Blend2 gain per CU

% 16) Define maximal losses per CU (default 1.0 -> 100 % allowed per CU)
s.maxCathodeLoss   = 1; % e.g. 0.5 -> limit cathode loss to 50% per CU
s.maxAnodeLoss     = 1; % e.g. 0.8 -> limit anode loss to 80% per CU
s.maxAnBlend1Loss  = 1; % e.g. 0.6 -> limit anode Blend1 loss to 60% per CU
s.maxAnBlend2Loss  = 1; % e.g. 0.6 -> limit anode Blend2 loss to 60% per CU

% 17) Decide, which figures should be plotted:
% plots a figure showing both electrodes and the corresponding DV including
% their scaling and shifting parameters
s.cellName              = 'P45B';   % cell name for plot titles
s.plotParamCU           = true;    % Plot parameter overview for each CU
s.plotDMAOverall        = true;    % Plot overall DMA results
% special figure for blend electrodes
s.plotBlend2Overall     = false;    % Plot overall Blend2 content
s.plotAllAcceptedCU     = false;   % Plot all accepted solutions (not just best)
% Labels for plots (only relevant for final plotting scripts; can be overwritten)
s.labelCathode          = 'Cathode';
s.labelAnode            = 'Anode';
s.labelAnodeBlend1      = 'An-blend1';
s.labelAnodeBlend2      = 'An-blend2';
s.labelCathodeBlend1    = 'Ca-blend1';
s.labelCathodeBlend2    = 'Ca-blend2';
s.labelChargeCarrierInv = 'Charge-carrier-inv'; % inventory; e.g. Li-inventory
% Show cathode curve in DMA plot (hide only for LFP, where cathode aging is not meaningful)
s.plotCathodeInDMA      = true;

% 18) Adjust path for both anode and cathode (see the next 15-20 lines);
% then your done ;) 

%% -------------- LOAD anode and cathode data --------------
% Load anode and cathode data, depending on charge/discharge direction
%
% fill in for charge direction
if strcmp(s.direction, 'charge')
    % path to Anode -> needed in case of charge
    s.blendAn1Data_path = "./InputData/Graphite/Gr_Lithiation_Kuecher.mat";
    % blendAn2 path; only needed when anode blend fitting is enabled
    s.blendAn2Data_path = "./InputData/Silicon/" + ...
        "SiReconstr_Lithiation_Kuecher_P45B_Anode_0C03.mat";
    % path to Cathode -> needed in case of charge    
    s.blendCa1Data_path = "./InputData/NCA/" + ...
        "GITT_P45b_Cat_NCA_JN_VS_Coin_1_GITT__Extracted_Continuous_pOCP.mat";
    % blendCa2 path; only needed when cathode blend fitting is enabled
    s.blendCa2Data_path = ".\InputData\NMC\MRe_M36_PE_coin_Cathode_Delithiation_0C02.mat";
%    
% fill in for discharge direction
elseif strcmp(s.direction, 'discharge')
    s.blendAn1Data_path = "./InputData/Graphite/Gr_Delithiation_Kuecher.mat";
    % blendAn2 path; only needed when anode blend fitting is enabled
    s.blendAn2Data_path = "./InputData/Silicon/" + ...
        "SiReconstr_Kuecher_P45B_Anode_Delithiation_0C03.mat";
    % path to Cathode -> needed in case of discharge
    s.blendCa1Data_path = "./InputData/NCA/P45B_Cathode_Delithiation_0C03.mat";
    % blendCa2 path; only needed when cathode blend fitting is enabled
    s.blendCa2Data_path = ".\InputData\NMC\MRe_M36_PE_coin_Cathode_Delithiation_0C02.mat";
else
    error('Invalid direction. Must be ''charge'' or ''discharge''.');
end

% ==== end of user settings ====

% If s.useAnodeBlend is false, treat the anode as purely graphite.
if ~s.useAnodeBlend
    s.blendAn2Data_path       = s.blendAn1Data_path;
    s.gammaAnBlend2_init      = 0;
end

% If s.useCathodeBlend is false, treat the cathode as single-component (Blend1 only)
if ~s.useCathodeBlend
    s.blendCa2Data_path       = s.blendCa1Data_path;
    s.gammaCaBlend2_init      = 0;
end

%% ------ Optional: Load settings from input -----------
% overwrites all settings done outside of the function
if exist("userSettingsOutside","var") && ...
        ~isempty(userSettingsOutside) && ...
        isstruct(userSettingsOutside)
    configFields = fieldnames(userSettingsOutside);
    for i = 1:numel(configFields)
        s.(configFields{i}) = userSettingsOutside.(configFields{i});
    end
end

%% --------------- PATH HANDLING ---------------
% Add path of all HelperFunctions
addpath(['.' filesep 'HelperFunctions']);
% Add PlottingFunctions and create save paths via path handling function
[savePath_model, savePathDMA, savePathBlend2Content] = handle_paths(s);

%% --- Extract Cell Identifier Key & Value from input -----------------
[cellIdentifier_keys, cellIdentifier_values] = ...
    extract_cell_identifier(s.tableFilter);

%% -------------- PREPARE OUTPUT STRUCT + REFERENCE DATA --------------
Data                    = struct();  % output
refData                 = [];        % reference capacities from first CU
LAM_anode_blend1_prev   = [];        % for tracking minimal Blend1 LAM
LAM_cathode_prev        = [];        % for tracking minimal cathode LAM
inhom_An_prev           = [];        % for tracking anode inhomogeneity
inhom_Ca_prev           = [];        % for tracking cathode inhomogeneity

% Calculate ROI (region of interest of all the used fitting methods
[lowestROI, highestROI] = calculate_ROI(s);

% We’ll track which EFC indices actually succeed:
CUsUsed = [];

%% -------------- MAIN LOOP OVER THE CHECKUPS --------------
% We will assume the i-th folder is "CUi" (e.g. "CU1", "CU2", ...),
% and that i runs from 1..length(s.nCUs).
numCUs = numel(s.nCUs);

if s.nCUs(1) > s.nCUs(end)
    % Vector is descending (reverse order)
    idxRange = numCUs:-1:1;
    fitReverse = true;
else
    % Vector is ascending (forward order)
    idxRange = 1:numCUs;

    fitReverse = false;
end
refDataAssigned = false; % required for reverse fitting
fitResults = struct();

% calculate half cell data
half_and_full_cell_data = struct();
half_and_full_cell_data = calculate_half_cell_data(half_and_full_cell_data, s);
dataStore = struct();

for i = idxRange
    fieldName   = sprintf('CU%d', i);
    if i == 1
        if s.allowFirstCycleInhomogeneity
            allowAnodeInhomogeneity = s.allowAnodeInhomogeneity;
            allowCathodeInhomogeneity = s.allowCathodeInhomogeneity;
        else
            allowAnodeInhomogeneity = false;
            allowCathodeInhomogeneity = false;
        end
    else
        allowAnodeInhomogeneity = s.allowAnodeInhomogeneity;  % restore setting
        allowCathodeInhomogeneity = s.allowCathodeInhomogeneity;
    end

    if exist('dataStore','var') && isfield(dataStore,'agingDataTable')... 
        && ~isempty(dataStore.agingDataTable)
        vararg_get_fullcell_data = {'agingDataTable', dataStore.agingDataTable};
    else
        vararg_get_fullcell_data = {};
    end


    [fullCellData, toggleMissingData, dataStore] = ...
        get_fullcell_data(s, s.pathAgingStudy, cellIdentifier_keys, ...
        cellIdentifier_values, s.nameTableColumnOCV, ...
        i, vararg_get_fullcell_data{:});

    % if no Data is available for a certain CU/RPT skip this CU/RPT
    if toggleMissingData
        continue;
    end

    % calculate full cell data
    half_and_full_cell_data = ...
        calculate_full_cell_data(fullCellData, half_and_full_cell_data, s);

    % Attempt to find up to s.reqAccepted solutions under s.rmseThreshold
    acceptedCount       = 0;
    totalTries          = 0;
    bestRMSE_fitRegion  = Inf;
    bestSolution        = struct();
    bestTryRMSE         = Inf;
    bestTrySolution     = struct();

    % We'll also store all accepted solutions if s.plotAllAcceptedCU
    allAcceptedSolutions = {};
    % Slightly dynamic upper bound for gamma(Blend2):
        if s.useAnodeBlend
            if s.gammaAnBlend2_upperBound < 0.0001
                warning(['gammaAnBlend2_upperBound is very low (%.4f). ' ...
                    'Please check.'], s.gammaAnBlend2_upperBound);
            end
        else
            s.gammaAnBlend2_upperBound = 0;  % pure Blend1
        end

    % ========== Start tries loop ==========
    if s.measureRuntimePerCU
        cuTimer = tic;
    end
    while (acceptedCount < s.reqAccepted) && (totalTries < s.maxTriesOverall)
        totalTries = totalTries + 1;

        if isempty(refData) || isempty(LAM_anode_blend1_prev)
            % First CU or no reference data: pass empty
            refArgs = {[], [], [], [], [], [], []};
        else
            % Later CUs: pass actual reference data
            refArgs = {refData, LAM_anode_prev, LAM_cathode_prev, ...
                LAM_anode_blend1_prev, LAM_anode_blend2_prev, ...
                inhom_An_prev, inhom_Ca_prev};
        end

        [meas_data, reconSOC, fcU_model, Q_DVA_meas, DVA_smooth_meas, ...
            Q_DVA_calc, DVA_smooth_calc, Q_ICA_meas, ICA_smooth_meas,...
            Q_ICA_calc, ICA_smooth_calc, capa_act, params, ~, ~, ...
            cathSOC, cathU_recon, anodeSOC, anodeU_recon] = ...
            dma_core(half_and_full_cell_data, s, ...
            refArgs{:}, allowAnodeInhomogeneity, allowCathodeInhomogeneity,...
            fitReverse);

        % Evaluate RMSE in the "fit region" (RMSE in the SOC region between
        % highest and lowest SOC used for OCV, DVA or ICA fitting) and
        % between 0 and 100% SOC
        rmse_fitRegion = calculate_RMSE(meas_data.OCV_cell, fcU_model, ...
            [lowestROI, highestROI]);
        rmse_fullRange = calculate_RMSE(meas_data.OCV_cell, fcU_model, [0, 1]);
        % RMSE between DVA measured and reconstructed
        rmse_DVA_full = calculate_RMSE(DVA_smooth_meas, DVA_smooth_calc, ...
            [0, 1]);

        fprintf('%s -> Try #%d: RMSE_fit = %.4f V | RMSE_full = %.4f V\n', ...
            fieldName, totalTries, rmse_fitRegion, rmse_fullRange);

        % Package solution attempt for later evaluation / storage
        candidateSol = store_solution_struct( ...
            meas_data, reconSOC, fcU_model, ...
            Q_DVA_meas, DVA_smooth_meas, ...
            Q_DVA_calc, DVA_smooth_calc, ...
            Q_ICA_meas, ICA_smooth_meas,...
            Q_ICA_calc, ICA_smooth_calc,...
            capa_act, params, ...
            s.algorithm, s.weightOCV, ...
            s.weightDVA, s.weightICA,...
            cathSOC, cathU_recon, ...
            anodeSOC, anodeU_recon,...
            rmse_fitRegion, rmse_fullRange, rmse_DVA_full);

        candidateSol.isAccepted = rmse_fitRegion < s.rmseThreshold;
        if candidateSol.isAccepted
            candidateSol.status = 'accepted';
        else
            candidateSol.status = 'rejected_above_threshold';
        end

        if rmse_fitRegion < bestTryRMSE
            bestTryRMSE     = rmse_fitRegion;
            bestTrySolution = candidateSol;
        end

        % Check acceptance
        if candidateSol.isAccepted
            acceptedCount = acceptedCount + 1;
            fprintf('  (ACCEPTED %d/%d)\n', acceptedCount, s.reqAccepted);

            if s.plotAllAcceptedCU
                allAcceptedSolutions{acceptedCount} = candidateSol; %#ok<AGROW>
            end

            % Update "best" if this solution shows the lowest RMSE within
            % SOC region used for fitting procedure
            if rmse_fitRegion < bestRMSE_fitRegion
                bestRMSE_fitRegion = rmse_fitRegion;
                bestSolution       = candidateSol;
            end
        else
            fprintf('    (REJECTED)\n');
        end
    end
    % ========== End tries loop ==========
    if s.measureRuntimePerCU
        cuRuntimeSeconds = toc(cuTimer);
        fprintf('%s runtime (single CU run): %.3f s\n', fieldName, cuRuntimeSeconds);
    end
    
    % case that no solution was below RMSE threshold, but number of maximum
    % tries was reached
    if isinf(bestTryRMSE)
        warning('%s: No optimisation after %d tries; skipping CU.', ...
            fieldName, totalTries);
        continue;
    end

    if isinf(bestRMSE_fitRegion)
        % No accepted solution -> fall back to best attempt above threshold
        finalSolution                = bestTrySolution;
        finalSolution.isAccepted     = false;
        finalSolution.status         = 'best_attempt_above_threshold';
        finalSolution.finalRMSE_fit  = finalSolution.rmse_fitRegion;
        finalSolution.finalRMSE_full = finalSolution.rmse_fullRange;

        warning('%s: No solution < %.4f V after %d tries; best %.4f V.', ...
            fieldName, s.rmseThreshold, totalTries, ...
            finalSolution.finalRMSE_fit);
        fprintf('%s -> Best attempt: RMSE_fit=%.4f V | RMSE_0_9=%.4f V\n', ...
            fieldName, finalSolution.finalRMSE_fit, ...
            finalSolution.finalRMSE_full);
    else
        % We have at least one accepted solution
        finalSolution                = bestSolution;
        finalSolution.isAccepted     = true;
        finalSolution.status         = 'accepted';
        finalSolution.finalRMSE_fit  = bestRMSE_fitRegion;
        finalSolution.finalRMSE_full = finalSolution.rmse_fullRange;

        fprintf(['%s -> Best accepted #%d: RMSE_fit=%.4f V', ...
            ' | RMSE_0_9=%.4f V\n'], fieldName, acceptedCount, ...
            finalSolution.finalRMSE_fit, finalSolution.finalRMSE_full);
    end

    % plot the current scaling and shifting of the OCPs for this CU
    if s.plotParamCU
        try
            plot_OCV_model_param_show(finalSolution, i, s.cellName, s);
            label = sprintf('%s_CU%d',s.cellName, i);
            save_figure(savePath_model, label, 'savepng', 1);
        catch ME
            warning('Plot parameter error for %s: %s', s.cellName, ME.message);
        end
    end

    % (Optional) plot the other accepted solutions (non-selected) if desired
    if s.plotAllAcceptedCU && ~isempty(allAcceptedSolutions)
        for idxSol = 1:numel(allAcceptedSolutions)
            thisSol = allAcceptedSolutions{idxSol};
            % Skip re-plotting the best solution
            if abs(thisSol.rmse_fitRegion - bestRMSE_fitRegion) < 1e-14
                continue;
            end

            if s.plotParamCU
                plot_OCV_model_param_show(thisSol, i, s.cellName, s);
                ax = gca;
                oldTitle = get(get(ax,'Title'),'String');
                set(get(ax,'Title'),'String',[oldTitle ' (NON-SELECTED)']);
            end
        end
    end

    fitResults.(fieldName) = finalSolution;
    % If this is the first successful CU, store reference data
    if ~refDataAssigned  % or if ~isfield(Data, 'CU1'), depending on naming
        if s.useAnodeBlend
            gamma_an_blend2_init = finalSolution.params(5);
        else
            gamma_an_blend2_init = 0;
        end
        if s.useCathodeBlend
            gamma_ca_blend2_init = finalSolution.params(6);
        else
            gamma_ca_blend2_init = 0;
        end
        capa_anode_init   = finalSolution.params(1) * finalSolution.capa_act;
        capa_cathode_init = finalSolution.params(3) * finalSolution.capa_act;
        capa_inventory_init = ( finalSolution.params(3) + ...
            finalSolution.params(4) - ...
            finalSolution.params(2) ) * finalSolution.capa_act;

        refData = struct('gamma_an_blend2_init',       gamma_an_blend2_init, ...
            'capa_anode_init',         capa_anode_init, ...
            'capa_cathode_init',       capa_cathode_init, ...
            'capa_inventory_init',     capa_inventory_init, ...
            'gamma_ca_blend2_init',    gamma_ca_blend2_init);
        refDataAssigned = true;   % required for reverse fitting
    end

    % Calculate final degradation metrics
    [LAM_anode, LAM_cathode, ~, LAM_anode_blend2, LAM_anode_blend1, ...
        LAM_cathode_blend2, LAM_cathode_blend1] = ...
        calculate_degradation_modes( ...
        finalSolution.params, ...
        finalSolution.capa_act, ...
        refData.capa_anode_init, ...
        refData.capa_cathode_init, ...
        refData.capa_inventory_init, ...
        refData.gamma_an_blend2_init, refData.gamma_ca_blend2_init, fitReverse);

    % Track for next iteration
    LAM_anode_prev        = LAM_anode;
    LAM_cathode_prev      = LAM_cathode;
    LAM_anode_blend1_prev = LAM_anode_blend1;
    LAM_anode_blend2_prev = LAM_anode_blend2;

    if numel(finalSolution.params) >= 7
        inhom_An_prev = finalSolution.params(7);
    else
        inhom_An_prev = 0;
    end
    if numel(finalSolution.params) >= 8
        inhom_Ca_prev = finalSolution.params(8);
    else
        inhom_Ca_prev = 0;
    end

    % Record the EFC corresponding to i-th CU
    CUsUsed(end+1) = s.nCUs(i); %#ok<AGROW>
end

% -------------------- added due to the reverse fitting --------------------
fitResults = orderfields(fitResults);

% Take CU1 as reference
refCU = 'CU1';
refSol = fitResults.(refCU);

gamma_an_blend2_init       = refSol.params(5);
gamma_ca_blend2_init       = refSol.params(6);
capa_anode_init         = refSol.params(1) * refSol.capa_act;
capa_cathode_init       = refSol.params(3) * refSol.capa_act;
capa_inventory_init     = ( refSol.params(3) + refSol.params(4) ...
    - refSol.params(2) ) * refSol.capa_act;

refData = struct('gamma_an_blend2_init',       gamma_an_blend2_init, ...
    'capa_anode_init',         capa_anode_init, ...
    'capa_cathode_init',       capa_cathode_init, ...
    'capa_inventory_init',     capa_inventory_init, ...
    'gamma_ca_blend2_init',    gamma_ca_blend2_init);

% Now: calculate LAM for ALL CUs using this reference
fields = fieldnames(fitResults);
for k = 1:numel(fields)
    fName = fields{k};
    sol = fitResults.(fName);
    [LAM_anode, LAM_cathode, LI, LAM_anode_blend2, LAM_anode_blend1, ...
        LAM_cathode_blend2, LAM_cathode_blend1] = ...
        calculate_degradation_modes(sol.params, sol.capa_act, ...
        refData.capa_anode_init, ...
        refData.capa_cathode_init, ...
        refData.capa_inventory_init, ...
        refData.gamma_an_blend2_init, refData.gamma_ca_blend2_init);

    % Store everything in Data
    [Data] = store_data_fields(Data, fName, sol, LAM_anode, LAM_cathode, ...
        LI, LAM_anode_blend2, LAM_anode_blend1, LAM_cathode_blend2, LAM_cathode_blend1);
    if isfield(sol,'isAccepted')
        Data.(fName).isAccepted = sol.isAccepted;
    end
    if isfield(sol,'status')
        Data.(fName).status     = sol.status;
    end
    Data.(fName).rmseThreshold = s.rmseThreshold;
end

% additionally store settings of the DMA (s)
Data.s = s;

% -------------- OVERALL DMA PLOTS --------------
if s.plotDMAOverall && ~isempty(fieldnames(Data))
    labelStruct = struct('labelCathode', s.labelCathode, ...
        'labelAnode', s.labelAnode, ...
        'labelAnodeBlend1', s.labelAnodeBlend1, ...
        'labelAnodeBlend2', s.labelAnodeBlend2, ...
        'labelCathodeBlend1', s.labelCathodeBlend1, ...
        'labelCathodeBlend2', s.labelCathodeBlend2, ...
        'labelChargeCarrierInv', s.labelChargeCarrierInv);
    plot_DMA(Data, 'useAnodeBlendModel', s.useAnodeBlend, ...
        'useCathodeBlendModel', s.useCathodeBlend, ...
        'CalendarOrCyclic', s.CalendarOrCyclic, 'EFC', CUsUsed, ...
        'labels', labelStruct, 'plotCathode', s.plotCathodeInDMA);

    label = sprintf('%s_DMA', s.cellName);
    save_figure(savePathDMA, label);

    label = sprintf('%s_DMA', s.cellName);
    save_figure(savePathDMA, label);
end

% -------------- OVERALL gamma PLOT --------------
if s.plotBlend2Overall && ~isempty(fieldnames(Data)) && s.useAnodeBlend
    plot_gamma(Data, CUsUsed, s.CalendarOrCyclic);
    if ~exist(savePathBlend2Content, 'dir')
        mkdir(savePathBlend2Content);  % Create the directory if necessary
    end

    label = sprintf('%s_SiContent', s.cellName);
    save_figure(savePathBlend2Content, label);
end

% ------ Saving parameters and fitting data ------
filenNameParameters = fullfile(s.pathSaveResults, ...
    sprintf('%s_parameters.mat', s.cellName));
save(filenNameParameters, "Data");

end % end of DMA_main
