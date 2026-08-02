%> regenerate_example_silicon_curves
%>
%> Regenerate the two bundled silicon OCP example curves under
%> input_data/silicon/ from the bundled blend and graphite curves with
%> generate_si_ocp. The two curves are the anode blend-2 input of the
%> minimal working example in main_dma, so this script documents how they
%> were made and keeps them reproducible.
%>
%> Recipe (kept exactly as the curves that shipped with 1.0.0 were made,
%> including the asymmetry between the two directions):
%>
%>   lithiation    blend    P45B_Anode_Lithiation_0C03.mat
%>                 graphite Gr_Lithiation_Rehm2026.mat
%>                 gammaSi  0.2446
%>                 filters  LOWESS on blend and on graphite
%>                 PAV      off
%>
%>   delithiation  blend    P45B_Anode_Delithiation_0C03.mat
%>                 graphite Gr_Delithiation_Rehm2026.mat
%>                 gammaSi  0.23   (Si2-peak value, see generate_si_ocp)
%>                 filters  off
%>                 PAV      off
%>
%> Provenance: the lithiation curve that shipped with 1.0.0 could not be
%> reproduced from the bundled inputs because the graphite reference lost
%> one off-grid sample (row 2342 at 0.695265 V) in the 2025-12-18 clean-up.
%> Regenerating against today's graphite moves the curve by at most
%> 6.16e-5 in normalized capacity, inside a 40 mV band around 0.695 V.
%>
%> Run from the project root.

%% Setup
addpath(fullfile('.', 'generate_si_ocp'));

blendDir    = fullfile('.', 'input_data', 'silicon-graphite');
graphiteDir = fullfile('.', 'input_data', 'graphite');
siliconDir  = fullfile('.', 'input_data', 'silicon');

%% Recipe, one row per direction
% direction | blend file | graphite file | gammaSi | input filters | target
recipe = { ...
    'lithiation'  , 'P45B_Anode_Lithiation_0C03.mat'  , 'Gr_Lithiation_Rehm2026.mat'  , ...
        0.2446, true , 'SiReconstr_Lithiation_Rehm2026_P45B_Anode_0C03.mat'   ; ...
    'delithiation', 'P45B_Anode_Delithiation_0C03.mat', 'Gr_Delithiation_Rehm2026.mat', ...
        0.23  , false, 'SiReconstr_Delithiation_Rehm2026_P45B_Anode_0C03.mat' };

%% Regenerate both curves
for k = 1:size(recipe, 1)
    direction    = recipe{k,1};
    blendPath    = fullfile(blendDir, recipe{k,2});
    graphitePath = fullfile(graphiteDir, recipe{k,3});
    gammaSi      = recipe{k,4};
    filterInputs = recipe{k,5};
    savePath     = fullfile(siliconDir, recipe{k,6});

    silicon = generate_si_ocp( ...
        'blendPath', blendPath, ...
        'graphitePath', graphitePath, ...
        'lithDirection', direction, ...
        'gammaSi', gammaSi, ...
        'filterBlend', filterInputs, ...
        'filterGraphite', filterInputs, ...
        'pavOutput', false, ...
        'interpolationReadyOutput', false, ...
        'plotFlag', false, ...
        'savePath', savePath);

    fprintf(['%-13s gammaSi %.4f, filters %d, PAV off -> %d samples, ' ...
             'U = [%.6g %.6g] V, q = [%.4g %.4g]\n'], ...
        direction, gammaSi, filterInputs, numel(silicon.voltage), ...
        min(silicon.voltage), max(silicon.voltage), ...
        min(silicon.normalizedCapacity), max(silicon.normalizedCapacity));
    fprintf('              written to %s\n', savePath);
end
