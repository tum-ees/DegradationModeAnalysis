%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author  : Mathias Rehm
%  E-mail  : mathias.rehm@tum.de
%
%  Script : test_generate_si_ocp_outputs
%  ------------------------------------------------------------------------
%  Tests for the optional second output of generate_si_ocp, the graphite
%  reference on the common voltage grid, and for the save behaviour of an
%  empty 'savePath'.
%
%  Run it as a script from any folder:
%     run('<repo>/generate_si_ocp/test_generate_si_ocp_outputs.m')
%  as a test suite:
%     runtests('<repo>/generate_si_ocp/test_generate_si_ocp_outputs.m')
%
%  Every section below is one test case and is self-contained. runtests
%  re-runs the shared block above the first section before each section and
%  does not carry variables from one section into the next.
%
%  Note on the reconstruction relation: qBlend = gamma*qSi + (1-gamma)*qGr
%  holds for the raw reconstruction only, that is with 'pavOutput' false and
%  at the samples that were not clipped to [0,1]. Under the default settings
%  the silicon curve departs from it on purpose, by up to 0.15 on the
%  bundled example data, which is exactly the clipping and the isotonic
%  monotonicity repair. It is therefore not asserted here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisDir  = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
addpath(thisDir);

workDir = fullfile(tempdir, 'test_generate_si_ocp_outputs');
if ~exist(workDir, 'dir')
    mkdir(workDir);
end

blendDir = fullfile(projRoot, 'input_data', 'silicon-graphite');
graphDir = fullfile(projRoot, 'input_data', 'graphite');
exampleCases = { ...
    'lithiation'  , fullfile(blendDir,'P45B_Anode_Lithiation_0C03.mat')  , ...
                    fullfile(graphDir,'Gr_Lithiation_Rehm2026.mat')      ; ...
    'delithiation', fullfile(blendDir,'P45B_Anode_Delithiation_0C03.mat'), ...
                    fullfile(graphDir,'Gr_Delithiation_Rehm2026.mat')    };
gammaSi = 0.245;

%% 02 second output on the P45B example data
for k = 1:size(exampleCases,1)
    direction = exampleCases{k,1};
    blendPath = exampleCases{k,2};
    graphitePath = exampleCases{k,3};
    if ~isfile(blendPath) || ~isfile(graphitePath)
        fprintf('SKIP  P45B %s (example data not found)\n', direction);
        continue
    end

    [silicon, graphite] = runPlain(blendPath, graphitePath, direction, gammaSi);

    for f = {'voltage','normalizedCapacity','meta'}
        assert(isfield(graphite, f{1}), 'graphiteStruct lacks %s.', f{1});
    end
    for f = {'source','direction','path','Vwindow'}
        assert(isfield(graphite.meta, f{1}), 'graphiteStruct.meta lacks %s.', f{1});
    end

    % the graphite curve lives on the common grid, which the default silicon
    % output shares
    assert(isequal(graphite.voltage, silicon.voltage), ...
        '%s: the graphite grid does not match the silicon output.', direction);
    assert(numel(graphite.normalizedCapacity) == numel(graphite.voltage), ...
        '%s: capacity and voltage lengths differ.', direction);
    assert(all(isfinite(graphite.normalizedCapacity)), ...
        '%s: graphite capacity contains non-finite values.', direction);

    % meta describes the aligned window the curve was built on
    assert(abs(graphite.meta.Vwindow(1) - min(graphite.voltage)) < 1e-12 && ...
           abs(graphite.meta.Vwindow(2) - max(graphite.voltage)) < 1e-12, ...
        '%s: Vwindow does not match the common grid.', direction);
    assert(strcmp(graphite.meta.direction, direction), ...
        '%s: graphite meta carries the wrong direction.', direction);

    fprintf('PASS  %-14s %d samples on the common grid, window [%.4g %.4g] V\n', ...
        direction, numel(graphite.voltage), graphite.meta.Vwindow);
end

%% 03 meta.source tells resolved references apart from an explicit path
direction = exampleCases{1,1};
blendPath = exampleCases{1,2};
graphitePath = exampleCases{1,3};
if isfile(blendPath) && isfile(graphitePath)
    % explicit file: no built-in reference was resolved
    [~, explicitGraphite] = runPlain(blendPath, graphitePath, direction, gammaSi);
    assert(strcmp(explicitGraphite.meta.source, 'custom'), ...
        'An explicit graphitePath reports meta.source "%s".', ...
        explicitGraphite.meta.source);
    assert(strcmp(explicitGraphite.meta.path, graphitePath), ...
        'meta.path does not repeat the path that was passed in.');

    % resolved from a built-in reference
    [~, resolvedGraphite] = generate_si_ocp('blendPath', blendPath, ...
        'graphiteSource', 'Rehm2026', 'lithDirection', direction, ...
        'gammaSi', gammaSi, 'filterBlend', false, 'filterGraphite', false, ...
        'pavOutput', false, 'interpolationReadyOutput', false, ...
        'plotFlag', false, 'savePath', "");
    assert(strcmp(resolvedGraphite.meta.source, 'Rehm2026'), ...
        'A resolved reference reports meta.source "%s".', ...
        resolvedGraphite.meta.source);
    assert(isfile(resolvedGraphite.meta.path), ...
        'meta.path of a resolved reference is not a file.');
    assert(isequal(explicitGraphite.normalizedCapacity, ...
                   resolvedGraphite.normalizedCapacity), ...
        'The same file gives different curves through the two path routes.');

    fprintf('PASS  meta.source: explicit path -> "custom", reference -> "Rehm2026"\n');
end

%% 04 the first output is unchanged by asking for the second
direction = exampleCases{1,1};
blendPath = exampleCases{1,2};
graphitePath = exampleCases{1,3};
if isfile(blendPath) && isfile(graphitePath)
    oneOut = generate_si_ocp('blendPath', blendPath, ...
        'graphitePath', graphitePath, 'lithDirection', direction, ...
        'gammaSi', gammaSi, 'filterBlend', false, 'filterGraphite', false, ...
        'pavOutput', true, 'interpolationReadyOutput', true, ...
        'plotFlag', false, 'savePath', "");
    [twoOut, graphite] = generate_si_ocp('blendPath', blendPath, ...
        'graphitePath', graphitePath, 'lithDirection', direction, ...
        'gammaSi', gammaSi, 'filterBlend', false, 'filterGraphite', false, ...
        'pavOutput', true, 'interpolationReadyOutput', true, ...
        'plotFlag', false, 'savePath', "");
    assert(isequal(oneOut, twoOut), ...
        'The silicon output differs between a one-output and a two-output call.');

    % interpolation-ready resamples the silicon table only, the graphite
    % curve stays on the common voltage grid
    assert(numel(twoOut.voltage) == 1001, ...
        'The interpolation-ready silicon table does not have 1001 points.');
    assert(~isequal(graphite.voltage, twoOut.voltage), ...
        'The graphite curve was resampled onto the silicon grid.');
    fprintf('PASS  one-output call identical, graphite keeps the voltage grid\n');
end

%% 05 savePath controls whether a file is written
direction = exampleCases{1,1};
blendPath = exampleCases{1,2};
graphitePath = exampleCases{1,3};
if isfile(blendPath) && isfile(graphitePath)
    saveDir = fullfile(workDir, 'save_probe');
    if exist(saveDir, 'dir')
        rmdir(saveDir, 's');
    end
    mkdir(saveDir);

    % empty savePath: nothing is written and no dialog is opened
    runPlain(blendPath, graphitePath, direction, gammaSi);
    written = dir(fullfile(saveDir, '*.mat'));
    assert(isempty(written), 'An empty savePath wrote %d file(s).', numel(written));

    % non-empty savePath: both structs land in the file
    target = fullfile(saveDir, 'si_curve.mat');
    generate_si_ocp('blendPath', blendPath, 'graphitePath', graphitePath, ...
        'lithDirection', direction, 'gammaSi', gammaSi, ...
        'filterBlend', false, 'filterGraphite', false, ...
        'pavOutput', false, 'interpolationReadyOutput', false, ...
        'plotFlag', false, 'savePath', target);
    assert(isfile(target), 'savePath did not produce %s.', target);
    stored = {whos('-file', target).name};
    for f = {'siliconStruct','graphiteStruct'}
        assert(any(strcmp(stored, f{1})), 'The saved file lacks %s.', f{1});
    end
    fprintf('PASS  empty savePath writes nothing, a set savePath stores 2 structs\n');
end

%% 06 cleanup
if exist(workDir, 'dir')
    rmdir(workDir, 's');
end

% ---- local helper functions ----------------------------------------------
% Declared without a section header on purpose: a "%%" line here would show
% up as an extra, empty test case under runtests.

function [silicon, graphite] = runPlain(blendPath, graphitePath, direction, gammaSi)
% Default output path, no filtering, no PAV, nothing saved.
    [silicon, graphite] = generate_si_ocp( ...
        'blendPath', blendPath, ...
        'graphitePath', graphitePath, ...
        'lithDirection', direction, ...
        'gammaSi', gammaSi, ...
        'filterBlend', false, 'filterGraphite', false, ...
        'pavOutput', false, 'interpolationReadyOutput', false, ...
        'plotFlag', false, ...
        'savePath', "");
end
