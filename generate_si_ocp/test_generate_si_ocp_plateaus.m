%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author  : Mathias Rehm
%  E-mail  : mathias.rehm@tum.de
%
%  Script : test_generate_si_ocp_plateaus
%  ------------------------------------------------------------------------
%  Regression tests for the plateau handling of the interpolation-ready
%  silicon OCP output ('interpolationReadyOutput', true).
%
%  Run it as a script from any folder:
%     run('<repo>/generate_si_ocp/test_generate_si_ocp_plateaus.m')
%  as a test suite:
%     runtests('<repo>/generate_si_ocp/test_generate_si_ocp_plateaus.m')
%  or in batch mode:
%     matlab -batch "run('<repo>/generate_si_ocp/test_generate_si_ocp_plateaus.m')"
%
%  Every section below is one test case and is self-contained. runtests
%  re-runs the shared block above the first section before each section and
%  does not carry variables from one section into the next, so no section
%  may rely on another one having run first.
%
%  The tests protect one property above all: the interpolation-ready output
%  must never silently lose voltage support. Collapsing the PAV plateaus
%  moves capacity coordinates around, and a shift that leaves the capacity
%  range used to be clamped back onto the range boundary. That re-created
%  exact ties, and the deduplication that followed dropped the tied samples
%  together with their voltages, which truncated the exported curve while
%  every assertion inside the generator stayed green.
%
%  Cases 1 to 3 are synthetic curves whose plateau levels sit closer
%  together than twice the tie-breaking shift, which is the configuration
%  that triggered the truncation. Case 4 is a smoke test on the measured
%  P45B data shipped with the repository. Cases 5 and 6 drive the plateau
%  levels down to a few ulps, where the shift is no longer representable
%  and the collapse has to pool the levels instead of shifting them apart.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisDir  = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
addpath(thisDir);

workDir = fullfile(tempdir, 'test_generate_si_ocp_plateaus');
if ~exist(workDir, 'dir')
    mkdir(workDir);
end

voltageSamples = (1:5).';

%% 02 rising boundary plateau
% Two plateaus 1e-6 apart, the upper one clipped at q = 1. The shifts of the
% two runs used to cross, the forward sweep pushed the samples past q = 1
% and the clamp collapsed them onto q = 1. The export then kept only the
% first of the tied samples and the curve ended at 3 V instead of 5 V.
qGraphite = linspace(0,1,5).';
qSilicon  = [0; 0.999999; 0.999999; 1; 1];
[blendPath, graphitePath] = writePair(workDir, 'rise', voltageSamples, ...
    qGraphite, qSilicon);

silicon = runGenerator(blendPath, graphitePath, 'lithiation', 0.5, workDir);
checkInterpolationReady(silicon, 1, 5, 'rising boundary plateau');

%% 03 falling boundary plateau
% Mirror image of case 1 on the qMin side of a delithiation curve. The
% values are dyadic so that the reconstruction q_Si = 2*q_blend - q_Gr is
% exact and the plateaus really are runs of bit-identical values. The
% falling direction is not structurally safer than the rising one: before
% the fix this construction lost the 4 V and 5 V samples the same way.
h = 2^-20;
qGraphite = [1; 0.75; 0.5; 0.25; 0];
qSilicon  = [1; h; h; 0; 0];
[blendPath, graphitePath] = writePair(workDir, 'fall', voltageSamples, ...
    qGraphite, qSilicon);

silicon = runGenerator(blendPath, graphitePath, 'delithiation', 0.5, workDir);
checkInterpolationReady(silicon, 1, 5, 'falling boundary plateau');

%% 04 interior near-coincident plateaus
% Levels 8e-6 apart, far away from both ends of the capacity range, so no
% clamping is involved. Here the shifts must stay small enough that the two
% plateaus keep their order and their levels.
voltageInterior = (1:9).';
qGraphite = linspace(0,1,9).';
qSilicon  = [0; 0.2; 0.4; 0.4; 0.400008; 0.400008; 0.7; 0.9; 1];
[blendPath, graphitePath] = writePair(workDir, 'interior', voltageInterior, ...
    qGraphite, qSilicon);

silicon = runGenerator(blendPath, graphitePath, 'lithiation', 0.5, workDir);
checkInterpolationReady(silicon, 1, 9, 'interior near-coincident plateaus');

% Both plateau levels must survive: reading the table back at the mid
% voltage of each plateau has to return its level to within one output grid
% cell plus the tie-breaking shift.
gridCell = (max(silicon.normalizedCapacity) - min(silicon.normalizedCapacity)) / 1000;
tolLevel = gridCell + 2e-5;
for plateau = [3.5 0.4; 5.5 0.400008].'
    qBack = interp1(silicon.voltage, silicon.normalizedCapacity, plateau(1));
    assert(abs(qBack - plateau(2)) <= tolLevel, ...
        'interior plateau at q = %.8g moved by %.3g (tolerance %.3g).', ...
        plateau(2), abs(qBack - plateau(2)), tolLevel);
end

%% 05 P45B example data
% The interpolation-ready output has to cover exactly the voltage window of
% the default output, in both directions.
blendDir = fullfile(projRoot, 'input_data', 'silicon-graphite');
graphDir = fullfile(projRoot, 'input_data', 'graphite');
realCases = { ...
    'lithiation'  , fullfile(blendDir,'P45B_Anode_Lithiation_0C03.mat')  , ...
                    fullfile(graphDir,'Gr_Lithiation_Rehm2026.mat')      ; ...
    'delithiation', fullfile(blendDir,'P45B_Anode_Delithiation_0C03.mat'), ...
                    fullfile(graphDir,'Gr_Delithiation_Rehm2026.mat')    };

for k = 1:size(realCases,1)
    direction = realCases{k,1};
    blendPath = realCases{k,2};
    graphitePath = realCases{k,3};
    if ~isfile(blendPath) || ~isfile(graphitePath)
        fprintf('SKIP  P45B %s (example data not found)\n', direction);
        continue
    end

    reference = runGeneratorPlain(blendPath, graphitePath, direction, 0.245, workDir);
    silicon = runGenerator(blendPath, graphitePath, direction, 0.245, workDir);
    checkInterpolationReady(silicon, min(reference.voltage), ...
        max(reference.voltage), ['P45B ' direction]);
end

%% 06 ulp-close plateau levels at the upper boundary
% Two plateaus whose levels are two ulps apart at q = 1. A quarter of that
% gap is half an ulp, so the shift rounded the endpoint back onto its own
% level, the pair tied, and the safety sweep pushed the last sample one ulp
% past the capacity range: the export ended at q = 1.0000000000000002.
ulpGap = 2^-52;
qGraphite = [0; 0.25; 0.5; 0.75; 1];
qSilicon  = [0; 1-ulpGap; 1-ulpGap; 1; 1];
[blendPath, graphitePath] = writePair(workDir, 'ulpmax', voltageSamples, ...
    qGraphite, qSilicon);

silicon = runGenerator(blendPath, graphitePath, 'lithiation', 0.5, workDir);
checkInterpolationReady(silicon, 1, 5, 'ulp-close levels at qMax');
checkCapacityRange(silicon, 'ulp-close levels at qMax');

%% 07 ulp-close plateau levels on the falling side
% Mirror of case 5 on a delithiation curve. This one never left the range,
% because the sweep pushes towards the range interior at that end, but it
% did make the sweep fire. After the pooling it does not.
ulpGap = 2^-52;
qGraphite = [1; 0.75; 0.5; 0.25; 0];
qSilicon  = [1; 1; 1-ulpGap; 1-ulpGap; 0];
[blendPath, graphitePath] = writePair(workDir, 'ulpmin', voltageSamples, ...
    qGraphite, qSilicon);

silicon = runGenerator(blendPath, graphitePath, 'delithiation', 0.5, workDir);
checkInterpolationReady(silicon, 1, 5, 'ulp-close levels, falling');
checkCapacityRange(silicon, 'ulp-close levels, falling');

%% 08 cleanup
if exist(workDir, 'dir')
    rmdir(workDir, 's');
end

% ---- local helper functions ----------------------------------------------
% Declared without a section header on purpose: a "%%" line here would show
% up as an extra, empty test case under runtests.

function [blendPath, graphitePath] = writePair(workDir, tag, voltage, ...
        qGraphite, qSilicon)
% Write a graphite/blend pair for gammaSi = 0.5, so that the generator
% reconstructs qSilicon as 2*qBlend - qGraphite.
    graphiteCurve.voltage = voltage;
    graphiteCurve.normalizedCapacity = qGraphite;
    blendCurve.voltage = voltage;
    blendCurve.normalizedCapacity = 0.5*qGraphite + 0.5*qSilicon;

    graphitePath = fullfile(workDir, ['gr_' tag '.mat']);
    blendPath    = fullfile(workDir, ['bl_' tag '.mat']);
    save(graphitePath, 'graphiteCurve');
    save(blendPath, 'blendCurve');
end

function silicon = runGenerator(blendPath, graphitePath, direction, gammaSi, workDir)
% Interpolation-ready output on top of the PAV result.
    silicon = generate_si_ocp('blendPath', blendPath, ...
        'graphitePath', graphitePath, ...
        'lithDirection', direction, ...
        'gammaSi', gammaSi, ...
        'filterBlend', false, 'filterGraphite', false, ...
        'pavOutput', true, 'interpolationReadyOutput', true, ...
        'plotFlag', false, ...
        'savePath', fullfile(workDir, 'si_interpolation_ready.mat'));
end

function silicon = runGeneratorPlain(blendPath, graphitePath, direction, gammaSi, workDir)
% Default output, used as the voltage-support reference.
    silicon = generate_si_ocp('blendPath', blendPath, ...
        'graphitePath', graphitePath, ...
        'lithDirection', direction, ...
        'gammaSi', gammaSi, ...
        'filterBlend', false, 'filterGraphite', false, ...
        'pavOutput', true, 'interpolationReadyOutput', false, ...
        'plotFlag', false, ...
        'savePath', fullfile(workDir, 'si_default.mat'));
end

function checkInterpolationReady(silicon, voltageLow, voltageHigh, tag)
% Assert the contract of the interpolation-ready output.
    voltage = silicon.voltage;
    capacity = silicon.normalizedCapacity;
    tol = 1e-9 * max(1, abs(voltageHigh - voltageLow));

    assert(numel(voltage) == 1001 && numel(capacity) == 1001, ...
        '%s: expected 1001 samples, got %d.', tag, numel(voltage));
    assert(all(isfinite(voltage)) && all(isfinite(capacity)), ...
        '%s: output contains non-finite values.', tag);
    assert(all(diff(capacity) > 0), ...
        '%s: capacity is not strictly increasing.', tag);
    assert(all(diff(voltage) > 0) || all(diff(voltage) < 0), ...
        '%s: voltage is not strictly monotone.', tag);
    assert(abs(min(voltage) - voltageLow) <= tol, ...
        '%s: lost the lower voltage support (%.12g instead of %.12g).', ...
        tag, min(voltage), voltageLow);
    assert(abs(max(voltage) - voltageHigh) <= tol, ...
        '%s: lost the upper voltage support (%.12g instead of %.12g).', ...
        tag, max(voltage), voltageHigh);

    fprintf('PASS  %-34s V = [%.6g %.6g], q = [%.6g %.6g]\n', tag, ...
        min(voltage), max(voltage), min(capacity), max(capacity));
end

function checkCapacityRange(silicon, tag)
% Assert that the exported capacity did not leave the normalized range.
    capacity = silicon.normalizedCapacity;
    assert(min(capacity) >= 0 && max(capacity) <= 1, ...
        '%s: capacity left [0,1] (%.17g .. %.17g).', ...
        tag, min(capacity), max(capacity));
end
