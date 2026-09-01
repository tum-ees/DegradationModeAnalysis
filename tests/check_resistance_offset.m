function ok = check_resistance_offset()
%> Author: Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2026-08-28
%
%CHECK_RESISTANCE_OFFSET Does the resistance offset do what it claims to do?
%
% The offset is only defensible if a handful of things hold. It has to be
% invisible to the objectives that differentiate the curve, it has to move
% the OCV level by exactly the ohmic drop and by nothing else, switching it
% off has to give back the model that existed before the slot was added,
% the capacity-normalised limit has to hand cells of different size the
% same room in volts, and the wiring through dma_core has to deliver all of
% that end to end: a drop injected into a synthetic measurement has to come
% back out of the fit, in both directions. Each check below is written so
% that it fails loudly if one of those stops being true.
%
% Lives in tests\ and locates the repository from its own position, so it
% runs from any folder and leaves the search path as it found it.
%
% OUTPUTS:
%   ok  - true when every check passed

thisDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(thisDir);
oldPath = path();
restorePath = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(repoRoot, fullfile(repoRoot, 'helper_functions'));

fprintf('\n=== resistance offset checks ===\n\n');
problems = strings(0, 1);

% ---------------------------------------------------------------- layout
% Every legacy length has one fixed meaning, so each is asserted as a full
% vector rather than by length alone: a padding that lands values in the
% wrong slots would keep the right length and change the physics. A length
% that no longer fits the layout has to be refused rather than silently
% reinterpreted.
layoutCases = { ...
    [1 2 3 4],           [1 2 3 4 0 0 0 0 0]; ...
    [1 2 3 4 5],         [1 2 3 4 5 0 0 0 0]; ...
    [1 2 3 4 5 6],       [1 2 3 4 0 0 5 6 0]; ...
    [1 2 3 4 5 6 7],     [1 2 3 4 5 0 6 7 0]; ...
    [1 2 3 4 5 6 7 8],   [1 2 3 4 5 6 7 8 0]; ...
    [1 2 3 4 5 6 7 8 9], [1 2 3 4 5 6 7 8 9]};
for k = 1:size(layoutCases, 1)
    got = expand_params_full(layoutCases{k, 1});
    if ~isequal(got(:).', layoutCases{k, 2})
        problems(end+1, 1) = sprintf( ...
            'expand_params_full(%d elements) gave %s, expected %s', ...
            numel(layoutCases{k, 1}), mat2str(got(:).'), ...
            mat2str(layoutCases{k, 2})); %#ok<AGROW>
    end
end
try
    expand_params_full(ones(1, 10));
    problems(end+1, 1) = "expand_params_full accepted 10 elements";
catch
    % expected
end
fprintf('layout      : every legacy length maps into its exact slots\n');

% ----------------------------------------------------------------- sign
% Charge lies above the OCV and discharge below it, so the same resistance
% has to come out with opposite signs in the two directions.
si = struct('rOffsetCurrent', 0.12, 'rOffsetSign', 1);
sj = struct('rOffsetCurrent', 0.12, 'rOffsetSign', -1);
uCharge    = resistance_offset_voltage(0.05, si);
uDischarge = resistance_offset_voltage(0.05, sj);
if abs(uCharge - 0.006) > 1e-12 || abs(uDischarge + 0.006) > 1e-12
    problems(end+1, 1) = sprintf( ...
        'ohmic drop is %.6f V in charge and %.6f V in discharge, expected +-0.006', ...
        uCharge, uDischarge);
end
% The reference curves carry their own drop, so a negative net resistance is
% a legitimate outcome and has to be mirrored exactly, not clipped.
if resistance_offset_voltage(-0.05, si) ~= -uCharge
    problems(end+1, 1) = ...
        "a negative resistance is not the exact mirror of a positive one";
end
% A solver bundle without a current may only pass while nothing is asked of
% it. The moment a non-zero resistance is requested it has to refuse.
if resistance_offset_voltage(0, struct()) ~= 0
    problems(end+1, 1) = "a zero offset did not return zero volts";
end
try
    resistance_offset_voltage(0.05, struct());
    problems(end+1, 1) = "a non-zero offset was accepted without a current";
catch
    % expected
end
fprintf('sign        : +%.1f mV in charge, %.1f mV in discharge at 0.12 A\n', ...
    1e3 * uCharge, 1e3 * uDischarge);

% ----------------------------------------------------------------- limit
% The limit is given per capacity so that the reachable voltage is the limit
% times the C-rate. Three cells of very different size, run at the same rate,
% must therefore end up with the same room in volts.
limit = 0.25;         % Ohm*Ah
cRate = 0.04;         % 1/h, the rate of every check-up in this study
rooms = zeros(1, 3);
capacities = [2.939, 1.868, 1.215];
for k = 1:3
    st = struct('resistanceOffsetLimit', limit, ...
        'allowNegativeResistanceOffset', false);
    b = resistance_offset_bounds(st, capacities(k));
    if b(1) ~= 0
        problems(end+1, 1) = sprintf( ...
            'with negatives off the lower bound is %g, expected 0', b(1)); %#ok<AGROW>
    end
    rooms(k) = b(2) * cRate * capacities(k);      % volts
end
if max(abs(rooms - rooms(1))) > 1e-12
    problems(end+1, 1) = sprintf( ...
        'the same limit gives %s V of room, expected one value', ...
        mat2str(rooms, 4));
end
% Switching negatives on has to open the lower half and nothing else.
stNeg = struct('resistanceOffsetLimit', limit, ...
    'allowNegativeResistanceOffset', true);
bNeg = resistance_offset_bounds(stNeg, capacities(1));
bPos = resistance_offset_bounds(struct('resistanceOffsetLimit', limit), ...
    capacities(1));
if bNeg(2) ~= bPos(2) || bNeg(1) ~= -bPos(2)
    problems(end+1, 1) = sprintf( ...
        'allowing negatives gave %s, expected the mirror of %s', ...
        mat2str(bNeg, 4), mat2str(bPos, 4));
end
fprintf(['limit       : 0.25 Ohm*Ah is %.2f mV for every cell from %.2f to ' ...
    '%.2f Ah\n'], 1e3 * rooms(1), min(capacities), max(capacities));

% ------------------------------------------------------------- objectives
% Build a small stand-in cell so the fit functions can be called directly.
solverInput = synthetic_solver_input();
q0 = solverInput.q0;
roiOCV = [0, 1];
X0 = [1.0, 0.0, 1.0, 0.0, 0, 0, 0, 0, 0];        % no offset
rTest = 0.05;                                     % 6 mV at 0.12 A
X1 = X0; X1(9) = rTest;

ocv0 = fit_ocv(X0, solverInput, roiOCV(1), roiOCV(2));
ocv1 = fit_ocv(X1, solverInput, roiOCV(1), roiOCV(2));
dva0 = fit_dva(X0, solverInput, q0, 0.1, 0.9);
dva1 = fit_dva(X1, solverInput, q0, 0.1, 0.9);
ica0 = fit_ica(X0, solverInput, q0, 0.1, 0.9);
ica1 = fit_ica(X1, solverInput, q0, 0.1, 0.9);

if dva0 ~= dva1
    problems(end+1, 1) = sprintf( ...
        'the DVA objective moved by %.3e when the offset was switched on', ...
        dva1 - dva0);
end
if ica0 ~= ica1
    problems(end+1, 1) = sprintf( ...
        'the ICA objective moved by %.3e when the offset was switched on', ...
        ica1 - ica0);
end
fprintf('derivatives : DVA and ICA are bit for bit unchanged by the offset\n');

% The OCV term has to move by exactly what a constant shift does to a mean
% squared error, and by nothing else: E[(d + u)^2] = E[d^2] + 2*u*E[d] + u^2.
% The tolerance is relative, because the two sides differ only by rounding
% order and an absolute threshold would break on other data or another BLAS.
uDrop = resistance_offset_voltage(rTest, solverInput);
resid = ocv_residual(X0, solverInput);
expected = ocv0 + 2 * uDrop * mean(resid) + uDrop^2;
if abs(ocv1 - expected) > 1e-12 * max(abs(ocv0), uDrop^2)
    problems(end+1, 1) = sprintf( ...
        'the OCV term is %.6e but a pure level shift predicts %.6e', ...
        ocv1, expected);
end
fprintf('level       : OCV term follows a pure level shift to %.1e\n', ...
    abs(ocv1 - expected));

% ---------------------------------------------------------------- off
% With the offset off the model has to be the one that existed before the
% slot was added, which is the length-8 vector.
if fit_ocv(X0(1:8), solverInput, roiOCV(1), roiOCV(2)) ~= ocv0
    problems(end+1, 1) = ...
        "an 8-element vector and a 9-element vector with a zero offset disagree";
end
fprintf('default off : an 8-slot vector reproduces the 9-slot result\n');

% ---------------------------------------------------------------- wiring
% The checks above exercise the helpers and fit_ocv in isolation. This one
% drives dma_core itself: a synthetic cell whose measurement carries a
% known ohmic drop, fitted with fmincon so the outcome is deterministic. It
% fails if the offset never reaches the solver (the mask), is not applied
% to the reconstruction (the lift), or is mapped to the wrong sign of the
% direction.
[hfc, st, ocvTrue, uTrue] = synthetic_core_input();

% Off, and without a direction field: a settings struct that predates the
% slot has to fit exactly as before and return a zero ninth slot.
stOff = st;
stOff.allowResistanceOffset = false;
hfc.fullCellU = ocvTrue;
[resid0, params0, capa0] = run_core(hfc, stOff);
if numel(params0) ~= 9 || params0(9) ~= 0
    problems(end+1, 1) = sprintf( ...
        'with the offset off, params came back as %s, expected slot 9 == 0', ...
        mat2str(params0(:).', 3));
end
if resid0 > 5e-5
    problems(end+1, 1) = sprintf( ...
        'the off-fit misses its own synthetic cell by %.2f mV', 1e3 * resid0);
end
if capa0 ~= hfc.capaAct
    problems(end+1, 1) = ...
        "dma_core did not return the capacity it was given";
end

% On, in charge: the measurement sits above the OCV by the drop, and the
% fitted offset has to give exactly that drop back.
stChg = st;
stChg.direction = 'charge';
hfc.fullCellU = ocvTrue + uTrue;
[residC, paramsC] = run_core(hfc, stChg);
dropC = paramsC(9) * stChg.pOCVCurrent;
if abs(dropC - uTrue) > 5e-5 || residC > 5e-5
    problems(end+1, 1) = sprintf( ...
        ['charge: %.2f mV were injected, %.2f mV came back, worst ' ...
        'residual %.2f mV'], 1e3 * uTrue, 1e3 * dropC, 1e3 * residC);
end

% On, in discharge: the same positive resistance now sits below the OCV.
stDis = st;
stDis.direction = 'discharge';
hfc.fullCellU = ocvTrue - uTrue;
[residD, paramsD] = run_core(hfc, stDis);
dropD = paramsD(9) * stDis.pOCVCurrent;
if abs(dropD - uTrue) > 5e-5 || residD > 5e-5
    problems(end+1, 1) = sprintf( ...
        ['discharge: %.2f mV were injected, %.2f mV came back, worst ' ...
        'residual %.2f mV'], 1e3 * uTrue, 1e3 * dropD, 1e3 * residD);
end
fprintf(['wiring      : dma_core returns an injected %.1f mV drop in both ' ...
    'directions\n'], 1e3 * uTrue);

% ---------------------------------------------------------------- verdict
fprintf('\n');
if isempty(problems)
    ok = true;
    fprintf('all checks passed\n\n');
else
    ok = false;
    fprintf(2, '%d problem(s):\n', numel(problems));
    fprintf(2, '  %s\n', problems);
    fprintf('\n');
end
end


function [worst, params, capaActOut] = run_core(hfc, st)
%RUN_CORE One dma_core fit, reduced to what the wiring check needs.
[~, reconSOC, uModel, ~, ~, ~, ~, ~, ~, ~, ~, capaActOut, params] = ...
    dma_core(hfc, st, [], [], [], [], [], [], [], false, false, false);
uOnMeas = interp1(reconSOC(:), uModel(:), hfc.fullCellSOC(:), 'linear', 'extrap');
worst = max(abs(uOnMeas - hfc.fullCellU(:)));
end


function [hfc, st, ocvTrue, uTrue] = synthetic_core_input()
%SYNTHETIC_CORE_INPUT A stand-in cell for driving dma_core end to end.
% The measurement is built from the very electrode curves dma_core is
% handed, so a perfect fit exists and any residual is wiring, not
% modelling. The electrode supports cover the full cell window under the
% true stretching, so no extrapolated zeros enter the comparison.
n = 400;
x = linspace(0, 1, n).';
anodeU   = 0.28 * exp(-12 * x) + 0.09;             % falls with lithiation
cathodeU = 3.00 + 0.20 * x + 0.25 * x.^8;          % rises towards the knee

trueAlphaAn  = 1.25;  trueBetaAn  = -0.05;
trueAlphaCat = 1.05;  trueBetaCat = -0.03;

q = linspace(0, 1, n).';
anodePot = interp1(trueAlphaAn  * x + trueBetaAn,  anodeU,   q, 'linear', 0);
cathPot  = interp1(trueAlphaCat * x + trueBetaCat, cathodeU, q, 'linear', 0);
ocvTrue = cathPot - anodePot;

current = 0.12;               % A
uTrue = 0.05 * current;       % 50 mOhm -> 6 mV

hfc = struct();
hfc.fullCellSOC          = q;
hfc.fullCellU            = ocvTrue;    % overwritten per case
hfc.anodeSOCSingle       = x;
hfc.anodeUSingle         = anodeU;
hfc.cathodeSOCSingle     = x;
hfc.cathodeUSingle       = cathodeU;
hfc.normCathodeSOC       = x;
hfc.normCathodeU         = cathodeU;
hfc.commonVoltageAnode   = [];
hfc.commonVoltageCathode = [];
hfc.qAnodeBlend1Interp   = [];
hfc.qAnodeBlend2Interp   = [];
hfc.qCathodeBlend1Interp = [];
hfc.qCathodeBlend2Interp = [];
hfc.q0                   = 1;
hfc.capaAct              = 2.0;

st = struct();
st.algorithm             = 'fmincon';
st.dataLength            = n;
st.smoothingPoints       = 5;
st.weightOCV             = 1;
st.weightDVA             = 0;
st.weightICA             = 0;
st.roiOCVMin             = 0;
st.roiOCVMax             = 1;
st.roiDVAMin             = 0.1;
st.roiDVAMax             = 0.9;
st.roiICAMin             = 0.1;
st.roiICAMax             = 0.9;
% Both electrodes are pinned at the true stretching, so the offset is the
% only degree of freedom. That is the point: the check isolates the wiring
% (mask, lift, sign) from the difficulty of a four-parameter fit, and the
% tolerances can be sharp because a perfect solution exists.
st.lowerBoundaries       = [1.25, -0.05, 1.05, -0.03];
st.upperBoundaries       = [1.25, -0.05, 1.05, -0.03];
st.gammaAnBlend2UpperBound = 1;
st.gammaCaBlend2UpperBound = 1;
st.useAnodeBlend         = false;
st.useCathodeBlend       = false;
st.maxInhomogeneity      = 0.3;
st.maxInhomogeneityDelta = 0.1;
st.inhomAnodeOffset      = 0;
st.inhomCathodeOffset    = 0;
st.maxAnodeGain    = 1;  st.maxAnodeLoss    = 1;
st.maxCathodeGain  = 1;  st.maxCathodeLoss  = 1;
st.maxAnBlend1Gain = 1;  st.maxAnBlend1Loss = 1;
st.maxAnBlend2Gain = 1;  st.maxAnBlend2Loss = 1;
st.allowResistanceOffset = true;
st.allowNegativeResistanceOffset = false;
st.resistanceOffsetLimit = 0.25;   % rMax = 0.125 Ohm, the true 0.05 inside
st.pOCVCurrent           = current;
% st.direction is set per case; the off-case leaves it out on purpose.
end


function solverInput = synthetic_solver_input()
%SYNTHETIC_SOLVER_INPUT A stand-in cell, so the checks need no measurement.
% The shapes only have to be monotone and smooth enough for the DVA to be
% finite. Nothing here is meant to resemble a real electrode.
n = 400;
soc = linspace(0, 1, n).';
cathodeU = 3.45 - 0.15 * soc - 0.25 * soc.^8;      % flat with an end drop
anodeU   = 0.28 * exp(-12 * soc) + 0.09;           % graphite-like decay

solverInput = struct();
solverInput.qCell                = soc.';
solverInput.ocvCell              = (cathodeU - anodeU).' + 0.004;
solverInput.anodeSOCSingle       = soc;
solverInput.anodeUSingle         = anodeU;
solverInput.cathodeSOCSingle     = soc;
solverInput.cathodeUSingle       = cathodeU;
solverInput.normCathodeSOC       = soc;
solverInput.normCathodeU         = cathodeU;
solverInput.useAnodeBlend        = false;
solverInput.useCathodeBlend      = false;
solverInput.qAnodeBlend1Interp   = [];
solverInput.qAnodeBlend2Interp   = [];
solverInput.qCathodeBlend1Interp = [];
solverInput.qCathodeBlend2Interp = [];
solverInput.inhomAnodeOffset     = 0;
solverInput.inhomCathodeOffset   = 0;
solverInput.q0                   = 1;
solverInput.rOffsetCurrent       = 0.12;
solverInput.rOffsetSign          = 1;
end


function resid = ocv_residual(X, solverInput)
%OCV_RESIDUAL The signed model minus measurement, the way fit_ocv squares it.
X = expand_params_full(X);
Q = solverInput.qCell(:);
anodePot = interp1(X(1) * solverInput.anodeSOCSingle + X(2), ...
    solverInput.anodeUSingle, Q, 'linear', 0);
cathPot = interp1(X(3) * solverInput.cathodeSOCSingle + X(4), ...
    solverInput.cathodeUSingle, Q, 'linear', 0);
resid = (cathPot - anodePot) - solverInput.ocvCell(:);
end
