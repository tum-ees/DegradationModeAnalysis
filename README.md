<div align="center">
  <h1>DegradationModeAnalysis</h1>

  Degradation mode analysis framework and tool to calculate silicon OCPs

  <p>
    <a href="https://www.mathworks.com/help/matlab/"><img src="https://img.shields.io/badge/Platform-MATLAB-blue.svg" alt="MATLAB"></a>
    <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="MIT License"></a>
    <a href="https://doi.org/10.1016/j.jpowsour.2026.239418"><img src="https://img.shields.io/badge/Paper-J.%20Power%20Sources-green.svg" alt="Journal of Power Sources"></a>
    <a href="https://doi.org/10.1039/D5EB00221D"><img src="https://img.shields.io/badge/Paper-EES%20Batteries-green.svg" alt="EES Batteries"></a>
    <a href="https://doi.org/10.5281/zenodo.17591931"><img src="https://img.shields.io/badge/DOI-10.5281%2Fzenodo.17591931-blue.svg" alt="DOI"></a>
  </p>

  <img src="doc/OCP_shift_over_SOC.gif" width="550">
</div>

<h2>Overview</h2>

This framework enables degradation mode analysis for lithium-ion and sodium-ion batteries.
For blend anodes, the silicon OCP can either be selected from literature data or reconstructed
from a measured blend OCP with the included silicon-OCP generation tool.

<br>

<div align="center"><img src="doc/flowChart2.jpg?raw=1" width="600" alt="Flow chart" /></div>

<h2>Installation</h2>
Clone the repository by running

```bash
git clone git@github.com:tum-ees/degradation-mode-analysis.git
```

<h2>Usage</h2>
Set the parameters and run the scripts in the MATLAB environment.

<h4>Silicon OCP generation (optional)</h4>

Optional silicon OCP generation for Si-Gr blends:
You can generate a cell-specific silicon OCP from a measured blend anode OCP
and a graphite OCP using `generate_si_ocp.m`.
This follows the algebraic reconstruction

<p><em>Q</em><sub>blend</sub> = <span style="font-style:italic;">gamma</span><sub>Si</sub>&middot;<em>Q</em><sub>Si</sub> + (1&minus;<span style="font-style:italic;">gamma</span><sub>Si</sub>)&middot;<em>Q</em><sub>Gr</sub></p>

and is robust even when <span style="font-style:italic;">gamma</span><sub>Si</sub>
is only roughly estimated. Filtering of the generated curve is available if you want
a strictly monotonic OCP. Use this only if you have a Si-Gr blend and want to avoid
mismatches from literature silicon OCPs.

For downstream interpolation routines that require unique, strictly increasing
capacity coordinates, enable the interpolation-ready output together with PAV:

```matlab
silicon = generate_si_ocp( ...
    'blendPath', blendPath, ...
    'savePath', savePath, ...
    'graphitePath', graphitePath, ...
    'lithDirection', 'lithiation', ...
    'gammaSi', 0.10, ...
    'pavOutput', true, ...
    'interpolationReadyOutput', true);
```

This optional mode uses shape-preserving PCHIP interpolation on a uniform
1001-point normalized-capacity grid and is disabled by default.

The default output changed in version 2.1.0. The aligned voltage window is
now sampled with exact boundary points, instead of dropping the first grid
sample and filling out-of-support samples with zero. This is an intentional
bugfix, so existing workflows do not reproduce their earlier output
bit-for-bit. Measured on the bundled P45B 0C03 example data with a silicon
fraction of 0.245 and the Rehm2026 graphite reference, the reconstructed
silicon curve is exported with two more samples (2238 to 2240 in lithiation
and 2333 to 2335 in delithiation) and shifts by up to 0.003 (lithiation) and
0.013 (delithiation) in normalized capacity.

<h4>Degradation mode analysis</h4>
Run the DMA by calling `[data, s] = main_dma(userSettingsOutside)`.
Set all defaults in `main_dma.m`. You can overwrite any subset of settings
from outside by passing a struct into `main_dma`; only the provided fields change.

* Data handling: all pOCV curves need to be stored in a table (see the minimal working example as reference).

* Resampling: use `s.dataLength` for resampling in SOC space and `s.smoothingPoints` for LOWESS smoothing of input curves.

* Cost function: combine OCV, DVA, and ICA via `s.weightOCV`, `s.weightDVA`, and `s.weightICA`. Focus the fit with `s.roiOCVMin` / `s.roiOCVMax`, `s.roiDVAMin` / `s.roiDVAMax`, and `s.roiICAMin` / `s.roiICAMax`.

* Solver and run control: choose `s.algorithm` such as `ga`, `particleswarm`, `patternsearch`, `GlobalSearch`, `fmincon`, or `lsqnonlin`. For non-deterministic algorithms use `s.rmseThreshold`, `s.reqAccepted`, and `s.maxTriesOverall`.

* Direction of pOCV: set `s.direction` to `'charge'` or `'discharge'`.

* Anode blend option: enable with `s.useAnodeBlend` and set `s.gammaAnBlend2Init` and `s.gammaAnBlend2UpperBound`.

* Cathode blend option: enable with `s.useCathodeBlend` and set `s.gammaCaBlend2Init` and `s.gammaCaBlend2UpperBound`; supply a second cathode OCP.

* Inhomogeneity: toggle `s.allowAnodeInhomogeneity` and `s.allowCathodeInhomogeneity`; limit with `s.maxInhomogeneity` and `s.maxInhomogeneityDelta`. Use `s.inhomAnodeOffset` and `s.inhomCathodeOffset` to define the fraction of maximum inhomogeneity already present at SOC = 0.

* Constraints and order: bound changes with `s.maxCathodeGain`, `s.maxAnodeGain`, `s.maxAnBlend1Gain`, `s.maxAnBlend2Gain`, `s.maxCathodeLoss`, `s.maxAnodeLoss`, `s.maxAnBlend1Loss`, and `s.maxAnBlend2Loss`. Control fitting order via the sort order of `s.nCUs`.

<h2>Content</h2>
Detailed documentation of the modules can be found below.
<br><br>

<details>
<summary><h4>Silicon OCP generation</h4></summary>

* Folder `generate_si_ocp`: all necessary scripts to generate the silicon OCP

* `generate_si_ocp.m`: entry script to perform the calculation in GUI or script mode

</details>

<details>
<summary><h4>Degradation mode analysis</h4></summary>

* `main_dma.m`: main entry point; all default settings are defined there

* `dma_core.m`: core routine that handles the main fitting workflow

* Folder `helper_functions`: required helper functions for import, preprocessing, fitting, plotting, and saving

</details>

<details>
<summary><h4>Input data</h4></summary>

* Folder `input_data`: literature OCPs and example data

* Subfolders: `graphite`, `silicon`, `LFP`, `NCA`, `NMC`, `test_data`

These OCPs originate from published sources. Add proper citations if you use them.
Check licenses and attribution requirements before redistribution.

</details>

<h2>Acknowledgments</h2>
We thank Johannes Natterer for providing a cyclic-aged P45B data set for testing the framework.
We also thank Maximilian Leitenstern for support in migration to GitHub.

<h2>Minimal workable example</h2>
In its current form, `main_dma.m` serves as the minimal working example.
It performs an anode-blend fitting for a cyclic-aged Molicel P45B cell.
The pOCV curves are stored in the table
`.\input_data\test_data\P45B_serial23_aging_data_table.mat`.
The OCP curves for the minimal working example are described in the accompanying publication.

<h2>Developers</h2>

* [Mathias Rehm](mailto:mathias.rehm@tum.de), Chair of Electrical Energy Storage Technology, School of Engineering and Design, Technical University of Munich, 80333 Munich, Germany

* [Josef Eizenhammer](mailto:josef.eizenhammer@tum.de), Chair of Electrical Energy Storage Technology, School of Engineering and Design, Technical University of Munich, 80333 Munich, Germany

* Moritz Guenthner (student research project)

* Can Korkmaz (student research project)

<h2>Citation</h2>

This framework is published alongside an open-source paper where the full method and code are described.
If you use this repository in any publication, please cite:

> M. Rehm et al., "How to determine the degradation modes of lithium-ion batteries with silicon-graphite blend electrodes,"
> *Journal of Power Sources*, 2026, DOI: [10.1016/j.jpowsour.2026.239418](https://doi.org/10.1016/j.jpowsour.2026.239418)

The framework is also applied and validated on commercial sodium-ion batteries in the following publication.
We appreciate citations of this work as well if your work involves sodium-ion cells:

> M. Rehm et al., "Aging of commercial sodium-ion batteries with layered oxides: how to measure and analyze it?,"
> *EES Batteries*, 2026, DOI: [10.1039/D5EB00221D](https://doi.org/10.1039/D5EB00221D)

To cite a specific version of the code for reproducibility, use the version-specific DOI from
[Zenodo](https://doi.org/10.5281/zenodo.17591931).
