# Aplysia Attractor Perturbation Analysis Toolbox

MATLAB toolbox for analysing population-level attractor dynamics after a
motor perturbation in *Aplysia* recordings.

A single user-edited settings script, a folder of small reusable
functions (one function per file, each documenting the inputs/outputs and
which window of the time series it operates on), and a top-level pipeline
script that only calls those functions in sequence.

## Top-level scripts

1. **`Config_UserSettings.m`** - the *only* place parameters should be
   edited (file paths, sampling rate, epoch timing, window sizes,
   thresholds). Every other script/function reads its parameters from the
   `cfg` struct this script defines; nothing is hard-coded elsewhere.
2. **`Run_Attractor_Analysis.m`** - the main pipeline. Loads `cfg`, then
   calls the analysis and plotting functions in sequence (load/clean/smooth
   -> PCA -> windowed PCA & participation ratio -> subspace alignment ->
   recurrence-density analysis -> RR-based attractor epoch detection),
   saving all figures to `cfg.FIGURES_DIR` and a results `.mat` file to
   `cfg.RESULTS_DIR`.

## Folders

- **`Functions/`** - one function per file, grouped by analysis stage:
  - `DataIO/` - loading the peaks matrix, deriving epoch timing from the protocol
  - `Preprocessing/` - neuron quality filtering, kernel-width estimation, Gaussian smoothing
  - `Dimensionality/` - PCA, windowed PCA, participation ratio (fixed & sliding)
  - `SubspaceAlignment/` - the alignment metric, K selection, epoch/K-sweep/sliding alignment
  - `Recurrence/` - recurrence density, cross-recurrence density, epsilon calibration, attractor epoch detection
  - `Plotting/` - one function per figure; each takes only the data it needs to draw and a save path
  - `Utils/` - small generic helpers (pass/fail check, ternary, run-length detection, synthetic subspace construction)
- **`Validation/`** - synthetic validation ("unit test") scripts for the two
  key metrics, checked against known theoretical values. These are
  development/QA scripts, not part of the analysis pipeline itself (they
  use no real recording data):
  - `Validate_SubspaceAlignment.m`
  - `Validate_RecurrenceDensity.m`
- **`Results/`**, **`Figures/`** - default output locations (see `Config_UserSettings.m`).

## How do I use this with my own data?

1. Make sure your recording is saved as a `.mat` file containing a
   variable `peaks` (nNeurons x nFrames, logical or numeric event matrix).
2. Open `Config_UserSettings.m` and set `DATA_FILE`, `RESULTS_DIR`,
   `FIGURES_DIR`, `recording_ID`, `protocol` (`'12min'` or `'20min'`), and
   `fs` to match your recording. Adjust epoch timing, window sizes, and
   thresholds if needed (all fields are documented in that script).
3. Run `Run_Attractor_Analysis.m`.
4. To add a new protocol duration, add a `case` to `fn_getEpochTiming.m`.
5. To add a new figure, write a new function in `Functions/Plotting/` that
   takes the data it needs plus a save path, and call it from
   `Run_Attractor_Analysis.m` - keep plotting code out of the analysis
   functions, and analysis code out of the plotting functions.

## Things to be aware of

- All functions assume `spike_conv` (Gaussian-smoothed activity) is
  `neurons x frames`, single precision, and that PCA bases (`V`, `mu`) are
  computed once and reused everywhere a shared reference subspace is
  needed (e.g. recurrence-density trajectories use the same 3-D basis
  throughout Section 5, calibrated from the Evoked epoch onward - see
  `Run_Attractor_Analysis.m`).
- Self-recurrence density for an epoch is always calibrated with that
  epoch's *own* epsilon (`fn_safeEpsilon`); cross-recurrence density
  between two epochs intentionally uses one shared epsilon (from the
  reference epoch) - see the comments in `Run_Attractor_Analysis.m`
  Section 5 for why these differ.
- `fn_recurrenceDensity` and `fn_crossRecurrenceDensity` are both
  O(N^2) in the number of trajectory points; always downsample long
  trajectories first (see `fn_getEpochTrajectory`, and `cfg.MAX_WIN_PTS` /
  `cfg.MAX_FULL_PTS`).

## Problems?

Check the docstring at the top of the relevant function first - every
function in `Functions/` documents its inputs, outputs, and (where
relevant) which time window it expects. If something still looks wrong,
re-run the matching script in `Validation/` to confirm the underlying
metric is behaving as expected before suspecting the pipeline wiring.
