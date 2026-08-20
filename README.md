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
3. **`Run_Cohort_Analysis.m`** - run this AFTER `Run_Attractor_Analysis.m`
   has been run for every animal/recording. Pools all
   `<recording_ID>_results.mat` files into one table for cohort-level
   statistics - see "Cohort-level analysis" below.

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
- **`CohortAnalysis/`** - the function that `Run_Cohort_Analysis.m` calls.
  Sits outside `Functions/`, like `Validation/`, because nothing in the
  single-recording pipeline calls it directly.
  - `fn_cohortAnalysis.m`
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

## Cohort-level analysis (across recordings)

Once you have run `Run_Attractor_Analysis.m` for every animal, run
**`Run_Cohort_Analysis.m`** to pool the results. Edit `recording_IDs` near the top first to list exactly which recordings
belong to be analysed in a cohort. This is deliberate rather than auto-pooling
everything in that folder, thus it keeps a record of exactly what went into
each cohort table. If you'd rather auto-discover every `*_results.mat` file instead, 
the script has a commented-out line for that.

It calls `fn_cohortAnalysis(RESULTS_DIR, recording_IDs)`, shows you
`cohort.summaryTable`, and saves two summary figures to
`RESULTS_DIR/CohortFigures/` (PR across epochs, one line per animal;
chance-corrected alignment, one bar per animal) via
`Functions/Plotting/fn_plotCohortPR.m` and
`fn_plotCohortAlignment.m`. See `help fn_cohortAnalysis` for what columns
you get and the menu of candidate tests.

This loads each recording's `.mat` file and pulls a fixed set of scalar metrics into one
row-per-recording table (`cohort.summaryTable`). 

### The two epoch definitions

Every metric below is computed twice, under two different ways of
deciding where "Evoked" ends and "Recovery" begins - the column name
tells you which one you're looking at:

- **Fixed epochs** (no suffix, e.g. `PR_norm_Evoked`, `align_ER`) - epoch
  boundaries come from protocol timing alone: fixed delays after the
  stimulus markers (`t_P9`, `t_C2`), the same for every recording of a
  given protocol duration. Simple and consistent across animals, but
  assumes every animal's population actually settles into
  Evoked/Recovery dynamics on the same fixed clock.
- **Attractor-defined / RR-defined epochs** (`_rr` suffix, e.g.
  `PR_norm_rr_AttractorEvoked`, `align_ER_rr`) - epoch boundaries come
  from the data itself: the point where recurrence density actually
  crosses `cfg.ONSET_RATIO`, i.e. when the population's own trajectory
  detectably "locks onto" a recurring state (`t_attractor_onset`) and
  when it detectably leaves it again (`t_attractor_return`). Adapts to
  each animal, but only as trustworthy as `onset_detected`/
  `return_detected` being `true` for that recording (see below).

Comparing the fixed vs `_rr` version of the same metric is itself a
useful check: if a result holds up under both, it's less likely to be an
artifact of how the epoch boundary was drawn.

### What each column means

**Subspace alignment** (Evoked subspace vs Recovery subspace)
| Column | Meaning |
|---|---|
| `align_ER` | Raw alignment score (0-1) between the Evoked and Recovery population subspaces (Elsayed & Cunningham method), fixed epochs. Higher = the two epochs occupy more similar directions in neural state space. |
| `chance_lvl` | The alignment score expected by chance for two *randomly oriented* subspaces of the same dimensionality (`= K/N`, i.e. `nDims_align / nN`). Not zero - random subspaces still overlap somewhat, more so in low-dimensional recordings. |
| `align_ER_corrected` | `align_ER` rescaled against `chance_lvl`: `(align_ER - chance_lvl) / (1 - chance_lvl)`. 0 = exactly chance, 1 = perfect alignment. **Use this one, not `align_ER`, when comparing across animals** - different animals can have different `chance_lvl` (different N or K), so raw `align_ER` values aren't directly comparable to each other, but the corrected value is. |
| `nDims_align` | K, the number of dimensions the alignment was computed in (chosen automatically from the Evoked epoch; same K used for both fixed and `_rr` alignment). |
| `align_ER_rr` | Same idea as `align_ER` (**raw**, not chance-corrected), but computed between the attractor-defined Evoked/Recovery epochs instead of the fixed-time ones. There is no chance-corrected counterpart for the `_rr` epochs in the current pipeline - if you need one, compute `(align_ER_rr - chance_lvl) / (1 - chance_lvl)` yourself before comparing `align_ER_rr` across animals, for the same reason `align_ER_corrected` exists for the fixed-epoch version. |

**Participation ratio** (dimensionality of population activity within an epoch, normalized by neuron count so it's comparable across recordings with different N)
| Column | Meaning |
|---|---|
| `PR_norm_Baseline` / `_Evoked` / `_Recovery` | Normalized participation ratio in each fixed epoch. Higher = activity is more spread across many dimensions (higher-dimensional); lower = more collapsed onto a few dominant directions. |
| `PR_norm_rr_Baseline` / `_AttractorEvoked` / `_AttractorRecovery` | Same metric, computed in the attractor-defined epochs instead. |

**Recurrence rate** (how often the population's trajectory revisits nearby states - a signature of being "on an attractor")
| Column | Meaning |
|---|---|
| `RQA_ev_RR` | Self-recurrence density within the Evoked epoch: how often Evoked-epoch states recur *within the Evoked epoch itself*, using an epsilon (distance threshold) calibrated on that epoch's own trajectory. |
| `RQA_re_RR` | Same idea, self-recurrence within the Recovery epoch, using Recovery's own calibrated epsilon. |
| `cross_recur_density` | Cross-recurrence: how often Evoked-epoch states recur *within the Recovery epoch*, using one shared epsilon (calibrated on Evoked). This is the one that speaks to whether Recovery actually returns to the same states visited during Evoked, rather than just being locally self-consistent. |
| `epsilon_rr` / `epsilon_re` | The calibrated distance thresholds themselves (Evoked-epoch and Recovery-epoch respectively) - useful for sanity-checking whether two recordings used comparable thresholds, less useful as an analysis variable in its own right. |

**Attractor onset/return timing** (attractor-defined epochs only - these are what define the `_rr` epoch boundaries above)
| Column | Meaning |
|---|---|
| `t_attractor_onset` | Time (s) the population first detectably locked onto a recurring state after the stimulus. |
| `t_attractor_return` | Time (s) it detectably left that state again. |
| `onset_detected` / `return_detected` | `true` if that lock-on/release was genuinely detected; `false` means detection failed and `t_attractor_onset`/`t_attractor_return` is a fallback value, **not a real timestamp** - check these flags before using the timing columns, and before trusting any `_rr` column for that recording, since the `_rr` epochs are built from these same timestamps. |

`recording_ID` and `protocol` are carried through as-is for reference/filtering.

By default it also saves the pooled table to
`<RESULTS_DIR>/CohortAnalysis_Summary.mat` (pass `'SaveTable', false` to
skip this). It is **not** called automatically by
`Run_Attractor_Analysis.m` or anywhere else - run it yourself once every
recording is done. **This is not a statistics tool**.

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
