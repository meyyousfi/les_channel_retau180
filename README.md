# Wall-resolved LES of turbulent channel flow at Re<sub>τ</sub> = 180 — mesh-convergence study

Part of **Project Hydroshark** (hydrodynamics of sharkskin denticles). This folder contains the
smooth-wall baseline: a wall-resolved Large Eddy Simulation (LES) of a periodic open channel at
Re<sub>τ</sub> = 180, run on three meshes (coarse, medium, fine) in OpenFOAM.

Each mesh starts from the converged solution of the next coarser mesh:

```
synthetic turbulent field ──► case_coarse (64³) ──► case_medium (96³) ──► case_fine (144³)
                                              mapFields            mapFields
```

---

## 1. Repository layout

```
les_channel_retau180/
├── README.md                     this file
├── case_coarse/                  64 × 64 × 64     (starts from synthetic initial field)
├── case_medium/                  96 × 96 × 96     (starts from case_coarse results)
├── case_fine/                    144 × 144 × 144  (starts from case_medium results)
│   ├── 0.orig/                   initial & boundary conditions (U, p, nut)
│   ├── constant/                 transportProperties, turbulenceProperties, fvOptions
│   ├── system/                   blockMeshDict, controlDict, fvSchemes, fvSolution, decomposeParDict
│   ├── a.pre-processing.sh       mesh + initial fields + domain decomposition
│   ├── q.Batch_AWS               SLURM job: runs pimpleFoam in parallel, then reconstructPar
│   ├── b.post-processing.sh      averaged profiles to CSV (ParaView pvpython)
│   ├── paraview_openfoam.py      the averaging script called by b.post-processing.sh
│   ├── w.clean                   resets the case (deletes all results!)
│   └── z.openfoam_output/        CSV output of the post-processing
├── tools/matlab/
│   ├── openfoam_grid_calc.m      computes the blockMesh wall-normal grading for a target y⁺
│   └── turbulent_field_generator.m   builds the synthetic initial U, p, nut fields
└── docs/LES_Case_Setup_Guide.md  detailed explanation of every file + student tasks
```

The three cases are **identical except for the mesh line in `system/blockMeshDict`,
the initial fields in `0.orig/` and the source case in `a.pre-processing.sh`**.
The cases must stay side by side in the same folder, because `case_medium` and
`case_fine` read their initial condition from `../case_coarse` and `../case_medium`.

---

## 2. Requirements

| Software | Version used | Needed for |
|---|---|---|
| OpenFOAM (ESI / openfoam.com) | **v2412** | everything (`blockMesh`, `mapFields`, `pimpleFoam`, …) |
| MPI (OpenMPI) | any recent | parallel run on 8 cores |
| SLURM | optional | `q.Batch_AWS` job script (can also run without it, see §5) |
| ParaView with `pvpython` | 5.x | post-processing only |
| MATLAB | optional | only to regenerate the synthetic initial field |

Other OpenFOAM v2xxx releases should work. The Foundation versions (openfoam.org, e.g. v10/v11)
use different file names (`momentumTransport`, `fvModels`) and will need changes.

---

## 3. Physical setup

| Quantity | Value |
|---|---|
| Domain (L<sub>x</sub> × h × L<sub>z</sub>) | 3π × 1 × π (streamwise × wall-normal × spanwise) |
| Boundaries | `bottomWall` no-slip; `topSurface` symmetry plane (half channel); x and z cyclic |
| Viscosity ν | 0.00556 (≈ 1/180) |
| Forcing | `meanVelocityForce`, bulk velocity U<sub>b</sub> = 15.67 |
| Target Re<sub>τ</sub> = u<sub>τ</sub>h/ν | 180 (u<sub>τ</sub> ≈ 1) |
| Solver | `pimpleFoam`, LES, WALE sub-grid model, Δ = cube-root volume |
| Schemes | `backward` in time, central (`Gauss linear`) convection |
| Time step / end time | Δt = 1.5·10⁻³, t<sub>end</sub> = 400, fields written every 6 time units |
| Statistics | `fieldAverage` of U and p (mean + variance), averaged from t = 100 to 400 in post-processing |

## 4. Meshes

All meshes are graded only in the wall-normal direction (cells clustered at `bottomWall`).

| Case | Cells (x × y × z) | Total cells | y-grading | Δy⁺ first cell | Δx⁺ | Δz⁺ |
|---|---|---|---|---|---|---|
| `case_coarse` | 64 × 64 × 64 | 0.26 M | 23.2 | 0.39 | 26.5 | 8.8 |
| `case_medium` | 96 × 96 × 96 | 0.88 M | 19.7 | 0.30 | 17.7 | 5.9 |
| `case_fine`   | 144 × 144 × 144 | 2.99 M | 19.7 | 0.20 | 11.8 | 3.9 |

The gradings follow from `tools/matlab/openfoam_grid_calc.m` with a first-cell-centre
target of y⁺ = 0.2 (coarse), 0.15 (medium) and 0.1 (fine).

---

## 5. How to run

Load OpenFOAM first (e.g. `source /opt/OpenFOAM/OpenFOAM-v2412/etc/bashrc`) and make sure
`echo $WM_PROJECT_DIR` prints a path. **Run the cases strictly in this order**, each one to
completion, because each mesh starts from the previous one.

### Step 1 — coarse mesh

```bash
cd case_coarse
./a.pre-processing.sh      # 0.orig -> 0, blockMesh, checkMesh, decomposePar
sbatch q.Batch_AWS         # pimpleFoam on 8 cores, then reconstructPar
```

### Step 2 — medium mesh (after the coarse run has finished)

```bash
cd ../case_medium
./a.pre-processing.sh      # blockMesh, mapFields from ../case_coarse (latest time), decomposePar
sbatch q.Batch_AWS
```

### Step 3 — fine mesh (after the medium run has finished)

```bash
cd ../case_fine
./a.pre-processing.sh      # blockMesh, mapFields from ../case_medium (latest time), decomposePar
sbatch q.Batch_AWS
```

`a.pre-processing.sh` stops with an error if the source case has no reconstructed results, so the
fine mesh can't accidentally be started from the wrong fields.

### Without SLURM (workstation)

Replace `sbatch q.Batch_AWS` with:

```bash
mpirun -np 8 pimpleFoam -parallel > log.pimpleFoam 2>&1
reconstructPar -newTimes
```

### Changing the number of cores

Change all three together: `numberOfSubdomains` and `coeffs n (… … …)` in
`system/decomposeParDict` (product of `n` = number of subdomains; keep 1 in y) and
`#SBATCH -n` in `q.Batch_AWS`.

### Before submitting

- Put your own e-mail and account in `q.Batch_AWS`.
- **Set a long enough wall time (`#SBATCH -t`).** If the job hits the limit, SLURM kills it
  before `reconstructPar` runs (see "Restarting" below). Run times on 8 cores for t = 0–400:

  | Case | Cells | Run time on 8 cores | `-t` in `q.Batch_AWS` |
  |---|---|---|---|
  | `case_coarse` | 0.26 M | fits in 8 h | `08:00:00` |
  | `case_medium` | 0.88 M | **more than 24 h** | `48:00:00` |
  | `case_fine`   | 2.99 M | not recorded; expect several times the medium run | set it yourself |

  If your cluster has a shorter maximum wall time, use more cores (see "Changing the
  number of cores") or run in several jobs (see below).
- **Restarting a job that was stopped early:** `controlDict` uses `startFrom latestTime`, so
  resubmitting `q.Batch_AWS` continues from the last time directory written in the
  `processor*` folders (every 6 time units). Do not rerun `a.pre-processing.sh`, because that
  would overwrite the decomposed fields. Once the run reaches t = 400, check that
  `reconstructPar` ran; if not, run `reconstructPar -newTimes` by hand.
- Do a short test first (e.g. set `endTime 1;` in `system/controlDict`) and check the Courant
  number and residuals in `out.o`.

---

## 6. Post-processing

After a run is reconstructed:

```bash
./b.post-processing.sh
# or, if pvpython is not on your PATH:
PV_PYTHON=/path/to/ParaView/bin/pvpython ./b.post-processing.sh
```

This writes `z.openfoam_output/openfoam_averaged.csv`: x–z plane averages as a function of y,
time-averaged over t = 100–400, for U, UMean, UPrime2Mean, p, pMean, pPrime2Mean and nut.
The fields, averaging window and averaging directions are set at the top of
`paraview_openfoam.py`.

---

## 7. Initial condition of the coarse mesh

`case_coarse/0.orig/initialVelocityProfile`, `initialPressureProfile` and
`initialViscosityProfile` hold one value per cell (262 144 entries) and are included in
`U`, `p` and `nut` with `#include`. They were generated by
`tools/matlab/turbulent_field_generator.m`: a law-of-the-wall mean profile with random,
smoothed fluctuations, a pressure field from a Poisson solve, and a mixing-length estimate of
ν<sub>sgs</sub>.

To regenerate them, set `nx = ny = nz = 64` and `d1_plus = 0.2` in the script, remove the
`return` after the call to `openfoam_grid_calc`, and run it. The script relies on helper
functions from the authors' MATLAB library (`Defs`, `PlotSpec`, `plot_1d`, `plot_2d`,
`compute_divergence`, `solve_poisson`, `compute_jacobian`, `compute_cell_volumes`,
`asymmetric_peak`, `set_contour_levels_wholefigure`), which are not included here.
The generated files are already in the repository, so none of this is needed to run the cases.

---

## 8. Cleaning a case

```bash
./w.clean     # deletes 0/, all time directories, processor*, mesh, logs and CSV output
```

`0.orig/` and all input dictionaries are kept, so the case can be rerun from Step 1.

---

## 9. Notes on the medium mesh

The original working folder of the 96³ case was overwritten when it was turned into the
144³ case, so `case_medium` was rebuilt from the fine case. Only the mesh line differs,
which is exactly how the medium case was originally created. The mesh line
`(96 96 96) simpleGrading (1 19.7 1)` is the output of
`openfoam_grid_calc.m` for the settings still present in `turbulent_field_generator.m`
(`nx = ny = nz = 96`, `d1_plus = 0.15`), and the original SLURM job name was `retau_180_96`.

## 10. Changes from the original working files

The physics and numerics are unchanged. The scripts were tidied for reuse:

- `a.pre-processing.sh` now copies `0.orig` to `0` (`restore0Dir`), runs `checkMesh`,
  and finds the latest time of the source case automatically.
- `q.Batch_AWS` runs `reconstructPar -newTimes` after the solver. This is needed both for
  post-processing and to map onto the next mesh.
- `b.post-processing.sh` finds `pvpython` on the PATH, or uses `$PV_PYTHON`.
- `w.clean` uses OpenFOAM's `cleanCase0`.
- `constant/fvOptions`: commented-out alternatives removed and line endings fixed (same
  `meanVelocityForce` source).
- Personal e-mail removed from the job scripts.
