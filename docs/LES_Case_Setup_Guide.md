# Setting Up the Wall-Resolved LES Channel Flow Case — A Simple Guide

This case is a **wall-resolved Large Eddy Simulation (LES) of turbulent channel flow**, solved
in OpenFOAM with `pimpleFoam`. This guide walks through the key files you need to
check before running it, in plain language.

---

## 1. The Big Picture

- **Solver:** `pimpleFoam` (incompressible, transient, PIMPLE algorithm).
- **Turbulence model:** LES with the **WALE** sub-grid model.
- **Geometry:** A single rectangular channel block, periodic in the streamwise (x)
  and spanwise (z) directions, with a solid wall at the bottom and a symmetry plane
  at the top (this models half of a channel, using symmetry). 
- **Flow driving:** A `meanVelocityForce` term that pushes the flow to keep a fixed
  bulk velocity (rather than a fixed pressure gradient).
- **Initial fields:** Pre-mapped turbulent profiles (from a precursor/previous
  simulation) rather than a uniform/zero start, so the flow starts "already turbulent".

Think of it as: a periodic box of fluid over a flat wall surface, 
being pushed at constant speed, with LES resolving the largest
turbulent eddies and modeling the small ones.

---

## 2. `system/blockMeshDict` — the mesh

```
hw = 1.0                 // channel half-height (we are modelling half channel)
lx = 3*pi                // streamwise length
lz = pi                  // spanwise width
blocks: hex (0 1 2 3 4 5 6 7) (64 64 64) simpleGrading (1 23.2 1)
```

**What to check:**
- **Cell counts (64 × 64 × 64):** This is deliberately a *coarse* wall-resolved LES
  mesh — enough to test the setup, but you would refine it (especially in y) for a
  production-quality run.
- **Grading `(1 23.2 1)`:** Only the *wall-normal* (y) direction is graded, with a
  ratio of 23.2. This clusters cells near the bottom wall — essential for a
  wall-resolved LES, since you need very fine cells at the wall to capture the
  viscous sublayer (How do we determine the sufficiently small grid size near the wall?). 
  The x and z directions have no grading (uniform cells),
  which is normal for periodic homogeneous directions.
- **Boundary patches:**
  - `bottomWall` → type `wall` (no-slip wall, the only physical wall)
  - `topSurface` → type `symmetryPlane` (mirrors the channel centerline)
  - `inlet`/`outlet` and `front`/`back` → type `cyclic` (periodic — flow wraps
    around, this is what lets the "channel" be infinite in x and z)

If you ever change the mesh resolution, make sure the grading still clusters cells
near `bottomWall` and that the domain size (`lx`, `lz`) stays large enough relative
to the channel half-height to capture the large turbulent structures.

---

## 3. `system/decomposeParDict` — parallel decomposition

```
numberOfSubdomains  8;
method              hierarchical;
coeffs { n (4 1 2); }
```

**What to check:**
- `numberOfSubdomains` (8) must match the number of MPI ranks you actually launch
  with (`mpirun -np 8 pimpleFoam -parallel`, as seen in `q.Batch_AWS`: `#SBATCH -n 8`).
- `n (4 1 2)` means: split into 4 pieces in x, 1 piece in y (don't split the
  wall-normal direction — good, since that's where all the small, expensive cells
  are and splitting it can hurt load balance and wall-modeling accuracy), and 2
  pieces in z. 4 × 1 × 2 = 8, which matches `numberOfSubdomains`. Always check that
  the product of `n` equals `numberOfSubdomains`.
- If you change core count, you must update **both** `numberOfSubdomains` and the
  `n` triplet (and the SLURM `-n` value) consistently.

---

## 4. `system/controlDict` — time control & solver

```
application   pimpleFoam
startFrom     latestTime
endTime       400
deltaT        1.5e-3
writeControl  timeStep
writeInterval 4000
purgeWrite    0
runTimeModifiable true
```

**What to check:**
- **`startFrom latestTime`**: the run will resume from whatever the most recent
  time directory is. For a fresh run this means the `0` (or `0.orig` → `0`)
  directory.
- **`endTime`/`deltaT`**: 400 time units at Δt = 1.5×10⁻³ ≈ 266,667 timesteps —
  a long LES run, as expected for gathering converged turbulence statistics.
- **`writeInterval 4000`** with **`purgeWrite 0`**: a full field write every 4000
  steps, and *all* of them are kept (purgeWrite 0 = never delete old ones). For a
  long run this can consume a lot of disk — consider `purgeWrite` > 0 if space is
  tight, unless you need every field snapshot.
- **`runTimeModifiable true`**: you can edit `controlDict` (and other dictionaries)
  while the simulation is running, and OpenFOAM will pick up the changes — handy
  for adjusting `endTime` or write settings on the fly.
- **`functions { fieldAverage1 {...} }`**: this activates on-the-fly averaging of
  `U` and `p` (mean and prime² / variance) — this is how you get the
  time-averaged turbulence statistics that LES studies need. Check that
  `writeControl writeTime` here matches how often you actually want averaged
  fields written (tied to the main `writeInterval` above).
- **Also check the Courant number implicitly**: with LES you generally want CFL < 1
  (often ≤ 0.5) for accuracy — verify this once running via `log.pimpleFoam`,
  since it isn't fixed directly in `controlDict` here (no adjustable time-stepping
  is set; `deltaT` is fixed).

---

## 5. `constant/fvOptions` — momentum source

```
momentumSource
{
    type          meanVelocityForce;
    selectionMode all;
    fields        (U);
    Ubar          (15.67 0 0);
}
```

**What to check:**
- This is what actually **drives the flow** in a periodic channel (since there's
  no real inlet/outlet to impose a pressure difference). It adds a uniform body
  force to `U` every timestep to force the volume-averaged streamwise velocity to
  equal `Ubar = 15.67 m/s`.

- **How was this bulk velocity value determined?** (See Pope 2000)

- **Consistency check:** `Ubar` here should match the `Ubar` value in
  `constant/transportProperties` (it does: both are `15.67 0 0`) — these are
  meant to be the same bulk velocity target used to set up the case (e.g. via the
  `channelFoam`-style Re_τ estimate), so don't let them drift apart if you tune one.

---

## 6. The `0` folder (`0.orig/U`, `p`, `nut`) — initial & boundary conditions

Rather than starting from zero/uniform fields, this case **maps in a pre-computed
turbulent flow field** as the initial condition:

```
#include  initialVelocityProfile;
internalField $velocityProfile;
```

(same pattern for `p` → `initialPressureProfile`, and `nut` → `initialViscosityProfile`)
(See the "turbulent_field_genrator.m" Matlab script)

**What to check:**
- These `initial*Profile` files are large lists of numbers — one value per cell —
  taken from a previous/precursor simulation (or synthetic turbulence generator).
  **Make sure the number of entries matches your mesh's cell count** (64×64×64 =
  262,144 cells here). If you ever change the mesh resolution in `blockMeshDict`,
  these mapped profiles will no longer match and you'll need to regenerate them
  (e.g., with `mapFields`, or a new precursor run) — this is the single most common
  mistake when reusing this kind of setup on a different mesh.
- **Boundary conditions**, consistent with the mesh patches:
  - `U`: `fixedValue uniform (0 0 0)` at `bottomWall` (no-slip), `symmetryPlane` at
    `topSurface`, `cyclic` elsewhere.
  - `p`: `zeroGradient` at the wall, `symmetryPlane` at the top, `cyclic` elsewhere.
  - `nut` (SGS turbulent viscosity): `zeroGradient` at the wall.
- `dimensions` in each file should match physical expectations: `U` → m/s,
  `p` → m²/s² (kinematic pressure, consistent with an incompressible solver where
  `p` is really `p/ρ`), `nut` → m²/s.

---

## 7. Other supporting files (quick check)

- **`constant/transportProperties`**: `nu = 0.00556` (kinematic viscosity) and
  `Ubar (15.67 0 0)` — together with the channel half-height (`hw = 1`), these
  define the friction Reynolds number of the case (consistent with "Re_τ = 180").
  If you change `nu`, remember the flow regime (Re_τ) changes too.
- **`constant/turbulenceProperties`**: `simulationType LES`, `LESModel WALE`. WALE
  is a good choice for wall-resolved LES because it naturally goes to zero at
  walls without extra damping functions — consistent with not needing a wall
  function on `nut`.
- **`system/fvSchemes` / `system/fvSolution`**: not detailed here, but for LES you
  should double check that:
  - Convection schemes are low-dissipation (e.g., a blended/central scheme) rather
    than a fully upwind scheme, since excess numerical dissipation can kill the
    resolved turbulence.
  - Time scheme is second order (e.g. `backward` or `CrankNicolson`) — LES needs
    good temporal accuracy.
- **`a.pre-processing.sh`**: runs `blockMesh` then `decomposePar` — run this first.
- **`q.Batch_AWS`**: the SLURM job script; launches `mpirun pimpleFoam -parallel`
  with `-n 8`, matching `decomposeParDict`.
- **`b.post-processing.sh`** / **`paraview_openfoam.py`**: post-processing via
  ParaView's `pvpython`, reading the `foam.foam` case marker file.
- **`w.clean`**: a cleanup script that removes processor directories, time
  directories, logs, and old `postProcessing` results — use this to reset the
  case before a fresh run (be careful, it's destructive).

---

## 8. Suggested Run Order

1. Check mesh cell count vs. size of the `initial*Profile` files (Section 6).
2. Check `numberOfSubdomains` matches your intended core count and SLURM script
   (Section 3).
3. Run `./a.pre-processing.sh` (builds mesh, decomposes for parallel run).
4. Submit `q.Batch_AWS` (or run `mpirun -np 8 pimpleFoam -parallel` directly).
5. Monitor `out.o` for Courant number and residuals and error.e for runtime errors.
6. Post-process with `./b.post-processing.sh` once complete.
7. Use `w.clean` if you need to wipe results and start over.

---

# Part 2 — Tasks for the Student

These tasks build on the case above. Work through them in order — don't skip
ahead to "fix the issue" before you've genuinely tried to diagnose it yourself.

## Task 1 — Run, Observe, Diagnose

**Do this:**
1. Run the pre-processing script (`./a.pre-processing.sh`) to build and decompose
   the mesh.
2. Start the simulation, but only let it run for a **short** time.
3. Once it stops, inspect the results. Useful things to check:
   - The solver log (residuals, Courant number, any warnings — scroll through
     the whole log, not just the last lines).
   - The flow field itself in ParaView — look at `U` over time.
   - Whether basic physical sanity checks hold: does the flow look like a
     channel flow? Is anything blowing up, oscillating wildly, symmetric when
     it shouldn't be, or flat/uniform when it should have structure?
4. Write down **what you observe** (be specific: what field, what location,
   what does it look like compared to what you'd expect).
5. Before asking anyone: brainstorm at least **2–3 possible explanations** for
   what you're seeing. For each one, ask yourself: "if this were the cause,
   what else would I expect to see?" — and check whether that's consistent
   with what you actually observe. This is how you narrow down a diagnosis
   instead of guessing.
6. Only after you've done steps 4–5, bring your observations and your
   explanations to your supervisor to confirm your diagnosis. 
   Bring evidence (log snippets, screenshots, field values).

## Task 2 — Compute Wall Shear Stress and Friction Velocity

This is a core LES/channel-flow diagnostic, and you'll need it both to validate
the case (Task 1) and to characterize the flow physically.

**Background:** 
What is the **wall shear stress**?
How is it calculated?
What is **friction velocity**?
How is it calculated?

**Do this:**
1. Implement a way to compute the wall-normal velocity gradient at
   `bottomWall:
   - using OpenFOAM's built-in **`wallShearStress`** function object (add it
     under `functions {}` in `controlDict`, alongside `fieldAverage1`)
2. From the wall shear stress, compute `u_τ`.
3. Sanity-check your result: given `nu = 0.00556` and the target friction
   Reynolds number implied by the case name (`retau_180`), does your computed
   `u_τ` give something close to Re_τ = u_τ · hw / ν ≈ 180? (Use the
   time-averaged field from `fieldAverage1`, not an instantaneous snapshot —
   instantaneous wall shear stress fluctuates a lot in a turbulent flow.)
4. Document your method and result (a short write-up: what you computed, how,
   and what value you got, compared to the expected Re_τ).

## Task 3 — Fix the Issue and Re-run

1. Based on the diagnosis confirmed with your supervisor in Task 1, correct
   the relevant file(s) in the case setup.
2. Re-run the pre-processing steps if you changed anything upstream of meshing
   or decomposition (e.g. `blockMeshDict`); otherwise you may be able to just
   restart the solver.
3. Repeat the short-run check from Task 1 to confirm the issue is actually
   resolved (don't jump straight to the full 400-time-unit run).
4. Once confirmed healthy, launch the full run (e.g. via `q.Batch_AWS`).
5. Once complete, recompute `u_τ` and Re_τ (Task 2) from the fixed run's
   time-averaged fields, and compare against the short faulty run — this
   comparison is good evidence for your write-up that the fix actually mattered.

## Task 4 — Periodically Check for Convergence

**Guiding questions — what does "convergence" even mean here?**
There are at least two distinct things people mean by "converged" in a transient
LES, and it's worth being precise about which one you're checking at any given
moment:
- Has the flow itself reached a **statistically stationary state** — i.e., it's
  no longer "developing" from its initial condition, and is just fluctuating
  around a stable mean?
- Have your **time-averaged statistics** (the mean/RMS profiles from
  `fieldAverage1`) stopped changing as you extend the averaging window — i.e.,
  the *averaging* has converged, which is a separate question from whether the
  underlying flow is stationary?

**Do this, periodically, as the run progresses:**
1. Monitor a simple global or near-wall quantity over time — for example, the volume-averaged
   turbulent kinetic energy. Plot it vs. simulation time. What should this curve
   look like once the flow is statistically stationary?
2. To check whether the *averaging* has converged (separately from the flow
   itself), compute your mean profiles (Task 5) using two different windows —
   e.g., the second half of the data collected so far vs. all of it — and compare.
   Do they agree within an acceptable tolerance, or are they still visibly
   different?
3. Keep a simple log/plot of your monitored quantity across checks so you can
   see the trend, not just a single snapshot (use Gnuplot or foamMonitor).

## Task 5 — Compute Flow Statistics

**Guiding questions — how long is "enough" before computing statistics?**
- Rather than thinking in raw simulation time units, think in terms of
  **flow-through times** (how long it takes fluid to traverse the streamwise
  domain length at the bulk velocity `Ubar`) or **eddy turnover times**. Given
  `lx` (from `blockMeshDict`) and `Ubar` (from `transportProperties`), how many
  flow-through times does the current `endTime` actually represent? Is that
  generous or marginal for converged turbulence statistics?
- Should you discard an initial transient period before you start averaging, or
  average from `t = 0`? What did Task 4 tell you about when the flow actually
  became statistically stationary?

**Guiding questions — what to compute, and what's "relevant" for channel flow:**
- Because the domain is **periodic in x and z**, what does that imply about how
  the true (statistical) flow should depend on those two directions? Which
  single coordinate should essentially all of your averaged profiles be plotted
  against?
- **Mean velocity profile:** what's the standard way of non-dimensionalizing
  this profile for channel flow, using the friction velocity from Task 2, so it
  can be compared across Reynolds numbers and studies (hint: this is the
  classic "law of the wall" form)?
- **Mean pressure:** since the flow is driven by `meanVelocityForce` rather than
  an imposed streamwise pressure gradient, what would you actually expect the
  mean pressure to do along x? Is a wall-normal mean-pressure *profile*
  meaningful/interesting here, or is the more relevant pressure-related quantity
  something else (think about what balances the imposed body force in the mean
  momentum equation)?
- **Turbulence statistics relevant to channel flow** — think about:
  - The Reynolds stress tensor components: which ones are non-zero by the
    symmetry of this flow, and which are related to each other or vanish?
  - RMS (root-mean-square) fluctuations of each velocity component individually
    — why might these differ from each other in magnitude, and in shape across
    the channel?
  - The turbulent kinetic energy profile (how is it built from the quantities
    above?).
  - The `fieldAverage1` function object already active in `controlDict` writes
    `UMean`, `UPrime2Mean`, `pMean`, `pPrime2Mean`. Which of the statistics
    above can you read directly from these, and which require you to do a bit
    more work (e.g., extracting individual Reynolds stress components from
    `UPrime2Mean`)?
- What non-dimensionalization should you use for the wall-normal coordinate
  itself when plotting all of these?

## Task 6 — Compare Against Literature and Discuss

**Do this:**
1. Search for established reference data for turbulent channel flow at
   Re_τ ≈ 180 — both DNS and, where available, experimental data. What makes a
   given dataset a fair comparison (matching Re_τ, matching how quantities are
   non-dimensionalized, similar domain size)? (the dns profiles are available in Github).
2. Overlay your mean velocity profile, RMS profiles, and Reynolds stresses with
   the reference data.
3. For each profile, note **where** the agreement is good and where it isn't —
   near the wall, in the log-law region, near the centerline? Is the pattern of
   disagreement consistent across different quantities, or does one particular
   statistic stand out?
4. Discuss **why** you might expect the differences you observe, considering
   in particular:
   - Mesh resolution — this is a deliberately **coarse** wall-resolved LES mesh
     (Section 1, point 2 above). Where would under-resolution typically show up
     most in the profiles?
   - Averaging time — does Task 4's convergence check suggest your statistics
     might still be under-converged?
