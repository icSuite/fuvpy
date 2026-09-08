# Repository guidance

## Purpose and project tracks

`fuvpy` reads and analyses global far-ultraviolet auroral images. Its main
tracks are:

- IMAGE/Polar file loading and camera geometry;
- B-spline dayglow and spherical-harmonic residual-background removal;
- auroral-boundary modelling; and
- plotting and interactive inspection.

Work on `background_uncertainty` is deliberately narrower: quantify
uncertainty introduced by the fitted background models without changing their
default fitted images or reopening the choice of background algorithm.

## Project memory

Project memory is stored in the repository-root `vault/`.
Before substantive work, read:

1. `vault/01_Project/Current State.md`
2. `vault/05_Handoff/Handoff - Latest.md`
3. `vault/01_Project/Project Brief.md`

Read `vault/02_Algorithm/Background Fit Uncertainty Architecture.md` before
changing the uncertainty calculation or public background-model interface.

Use this source-of-truth order:

1. live code, Git state, tests, and regenerated results;
2. `Current State.md` and `Handoff - Latest.md`;
3. decision and algorithm notes;
4. dated session history.

If memory conflicts with live evidence, follow the repository and update the
live notes when the task changes project understanding.

## Scientific boundaries

- Distinguish detector/count noise, conditional fitted-background uncertainty,
  and background-model discrepancy. They are not interchangeable.
- `dgsigma` is a residual-derived spread model. It is not currently the
  posterior standard deviation of `dgmodel`.
- Tukey and spatial weights are not automatically inverse measurement
  variances. Do not call the inverse regularized normal matrix a posterior
  covariance without stating and checking the required likelihood and prior
  assumptions.
- The final weights depend on the data. Any analytic covariance that freezes
  them is conditional on the converged IRLS solution and should be checked
  against repeated fits.
- The WIC correction uses sequential BS and SH fits. Their errors are
  dependent; do not add two separately estimated variances as if the stages
  were independent.
- Use complete image sequences for scientific validation of temporal spline
  fits. A few extracted frames do not reproduce a full-orbit fit.
- Background-fit uncertainty can quantify sensitivity around the selected
  model. It does not include calibration, flat-field, atmospheric-response,
  or wrong-model systematic error unless those are varied explicitly.

## Working rules

- Preserve unrelated work and public behavior. The default background output
  must remain numerically unchanged while uncertainty support is prototyped.
- Keep the main scientific calculation visible. Prefer direct NumPy/SciPy
  operations and ordinary dictionaries over framework or class machinery.
- Use short helpers for distinct calculations such as the weighted solve,
  covariance calculation, or prediction variance. Do not hide the fitting
  sequence behind a generic solver framework.
- Make units, coefficient ordering, design-matrix construction, regularization,
  convergence, and uncertainty interpretation explicit.
- Do not overwrite tracked examples or generated figures during exploratory
  runs. Use a scratch directory and begin with the bundled WIC sequence or a
  copied complete orbit.
- Do not run a large external IMAGE corpus, download data, install packages,
  push branches, or alter upstream remotes unless the user explicitly places
  that action in scope.

## Verification

No automated test suite is present at the initial `background_uncertainty`
checkpoint. Add focused tests before changing public behavior.

Safe syntax verification without writing bytecode:

```bash
python -c "import ast, pathlib; files=list(pathlib.Path('fuvpy').rglob('*.py')); [ast.parse(p.read_text(), filename=str(p)) for p in files]; print(f'parsed {len(files)} files')"
```

For documentation-only changes, run:

```bash
git diff --check
```

Once focused tests exist, run only the background-uncertainty tests first.
Full example processing may require `polplot`, old environment versions, and
external IMAGE files.

## Automatic memory checkpoints

Project-memory maintenance is required after a durable decision, verified
implementation, scientific result, changed blocker, or changed next action.

At a meaningful checkpoint:

1. create `vault/04_Sessions/YYYY-MM-DD.md` only when historical detail is
   worth retaining;
2. rewrite `vault/01_Project/Current State.md` with verified state;
3. append only durable choices to the decision log;
4. replace obsolete content in `vault/05_Handoff/Handoff - Latest.md`;
5. update algorithm notes when interpretation or architecture changes; and
6. refresh the handoff's `Portfolio impact` section.

Do not copy raw logs or long conversation history into the vault. Project
agents write project memory; the central portfolio agent owns the second
brain.
