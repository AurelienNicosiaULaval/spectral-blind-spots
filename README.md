# Spectral Blind Spots in Circular Imputation Diagnostics

Aurélien Nicosia, Université Laval

Code, data and recorded results for the circular diagnostic simulations and
the Texas wind direction example.

## Run the analyses

Run the commands from this directory in a working copy. The original run used
R 4.5.0 and ggplot2 4.0.2. The environment files record the versions used for
the original and subsequent runs. Install ggplot2 in R if needed.

```r
install.packages("ggplot2")
```

Rebuild the original tables and figures and check 44 datasets against the
440 recorded outcomes.

```sh
Rscript --vanilla reproduce.R --verify
```

To rerun the additional simulations and the wind analysis as well, use the
following command. The rotation and sensitivity scripts save checkpoints
in `renforcement/` during execution. Use a fresh copy for a fresh run.

```sh
Rscript --vanilla reproduire_revision.R --verify
```

Recalculate all 5,600 original datasets in a separate output directory
and compare them with the recorded results.

```sh
Rscript --vanilla benchmark.R --output=resultats_recalcules
Rscript --vanilla compare_full_rerun.R resultats_recalcules
```

The full run uses the seeds and settings in [the simulation plan](plan_fixe.json).
The comparison checks the numerical results and excludes elapsed time.
The small paired example can be run separately with `Rscript --vanilla toy.R`.

## Files

- `diagnostics_core.R` defines the generators and diagnostics.
- `resultats_simulation/` holds the original outcomes and summaries.
- `complement/` and `renforcement/` hold the additional simulation results.
  The [extension plan](plan_renforcement_20261102.md) records the design and seeds.
- `donnees/` holds the hourly wind directions. See the [data notice](donnees/PROVENANCE.md)
  for the data source and citation. The [application plan](plan_application.md)
  describes the chronological training and evaluation split.
- `application/` holds the fitted parameters, scores and fields for the wind
  example. The plotting scripts write figures in `figures/` and `application/`.
- `controles_revision/` holds the additional checks. The optional Python checks
  use NumPy and SciPy. Run `python3 controle_renforcement.py` to check the extension
  against its recorded outputs.

The SHA-256 manifest describes the files before execution. The scripts write
tables, figures and checks in the working copy.
