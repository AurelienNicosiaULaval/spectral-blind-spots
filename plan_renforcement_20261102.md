# Targeted post-review extension, 2 November 2026

This extension was specified after the original study, its first revision and an exploratory inspection of the proposed circula construction and the wind example. The original 5,600-dataset plan and its outcomes remain unchanged. These choices are fixed before the new simulation runs; no condition will be selected or replaced based on its rejection rates.

## Rotation shifts

Let g_rho be the wrapped Cauchy density. With uniform X, reference and alternative conditional laws are g_rho(theta-X-delta) and g_rho(theta-X+delta). The null draws both responses from the reference. For rho in {0.5,0.8,0.9}, fix delta=0.5; also use an axial control g_rho(2*(theta-X-delta)) versus the opposite shift with rho=0.8. Each of the eight conditions has 500 datasets, n=240 and 199 shared sign vectors. Five existing projected configurations, full Fourier orders 8 and 12, the original four-band maximum field, and kernel discrepancies at scales 0.3,0.5,1 and their equal-weight kernel mixture are compared. The axial case is a mathematical control, not a claim of empirical prevalence.

For the nine nonaxial rho/delta combinations rho={0.5,0.8,0.9}, delta={0.2,0.5,1}, compare the squared-L2 even-sine projection ratio from the Fourier series and angular quadrature at grid sizes 256 and 512, through even order 128. Also check the axial ratio. Tolerance 1e-8 for ratios and 1e-10 for numerical integrals. Undefined ratios under fidelity remain NA.

## Sensitivity and calibration

Replay all 800 original datasets for phase orders 4 and 6 at n=240 and 720, reproducing the original predictor generation, response draws and 199 signs from recorded seeds. Preserve the five projected diagnostics, Fourier-8, the four-band field and scale-0.5 MMD for comparisons against the recorded outcomes. Add each of the four individual bandwidth combinations, scales 0.3 and 1, and their equal-weight kernel mixture. New comparisons on these data are post-review analyses. Summaries include paired differences, Wilson intervals, Monte Carlo error and all unsuccessful outcomes.

Run exactly 2,000 independent datasets from the existing conditional von Mises fidelity law at n=240 for Fourier-12. Report the prior 10/100 and 33/500 separately, the new run, and the pooled count with the same design. An unexpected size estimate triggers investigation, not repeated sampling until it approaches 0.05.

## Randomness and execution

RNGkind Mersenne-Twister, Inversion, Rejection. Rotation seed=2026110200+10000*condition+rep; new calibration seed=2026211000+rep. Rotation seeds do not overlap the original dataset seeds. New calibration seeds differ from both earlier same-law calibration runs. A later identifier audit found that 200 calibration seed numbers also occur in an original phase-alternative condition: these use different generators and are not extra independent replications of that phase study. The calibration counts are pooled only with the two earlier fidelity experiments, whose seeds are disjoint. No simulation was replaced after this identifier check. Parallelize at most two rotation conditions at a time with mc.set.seed=FALSE; all task seeds are explicit. Existing data seeds and draw order are retained for sensitivity. Store condition-level checkpoints and raw outcomes; do not silently drop failures.

## Wind comparison

Keep 600 training and 276 evaluation pairs and all archived draws of the original model. Fit a von Mises law to Y-X only on training pairs; produce 50 candidate draws per evaluation predictor with seed 2026110250. Compare both models on the same 32x64 grid, support threshold 20 and main bandwidths (0.70,0.30). Use a shared color scale. Report average held-out log densities and daily differences descriptively; model comparison is post hoc and temporal dependence remains.

Extend leave-one-day-out bandwidth CV on training pairs to hx={0.0875,0.175,0.35,0.70,1.20}, ht={0.0375,0.075,0.15,0.30,0.60}. Retain the full score grid and a boundary optimum. Evaluate both model fields at the chosen common bandwidths; no bandwidth is chosen from evaluation performance. Projected displays use 24x24 fixed bins on [-1,1]^2, observed probabilities minus candidate probabilities weighted 1/50 per draw, with the displayed difference consistently candidate minus observed. They are descriptive with no p-values.
