# CCST-structural-stability-recovery-env
Stability Geometry of Branch Structures in Open Quantum Systems (Extended): Recovery Geometry and Environmental Coupling

## A.1 Purpose and scope (honesty statement)

ccst_env_recovery_sweep_vnext.py is a toy-model figure generator for a reduced, bounded proxy state η(t) ∈ [0, 1] interpreted as a structural coherence / stability indicator used in the manuscript’s stability-geometry discussion. The script is designed to produce paper-ready, reproducible visualizations that cleanly separate and illustrate:

Baseline environment coupling (deterministic drift backbone)

Recovery geometry (restoring drift; symmetric quad or direction-asymmetric asym)

Environmental fluctuations (fast white noise or correlated ou noise), acting on either:

target=lambda (effective coupling-rate channel), or

target=Phi (environmental drive channel)

Honesty note (what this is not).
This script is not a microscopic branching model. It does not define probabilities, Born weights, or any selection rule. It is intentionally a visual regime toy: it isolates a deterministic drift backbone, then adds (i) recovery geometry and (ii) environment variability as fluctuations acting on a specified target channel. Its value is interpretability and auditability, not microphysical fidelity.

## A.2 Model definition

We simulate an ensemble of N trajectories of a bounded proxy state η_b(t) ∈ [0, 1] for branch label b.

### A.2.1 Baseline deterministic backbone (environment-only)

Each trajectory draws a heterogeneous coupling parameter λ_b once at initialization:

λ_b ~ Uniform(λ_min, λ_max)

The environment-only deterministic dynamics are:

dη_b/dt = − λ_b η_b − γ Φ(η_b)

Φ(η) = Φ0 + α (1 − η)

Here, the “Φ-channel” is an environment-drive term with a simple state-dependence. The backbone is deterministic once λ_b is sampled.

### A.2.2 Recovery geometry (optional)

Recovery adds a restoring drift toward a target level η*:

dη_b/dt = − λ_b η_b − γ Φ(η_b) + k_side(η_b) (η* − η_b)

Supported recovery types:

recovery=none: k_side(η) = 0

recovery=quad: k_side(η) = k (symmetric recovery strength)

recovery=asym (direction-asymmetric):

k_side(η) = k_minus if η < η*
k_plus if η ≥ η*

This implements direction-asymmetric recovery while keeping the same basin center η*.

Bounding. After each update, η is clipped to [0, 1] to maintain the proxy interpretation.

## A.3 Environmental fluctuations (optional)

Fluctuations can act on λ or Φ.

### A.3.1 Fluctuation targets

target=lambda: λ_b → λ_b + δλ_b(t)
This enters the η-dynamics as an extra drift contribution: − δλ_b(t) · η_b

target=Phi: Φ → Φ + δΦ_b(t)
This enters the η-dynamics as an extra drift contribution: − γ · δΦ_b(t)

### A.3.2 Fluctuation models

Two fluctuation models are supported:

(i) white noise (Euler–Maruyama increment)
At each timestep, a Gaussian increment is applied (Euler–Maruyama):

target=lambda: η ← η − (σ_env · η) √(Δt) Z

target=Phi : η ← η − (γ σ_env) √(Δt) Z

where Z ~ Normal(0, 1). In this mode, σ_env is a diffusion amplitude (per √time).

(ii) ou noise (correlated fluctuations)
An OU state X(t) is evolved with correlation time τ and stationary standard deviation σ_env:

dX = −(1/τ) X dt + √(2/τ) σ_env dW

The script uses an exact OU discretization:

a = exp(−Δt/τ)

X_{t+Δt} = a X_t + σ_env √(1 − a²) Z, Z ~ Normal(0, 1)

After updating X, it is inserted into the drift channel:

target=lambda: add drift − (X · η)

target=Phi : add drift − γ X

Bounding. After each update (drift + noise), η is clipped to [0, 1].

## A.4 Diagnostic: structural survival fraction R_struct

For a threshold η_crit, a trajectory survives on [0, T] iff:

min_{t ∈ [0, T]} η(t) ≥ η_crit

The structural survival fraction is:

R_struct(η_crit) = (1/N) Σ_{b=1..N} 1[ min_{t∈[0,T]} η_b(t) ≥ η_crit ]

All main figures plot R_struct (typically in %) as a function of η_crit.

Default threshold grid (if not specified):

10 evenly spaced points in [0.50, 0.95]

## A.5 Figure types produced by the script
### A.5.1 Overlay (three-stage) plot

A single plot with three curves:

env-only: no recovery, no fluctuations

+ recovery: recovery enabled, no fluctuations

+ recovery + fluctuations: recovery enabled and fluctuations enabled

This isolates the incremental effect of (i) recovery geometry and (ii) environmental variability.

### A.5.2 Grid plot (paper-ready)

A panel grid with:

columns = σ_env values (environment fluctuation strength)

rows = recovery strength values

for quad: rows are literal k values

for asym: rows are multipliers s, applied as (k_minus, k_plus) ← s·(k_minus, k_plus)

s = 0 effectively turns recovery off

The grid is paper-ready: panels share the same y-axis limits (shared scale), making small differences visually comparable. Per-panel legends are removed to reduce clutter; outer labels provide readability.

Note. If recovery=none, sweeping rows by k-list is redundant (rows will be identical up to Monte Carlo noise).

## A.6 Reproducibility: RNG, seeds, and CRN

The script supports stable runs via --seed.

### A.6.1 CRN (common random numbers)

By default, CRN is ON (crn=True): the same underlying randomness is reused across conditions to reduce Monte Carlo variance in differences between curves (stabilizing comparisons). Disable via:

--no-crn

Turning CRN off is useful as a robustness check: qualitative patterns should persist without variance reduction.

## A.7 Command-line usage (complete)
### A.7.1 Modes

--overlay : generate a 3-stage overlay plot

--grid : generate a paper-ready grid plot

If neither is specified, the default behavior is to generate an overlay.

### A.7.2 Presets

--preset <name> : run one named preset configuration (overlay or grid)

--run-all-presets : generate one figure per preset (“one of each kind” final test set)

--list-presets : list available preset names

### A.7.3 Core parameters

--N : ensemble size (number of trajectories)

--T : time horizon T

--dt : timestep Δt

--eta0 : initial condition η(0)

### A.7.4 Baseline coupling and drive

--lambda-min, --lambda-max : sampling range for λ_b

--gamma : strength of the Φ channel

--Phi0 : base drive level

--alpha : state dependence in Φ(η) = Φ0 + α(1 − η)

### A.7.5 Recovery geometry

--recovery {none,quad,asym}

--eta-star sets η*

For quad: --k sets k

For asym: --k-minus, --k-plus set (k_minus, k_plus)

### A.7.6 Fluctuations

--env-model {none,white,ou} (default is white when not using presets)

--fluct-target {lambda,Phi}

--env-sigma : fluctuation amplitude (fluctuations are effectively off if 0)

--tau : OU correlation time (only used if env-model is ou)

### A.7.7 Threshold grid

--eta-crit-list : comma-separated list; if omitted, uses 10 points in [0.50, 0.95]

### A.7.8 Grid sweep lists

--env-sigma-list : comma-separated values for grid columns

--k-list : row values

for quad: literal k values

for asym: multipliers s applied as (k_minus, k_plus) ← s·(k_minus, k_plus); s=0 disables recovery

Grid cosmetics:

--grid-annotate {none,int,2dp} : annotate values inside panels (normally none for paper)

--grid-draft : larger draft styling (debugging)

RNG control:

--seed

--no-crn

## A.8 Presets included

Presets are defined in preset_table() and are intended to provide “paper-meaningful” configurations. At the time of writing, presets include:

Overlay presets:

overlay_white_lambda_quad

overlay_white_lambda_asym

overlay_ou_phi_quad

overlay_ou_phi_asym

Grid presets:

grid_white_lambda_quad

grid_white_lambda_asym

grid_ou_phi_quad

grid_ou_phi_asym

Example:

python3 ccst_env_recovery_sweep_vnext.py --preset grid_ou_phi_asym --outdir figs

python3 ccst_env_recovery_sweep_vnext.py --run-all-presets --outdir figs

## A.9 Practical parameter intuition (how to explore)

Common ways to decrease survival (lower R_struct):

increase --lambda-max (stronger coupling)

increase --gamma or --Phi0 (stronger Φ channel)

increase --env-sigma (stronger fluctuations)

increase --T (more time to degrade)

decrease --eta0 (worse initial condition)

Common ways to increase survival (higher R_struct):

increase recovery strength (--k or --k-minus/--k-plus)

choose η* in the high-coherence region (stable basin center)

Asymmetry guidance:

If the narrative is “fast recovery from degraded states, slower erosion otherwise,” typically choose k_minus > k_plus.

## A.10 Troubleshooting

Curves look too similar:

increase --env-sigma

increase --lambda-max

increase --T

widen the recovery sweep range in --k-list

increase --N

For clean comparisons: keep CRN on

For robustness: turn CRN off (--no-crn) and confirm qualitative stability

Numerical stability:

If artifacts appear, reduce --dt and verify that results are stable under timestep refinement.
