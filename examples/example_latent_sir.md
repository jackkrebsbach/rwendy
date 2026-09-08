# SIR with only infected observed

Run `examples/example_latent_sir.R` from the repository root. It simulates the
three compartments in people, observes only infected prevalence, and knows the
constant population `N = 10000`. The example uses 64 observations over days
0–40, covering growth, the peak and decline, with additive Gaussian noise at
1% of infected-trajectory RMS. Transmission and recovery rates are constant.
There is no unknown reporting fraction, incidence observation model, birth,
death, migration or reinfection.

## Why the model is identifiable

The physical model is

\[
\dot S=-\beta SI/N,\qquad
\dot I=\beta SI/N-\gamma I,\qquad
\dot R=\gamma I,\qquad S+I+R=N.
\]

Let `i=I/N`, `s=S/N`, and `q=i'/i` on an interval with positive infection.
Then

\[
q=\beta s-\gamma,\qquad
\dot q=-\beta i(q+\gamma),\qquad
-\dot q/i=\beta q+\beta\gamma.
\]

With an exact observed infection curve and nonconstant `q`, this last relation
uniquely determines its slope `beta` and intercept `beta*gamma`. For `beta>0`,
it therefore determines `gamma` too. Both hidden states follow uniquely:

\[
S=N(q+\gamma)/\beta,\qquad R=N-S-I.
\]

Their initial values follow by continuity, while `I(0)` is itself the observed
output's initial value. Thus both rates and all initial states are generically
globally structurally identifiable, subject to the one population constraint.
The fit does **not** assume `R(0)=0` or provide true susceptible/recovered initial
values. The susceptible starting guess of 70% is an optimizer initialization.

This input-output argument agrees with the infected-only SIR result in
[Example 3.2 of Heitzman-Breen et al., *A Practical Identifiability Criterion Leveraging
Weak-Form Parameter Estimation*](https://link.springer.com/article/10.1007/s11538-026-01639-x).
Their transmission coefficient uses mass-action count units; here the
coefficient multiplying `SI` is `beta/N`.

The disease-free case, zero transmission/susceptibles, or an effectively
exponential short observation window do not supply the same information.
Structural identifiability concerns ideal curves; it does not guarantee
precision from every finite noisy dataset.

## How conservation enters the solver

The simulation starts with the full three-state solution. The observation
matrix has only its infected column filled. For fitting, conservation eliminates
recovered exactly, leaving `(i,s)` with `s` as the one independent latent state.
The solver automatically detects this entirely NA column. After fitting,
recovered is reconstructed as `N*(1-i-s)` and all three states are plotted.
No synthetic susceptible/recovered observations or soft population penalty are
introduced. The example checks positive rates and nonnegative reconstructed
states without clipping.

Both components use population-fraction ODE scales of one. `gp_jitter=1e-5`
is an explicit numerical stabilization for the smooth latent trajectory's
stage-3 coordinate map. The default physical stationarity tolerance remains
`1e-4`; the latent GP quadratic remains disabled. This is a numerical setting,
not an identifiability assumption or a prior that fixes unknown initial states.

## Independent numerical information check

The example integrates analytic sensitivities of the **infected output only**
with respect to four independent unknowns: `(beta, gamma, s0, i0)`. The remaining
initial fraction is `r0=1-s0-i0`. It reports the rank and column-scaled condition
number, plus local noise scales from the inverse data-information matrix.
These contain no GP-prior rows and are not confidence intervals for WENDyGP.
The weak residual's joint rank is not used as proof of identifiability.

The seeded run passed physical stationarity, observation interpolation,
quadrature, grid resolution and physical admissibility. It gave:

| Quantity | True | Estimated |
| --- | ---: | ---: |
| Transmission rate beta | 0.400 | 0.3893 |
| Recovery rate gamma | 0.120 | 0.1231 |
| Initial susceptible | 8500 | 8772.4 |
| Initial infected | 500 | 487.4 |
| Initial recovered | 1000 | 740.2 |

The infected-only sensitivity rank was **4/4**, with column-scaled condition
number **49.48**. Initial recovered was less precise than the rates or infected
initial value; its local noise scale was about 147 people. Population
conservation error was below `2e-12` people. The physical scaled gradient was
`5.56e-5`; native nlminb reported singular convergence, which the solver retains
separately from the passing physical stationarity check.

These are one-run diagnostics, not a claim that every parameter or hidden
initial condition is estimated exactly, or that all noisy designs are equally
informative. The example stops if its physical or accuracy checks fail.
