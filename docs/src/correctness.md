# Sampler Correctness

Every sampler in this package is statistically tested against the definition
of the process it claims to sample. This page says what is tested, how, and
what the current results are, so that you can judge how much to trust each
sampler — and rerun the evidence yourself.

## The trusted reference

`FirstReactionMethod` is the reference sampler. It is a direct transcription
of the definition of a competing-clocks process: at every query it draws each
enabled clock's firing time from its distribution, conditioned on survival to
the current time, and returns the minimum. It retains no state between
queries, so there is nothing to cache incorrectly. It is the slowest sampler
here and the easiest to believe; every other sampler earns trust by agreeing
with it.

## The test battery

Each sampler runs under each test condition through three statistical tests.

**Doob–Meyer step cumulant.** For each simulated step, the integrated hazard
of the waiting time, transformed as ``U = 1 - e^{-H}``, must be uniform on
``(0,1)`` if the sampler draws from the correct waiting-time distribution.
The recorded step cumulants are pooled and tested for uniformity
(Kolmogorov–Smirnov). This checks the *when* of each event.

**Mark calibration.** Conditional on the firing time, the probability that a
particular clock is the one that fires is the ratio of its hazard to the
total hazard. Comparing observed clock frequencies against these conditional
probabilities checks the *which* of each event.

**Two-sample comparison against the reference.** A fixed history is condensed
into its final enabled set; many replicas of the sampler under test and of
`FirstReactionMethod` each draw the next event from that state, and the two
firing-time samples are compared with the k-sample Anderson–Darling test.
The query respects the [`next`](@ref) contract — every replica is asked at
the time of its most recent state change — so the comparison measures the
distribution of the next event given the history, which is exactly the
quantity the samplers are supposed to agree on.

## The matrix

The battery runs over a matrix of samplers and conditions on a travel model
(a walker on a graph whose transition clocks compete):

* **Distribution family**: all-exponential; a Weibull mix with both
  decreasing- and increasing-hazard clocks; and exponential with shifted
  (delayed) enabling times.
* **Memory**: `forget` (a re-enabled clock starts fresh) and `remember`
  (a re-enabled clock's enabling time is shifted left by its accumulated
  waiting, so it keeps its age).
* **Graph**: a cycle and a complete graph, which change how many clocks
  compete and how often clocks are enabled and disabled.

Samplers that accept only exponential distributions (`DirectMethod` in its
four variants, `RejectionMethod`, `PartialPropensityMethod`) run only the
exponential conditions, by design. `MultiSampler` compositions — clocks
split across sub-samplers, including a mixed composition — are tested as
samplers in their own right, because the superposition re-queries losing
sub-samplers at times they did not choose, and that deserves direct evidence
rather than argument.

Because the matrix runs hundreds of simultaneous tests, some raw p-values
below 0.05 are expected by chance. Two disciplines control this:

1. **False-discovery-rate control.** Benjamini–Hochberg correction is applied
   across the entire matrix; only corrected flags count.
2. **A built-in null control.** `FirstReactionMethod` runs as its own subject
   — compared against itself through the identical machinery. Its raw flag
   rate is the noise floor; a sampler is only suspect if it exceeds what the
   reference produces against itself.

## Current results

Run of 2026-07-05: 12 samplers and compositions × up to 12 conditions =
96 cells × 3 tests = **288 results**, ``n = 1000`` replications per cell,
base seed 20260705.

* **Zero of 288 results survive Benjamini–Hochberg correction** at
  ``\alpha = 0.05``; the smallest corrected p-value is 0.67.
* Twenty raw p-values fall below 0.05 — consistent with the ~14 expected
  from 288 simultaneous tests — and the null control (`FirstReactionMethod`
  against itself) shows the highest per-sampler raw flag rate, 0.139. No
  sampler under test meaningfully exceeds the reference's own noise floor.
* No sampler shows a statistically supported defect.

The verdict matrix from that run, one `DM/AD/MK` glyph triple per cell
(Doob–Meyer / two-sample Anderson–Darling / mark calibration;
`.` p > 0.10, `~` 0.05 ≤ p ≤ 0.10, `X` raw p < 0.05 before correction,
`-` condition excluded for that sampler):

| sampler | exp/f/cyc | exp/f/cmp | exp/r/cyc | exp/r/cmp | wbl/f/cyc | wbl/f/cmp | wbl/r/cyc | wbl/r/cmp | shf/f/cyc | shf/f/cmp | shf/r/cyc | shf/r/cmp |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| FirstReactionMethod (null control) | `X..` | `...` | `X~.` | `...` | `...` | `X..` | `.X.` | `...` | `..X` | `...` | `...` | `.~.` |
| NextReactionMethod | `X..` | `...` | `...` | `...` | `...` | `...` | `~..` | `...` | `...` | `.X.` | `...` | `...` |
| FirstToFireMethod | `..~` | `...` | `.~.` | `.~.` | `...` | `...` | `...` | `...` | `...` | `...` | `...` | `...` |
| DirectMethod(:keep,:array) | `...` | `.X.` | `...` | `...` | `-` | `-` | `-` | `-` | `-` | `-` | `-` | `-` |
| DirectMethod(:keep,:tree) | `...` | `..X` | `...` | `...` | `-` | `-` | `-` | `-` | `-` | `-` | `-` | `-` |
| DirectMethod(:remove,:array) | `...` | `...` | `...` | `...` | `-` | `-` | `-` | `-` | `-` | `-` | `-` | `-` |
| DirectMethod(:remove,:tree) | `...` | `...` | `...` | `..~` | `-` | `-` | `-` | `-` | `-` | `-` | `-` | `-` |
| RejectionMethod (RSSA) | `...` | `..X` | `...` | `...` | `-` | `-` | `-` | `-` | `-` | `-` | `-` | `-` |
| PartialPropensityMethod (PSSACR) | `...` | `.X.` | `...` | `X~.` | `-` | `-` | `-` | `-` | `-` | `-` | `-` | `-` |
| MultiSampler(all⇒CNR) | `...` | `...` | `~..` | `...` | `...` | `..X` | `...` | `X..` | `...` | `...` | `...` | `...` |
| MultiSampler(even⇒CNR, odd⇒CNR) | `...` | `...` | `..X` | `..X` | `...` | `...` | `...` | `X~.` | `X..` | `...` | `X..` | `~..` |
| MultiSampler(even⇒CNR, odd⇒FirstToFire) | `..~` | `...` | `...` | `...` | `...` | `...` | `...` | `...` | `..X` | `...` | `...` | `~..` |

Column labels abbreviate distribution family (`exp`onential, Wei`b`u`l`l,
`sh`i`f`ted) / memory (`f`orget, `r`emember) / graph (`cyc`le, `c`o`mp`lete).
Every `X` above is a pre-correction value that did not survive
false-discovery-rate control and did not exceed the null control's rate.

Two honest caveats. Under all-exponential conditions the two-sample test has
no discriminating power against the reference: inside the `next` contract,
the next-reaction samplers and the reference apply the identical truncated
inversion to identical uniforms, so those columns verify pathwise equivalence
rather than provide independent statistical evidence — power comes from the
Weibull and shifted conditions. And the matrix currently fixes the model
size and history length; population size, step count, and further graph
topologies are future axes.

## Reproducing the evidence

The full matrix, with per-cell raw and corrected p-values in CSV and a
markdown report, regenerates with:

```
julia --project=test test/gauntlet/run_matrix.jl 1000 20260705
```

which writes `results/gauntlet/matrix_report.md` and
`matrix_results.csv`. Any single cell reproduces with, for example:

```
julia --project=test test/gauntlet/run_matrix.jl cell "NextReactionMethod" "weibull/remember/complete" 1000 20260705
```

A fast smoke test of the gauntlet machinery runs in the ordinary test suite
so the harness cannot rot; the statistical runs themselves are deliberately
outside CI, because a p-value band is a diagnostic to be read, not a pass
bit to be gamed. The methodology — including why the two-sample query is
structured around the `next` contract — is documented in
`test/gauntlet/runner.jl`.
