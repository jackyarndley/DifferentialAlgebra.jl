# Interval kernel experiment

Run `julia --startup-file=no --project=benchmark benchmark/interval_kernels.jl`
from the repository root. Add `--profile` for CPU samples and allocation stacks.
The shared measurement harness records first-call compilation separately from
warmed runtime. The workloads use six variables, order three, and both zero and
32 MiB multiplication-table budgets. These short kernels request 1,000 samples
with a one-second trial budget to reduce clock and scheduling variation.

## Scope and rounding

The private Float64 kernels specialize two operations in the optional interval
extension. They add no public type, dependency or configuration switch.

- Normalized coordinates `{0}` and `[-1,1]` have exact positive integer powers.
  Returning the appropriate decorated zero, odd power or absolute value avoids
  the general scalar-power implementation. Other coordinates and BigFloat use
  `IntervalArithmetic.pown`.
- Discarded products on normalized unit boxes multiply coefficient intervals
  with `BareInterval` arithmetic. Multiplication by the monomial range `{0}`,
  `[0,1]` or `[-1,1]` then only selects endpoints. Both coefficient products and
  accumulation retain IntervalArithmetic's outward-rounded primitives.

The second kernel first checks for finite, nonempty, guaranteed Float64
coefficients with supported decorations. The decorated implementation remains
the fallback and a direct benchmark reference. General boxes retain their
existing power-based bound. Metadata is accumulated from participating inputs
and every intermediate result, including an overflowing product subsequently
multiplied by a fixed zero coordinate. This is essential even when the final
endpoints are finite. Global rounding, tolerance and cutoff are unchanged.

## Profile findings

The baseline profile at `f00e583` placed discarded-term bounding in 109 of 127
active dense-product CPU samples. The uncertainty expression allocated 1,120
objects; 828 allocation stacks originated in scalar powers inside `power_cache`.
The exact normalized powers remove all 828 of those allocations. The remaining
allocations primarily hold coefficient buffers, power matrices and model values.
CPU samples identify a hotspot; they are not precise elapsed-time attribution.

## Verification and measurements

The focused tests use exact rational ordered-product bounds, including an exact
rational check on rounded subnormal and overflowing products. They also compare
endpoints, decorations and guarantee flags with the decorated implementation,
covering fixed axes, degraded zeros, NG fallbacks and both table budgets.

The baseline for paired measurements is `cec303d`, which adds only the profiling
script to `f00e583`. The candidate kernel revision is `daec28a`. Use the current
script and matching dependency environments with each package revision. Full
precision endpoint diagnostics allow the resulting enclosures and remainders to
be compared as well as their widths.

Current output is in [measurements/interval_kernels.txt](measurements/interval_kernels.txt).

Paired runs used driver `c81d2fb`, IntervalArithmetic 1.0.12 and BenchmarkTools
1.8.0 on Windows, on an Intel Core i7-12700. Fresh processes ran sequentially
after validation, with one Julia thread pinned to logical CPU 2 and AboveNormal
process priority. Local timings still vary; the values below are the warmed
minima recorded by the shared harness.

| Julia | Table budget | Workload | Before, μs | After, μs | Allocations before → after |
|:--|:--|:--|--:|--:|--:|
| 1.13.2 | 0 | Normalized powers | 6.08 | 0.300 | 182 → 2 |
| 1.13.2 | 32 MiB | Normalized powers | 4.85 | 0.256 | 182 → 2 |
| 1.13.2 | 0 | Dense model product | 318.4 | 269.8 | 186 → 6 |
| 1.13.2 | 32 MiB | Dense model product | 268.6 | 245.1 | 186 → 6 |
| 1.13.2 | 0 | Uncertainty expression | 81.2 | 62.7 | 1,120 → 292 |
| 1.13.2 | 32 MiB | Uncertainty expression | 67.1 | 51.4 | 1,120 → 292 |
| 1.10.12 | 0 | Normalized powers | 4.53 | 0.369 | 361 → 1 |
| 1.10.12 | 32 MiB | Normalized powers | 4.67 | 0.245 | 361 → 1 |
| 1.10.12 | 0 | Dense model product | 453.6 | 395.2 | 363 → 3 |
| 1.10.12 | 32 MiB | Dense model product | 289.4 | 251.0 | 363 → 3 |
| 1.10.12 | 0 | Uncertainty expression | 73.1 | 51.5 | 1,865 → 209 |
| 1.10.12 | 32 MiB | Uncertainty expression | 92.6 | 52.2 | 1,865 → 209 |

The same-process discarded-product comparison isolates that kernel from the
power-cache change. Decorated reference → private kernel timings were
230.6 → 161.7 μs and 280.2 → 211.0 μs on Julia 1.13.2, and 389.2 → 238.1 μs
and 242.2 → 212.3 μs on Julia 1.10.12, for the two respective table budgets.
Both implementations allocate zero bytes in this workload.

On Julia 1.13.2, dense model multiplication dropped from 14,807 to 3,287 bytes;
the uncertainty expression dropped from about 113.7 to 60.7 kB. Julia 1.10.12
dropped from 17,696 to 3,296 bytes and from 127,504 to 61,264 bytes respectively.
Every printed enclosure and remainder endpoint, decoration and guarantee flag
matched exactly between revisions, for both budgets and both Julia versions.

First-call costs are separate from those warmed timings. The dense product is
measured before its subkernel benchmarks; subsequent workloads reuse compiled
methods. Package loading and harness compilation are excluded.

| Julia | Workload | First call before → after, ms | Compilation before → after, ms |
|:--|:--|--:|--:|
| 1.13.2 | Dense model product | 3.86 → 6.70 | 3.50 → 6.29 |
| 1.13.2 | Uncertainty expression | 6.17 → 7.49 | 5.87 → 3.79 |
| 1.10.12 | Dense model product | 775.43 → 573.60 | unavailable |
| 1.10.12 | Uncertainty expression | 637.72 → 646.13 | unavailable |

## Validation

- Bounds-checked `Pkg.test()` passed on Julia 1.10.12: 98,503 checks across
  33 subjects, including the new kernel tests.
- Bounds-checked `Pkg.test()` passed on Julia 1.13.2: 96,047 existing checks;
  the 2,456 new kernel checks also passed in a separate bounds-checked process.
- Aqua passed on both versions with all optional integrations loaded, without
  exceptions. Runic checks passed on both versions.
- The Julia 1.13.2 documentation build executed all 30 standalone Literate
  tutorials, passed 10 Literate rendering/catalog checks, and passed 6,965
  built-catalog, local-link and asset checks.

There were no remaining test, documentation or formatting failures. These
results measure DA-specific kernels over IntervalArithmetic's outward-rounded
primitives.
