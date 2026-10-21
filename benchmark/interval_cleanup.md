# Interval cleanup measurements

Run from the repository root:

```sh
julia --startup-file=no --project=benchmark benchmark/interval_cleanup.jl --bigfloat
julia --startup-file=no --project=benchmark benchmark/interval_kernels.jl
```

The baseline is `d028f67`, whose numerical source matches `7efc78c`. The candidate
is `5b7ac32`; `62215c5` separates the preceding code moves from calculation changes.
The scripts and shared harness are identical across the measured revisions.
Full output is in [measurements/interval_cleanup.txt](measurements/interval_cleanup.txt).

## Scope

The optional extension continues to use IntervalArithmetic for scalar arithmetic,
outward rounding, decorations and guarantee flags. There is no new interval type,
dependency or public API. Generic polynomial arithmetic is unchanged.

- Polynomial range and discarded-product kernels now live together. Scalar
  validation lives beside the coefficient hooks.
- The decorated and finite Float64 accumulators share discarded-product traversal.
  Each retains its existing scalar arithmetic and intermediate-decoration rules.
- Live models and compiled snapshots share enclosure calculations. Public queries
  retain their independent storage and different algebra-lifetime requirements.
- Scalar addition and subtraction avoid constructing a constant Taylor model.
  BigFloat endpoint precision, owned zeros, domain metadata and range certificates
  retain the previous semantics. Multiplication and division keep their existing
  model paths.

The shared accumulator specializes on its existing state: a decorated interval,
or a bare Float64 interval paired with its accumulated decoration. This removes
one duplicate monomial traversal without introducing a backend interface.

## Runtime and allocations

Fresh processes ran sequentially on Windows with an Intel Core i7-12700, one Julia
thread pinned to logical CPU 2 and AboveNormal priority. Dependencies were
IntervalArithmetic 1.0.12 and BenchmarkTools 1.8.0. BigFloat used 256-bit precision.
Both zero and 32 MiB multiplication-table budgets were measured. Timings below
are warmed minima from the shared harness; local clock and scheduling variation
remain visible, particularly in the dense-product measurements.

The scalar case is `m + 1//3`, where `m` is the six-coordinate uncertainty model
constructed by the script. The raw output also covers both operand orders and
subtraction. Construction evaluates the original expression using model inputs.

| Julia | Endpoints | Budget | Workload | Before → after, μs | Allocations |
|:--|:--|:--|:--|--:|--:|
| 1.13.2 | Float64 | 0 | Construction | 66.6 → 47.0 | 292 → 233 |
| 1.13.2 | Float64 | 32 MiB | Construction | 66.0 → 47.5 | 292 → 233 |
| 1.13.2 | BigFloat | 0 | Construction | 1663.4 → 1164.7 | 49856 → 42738 |
| 1.13.2 | BigFloat | 32 MiB | Construction | 1815.7 → 1140.5 | 49856 → 42738 |
| 1.13.2 | Float64 | 0 | Scalar sum | 7.18 → 4.67 | 10 → 5 |
| 1.13.2 | Float64 | 32 MiB | Scalar sum | 6.78 → 4.80 | 10 → 5 |
| 1.13.2 | BigFloat | 0 | Scalar sum | 170.5 → 105.2 | 5529 → 3866 |
| 1.13.2 | BigFloat | 32 MiB | Scalar sum | 190.8 → 107.6 | 5529 → 3866 |
| 1.10.12 | Float64 | 0 | Construction | 60.1 → 52.2 | 209 → 171 |
| 1.10.12 | Float64 | 32 MiB | Construction | 69.0 → 50.8 | 209 → 171 |
| 1.10.12 | BigFloat | 0 | Construction | 2154.7 → 1816.7 | 76478 → 65911 |
| 1.10.12 | BigFloat | 32 MiB | Construction | 2361.5 → 1868.8 | 76478 → 65911 |
| 1.10.12 | Float64 | 0 | Scalar sum | 5.52 → 6.03 | 6 → 3 |
| 1.10.12 | Float64 | 32 MiB | Scalar sum | 7.05 → 4.84 | 6 → 3 |
| 1.10.12 | BigFloat | 0 | Scalar sum | 225.3 → 177.6 | 8532 → 5984 |
| 1.10.12 | BigFloat | 32 MiB | Scalar sum | 230.8 → 175.1 | 8532 → 5984 |

Scalar-sum bytes fell from 2560 to 2224 on Julia 1.13.2 with Float64, and from
297768 to 210088 with BigFloat. Julia 1.10.12 fell from 2512 to 2192 and from
370064 to 261000 bytes respectively. Construction also allocated less: roughly
60.7 → 57.9 kB and 2.69 → 2.33 MB on Julia 1.13.2, and 61.3 → 58.5 kB and
3.34 → 2.89 MB on Julia 1.10.12. These are decimal units.

Live/snapshot queries retain their allocation counts and bytes. The Float64
full-domain queries allocate 672 bytes: two allocations on Julia 1.13.2 and one
on Julia 1.10.12. BigFloat full-domain queries retain 104776 bytes / 2180
allocations and 128720 bytes / 3172 allocations respectively. Subbox queries
also retain their counts and bytes.

Dense multiplication is measured separately to check the shared traversal:

| Julia | Budget | Before → after, μs | Bytes | Allocations |
|:--|:--|--:|--:|--:|
| 1.13.2 | 0 | 204.6 → 227.6 | 3287 → 3287 | 6 → 6 |
| 1.13.2 | 32 MiB | 213.3 → 259.8 | 3287 → 3287 | 6 → 6 |
| 1.10.12 | 0 | 283.8 → 260.8 | 3296 → 3296 | 3 → 3 |
| 1.10.12 | 32 MiB | 273.0 → 250.6 | 3296 → 3296 | 3 → 3 |

The dense-product timings increased on Julia 1.13.2 and decreased on Julia 1.10.12
in these paired runs. The shared traversal is a maintenance change; these data
do not establish a speedup. Both discarded-product accumulators still allocate
zero bytes. All printed range and remainder endpoints, decorations and guarantee
flags match the baseline exactly, including BigFloat and both table budgets.

## First-call costs

Package loading and harness compilation are excluded. Setup constructs a model
before the cleanup script's timers, so scalar arithmetic is already specialized;
the Julia 1.13.2 scalar rows report zero additional compilation. Later table-budget
rows also reuse earlier specializations. These are script-order costs, not isolated
function compilation measurements.

| Julia | Workload | First call before → after, ms | Compilation, ms |
|:--|:--|--:|--:|
| 1.13.2 | Float64 construction | 6.00 → 3.74 | 5.81 → 3.61 |
| 1.13.2 | BigFloat construction | 7.66 → 5.28 | 5.52 → 3.91 |
| 1.13.2 | Float64 live query | 4.41 → 3.09 | 4.39 → 3.07 |
| 1.13.2 | Float64 snapshot query | 37.84 → 26.17 | 37.82 → 26.16 |
| 1.13.2 | Dense model product | 4.41 → 4.39 | 4.11 → 4.09 |
| 1.13.2 | Direct discarded bound | 3.60 → 7.77 | 3.36 → 7.51 |
| 1.10.12 | Float64 construction | 4.18 → 4.09 | unavailable |
| 1.10.12 | BigFloat construction | 6.70 → 6.68 | unavailable |
| 1.10.12 | Float64 snapshot query | 37.17 → 34.07 | unavailable |
| 1.10.12 | Dense model product | 468.64 → 529.64 | unavailable |

Julia 1.10 does not expose compilation time through this harness. The direct-bound
entry incurs more first-call compilation on Julia 1.13.2; the shared traversal
does not improve every compilation path.

## Verification

At candidate `5b7ac32`, on Windows:

- Bounds-checked `Pkg.test()` passed on Julia 1.10.12 and 1.13.2: 102527 checks
  across 33 sequential subjects on each version. This includes fresh-process
  extension loading, independent rational bounds, both multiplication paths,
  precision, aliasing, snapshot ownership and interval metadata checks.
- Aqua passed on both versions, including a separate run with IntervalArithmetic,
  ForwardDiff, DiffEqBase, OrdinaryDiffEqCore and GenericLinearAlgebra loaded.
  No exceptions were added.
- Runic 1.11.1 passed on both versions for source, tests, extensions, examples,
  benchmarks and documentation scripts. Added Julia lines fit within 92 columns.
- The Julia 1.13.2 documentation build executed all 30 standalone tutorials.
  All 10 catalog/rendering checks and 6965 built catalog/link/asset checks passed.

The 4024 new scalar checks compare retained coefficients, remainders and range
certificates with binary model arithmetic. They include fixed and degraded axes,
point and positive-width scalars, independent output storage, and stored/current
BigFloat precision changes. Existing independent numerical oracles are retained.
Expected warnings arise from deliberate invalid-interval and rounding-policy
tests. There are no remaining validation failures. Ubuntu CI configuration is
unchanged; these recorded runs were local Windows checks.
