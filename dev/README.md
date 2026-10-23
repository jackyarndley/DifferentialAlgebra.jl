# Function inventory

Generate the complete source inventory without loading the package or its optional
dependencies:

```sh
julia --startup-file=no dev/inventory.jl dev/functions.tsv
```

`functions.tsv` lists every function definition parsed from Julia files in `src/`,
`ext/`, `examples/`, `test/`, `docs/`, `benchmark/` and `dev/`. It records source
locations, signatures, lexical context and a hash of each body with line numbers
removed. Nested functions, callable objects, anonymous functions and empty generic
declarations are included. Generated build output is excluded.

An `@eval` loop is recorded as a **template**, with its loop in the context column;
the script does not execute code or count Julia's generated keyword/default-argument
wrappers. Definitions inside other quoted expressions are also marked as templates.
Anonymous functions use the nearest line recorded by Julia's parser.

Repeated names usually represent multiple dispatch. Equal body hashes identify
review candidates, not automatic deletions: constructors, forwarding methods and
scalar-type specializations can intentionally have identical bodies. Numerical
paths with different rounding, filtering or conditioning contracts must remain
distinct unless a shared implementation preserves those contracts.
