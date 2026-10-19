# proxyC benchmarks

Performance benchmarks for proxyC's `simil()`, `dist()` and matrix helpers (`crossprod()`, `colSds()`, `colZeros()`, ...), measured with [`bench::mark()`](https://bench.r-lib.org/) and tracked over time on [bencher.dev](https://bencher.dev/perf/proxyc).

Every push to `master` is benchmarked on three testbeds (single-threaded, OpenBLAS, reference BLAS), and every pull request that touches the package code is compared against `master`, with one bencher comment per testbed on the PR. There are three testbeds because a change can be faster on one BLAS/threading setup and slower on another: [#74] sped up dense input with a single-threaded BLAS but slowed it down with a multithreaded one.

## Prerequisites

- **R packages**: `bench`, `callr` and `jsonlite`, on top of proxyC's own dependencies.

  ```r
  install.packages(c("bench", "callr", "jsonlite"))
  ```

- **Linux** for the memory measure: each case reads its peak memory from `/proc/self/status`. On macOS or Windows the cases still run, but the JSON has no `memory` entries.
- **An installed proxyC**: every case runs in a fresh R process started with `callr`, which loads proxyC from the library, so `devtools::load_all()` has no effect. Install the working copy from the repository root before each run with `R CMD INSTALL .`. The log header prints the version and path of the proxyC that is used.
- **TBB** (optional): to benchmark proxyC's parallel code, install TBB before building proxyC (`libtbb-dev` on Debian/Ubuntu, as in CI). The log header shows `TBB: TRUE` when it is enabled.

## Running locally

From the repository root:

```bash
# install the working copy; the benchmarks always use the installed proxyC
R CMD INSTALL .

# list the cases without running them (doesn't need proxyC)
Rscript benchmarks/run_benchmarks.R --list
Rscript benchmarks/run_benchmarks.R --list --tier pair --filter jaccard

# quick smoke test of two cases
Rscript benchmarks/run_benchmarks.R --filter 'simil_cosine_(sparse|dense)_m2$' --iterations 1

# one thread vs all threads
Rscript benchmarks/run_benchmarks.R --filter 'simil_cosine_(sparse|dense)_m2$' --threads 1 --output benchmarks/results-t1.json
Rscript benchmarks/run_benchmarks.R --filter 'simil_cosine_(sparse|dense)_m2$' --threads -1 --output benchmarks/results-tall.json

# all cases, as one CI testbed does
Rscript benchmarks/run_benchmarks.R
```

The run starts with a header (proxyC version and path, BLAS, whether TBB is enabled, thread settings), then prints one line per case with the median latency, the ratio of CPU time to latency and the peak memory:

```text
[  1/2] simil_cosine_dense_m2                        2506.7 ms  cpu/latency  0.99  peak    632.8 MB
[  2/2] simil_cosine_sparse_m2                       2616.0 ms  cpu/latency  0.99  peak    479.7 MB

0 of 2 cases failed; results written to benchmarks/results.json
```

Files matching `benchmarks/results*.json` are ignored by git.

### CLI options

Options take `--key value` or `--key=value`; unknown options are an error. `--tier` and `--filter` can be combined.

| Option | Default | Description |
|---|---|---|
| `--iterations N` | `3` | Timed calls per case. |
| `--tier LIST` | `all` | Comma-separated tiers to run, e.g. `large,pair` (see [Tiers](#tiers)). |
| `--filter REGEX` | none | Only run cases whose name matches this regular expression. |
| `--threads N` | `-1` | Sets `options(proxyC.threads)` for each case; `-1` uses all threads. `1` also sets `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1` for the case, so a multithreaded BLAS runs on one thread too. |
| `--timeout SECONDS` | `300` | Wall-clock limit per case; a case that runs longer is killed and counts as failed. |
| `--output FILE` | `benchmarks/results.json` | Where to write the results. The default is next to `run_benchmarks.R`, wherever you run it from. |
| `--list` | off | Print the selected cases and their count, then exit without running anything. |

## Output format

Results are written as [Bencher Metric Format (BMF)](https://bencher.dev/docs/reference/bencher-metric-format/) JSON, one entry per case:

```json
{
  "simil_cosine_dense_m2": {
    "latency": {
      "value": 2506651133.06604,
      "lower_value": 2480675253.08579,
      "upper_value": 2517504118.1501
    },
    "cpu_time": {
      "value": 2487666666.66667
    },
    "memory": {
      "value": 663523328
    }
  }
}
```

| Measure | Unit | Value |
|---|---|---|
| `latency` | nanoseconds | Wall time per call: `value` is the median, `lower_value` and `upper_value` the minimum and maximum of the raw `bench::mark()` iterations (no GC filtering). |
| `cpu_time` | nanoseconds | User + system CPU time of the process during the timed calls, divided by the number of iterations. Includes the time spent in TBB and BLAS worker threads. |
| `memory` | bytes | Peak resident set size (RSS) during the timed calls, minus the RSS just before them. |

`latency` and `memory` are built-in bencher measures; `cpu_time` is a custom measure that bencher creates on the first upload.

### How a case is measured

Each case runs in a fresh R process (`callr::r_bg()`), which:

1. reads the case's input data and builds the arguments (see [Data](#data));
2. warms up `bench` with a dummy call and runs `gc()`;
3. records `VmRSS` from `/proc/self/status` and writes `5` to `/proc/self/clear_refs`, which resets the peak (`VmHWM`) to the current RSS;
4. runs `bench::mark()` for `--iterations` calls, with `proc.time()` taken before and after;
5. reads `VmHWM` and returns the raw timings, the CPU time and `VmHWM - VmRSS`.

A fresh process per case means peaks don't carry over between cases and a crash only takes down that case. The input data is loaded before the baseline, so `memory` is what the calls allocate on top of their inputs, including copies such as proxyC converting dense input to a sparse matrix. Unlike the `mem_alloc` column of `bench::mark()`, which only counts memory allocated by R, the RSS peak also includes what proxyC allocates in C++ (Armadillo, TBB, BLAS). If `/proc/self/clear_refs` can't be written, `memory` is left out for that case and `latency` and `cpu_time` are still reported.

### Reading `cpu_time`

`cpu_time / latency` is roughly the number of cores kept busy during a call (the console prints it as `cpu/latency`):

- about 1: one thread does the work, as expected on the `single-thread` testbed;
- close to the runner's core count: the work is spread over all cores;
- `cpu_time` going up while `latency` stays the same or gets worse: threads are competing or spinning rather than working, e.g. TBB workers that each start a multithreaded BLAS (the [#74] regression).

## Benchmark cases

All cases come from one table built by `build_grid()` in `lib/grid.R`, and one generic runner executes them; `--list` prints the table. There are currently 124 cases: 44 on `large`, 16 on `xsparse` and 64 on `pair`.

### Tiers

`TIERS` defines the data shapes. `nf` is the number of features, `nx` and `ny` the number of items in `x` and `y` that are compared with each other.

| Tier | `nf` | `nx` | `ny` | Density | Inputs | Runs |
|---|---|---|---|---|---|---|
| `large` | 50,000 | 100 | 200 | 0.9 | sparse, dense | `LINEAR` methods; `rank`, `min_simil` and `mask` variants; `IMPLIED` functions |
| `xsparse` | 50,000 | 2,000 | none (`y = NULL`) | 0.01 | sparse | `LINEAR` methods; `IMPLIED` functions |
| `pair` | 10,000 | 100 | 100 | 0.3 | sparse, dense | `PAIR` methods; `diag` variant |

`large` uses the shapes from [#74] with half the features, to keep CI time down. `xsparse` is the typical very sparse use case, e.g. a document-feature matrix compared with itself. `pair` is smaller because the methods computed by `cpp_pair()` are much slower than those computed by `cpp_linear()`.

With `margin = 1` the items are rows (`x` is `nx × nf`), with `margin = 2` they are columns (`x` is `nf × nx`), so both margins compare the same number of items.

### Methods and functions

| Object | Contents | Used for |
|---|---|---|
| `LINEAR` | `simil()`: cosine, correlation, dice, edice; `dist()`: euclidean | Methods computed by `cpp_linear()`; run on `large` and `xsparse`. |
| `PAIR` | `simil()`: jaccard, ejaccard, fjaccard, hamann, faith, simple matching; `dist()`: chisquared, kullback, jeffreys, jensen, manhattan, maximum, canberra, minkowski, hamming | Methods computed by `cpp_pair()`; run on `pair`. |
| `BINARY` | jaccard, hamann, faith, simple matching, hamming | Get 0/1 data. |
| `SMOOTHED` | chisquared, kullback, jeffreys, jensen | Called with `smooth = 1`. |
| `IMPLIED` | `crossprod()` (2), `tcrossprod()` (1), `colSds()` (2), `rowSds()` (1), `colZeros()` (2), `rowZeros()` (1) | Functions whose margin follows from their name; run on `large` and `xsparse`, at that margin only. |

Together, `LINEAR` and `PAIR` cover every method of `simil()` and `dist()`.

- **Why `smooth = 1`**: with the default `smooth = 0`, chisquared, kullback, jeffreys and jensen return `NaN` as soon as either vector contains a zero, without computing anything (`src/pair.cpp`). Every vector in the benchmark data contains zeros, so these cases would only time that early return.
- **Why 0/1 data**: these methods are meant for presence/absence data. `simil()` converts the input to logical for jaccard, hamann, faith and simple matching, and hamming counts unequal entries, which on continuous values is nearly all of them. Their cases get the same random sparsity pattern with every non-zero set to 1.

### Variants

`VARIANTS` holds extra arguments that are added to a case:

| Variant | Arguments | Cases |
|---|---|---|
| `base` | none | all other cases |
| `rank` | `rank = 10` | `simil(method = "cosine")` on `large` |
| `min_simil` | `min_simil = 0.5` | `simil(method = "cosine")` on `large` |
| `mask` | `mask = proxyC::mask(lx, ly)` | `simil(method = "cosine")` on `large`; the mask is built from random labels (5 groups) before timing |
| `diag` | `diag = TRUE` | `simil(method = "cosine")` on `pair`, because `diag = TRUE` sends cosine through `cpp_pair()` |

### Naming

Case names follow `fn_method[_variant]_input_mN`, with spaces replaced by `_`, the method left out for `IMPLIED` functions and the variant left out for `base`:

- `simil_cosine_dense_m1`
- `simil_cosine_rank_sparse_m2`
- `simil_simple_matching_dense_m1`
- `colSds_sparse_m2`

When a name is already taken by an earlier tier (in the order of `TIERS`), the tier is appended, e.g. `simil_cosine_sparse_m2_xsparse`. These names are the benchmark names on bencher, so renaming a case starts a new history there.

### Data

`make_data()` in `lib/data.R` generates the input with `Matrix::rsparsematrix()`, with non-zero values drawn uniformly from [0.1, 1] (positive, so kullback, chisquared and jensen are defined and fjaccard stays within [0, 1]), and converts it to a dense matrix for `dense` cases. The seed is fixed at 42, so every commit and testbed benchmarks the same data.

Cases with the same tier, input, margin and binary flag (`data_key`) share their data. The runner generates each group's data once, saves it as an uncompressed `.rds` in `tempdir()`, runs all cases of the group against it and deletes it afterwards.

Each case is called with `x`, `y`, `margin`, `method` and `use_nan = FALSE` plus its variant and smoothing arguments, keeping only the arguments that the function has (`formals()`). So `colSds()` gets just `x` and `crossprod()` gets `x` and `y`, without a per-function list. `use_nan = FALSE` skips proxyC's check for zero-variance or all-zero vectors and its warning.

### Adding a case

- A new method on an existing tier: add it to `LINEAR` or `PAIR`, and to `BINARY` or `SMOOTHED` if it needs 0/1 data or smoothing.
- An existing method or function on another tier or variant: add a `spec(fn, method, tier, variant)` row to the `rbind()` in `build_grid()`. For a function without a method, use `NA` as method, and add it to `IMPLIED` if its margin is fixed.
- New arguments: add an entry to `VARIANTS` and a `spec()` row that uses it. Arguments the function has are passed as they are; an argument that has to be built from the case's data (like `mask`) needs a line in `measure_case()`.
- A new data shape: add an entry to `TIERS` and `spec()` rows that use it.

For example, to benchmark `drop0 = TRUE` for cosine on the `large` tier:

```r
VARIANTS <- list(..., drop0 = list(drop0 = TRUE))

build_grid <- function() {
    grid <- rbind(
        ...,
        spec("simil", "cosine", "large", "drop0")
    )
    ...
```

Then check the names with `--list` and run the new cases with `--filter`. A PR that adds a case has no `master` history to compare it with yet; it gets one once the change is merged.

## Failed cases

If a case errors, crashes or exceeds `--timeout`, the runner prints a GitHub Actions warning annotation and leaves the case out of the JSON, so bencher shows a missing data point rather than a wrong number:

```text
::warning::FAILED dist_kullback_dense_m1: <error message>
```

The JSON is always written. The runner exits with status 1 when more than 20 % of the cases failed, which fails the CI job, so that testbed's results are not uploaded. Fewer failures only produce the warnings.

The output of the child process is discarded, so the warning is all you get. To dig further, reproduce the case in an R session:

```r
source("benchmarks/lib/grid.R")
source("benchmarks/lib/data.R")
d <- make_data("pair", "dense", margin = 1, binary = FALSE)
proxyC::dist(d$x, d$y, margin = 1, method = "kullback", smooth = 1, use_nan = FALSE)
```

## Testbeds

Each testbed is a separate CI job on a GitHub-hosted `ubuntu-latest` runner, all with proxyC built with TBB:

| Testbed | BLAS | `--threads` | `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` | What it shows |
|---|---|---|---|---|
| `single-thread` | reference BLAS | `1` | `1` | Changes in the algorithms themselves, with the least noise. |
| `openblas` | OpenBLAS | `-1` | unset | How proxyC's TBB threads interact with a multithreaded BLAS. |
| `refblas` | reference BLAS | `-1` | unset | How proxyC's own TBB parallelism scales. |

The composite action `.github/actions/bench-run/action.yaml`, shared by `bench-base.yaml` and `bench-pr-run.yaml`, sets them up:

1. BLAS: apt installs `libtbb-dev`, `libblas-dev` and `liblapack-dev`. The `openblas` testbed also installs `libopenblas-dev`, which has a higher `update-alternatives` priority, so `libblas.so.3` and `liblapack.so.3` then point at OpenBLAS. The other testbeds pin the reference BLAS and LAPACK with `update-alternatives --set`, in case the runner image ships OpenBLAS. The resulting alternatives are printed in the log.
2. R and proxyC: `r-lib/actions/setup-r` installs R, and `r-lib/actions/setup-r-dependencies` installs proxyC from the checkout (`local::.`) together with `bench`, `callr` and `jsonlite`. This runs after apt, so proxyC's `configure` finds TBB.
3. Threads: `single-thread` exports `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1` and runs with `--threads 1`; the others run with `--threads -1`.
4. Check testbed: prints the testbed, proxyC version and build, BLAS and LAPACK from `sessionInfo()`, whether TBB is enabled, the maximum number of threads, `parallel::detectCores()` and the thread variables. It fails the job unless TBB is enabled and OpenBLAS is loaded on the `openblas` testbed and only there.
5. Run: `Rscript benchmarks/run_benchmarks.R --iterations 3 --threads $THREADS --output results-<testbed>.json`.

To get close to a testbed locally, pick the BLAS through your system (e.g. `sudo update-alternatives --config libblas.so.3-x86_64-linux-gnu` on Debian/Ubuntu, or `FLEXIBLAS=NETLIB` for the reference BLAS on Fedora) and check the `BLAS:` line in the log header. `--threads 1` matches the thread settings of `single-thread`.

## CI workflows

The setup follows the benchmarks of [anndataR](https://github.com/scverse/anndataR): three workflows, so that pull requests from forks can be benchmarked without their code getting access to secrets.

### `bench-base.yaml`: Benchmarks (base)

- Trigger: push to `master`, or by hand via `workflow_dispatch`.
- Runs the three testbeds as a matrix (`fail-fast: false`, 60 minutes per job) and uploads each result to bencher on branch `master`.
- Defines the thresholds: a t-test on `latency`, `memory` and `cpu_time` with an upper boundary of 0.99, over the last 5 to 64 reports (`--threshold-min-sample-size 5`, `--threshold-max-sample-size 64`). `--thresholds-reset` removes any other thresholds for that branch and testbed.
- `--err` fails the job when a threshold raises an alert.
- A newer push to `master` cancels a run that is still in progress.

### `bench-pr-run.yaml`: Benchmarks (PR run)

- Trigger: a pull request is opened, reopened or updated, and it changes `R/`, `src/`, `configure`, `configure.ac`, `DESCRIPTION`, `benchmarks/` or the benchmark workflows and action.
- Runs the same three testbeds on the PR's code and uploads `results-<testbed>.json` as artifact `results-<testbed>`. A separate small job uploads the PR event payload as artifact `event`.
- Fork-safe: it has a read-only token (`contents: read`) and no secrets, and the checkout doesn't keep credentials.
- A new push to the PR cancels the run that is still in progress.

### `bench-pr-track.yaml`: Benchmarks (PR track)

- Trigger: `workflow_run`, when Benchmarks (PR run) completes. It is skipped when that run was cancelled or wasn't for a pull request.
- Runs the version of the workflow on `master`, with access to `BENCHER_API_TOKEN` and a token that can comment on PRs. It never checks out or runs PR code: the artifacts are untrusted data. The PR event is validated (PR number, SHAs, refs), cross-checked against the triggering run and the live PR, and its values reach the shell only through environment variables.
- Uploads each testbed's results to bencher on the PR's branch, with the base branch at the PR's base commit as start point and its thresholds cloned. Branches from forks are prefixed with the fork owner (`owner/branch`), so they can't collide with branches of this repository.
- Each testbed gets its own PR comment. If one testbed raises an alert, the others are still uploaded and the job fails at the end. If a testbed of the PR run failed, the ones that produced results are still tracked.

### Only on `master`

GitHub only runs `workflow_run` workflows from the default branch, and only offers "Run workflow" for `workflow_dispatch` once the workflow is on the default branch. Until these files are on `master`, PRs are benchmarked by Benchmarks (PR run) but nothing reaches bencher, and Benchmarks (base) can't be started by hand. For the same reason, changes to `bench-pr-track.yaml` only take effect after they are merged.

### Secrets

| Secret | Description |
|---|---|
| `BENCHER_API_TOKEN` | API token of a [bencher.dev](https://bencher.dev) user with access to the `proxyc` project. Used by `bench-base.yaml` and `bench-pr-track.yaml`. |

`GITHUB_TOKEN` is provided automatically.

## First-time setup

- [ ] The bencher project `proxyc` exists: <https://bencher.dev/perf/proxyc>.
- [ ] The repository secret `BENCHER_API_TOKEN` exists (Settings → Secrets and variables → Actions).
- [ ] The workflows are merged into `master` (see [Only on `master`](#only-on-master)).
- [ ] `master` has a report on all three testbeds, so PRs have a start point. The merge itself triggers Benchmarks (base) through the push to `master`; if that run didn't complete on every testbed, run it once by hand: Actions → Benchmarks (base) → Run workflow on `master`, or `gh workflow run bench-base.yaml --ref master`.
- [ ] Collaborators are invited to the bencher organization, so they can manage thresholds and alerts.

Thresholds only raise alerts once `master` has at least 5 reports per testbed. Until then, PR results are uploaded and commented on without alerts. To get there sooner, run Benchmarks (base) a few more times by hand.

## Interpreting alerts

bencher compares each PR result with the recent `master` results of the same benchmark, testbed and measure, and raises an alert when it lies above the upper boundary of a t-test at 0.99. Only upper boundaries are set, so an alert means slower, more CPU time or more memory; improvements never raise one. An alert is listed in the PR comment and fails Benchmarks (PR track) (`--err`).

- `memory` is nearly deterministic (fixed data, fresh process per case), so a memory alert is almost always real.
- `single-thread` is the most stable testbed. Latency on `openblas` and `refblas` depends more on the shared runner and can vary by 10–20 % between runs.
- An alert on all three testbeds points at the algorithm. An alert only on `openblas` and/or `refblas`, often with a much larger `cpu_time`, points at threading, as in [#74].
- A case missing from the report failed: look for `::warning::FAILED` in the log of Benchmarks (PR run).

If an alert looks like noise, re-run all jobs of Benchmarks (PR run); Benchmarks (PR track) runs again once it completes. Alerts can also be dismissed on bencher. If a regression is expected (e.g. a slower but correct fix), merge anyway: the `master` reports after the merge become the new baseline.

## Runtime

Each case starts a fresh R process, loads its data and runs 3 iterations. Each testbed is one job with a 60-minute timeout, the three testbeds run in parallel, and each case has a 300 s timeout. If CI gets too slow, these are the levers:

- Fewer iterations: the `iterations` input in `bench-base.yaml` and `bench-pr-run.yaml` (keep both the same). `latency` still gets a min/median/max with fewer iterations, but the median gets noisier.
- Smaller data: `nf`, `nx` and `ny` in `TIERS`, e.g. `nf` of `pair`.
- Lower density, e.g. of `large`.
- Fewer inputs: drop `dense` from a tier's `inputs`, or move the slowest methods to a tier of their own with `inputs = "sparse"`.

Changing data sizes, density or the seed changes every affected number, so do it in a PR of its own. Expect a step in the history on bencher, and alerts on that PR if the numbers go up. After the merge, PRs are compared with a mix of old and new `master` reports until the old ones drop out of the 64-report window.

## File structure

```text
benchmarks/
├── README.md            # this file
├── run_benchmarks.R     # CLI entry point: selects cases, runs them, writes the JSON
├── .gitignore           # ignores results*.json
└── lib/
    ├── cli.R            # parse_args()
    ├── grid.R           # TIERS, VARIANTS, method lists, build_grid()
    ├── data.R           # make_data(), ensure_data()
    ├── measure.R        # measure_case() (runs in the child process), run_case()
    └── bmf.R            # to_bmf(), write_bmf_json()

.github/
├── actions/bench-run/
│   └── action.yaml          # testbed setup and benchmark run, shared by base and PR run
└── workflows/
    ├── bench-base.yaml      # Benchmarks (base)
    ├── bench-pr-run.yaml    # Benchmarks (PR run)
    └── bench-pr-track.yaml  # Benchmarks (PR track)
```

`benchmarks/` is listed in `.Rbuildignore`, so it is not part of the package build.

[#74]: https://github.com/koheiw/proxyC/pull/74
