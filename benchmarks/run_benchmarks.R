# continuous benchmarks for proxyC, written in bencher metric format (see README.md);
# always measures the installed proxyC, so run `R CMD INSTALL .` first
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
for (f in list.files(file.path(dirname(script), "lib"), "\\.R$", full.names = TRUE))
    source(f)

opts <- parse_args(list(iterations = 3L, tier = "all", filter = "", threads = -1L, timeout = 300,
                        output = file.path(dirname(script), "results.json"), list = FALSE))

# select cases
grid <- build_grid()
if (opts$tier != "all")
    grid <- grid[grid$tier %in% strsplit(opts$tier, ",")[[1]], ]
if (nzchar(opts$filter))
    grid <- grid[grepl(opts$filter, grid$name), ]
if (opts$list) {
    options(width = 200)
    print(grid[c("name", "fn", "method", "variant", "tier", "input", "margin")], row.names = FALSE)
    cat(nrow(grid), "cases\n")
    quit(status = 0)
}
if (nrow(grid) == 0)
    stop("no cases match --tier/--filter", call. = FALSE)

# self-describing log, which also fails fast when proxyC is not installed
cat(sprintf("proxyC %s (%s)\nBLAS: %s\nTBB: %s, max threads: %d, proxyC.threads: %d, iterations: %d\n\n",
            packageVersion("proxyC"), find.package("proxyC"), sessionInfo()$BLAS,
            proxyC:::cpp_tbb_enabled(), proxyC:::cpp_get_max_thread(), opts$threads, opts$iterations))

# one data file per data_key, one child process per case
results <- list()
for (key in unique(grid$data_key)) {
    rows <- grid[grid$data_key == key, ]
    data_file <- ensure_data(rows[1, ])
    for (i in seq_len(nrow(rows))) {
        res <- run_case(rows[i, ], data_file, opts$iterations, opts$threads, opts$timeout)
        results[rows$name[i]] <- list(if (!is.null(res)) to_bmf(res))
        if (!is.null(res))
            cat(sprintf("[%3d/%d] %-40s %10.1f ms  cpu/latency %5.2f  peak %8.1f MB\n",
                        length(results), nrow(grid), rows$name[i], stats::median(res$times_ns) / 1e6,
                        res$cpu_ns / stats::median(res$times_ns),
                        if (is.null(res$peak_bytes)) NA else res$peak_bytes / 2^20))
    }
    unlink(data_file)
}

# fail when nothing or more than 20% of the cases failed; the JSON is written either way
write_bmf_json(results, opts$output)
failed <- sum(vapply(results, is.null, logical(1)))
cat(sprintf("\n%d of %d cases failed; results written to %s\n", failed, length(results), opts$output))
if (failed == length(results) || failed > 0.2 * length(results))
    quit(status = 1)
