# body of a callr child: self-contained, only `::` calls and local helpers
measure_case <- function(fn, method, margin, args, data_file, iterations, threads) {
    options(proxyC.threads = threads)
    status_kb <- function(key) {
        line <- grep(paste0("^", key, ":"), readLines("/proc/self/status"), value = TRUE)
        as.numeric(gsub("[^0-9]", "", line))
    }
    d <- readRDS(data_file)
    f <- getExportedValue("proxyC", fn)

    # candidate arguments, filtered by the formals (crossprod and colSds drop margin, method, use_nan)
    cand <- c(list(x = d$x, y = d$y, margin = margin, method = method, use_nan = FALSE), args)
    if (isTRUE(cand[["mask"]]))
        cand$mask <- proxyC::mask(d$lx, d$ly)
    a <- cand[names(cand) %in% names(formals(f))]

    # warm up bench (loads its dependencies) and reset VmHWM to the current RSS, so CPU time and
    # peak cover only the timed calls
    bench::mark(NULL, iterations = 1L, check = FALSE, memory = FALSE)
    invisible(gc(reset = TRUE))
    # peak memory needs /proc (Linux); elsewhere only latency and cpu_time are reported
    reset <- tryCatch({
        rss0 <- status_kb("VmRSS")
        writeLines("5", "/proc/self/clear_refs")
        TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    p0 <- proc.time()
    bm <- bench::mark(do.call(f, a), iterations = iterations, check = FALSE,
                      memory = FALSE, filter_gc = FALSE)
    p1 <- proc.time() - p0

    # user + sys time of the process includes TBB and BLAS worker threads
    list(times_ns = as.numeric(bm$time[[1]]) * 1e9,
         cpu_ns = sum(p1[c("user.self", "sys.self")]) / iterations * 1e9,
         peak_bytes = if (reset) max(status_kb("VmHWM") - rss0, 1) * 1024)
}

# one fresh process per case: isolates crashes and keeps peaks from bleeding between cases
run_case <- function(row, data_file, iterations, threads, timeout = 300) {
    env <- c(callr::rcmd_safe_env(),
             if (threads == 1) c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1"))
    args <- list(fn = row$fn, method = row$method, margin = row$margin, args = row$args[[1]],
                 data_file = data_file, iterations = iterations, threads = threads)
    # wait() polls with a timeout; callr::r(timeout =) relies on the process start time, which
    # processx can misread on hosts whose boot-time clock drifted (e.g. after suspend)
    tryCatch({
        p <- callr::r_bg(measure_case, args = args, libpath = .libPaths(), env = env,
                         stdout = NULL, stderr = NULL)
        p$wait(timeout * 1000)
        if (p$is_alive()) {
            p$kill()
            stop(sprintf("timed out after %g s", timeout), call. = FALSE)
        }
        p$get_result()
    }, error = function(e) {
        msg <- conditionMessage(if (inherits(e$parent, "condition")) e$parent else e)
        cat(sprintf("::warning::FAILED %s: %s\n", row$name, gsub("\\s+", " ", trimws(msg))))
        NULL
    })
}
