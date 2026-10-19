# bencher metric format; memory only when the peak could be measured
to_bmf <- function(res) {
    out <- list(latency = list(value = stats::median(res$times_ns),
                               lower_value = min(res$times_ns),
                               upper_value = max(res$times_ns)),
                cpu_time = list(value = res$cpu_ns))
    if (!is.null(res$peak_bytes))
        out$memory <- list(value = res$peak_bytes)
    out
}

# failed cases (NULL) are left out, never reported as a number
write_bmf_json <- function(results, path) {
    json <- jsonlite::toJSON(Filter(Negate(is.null), results), auto_unbox = TRUE,
                             digits = NA, pretty = TRUE)
    writeLines(json, path)
}
