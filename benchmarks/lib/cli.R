# parse `--key value`, `--key=value` and bare `--flag` (logical defaults) into a typed list
parse_args <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
    opts <- defaults
    i <- 1L
    while (i <= length(args)) {
        if (!startsWith(args[i], "--"))
            stop("unexpected argument '", args[i], "'", call. = FALSE)
        key <- sub("^--([^=]*).*$", "\\1", args[i])
        if (!key %in% names(defaults))
            stop("unknown option --", key, call. = FALSE)
        if (grepl("=", args[i])) {
            value <- sub("^[^=]*=", "", args[i])
        } else if (is.logical(defaults[[key]])) {
            value <- "TRUE"
        } else {
            i <- i + 1L
            value <- args[i]
        }
        suppressWarnings(storage.mode(value) <- storage.mode(defaults[[key]]))
        if (is.na(value))
            stop("invalid value for --", key, call. = FALSE)
        opts[[key]] <- value
        i <- i + 1L
    }
    opts
}
