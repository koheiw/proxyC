# fixed seed, so every commit and testbed benchmarks identical data
make_data <- function(tier, input, margin, binary, seed = 42L) {
    tr <- TIERS[[tier]]
    set.seed(seed)
    gen <- function(n) {
        if (is.na(n))
            return(NULL)
        dims <- if (margin == 1) c(n, tr$nf) else c(tr$nf, n)
        sp <- Matrix::rsparsematrix(dims[1], dims[2], tr$density,
                                    rand.x = function(k) stats::runif(k, 0.1, 1))
        if (binary)
            sp@x[] <- 1
        if (input == "dense") as.matrix(sp) else sp
    }
    x <- gen(tr$nx)
    y <- gen(tr$ny)
    # group labels for the mask variant; the child turns them into a mask before timing
    list(x = x, y = y, lx = sample(letters[1:5], tr$nx, TRUE),
         ly = if (!is.na(tr$ny)) sample(letters[1:5], tr$ny, TRUE))
}

# one uncompressed .rds per data_key, shared by all cases of that group
ensure_data <- function(row, dir = tempdir()) {
    path <- file.path(dir, paste0(gsub("[^[:alnum:]]+", "_", row$data_key), ".rds"))
    if (!file.exists(path))
        saveRDS(make_data(row$tier, row$input, row$margin, row$binary), path, compress = FALSE)
    path
}
