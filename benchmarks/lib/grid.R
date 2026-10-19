TIERS <- list(   # nf = features; nx/ny = items compared (ny NA -> y = NULL)
    large   = list(nf = 5e4, nx = 100,  ny = 200, density = 0.9,  inputs = c("sparse", "dense")),  # half of #74's shapes, for CI time
    xsparse = list(nf = 5e4, nx = 2000, ny = NA,  density = 0.01, inputs = "sparse"),              # real sparse use case
    pair    = list(nf = 1e4, nx = 100,  ny = 100, density = 0.3,  inputs = c("sparse", "dense"))
)
VARIANTS <- list(base = list(), diag = list(diag = TRUE), rank = list(rank = 10L),
                 min_simil = list(min_simil = 0.5), mask = list(mask = TRUE))  # mask built in child from labels
LINEAR   <- list(simil = c("cosine", "correlation", "dice", "edice"), dist = "euclidean")
PAIR     <- list(simil = c("jaccard", "ejaccard", "fjaccard", "hamann", "faith", "simple matching"),
                 dist  = c("chisquared", "kullback", "jeffreys", "jensen", "manhattan", "maximum",
                           "canberra", "minkowski", "hamming"))
BINARY   <- c("jaccard", "hamann", "faith", "simple matching", "hamming")   # 0/1 data
SMOOTHED <- c("chisquared", "kullback", "jeffreys", "jensen")  # smooth = 0 returns NaN early on any zero
IMPLIED  <- c(crossprod = 2, tcrossprod = 1, colSds = 2, rowSds = 1, colZeros = 2, rowZeros = 1)  # margin from fn name

spec <- function(fn, method, tier, variant = "base") {
    expand.grid(fn = fn, method = method, tier = tier, variant = variant, stringsAsFactors = FALSE)
}

build_grid <- function() {
    # which functions and methods run on which tier and variant
    grid <- rbind(
        spec("simil", LINEAR$simil, c("large", "xsparse")),
        spec("dist", LINEAR$dist, c("large", "xsparse")),
        spec("simil", PAIR$simil, "pair"),
        spec("dist", PAIR$dist, "pair"),
        spec("simil", "cosine", "large", c("rank", "min_simil", "mask")),
        spec("simil", "cosine", "pair", "diag"),
        spec(names(IMPLIED), NA, c("large", "xsparse"))
    )
    # cross with input x margin, keeping inputs the tier allows and the margin a function implies
    grid <- merge(grid, expand.grid(input = c("sparse", "dense"), margin = 1:2,
                                    stringsAsFactors = FALSE), by = NULL)
    allowed <- mapply(function(tier, input) input %in% TIERS[[tier]]$inputs, grid$tier, grid$input)
    implied <- IMPLIED[grid$fn]
    grid <- grid[allowed & (is.na(implied) | grid$margin == implied), ]

    # fn_method[_variant]_input_mN, with a tier suffix only on names that clash with an earlier tier
    grid <- grid[order(match(grid$tier, names(TIERS))), ]
    parts <- cbind(grid$fn, grid$method, ifelse(grid$variant == "base", NA, grid$variant),
                   grid$input, paste0("m", grid$margin))
    grid$name <- gsub(" ", "_", apply(parts, 1, function(p) paste(p[!is.na(p)], collapse = "_")))
    clash <- duplicated(grid$name)
    grid$name[clash] <- paste0(grid$name[clash], "_", grid$tier[clash])
    stopifnot(!anyDuplicated(grid$name))

    grid$binary <- grid$method %in% BINARY
    grid$args <- lapply(seq_len(nrow(grid)), function(i)
        c(VARIANTS[[grid$variant[i]]], if (grid$method[i] %in% SMOOTHED) list(smooth = 1)))
    grid$data_key <- paste(grid$tier, grid$input, paste0("m", grid$margin),
                           ifelse(grid$binary, "binary", "real"), sep = "|")
    grid <- grid[order(grid$data_key, grid$name), ]
    rownames(grid) <- NULL
    grid
}
