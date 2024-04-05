require(proxyC)
require(testthat)
Sys.setenv("RCPP_PARALLEL_NUM_THREADS" = 2)
test_check("proxyC")
