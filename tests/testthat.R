require(proxyC)
require(testthat)

options(proxyC.threads = 2)
test_check("proxyC")
options(proxyC.threads = NULL)
