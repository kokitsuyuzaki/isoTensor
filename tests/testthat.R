library("testthat")
library("isoTensor")

options(testthat.use_colours = FALSE)

test_file("testthat/test_isoToyModel.R")
test_file("testthat/test_isoTensor.R")
test_file("testthat/test_isoTensor_tensor.R")
test_file("testthat/test_isoTensor_lbfgs.R")
