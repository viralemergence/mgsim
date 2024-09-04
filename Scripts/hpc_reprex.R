library(doFuture)
library(future.batchtools)
registerDoFuture()
plan(list(batchtools_torque, multicore))
y <- foreach(x = 1:4, y = 1:10) %dopar% {
  z <- x + y
  sqrt(z)
}
