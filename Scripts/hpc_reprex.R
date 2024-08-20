library(doParallel) 
library(foreach)

# Detect the number of available cores 
numCores <- detectCores()

# Create a cluster with the number of cores detected 
cl <- makeCluster(numCores)

# Register the parallel backend 
registerDoParallel(cl)

# Define a function to run in parallel 
myFunction <- function(x) { return(x^2) }

# Use foreach to run the function in parallel 
results <- foreach(i = 1:10, .combine = 'c') %dopar% { myFunction(i) }

 # Print the results print(results)

 # Stop the cluster 
stopCluster(cl)

