cl <- makeCluster(type = "MPI")
registerDoParallel(cl)

print(parallel::detectCores())

