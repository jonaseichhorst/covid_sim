# set up simulation
N_pop_size <- 10000
N_periods <- 300
N_initial_patients <- 5

run_simulation <- function(iterations){
  results_df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("new_i_period", "new_d_period", "new_r_period", "total_s_period"))  
  for (i in 1:iterations){
    x <- sim_city_epidemic()
    results_df <- rbind(results_df, x)
  }
  return(results_df)
}

# create cluster object
numCores <- detectCores()
target_trials <- 1
trials_per_core <- target_trials / numCores

cl <- makeCluster(numCores)
clusterExport(cl, c("sim_city_epidemic", "run_simulation", "calculate_ifr", "setup_population"))
clusterExport(cl, "trials_per_core")

start <- proc.time()
cluster_results <- clusterEvalQ(cl, run_simulation(trials_per_core))

# merge the different cluster results back into one frame
output_df <- cluster_results[[1]]
for (i in 2:numCores){
  output_df <- rbind(output_df, cluster_results[[i]])
}

rownames(output_df) <- seq(length=nrow(output_df))
end <- proc.time()
print(end - start) 

# close
stopCluster(cl)

write.csv(output_df,"sim_results_k_0_2_mu_2_5_aggressive_distancing_2.csv")