library("rethinking")

output <- read.csv("sim_results_k_0_3.csv")

output$X <- output$X %% 300

i_cumsum <- c(0)
d_cumsum <- c(0)
r_cumsum <- c(0)

for(i in 1:(length(output$X)/300)){
  i_cumsum <- c(i_cumsum, cumsum(output$new_i_period[((i-1)*300 + 1):(((i-1)*300)+300)]))
  d_cumsum <- c(d_cumsum, cumsum(output$new_d_period[((i-1)*300 + 1):(((i-1)*300)+300)]))
  r_cumsum <- c(r_cumsum, cumsum(output$new_r_period[((i-1)*300 + 1):(((i-1)*300)+300)]))
}

output$i_cumsum <- i_cumsum[2:length(i_cumsum)] # drop the first element that was 0 and added for coding convenience
output$d_cumsum <- d_cumsum[2:length(d_cumsum)] # drop the first element that was 0 and added for coding convenience
output$r_cumsum <- r_cumsum[2:length(r_cumsum)] # drop the first element that was 0 and added for coding convenience

sim_results <- output[seq(300, length(output$X)+1, by=300),]
plot(sim_results$d_cumsum, main="Total Deaths By Simulation", xlab="Simulation", ylab="Deaths", col="#1b4c9b")
dens(sim_results$d_cumsum, main="Density Distribution")
precis(sim_results$d_cumsum)

plot(sim_results$i_cumsum, main="Total Infected By Simulation", xlab="Simulation", ylab="Infected", col="#1b4c9b")
dens(sim_results$i_cumsum, main="Density Distribution")
precis(sim_results$i_cumsum)

HPDI(sim_results$i_cumsum)

sim_results$cfr <- sim_results$d_cumsum / sim_results$i_cumsum
dens(sim_results$cfr, main="Infection Fatality Rate Distribution", xlab="IFR [%]", ylab="Density")
precis(cfr, prob=0.97, digits=0)
mean(cfr)
