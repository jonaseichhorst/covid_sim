library("rethinking")

#input_df <- read.csv("sim_results_k_0_2_mu_2_5.csv")

input_df <- output_df
input_df$X <- seq(1, length(input_df[,1]))

i_cumsum <- c(0)
d_cumsum <- c(0)
r_cumsum <- c(0)

for(i in 1:(length(input_df$X)/300)){
  i_cumsum <- c(i_cumsum, cumsum(input_df$new_i_period[((i-1)*300 + 1):(((i-1)*300)+300)]))
  d_cumsum <- c(d_cumsum, cumsum(input_df$new_d_period[((i-1)*300 + 1):(((i-1)*300)+300)]))
  r_cumsum <- c(r_cumsum, cumsum(input_df$new_r_period[((i-1)*300 + 1):(((i-1)*300)+300)]))
}

input_df$i_cumsum <- i_cumsum[2:length(i_cumsum)] # drop the first element that was 0 and added for coding convenience
input_df$d_cumsum <- d_cumsum[2:length(d_cumsum)] # drop the first element that was 0 and added for coding convenience
input_df$r_cumsum <- r_cumsum[2:length(r_cumsum)] # drop the first element that was 0 and added for coding convenience

par(mfrow=c(2,2))
sim_results <- input_df[seq(300, length(input_df$X)+1, by=300),]
plot(sim_results$d_cumsum, main="Total Deaths By Simulation", xlab="Simulation", ylab="Deaths", col="#1b4c9b")
dens(sim_results$d_cumsum, main="Density Distribution")
precis(sim_results$d_cumsum)

plot(sim_results$i_cumsum, main="Total Infected By Simulation", xlab="Simulation", ylab="Infected", col="#1b4c9b")
dens(sim_results$i_cumsum, main="Density Distribution")
precis(sim_results$i_cumsum)

HPDI(sim_results$i_cumsum)

sim_results$cfr <- sim_results$d_cumsum / sim_results$i_cumsum
dens(sim_results$cfr, main="Infection Fatality Rate Distribution", xlab="IFR [%]", ylab="Density")
precis(sim_results$cfr, prob=0.97, digits=0)
mean(sim_results$cfr)
sim_results
input_df[300,]
input_df[301,]
input_df[input_df$X == 0, ]
plot(input_df$X, input_df$i_cumsum)
