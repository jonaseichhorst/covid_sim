library("rethinking")

input_df <- read.csv("./results/sim_results-0.2-0-9999-2.5-0.001-40-younger-10000-.csv")
#output_df <- t
#input_df <- output_df
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

sim_results <- input_df[seq(300, length(input_df$X)+1, by=300),]
sim_results$cfr <- sim_results$d_cumsum / sim_results$i_cumsum

mean(sim_results$cfr)
sim_results
plot(input_df$X, input_df$i_cumsum, main="Infections by Simulation", xlab="", xaxt = "n", ylab="Infections", col="#1b4c9b", pch=1)
dens(sim_results$cfr, main="Implied IFR Distribution", xlab="Implied IFR [%]", ylab="Density")

#x0.15_0.25_10_2.5_older <- sim_results
#x0.15_0.25_10_2.5_younger <- sim_results
#x0.15_0.25_10_2.5_US <- sim_results
#x0.15_0.25_10_2.5_very_old <- sim_results

#b_very_old <- sim_results
#b_older <- sim_results
#b_US <- sim_results
#b_younger <- sim_results
color_vals <- c("#2176FF", "#A30B37", "#FF9505", "#391463")
par(mfrow = c(2,1))

dens(b_younger$cfr, main="Baseline Model", xlab="", ylab="Density", xaxt='n', xlim=c(0, 0.15), col=color_vals[1], lty=2, lwd=2)
dens(b_US$cfr, add = TRUE, col=color_vals[2], lty=1, lwd=2)
dens(b_older$cfr, add = TRUE, col=color_vals[3], lty=3, lwd=2)
dens(b_very_old$cfr, add = TRUE, col=color_vals[4], lty=6, lwd=2)
legend(0.125, 185, c("Younger", "US Baseline", "Older", "Very Old"), bg="white", lty = c(2, 1, 3, 6), col=color_vals, text.col = color_vals)

dens(x0.15_0.25_10_2.5_younger$cfr, main="k: 0.15, social distancing: 0.25, group limit: 10",xlab="Implied IFR [%]", ylab="Density", xlim=c(0, 0.15), col=color_vals[1], lty=2, lwd=2)
dens(x0.15_0.25_10_2.5_US$cfr, add = TRUE, col=color_vals[2], lty=1, lwd=2)
dens(x0.15_0.25_10_2.5_older$cfr, add = TRUE, col=color_vals[3], lty=3, lwd=2)
dens(x0.15_0.25_10_2.5_very_old$cfr, add = TRUE, col=color_vals[4], lty=6, lwd=2)
legend(0.125, 48, c("Younger", "US Baseline", "Older", "Very Old"), bg="white", lty = c(2, 1, 3, 6), col=color_vals, text.col = color_vals)


#legend(1, 95, legend=c("Line 1", "Line 2"),
#       col=c("red", "blue"), lty=1:2, cex=0.8,
#       title="Line types", text.font=4, bg='lightblue')