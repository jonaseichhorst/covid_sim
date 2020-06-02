#0: set up simulation with p=0 and a randomly picked number of sick patients
### Aferwards, for each period:
#1: determine the number p_inf of potenially infectuous patients (have been infected within 14 days)
#2: determine the total sum s_inf of secondary infections across all infectuous primary infections
#3: Randomly set s_inf people from susceptible population to become infected in p+1 adjusted for the susceptible
#   portion of the total active (S + R) portion of the population
#4: For new patients, determine their mortality risk and set their period of death or recovery

# Target Output:
# a data frame with n=population size rows and 4 columns:
# 1) is person still susceptible
# 2) what period (if any) did that person get infected
# 3) what period (if any) did that person die
# 4) what period (if any) did that person recover


library("rethinking")
library("parallel")
library("collections")

set.seed(101)

calculate_ifr <- function(){
  # TODO: Further define the correct mortality by age bracket
  #likelihood <- c(0.01, 0.01, 0.0001, 0.04, 0.20, 0.86, 3.50, 10.3, 27.06, 40.98, 17.06)
  mortality_likelihood <- c(0.01, 0.01, 0.0001, 0.04, 0.20, 0.86, 3.50, 10.3, 27.06, 40.98, 60.06) # last value cleansed
  return(mortality_likelihood / sum(mortality_likelihood))
}

setup_population <- function( N_pop_size ){
  s <- rep(1,N_pop_size) # All are susceptible
  i <- rep(-1,N_pop_size) # nobody infected & contagious
  d <- rep(-1,N_pop_size) # nobody died
  r <- rep(-1,N_pop_size) # nobody recovered
  
  # age distribution:
  # work with US demographic data in 5 year increments, then simply randomize a bit to have nicer distribution
  a <- sample(0:17, size = N_pop_size, replace = TRUE, prob=c(19.81, 20.2, 20.88, 21.09, 21.87, 23.56, 22.13, 21.56, 19.72, 20.74, 20.89, 21.94, 20.33, 17.08, 13.4, 9.26, 6.13, 6.55))
  a <- a*5+sample(0:4, N_pop_size, replace = TRUE)
  return(data.frame(s, i, d, r, a))
}

sim_city_epidemic <- function( N_pop_size = 10000 , N_periods = 300 , 
                               N_initial_patients = 5 ,
                               agg_distancing = c(0.05, 20), mod_distancing = c(0.02, 50) ) {
  
  ifr <- calculate_ifr()
  pop <- setup_population( N_pop_size )
  
  for ( period in 0:N_periods ) {
    
    # determine how many patients will be infected in this period
    if(period < 8){ # seed the simulation first with a constant stream of inbound patients
      new_patients_count <- N_initial_patients
    }
    else{
      
      # patients between 4-7 days post infection are most dangerous
      infectious_patients <- sum((period - 3 > i) & (period - 8 < i ))
      
      # TODO: normalize for difference in infectiousness by day:
      # currently flat infectiousness during peak days assumed (day 4-7)
      # Mu & K default (2.5 & 0.1) set in accordance with Endo et al: https://wellcomeopenresearch.org/articles/5-67
      secondary_infections <- rnbinom(infectious_patients / 4, 0.1, mu=2.5)
      
      # active cases are all infected cases minus the recovered and/or deceased
      active_cases <- length(i[i >= 0]) - length(r[r <= period & r > 0]) - length(d[d <= period & d > 0])
      
      
      # apply distancing measures at thresholds to reduce super spreader events
      if (active_cases / length(s[s>0]) > mod_distancing[1]){
         secondary_infections[secondary_infections > mod_distancing[2]] <- 0
      }
      if (active_cases / length(s[s>0]) > agg_distancing[1]){
        secondary_infections[secondary_infections > agg_distancing[2]] <- 0
      }
      
      # this is the total number of newly infected people inthis round
      new_patients_count <- sum(secondary_infections)
      
      # instead of assuming that infectious patients will only interact with susceptible patients
      # assume they deal with any healthy patients which means either susceptible or recovered, hence, the
      # ratio of susceptible people of the exposed people over time reduces as more people who are recovered join
      new_patients_count <- new_patients_count * (length(s[s == 1]) / ( length(s[s == 1]) + length(r[r > 0 & r < period])))
      cat("\nRound: ", period, " \tInfectious Patients: ", infectious_patients, 
          " \tNew Patients: ", new_patients_count, "\tActive Cases: ", active_cases)
      
    }
    
    # if we hit 100% penetration, just infect remaining potential patients
    if (sum(s > 0) < new_patients_count)
      new_patients_count = sum(s > 0)
    
    # abort simulation if everybody has been infected
    if(sum(s > 0) == 0){
      print("Entire population infected")
      break
    }
    
    # randomly pick patients from the susceptible pool - by using prob = s we can ensure we only pick
    # susceptible patients as the others have s:0, i.e. will have a 0 probability
    new_patients_idx <- sample(1:N_pop_size, new_patients_count, prob=s)
    
    # set the susceptible patients now as infected and record their infection period
    s[new_patients_idx] <- 0
    i[new_patients_idx] <- period
    
    # determine patients age and randomly simulate whether the patient dies or recovers
    # based on their age bracket risk (mortality risk is given by decade)
    sick_ages <- a[new_patients_idx]
    ifr_scores <- ifr[round(a[new_patients_idx]/10, digits=0)+1]
    mortality <- rbinom(new_patients_count, 1, ifr_scores)
    
    # set the dates for the recovery or decease for the people whose fate we now know
    recovery <- rep(1, new_patients_count) - mortality
    d[new_patients_idx] <- mortality * (period + 21) # patients die after 21 days
    r[new_patients_idx] <- recovery * (period + 21) # patients recover within 14 days
  }
  
  # summarize the events by day and return 
  new_i_period <- new_d_period <- new_r_period <- total_s_period <- 0
  for (p in 1:N_periods){
    new_i_period <- c(new_i_period, length(i[i==p]))
    new_d_period <- c(new_d_period, length(d[d==p]))
    new_r_period <- c(new_r_period, length(r[r==p]))
    total_s_period <- c(total_s_period, N_pop_size - length(i[(i > -1) & (i < p)]))
  }
  
  sim_result_by_day <- data.frame(new_i_period, new_d_period, new_r_period, total_s_period)
  return(sim_result_by_day)
}

t <- sim_city_epidemic()
# set up simulation
N_pop_size <- 10000
N_periods <- 300
N_initial_patients <- 5

run_simulation <- function(iterations){
  results_df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("new_i_period", "new_d_period", "new_r_period", "total_s_period"))  
  for (i in 1:iterations){
    x <- sim_city_epidemic()
    results_df <- rbind(results_df, x[2:length(x$new_i_period),])
  }
  return(results_df)
}

# create cluster object
numCores <- detectCores()
target_trials <- 1
trials_per_core <- target_trials / numCores

cl <- makeCluster(numCores)
clusterExport(cl, "sim_city_epidemic")
clusterExport(cl, "run_simulation")
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

#write.csv(output_df,"sim_results_k_0_2_mu_2_5_aggressive_distancing_2.csv")
