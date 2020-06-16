library("rethinking")
library("parallel")
library("collections")

set.seed(101)

###############################
# return the ifr by age bracket
calculate_ifr <- function(){
  mortality_likelihood <- c(0.01, 0.0001, 0.04, 0.20, 0.86, 3.50, 10.3, 27.06, 40.98, 17.06)
  return(mortality_likelihood / sum(mortality_likelihood))
}

##############################
# population is connected by:
# 1) randomly determining the number of links each person will have
# 2) assigning those links per individual with clustering of similarly aged participants
setup_population_network <- function( ages ){
  
  if(rebuild_pop_network == FALSE){
    p_pairs <- read.csv("./params/pop_network_10000.csv")
    return(p_pairs)
  }
  
  pop_count <- length(ages)
  contacts <- rnbinom(pop_count, 2, 0.05)+1 #everybody must have at least one contact
  
  population_idx <- seq(pop_count)
  p_pairs <- c()
  
  for (i in 1:pop_count){
    # we cluster by calculating squared (age) distance of people to provide heavier weightage for similarly aged
    prob <- abs(ages-ages[i])^2
    prob <- max(prob) - prob
    
    new_links <- sample(population_idx, contacts[i], prob=prob, replace=TRUE)
    new_pairs <- cbind(rep(i, length(new_links)), new_links)
    p_pairs <- rbind(p_pairs, new_pairs)
  }
  return(p_pairs)
}

###################################
setup_population <- function( N_pop_size = 10000 , age_range = "US"){
  s <- rep(1,N_pop_size) # All are susceptible
  i <- rep(-1,N_pop_size) # nobody infected & contagious
  d <- rep(-1,N_pop_size) # nobody died
  r <- rep(-1,N_pop_size) # nobody recovered
  a <- setup_population_ages( N_pop_size, age_range )
  return(data.frame(s, i, d, r, a))
}

###################################
setup_population_ages <- function( N_pop_size = 10000, age_range = "US" ){
  # age distribution:
  # work with US demographic data in 5 year increments, then simply randomize a bit to have nicer distribution
  US_age_range <- c(19.81, 20.2, 20.88, 21.09, 21.87, 23.56, 22.13, 21.56, 19.72, 20.74, 20.89, 21.94, 20.33, 17.08, 13.4, 9.26, 6.13, 6.55)
  
  # double the proportion of younger half of population
  younger_age_range <- US_age_range
  younger_age_range[1:(length(US_age_range)/2)] <- US_age_range[1:(length(US_age_range)/2)]*2
  
  # double the proportion of older half of population
  older_age_range <- US_age_range
  older_age_range[(length(US_age_range)/2):length(US_age_range)] <- US_age_range[(length(US_age_range)/2):length(US_age_range)]*2
  older_age_range[(length(US_age_range)/4*3):length(US_age_range)] <- US_age_range[(length(US_age_range)/4*3):length(US_age_range)]*2
  
  # double the proportion of older half of population with 3x the top 10 percentile
  very_old_age_range <- US_age_range
  very_old_age_range[(length(US_age_range)/2):length(US_age_range)] <- US_age_range[(length(US_age_range)/2):length(US_age_range)]*2
  very_old_age_range[(length(US_age_range)/5*4):length(US_age_range)] <- US_age_range[(length(US_age_range)/5*4):length(US_age_range)]*3
  
  if(age_range == "US")
    prob <- US_age_range  
  if(age_range == "younger")
    prob <- younger_age_range
  if(age_range == "older")
    prob <- older_age_range
  if(age_range == "very_old")
    prob <- very_old_age_range
  
  a <- sample(0:17, size = N_pop_size, replace = TRUE, prob=prob)  
  a <- a*5+sample(0:4, N_pop_size, replace = TRUE)
  return(a)  
}


#######################################
get_infectious_p <- function( period, pop ){
  # check for any patient that was infected 3-8 days ago
  return(which(pop$i %in% pop$i[(pop$i >= 0 & pop$i > (period - 8) & pop$i < (period - 3))]))
}

#######################################
get_susceptible_p <- function( period, pop ){
  # any susceptible patient
  return(which(pop$s %in% pop$s[pop$s == 1]))
}

#######################################
get_active_p <- function( period, pop ){
  # any patient who's been infected but hasn't recovered or died yet
  return(which(pop$i %in% pop$i[(pop$i >= 0 & (period < pop$r) | (period < pop$d))]))
}

#######################################
get_last_period_transmissions <- function( period, pop){
  return(which(pop$i %in% pop$i[pop$i == (period - 1) & period > 0]))
}


#######################################
calc_secondary_infect <- function( mu=2.5, k=0.2, period , social_distancing_status, pop , 
                                     N_initial_p_perc_pop = 0.001, seeding = FALSE){
  
  # patients between 4-7 days post infection are most dangerous
  n_infectious_patients <- length(get_infectious_p(period , pop))
  
  # active cases are all infected cases minus the recovered and/or deceased
  active_cases <- length(pop$i[pop$i >= 0]) - length(pop$r[pop$r <= period & pop$r > 0]) - length(pop$d[pop$d <= period & pop$d > 0])
  
  # determine how many patients will be infected in this period
  if(seeding){ # seed the simulation first with a constant stream of inbound patients  
    
    new_patients_count <- as.numeric(round(N_initial_p_perc_pop * nrow(pop), digits=0))
    
    # ensure we always have at least one new infection
    if (new_patients_count == 0)
      new_patients_count <- 1
    
    secondary_infections <- rep(1, new_patients_count)
  }
  else{
    
    # currently flat infectiousness during peak days assumed (day 4-7)
    # Endo et al estimate Mu & K as 2.5 and 0.1 but the default value for k is set
    # slightly more conservative at 0.2
    # https://wellcomeopenresearch.org/articles/5-67
    secondary_infections <- rnbinom(n_infectious_patients / 4, k, mu=mu)
    
    # apply social distancing factor to reduce number of infections
    
    secondary_infections <- secondary_infections * rbinom(length(secondary_infections), 1, (1 - social_distancing_status[period, 1]))

    if (social_distancing_status[period, 2] > 0)
    {
        secondary_infections[secondary_infections > social_distancing_status[period, 2]] <- 0
    }
  }
  return(secondary_infections)
}

################################
manage_social_distancing <- function(new_patients_count, pop, period, social_distancing_status, 
                                     sd_ap_threshold = 0.01, sd_duration = 21, sd_reduction = 0.3, sd_max_group_size = 15){
  ap_count <- length(get_active_p(period, pop))
  if (social_distancing_status[period] == 0 & ap_count > sd_ap_threshold * nrow(pop)){
    end_period <- min((period + sd_duration), N_periods)
    social_distancing_status[period:end_period,1] <- sd_reduction
    social_distancing_status[period:end_period,2] <- sd_max_group_size
  }
  return(social_distancing_status)
}

#################################
identify_new_patients <- function(sec_infections, pop, period, p_pairs, directed = TRUE){
  
  if (directed){
    if(sum(sec_infections) > sum(pop$s)){
      print("exceeded")
      return(c())
    }
    
    # 1. if a patient causes secondary infections
    # 2.  find all connected people
    # 3.  infect the maximum possible number of secondary infections up to the expected maximum secondary infections
    infectious_p <- get_infectious_p(period, pop)
    new_patients_idx <- c()

    for (i in 1:length(sec_infections)){
      if(sec_infections[i] > 0)
      {
        connected_p <- p_pairs[which(p_pairs[,1]  == infectious_p[i]),2]
        potential_p <- connected_p[pop$s[connected_p] == 1] # anybody who's still susceptible is a potential patient

        # if we're running out of targets, limit to the amount of targets left
        if(sec_infections[i] > length(potential_p)){
          sec_infections[i] = length(potential_p)
        }
        
        # randomly pick the number of patients to infect from the connected patients
        new_patients_idx <- c(new_patients_idx, sample(potential_p, sec_infections[i]))
      }
    }
  }
  # randomly pick new patients to infect
  else{
    new_patients_idx <- sample(1:length(pop[,1]), sum(sec_infections), prob=pop$s)
  }
    
  
  return(new_patients_idx)
}

################################
sim_city_epidemic <- function( p_size = 10000, age_range = "US", 
                               N_periods = 300, N_seed = 5, N_initial_p_perc_pop = 0.001,
                               mu = 2.5, k = 0.2,
                               sd_ap_threshold = 0.01, sd_duration = 21, sd_reduction = 0.3, sd_max_group_size = 15) {

  ifr <- calculate_ifr( )
  pop <- setup_population( p_size , age_range)
  p_pairs <- setup_population_network(pop$a)
  social_distancing_status <- matrix(rep(0,600), ncol = 2)
  
  for ( period in 1:N_periods ) {
  
    # abort simulation if everybody has been infected
    if(sum(pop$s > 0) == 0){
      print("Entire population infected")
      break
    }
    
    seed <- ifelse(period <= N_seed, TRUE, FALSE)

    secondary_infections_idx <- calc_secondary_infect(mu = mu, k = k, period = period, N_initial_p_perc_pop = N_initial_p_perc_pop,
                                                      social_distancing_status = social_distancing_status, 
                                                      pop = pop, seeding = seed)
    
    new_patients_count <- sum(secondary_infections_idx)
    social_distancing_status <- manage_social_distancing(new_patients_count, pop, period, social_distancing_status, 
                                                         sd_ap_threshold, sd_duration, sd_reduction, sd_max_group_size)
    
    if (new_patients_count > 0){
      new_patients_idx <- identify_new_patients(sec_infections = secondary_infections_idx, pop = pop, p_pairs = p_pairs, period = period, directed=!seed)
      
      # set the susceptible patients now as infected and record their infection period
      pop$s[new_patients_idx] <- 0
      pop$i[new_patients_idx] <- period
      
      # determine patients age and randomly simulate whether the patient dies or recovers
      # based on their age bracket risk (mortality risk is given by decade)
      sick_ages <- pop$a[new_patients_idx]
      ifr_scores <- ifr[round(pop$a[new_patients_idx]/10, digits=0)+1]
      mortality <- rbinom(length(new_patients_idx), 1, ifr_scores)
      
      # set the dates for the recovery or decease for the people whose fate we now know
      recovery <- rep(1, length(new_patients_idx)) - mortality
      
      pop$d[new_patients_idx] <- mortality * (period + 18) # patients die after 18 days
      pop$r[new_patients_idx] <- recovery * (period + 22) # patients recover within 22 days
    }
  }
  
  # Dataframe summarizing
  {
  # summarize the events by day and return
  new_i_period <- new_d_period <- new_r_period <- total_s_period <- 0
  for (p in 1:N_periods){
    new_i_period <- c(new_i_period, length(pop$i[pop$i==p]))
    new_d_period <- c(new_d_period, length(pop$d[pop$d==p]))
    new_r_period <- c(new_r_period, length(pop$r[pop$r==p]))
    total_s_period <- c(total_s_period, length(pop$s) - length(pop$i[(pop$i > -1) & (pop$i < p)]))
  }
  
  sim_result_by_day <- data.frame(seq(1, length(new_i_period) - 1), new_i_period[-1], new_d_period[-1], new_r_period[-1], total_s_period[-1], social_distancing_status[,1], social_distancing_status[,2])
  colnames(sim_result_by_day) <- c("X", "new_i_period", "new_d_period", "new_r_period", "total_s_period", "sd_reduction", "sd_max_group_size")
  }

  return(sim_result_by_day)
}

run_simulation <- function(iterations = 1, row){
  results_df <- data.frame(matrix(ncol = 6, nrow = 0))
  for (i in 1:iterations){
    x <- sim_city_epidemic(p_size = sim_params[row,]$pop, age_range = sim_params[row,]$age,
                           k = sim_params[row,]$k_vals, 
                           mu = sim_params[row,]$mu_vals,
                           N_initial_p_perc_pop = sim_params[row,]$seed_inf_p,
                           sd_reduction = sim_params[row,]$sd,
                           sd_max_group_size = sim_params[row,]$gs
                           )
    results_df <- rbind(results_df, x)
  }
  return(results_df)
}

# create cluster object
numCores <- detectCores()

cl <- makeCluster(numCores)
clusterExport(cl, c("sim_city_epidemic", "run_simulation", 
                    "calculate_ifr", "setup_population",
                    "calc_secondary_infect", "get_infectious_p", "get_active_p", "manage_social_distancing", 
                    "identify_new_patients", "setup_population_network", "setup_population_ages"))

sim_params <- read.csv("./params/sim_params.csv")

for (r in 1:nrow(sim_params)){
  start <- proc.time()
  
  rebuild_pop_network <- sim_params[r,]$rebuild_pop_network
  trials_per_core <- sim_params[r,]$sims / numCores
  clusterExport(cl, c("trials_per_core", "N_pop_size", "sim_params", "rebuild_pop_network", "r", "N_periods"), envir=environment())
  cluster_results <- clusterEvalQ(cl, run_simulation(trials_per_core, r))
  
  # merge the different cluster results back into one frame
  output_df <- cluster_results[[1]]
  for (i in 2:numCores){
    output_df <- rbind(output_df, cluster_results[[i]])
  }
  
  rownames(output_df) <- seq(length=nrow(output_df))

  file_name <- paste("./results/sim_results", 
                     sim_params[r,1], sim_params[r,2], sim_params[r,3], sim_params[r,4], 
                     sim_params[r,5], sim_params[r,6], sim_params[r,7], sim_params[r,8],  
                     ".csv", sep = "-")
  write.csv(output_df, file_name)
  end <- proc.time()
  print(end - start)
  
}

# close
stopCluster(cl)