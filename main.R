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

# return the ifr by age bracket
calculate_ifr <- function(){
  # TODO: Further define the correct mortality by age bracket
  #likelihood <- c(0.01, 0.01, 0.0001, 0.04, 0.20, 0.86, 3.50, 10.3, 27.06, 40.98, 17.06)
  mortality_likelihood <- c(0.01, 0.01, 0.0001, 0.04, 0.20, 0.86, 3.50, 10.3, 27.06, 40.98, 60.06) # last value cleansed
  return(mortality_likelihood / sum(mortality_likelihood))
}


create_population_network <- function(population = 10000){
  start <- proc.time()
  contacts <- rnbinom(population, 2, 0.05)+1 #everybody must have at least one contact
  
  population_idx <- seq(population)
  p_pairs <- c()
  
  for (i in 1:population){
    new_links <- sample(population_idx, contacts[i], prob=contacts, replace=TRUE)
    new_pairs <- cbind(rep(i, length(new_links)), new_links)
    p_pairs <- rbind(p_pairs, new_pairs)
  }
  end <- proc.time()
  print(end - start) 
  return(p_pairs)
}

#x <- create_population_network(5000)


# initialitze data frame
# set up age distribution of population 
setup_population <- function( N_pop_size = 10000 ){
  s <- rep(1,N_pop_size) # All are susceptible
  i <- rep(-1,N_pop_size) # nobody infected & contagious
  d <- rep(-1,N_pop_size) # nobody died
  r <- rep(-1,N_pop_size) # nobody recovered
  a <- setup_population_ages( N_pop_size )
  return(data.frame(s, i, d, r, a))
}

setup_population_ages <- function( N_pop_size = 10000 ){
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
  
  a <- sample(0:17, size = N_pop_size, replace = TRUE, prob=very_old_age_range)
  a <- a*5+sample(0:4, N_pop_size, replace = TRUE)
  return(a)  
}

get_infectious_p <- function( period, pop ){
  # check for any patient that was infected 3-8 days ago
  i <- which(pop$i %in% pop$i[(pop$i >= 0 & pop$i > (period - 8) & pop$i < (period - 3))])
  return(i)
}

calc_secondary_infect <- function( period , pop , 
                                     N_initial_p_perc_pop = 0.001, seeding = FALSE, 
                                     agg_distancing = c(0.05, 10), mod_distancing = c(0.02, 50)){
  
  # patients between 4-7 days post infection are most dangerous
  n_infectious_patients <- length(get_infectious_p(period , pop))
  
  # active cases are all infected cases minus the recovered and/or deceased
  active_cases <- length(pop$i[pop$i >= 0]) - length(pop$r[pop$r <= period & pop$r > 0]) - length(pop$d[pop$d <= period & pop$d > 0])
  
  # determine how many patients will be infected in this period
  if(seeding){ # seed the simulation first with a constant stream of inbound patients  
    
    new_patients_count <- as.numeric(round(N_initial_p_perc_pop * nrow(pop), digits=0))
    secondary_infections <- rep(1, new_patients_count)
  }
  else{
    # TODO: normalize for difference in infectiousness by day:
    # currently flat infectiousness during peak days assumed (day 4-7)
    # Mu & K default (2.5 & 0.1) set in accordance with Endo et al: https://wellcomeopenresearch.org/articles/5-67
    secondary_infections <- rnbinom(n_infectious_patients / 4, 0.2, mu=2.5)
    
    
    # apply distancing measures at thresholds to reduce super spreader events
    if (active_cases / length(pop$s[pop$s>0]) > mod_distancing[1]){
      secondary_infections[secondary_infections > mod_distancing[2]] <- 0
    }
    if (active_cases / length(pop$s[pop$s>0]) > agg_distancing[1]){
      secondary_infections[secondary_infections > agg_distancing[2]] <- 0
    }
  }
  return(secondary_infections)
}

identify_new_patients <- function(sec_infections, pop, period, p_pairs, directed = TRUE){
  # randomly pick patients from the susceptible pool - by using prob = s we can ensure we only pick
  # susceptible patients as the others have s:0, i.e. will have a 0 probability
  
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
        
        new_patients_idx <- c(new_patients_idx, sample(potential_p, sec_infections[i]))
      }
    }
  }
  else{ # randomly pick new patients to infect
    new_patients_idx <- sample(1:length(pop[,1]), sum(sec_infections), prob=pop$s)
  }
    
  
  return(new_patients_idx)
}


sim_city_epidemic <- function( p_size = 10000, p_pairs, N_periods = 300, N_seed = 5) {
  
  ifr <- calculate_ifr( )
  pop <- setup_population( p_size )
  
  for ( period in 1:N_periods ) {
  
    # abort simulation if everybody has been infected
    if(sum(pop$s > 0) == 0){
      print("Entire population infected")
      break
    }
    
    seed <- ifelse(period <= N_seed, TRUE, FALSE)

    secondary_infections <- calc_secondary_infect(period = period, pop = pop, seeding = seed)
    new_patients_count <- sum(secondary_infections)
    
    if (new_patients_count > 0){
      new_patients_idx <- identify_new_patients(sec_infections = secondary_infections, pop = pop, p_pairs = p_pairs, period = period, directed=!seed)
      
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
      
      pop$d[new_patients_idx] <- mortality * (period + 21) # patients die after 21 days
      pop$r[new_patients_idx] <- recovery * (period + 21) # patients recover within 14 days
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
  
  sim_result_by_day <- data.frame(seq(1, length(new_i_period) - 1), new_i_period[-1], new_d_period[-1], new_r_period[-1], total_s_period[-1])
  names(sim_result_by_day)[1] <- "X"
  }

  return(sim_result_by_day)
  
}

#p_pairs <- read.csv("10000_population_p_pairs.csv")
#p_pairs$X <- NULL

#p_pairs <- create_population_network(1000)
#t <- sim_city_epidemic()


# set up simulation
N_pop_size <- 10000

run_simulation <- function(iterations = 1){
  p_pairs <- create_population_network(N_pop_size)
  results_df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("new_i_period", "new_d_period", "new_r_period", "total_s_period"))  
  for (i in 1:iterations){
    x <- sim_city_epidemic(p_size = N_pop_size, p_pairs = p_pairs)
    results_df <- rbind(results_df, x)
  }
  return(results_df)
}

# create cluster object
numCores <- detectCores()
target_trials <- 0
trials_per_core <- target_trials / numCores

cl <- makeCluster(numCores)
clusterExport(cl, c("sim_city_epidemic", "run_simulation", 
                    "calculate_ifr", "setup_population",
                    "calc_secondary_infect", "get_infectious_p",
                    "identify_new_patients", "create_population_network", "setup_population_ages"))

clusterExport(cl, c("trials_per_core", "N_pop_size"))

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

#write.csv(output_df,"sim_results_k_0_2_mu_2_5_with_dist_02_50_05_10_very_old.csv")
