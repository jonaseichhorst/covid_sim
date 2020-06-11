library("rethinking")

set.seed(101)

population <- 25000
contacts <- rnbinom(population, 2, 0.05)+1 #everybody must have at least one contact
#dens(contacts)

#population <- 50
#contacts <- rnbinom(population, 1, 0.1)+1 #everybody must have at least one contact
#dens(contacts)

#population <- 100
# every person has a different number of contacts modeled through this negative
# binomial distribution with a mean of about 40 x 2 = 80 (~1/2 of Dunbar's number)
#contacts <- rnbinom(population, 1, 0.1)

start <- proc.time()

population_idx <- seq(population)
p_pairs <- c()

for (i in 1:population){
  new_links <- sample(population_idx, contacts[i], prob=contacts, replace=TRUE)
  new_pairs <- cbind(rep(i, length(new_links)), new_links)
  p_pairs <- rbind(p_pairs, new_pairs)
}
end <- proc.time()

print(end - start) 

#write.csv(p_pairs,"10000_population_pairs.csv")

# starting value:  1,000:   2.24s ->   0.17s
# starting value: 10,000: 249.89s ->   7.87s
# starting value: 25,000: ??????? ->  58.15s
# starting value: 50,000: ??????? -> 268.24s