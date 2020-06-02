library("rethinking")

set.seed(101)

#population <- 10000
#contacts <- rnbinom(population, 2, 0.05)+1 #everybody must have at least one contact
dens(contacts)

population <- 50
contacts <- rnbinom(population, 1, 0.1)+1 #everybody must have at least one contact
dens(contacts)

#population <- 100
# every person has a different number of contacts modeled through this negative
# binomial distribution with a mean of about 40 x 2 = 80 (~1/2 of Dunbar's number)
#contacts <- rnbinom(population, 1, 0.1)

population_idx <- seq(population)
pairs <- c()


for (i in 1:population){
  new_links <- sample(population_idx, contacts[i], prob=contacts, replace=TRUE)
  for (l in 1:length(new_links)){
    pairs <- rbind(pairs, c(i, new_links[l]))
  }
}

table(pairs[,2])
