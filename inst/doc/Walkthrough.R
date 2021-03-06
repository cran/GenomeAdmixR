## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width = 6)

## -----------------------------------------------------------------------------
library(GenomeAdmixR)
library(ggplot2)
packageVersion("GenomeAdmixR")

## ----create wildpop-----------------------------------------------------------
wildpop <-  simulate_admixture(
                    module = ancestry_module(number_of_founders = 4,
                                             morgan = 1),
                    pop_size = 100,
                    total_runtime = 1000)


## -----------------------------------------------------------------------------
wildpop

## -----------------------------------------------------------------------------
wildpop[[1]]

## ----create isofemale---------------------------------------------------------
isofemale <- create_iso_female(
                  module = ancestry_module(input_population = wildpop,
                                           morgan = 1),
                  n = 1,
                  inbreeding_pop_size = 1000,
                   run_time = 200)

## -----------------------------------------------------------------------------
plot(isofemale[[1]])

## ----toyexample---------------------------------------------------------------
wildpop <-  simulate_admixture(
                  module = ancestry_module(number_of_founders = 4,
                                           morgan = 1),
                  pop_size = 100,
                  total_runtime = 100)

isofemale <- create_iso_female(
                  module = ancestry_module(input_population = wildpop,
                                           morgan = 1),
                  n = 1,
                  inbreeding_pop_size = 100,
                  run_time = 10000)
plot(wildpop$population[[1]])
plot(isofemale[[1]])

## ----plot selection of chromosome---------------------------------------------
plot_chromosome(isofemale[[1]]$chromosome1, xmin = 0.0, xmax = 0.5)

## ----create two populations---------------------------------------------------
both_populations <-
  simulate_admixture(
    module = ancestry_module(morgan = 1),
    migration = migration_settings(population_size = c(100, 100),
                        migration_rate = 0.0,
                        initial_frequencies =
                                 list(c(rep(1, 20), rep(0, 20)),
                                      c(rep(0, 20), rep(1, 20)))
                      ),
    total_runtime = 1000)

population_1 <- both_populations$population_1

population_2 <- both_populations$population_2

## ----draw two isofemales------------------------------------------------------
isofemales <- create_iso_female(
                  module = ancestry_module(input_population = population_1,
                                           morgan = 1),
                  n = 2,
                  inbreeding_pop_size = 100,
                  run_time = 10000)

plot_chromosome(isofemales[[1]]$chromosome1, 0, 1)
plot_chromosome(isofemales[[2]]$chromosome1, 0, 1)

## ----seed mixed population----------------------------------------------------
mixed_population <-
  simulate_admixture(
        module = ancestry_module(input_population = list(isofemales[[1]],
                                                         isofemales[[2]]),
                                 morgan = 1),
        pop_size = 100, 
        total_runtime = 100)

## ----plot mixed_population----------------------------------------------------
plot(mixed_population$population[[1]])

## ----calc FST-----------------------------------------------------------------
fst <- calculate_fst(population_1,
                 population_2,
                 sampled_individuals = 10,
                 number_of_markers = 100,
                 random_markers = TRUE)
fst

## -----------------------------------------------------------------------------
  ld_results <- calculate_ld(wildpop,
                             markers = 10)

  plot(ld_results$ld_matrix~ld_results$dist_matrix,
       xlab = "Genetic Distance (Morgan)",
       ylab = "LD",
       pch = 16)
  plot(ld_results$rsq_matrix~ld_results$dist_matrix,
       xlab = "Genetic Distance (Morgan)",
       ylab = "r_sq",
       pch = 16)
  plot(ld_results$ld_matrix~ld_results$rsq_matrix,
       xlab = "r_sq",
       ylab = "LD",
       pch = 16)

## ----no LD--------------------------------------------------------------------
no_ld_pop <- simulate_admixture(
                  module = ancestry_module(number_of_founders = 4,
                                           morgan = 1),
                  pop_size = 100,
                  total_runtime = 1000)

ld_results <- calculate_ld(no_ld_pop,
                           sampled_individuals = 10,
                           markers = 10)

plot(ld_results$ld_matrix~ld_results$dist_matrix,
     pch = 16,
     xlab = "Distance",
     ylab = "LD",
     xlim = c(0, 1),
     ylim = c(0, 1))

## ----strong LD----------------------------------------------------------------
strong_ld_pop <- simulate_admixture(
                  module = ancestry_module(number_of_founders = 4,
                                           morgan = 1),
                  pop_size = 100,
                  total_runtime = 10)

ld_results <- calculate_ld(strong_ld_pop,
                           sampled_individuals = 10,
                           markers = 10)

plot(ld_results$ld_matrix~ld_results$dist_matrix,
     pch = 16,
     xlab = "Distance",
     ylab = "LD",
     xlim = c(0, 1),
     ylim = c(0, 1))

## ----heterozygote selection---------------------------------------------------
s <- 0.1
selection_matrix <- matrix(nrow = 1, ncol = 5)
selection_matrix[1, ] <- c(0.5,
                           1.0, 1.0 + s, 1.0,
                           0)

markers <- 0.5

selected_pop <- simulate_admixture(
                    module = ancestry_module(number_of_founders = 2,
                                             morgan = 1,
                                             markers = markers),
                    pop_size = 1000,  
                    total_runtime = 100,
                    select_matrix = selection_matrix)

plot_over_time(selected_pop$frequencies, markers)

## ----homozygote selection-----------------------------------------------------
s <- 0.1
selection_matrix[1, ] <- c(0.5,
                           1.0, 1.0, 1.0 + s,
                           0)

selected_pop <- simulate_admixture(
                    module = ancestry_module(number_of_founders = 10,
                                             morgan = 1,
                                             markers = markers),
                    pop_size = 1000,  
                    total_runtime = 300,
                    select_matrix = selection_matrix)

plot_over_time(selected_pop$frequencies, markers)

