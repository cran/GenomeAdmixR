## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width = 6)

## ----load libraries-----------------------------------------------------------
library(GenomeAdmixR)
library(ggplot2)
packageVersion("GenomeAdmixR")

## ----create populations-------------------------------------------------------
pops <-
  simulate_admixture(
        module = ancestry_module(),
        migration = migration_settings(migration_rate = 0,
                                         population_size = c(100, 100),
                                         initial_frequencies =
                                 list(c(1, 1, 1, 1, 0, 0, 0, 0),
                                      c(0, 0, 0, 0, 1, 1, 1, 1))),
     total_runtime = 1000)

pop_1 <- pops$population_1
pop_2 <- pops$population_2

## ----calculate LD-------------------------------------------------------------
mean(calculate_ld(pop = pop_1,
                  sampled_individuals = 10,
                  markers = 30)$ld_matrix,
     na.rm = TRUE)
mean(calculate_ld(pop = pop_2,
                  sampled_individuals = 10,
                  markers = 30)$ld_matrix,
     na.rm = TRUE)

## ----create isofemale lines---------------------------------------------------
iso_females_pop_1 <- create_iso_female(
                          module = ancestry_module(input_population = pop_1),
                                       n = 2,
                                       inbreeding_pop_size = 100)
iso_females_pop_2 <- create_iso_female(
                          module = ancestry_module(input_population = pop_2),
                                       n = 2,
                                       inbreeding_pop_size = 100)

## ----create from individuals--------------------------------------------------
pop_1_1 <- simulate_admixture(
               module = ancestry_module(input_population = iso_females_pop_1),
                              pop_size = 1000,
                              total_runtime = 1000)

pop_1_2 <- simulate_admixture(
                    module = ancestry_module(input_population =
                                               list(iso_females_pop_1[[1]],
                                                    iso_females_pop_2[[1]])),
                               pop_size = 1000,
                               total_runtime = 1000)

pop_2_2 <- simulate_admixture(
                module = ancestry_module(input_population = iso_females_pop_2),
                              pop_size = 1000,
                              total_runtime = 1000)

## ----FST calculation----------------------------------------------------------
f1 <- calculate_fst(pop_1_1, pop_1_2,
                    sampled_individuals = 10, number_of_markers = 100)
# this one should be highest
f2 <- calculate_fst(pop_1_1, pop_2_2,
                    sampled_individuals = 10, number_of_markers = 100)
f3 <- calculate_fst(pop_1_2, pop_2_2,
                    sampled_individuals = 10, number_of_markers = 100)
f1
f2
f3

