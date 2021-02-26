## ----setup, include=FALSE-----------------------------------------------------
library(GenomeAdmixR)
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width = 6)
print_mat <- function(mat) {
  n <- nrow(mat)
  c("\\begin{bmatrix}",
    paste0(sapply(seq_len(n - 1),
                  function(i) paste0(mat[i, ], collapse = " & ")),
           " \\\\"),
    paste0(mat[n, ], collapse = " & "),
    "\\end{bmatrix}")
}

## ----creating fake data 2-----------------------------------------------------
  fake_input_data <-
    create_artificial_genomeadmixr_data(number_of_individuals = 10,
                                        marker_locations = 1:10)

  fake_input_data$markers
  fake_input_data$genomes


## ----simulate admix data------------------------------------------------------
simulated_pop <- simulate_admixture(
                    module = sequence_module(molecular_data = fake_input_data),
                    pop_size = 100,
                    total_runtime = 100)

## ----creating fake data-------------------------------------------------------
  chosen_markers <- 1:100
   fake_input_data1 <-
    create_artificial_genomeadmixr_data(number_of_individuals = 100,
                                        marker_locations = chosen_markers,
                                        used_nucleotides = 1:2)
  fake_input_data2 <-
    create_artificial_genomeadmixr_data(number_of_individuals = 100,
                                        marker_locations = chosen_markers,
                                        used_nucleotides = 3:4)
  
  combined_data <- combine_input_data(input_data_list = list(fake_input_data1,
                                                             fake_input_data2),
                                      frequencies = c(0.5, 0.5),
                                      pop_size = 1000)

## ----simulate-----------------------------------------------------------------
simulated_pop <- simulate_admixture(
                      module = sequence_module(molecular_data = combined_data,
                                               markers = 1:100,
                                               morgan = 1),
                      pop_size = 1000,
                      total_runtime = 100)

## ----visualise----------------------------------------------------------------
plot_over_time(simulated_pop$frequencies, focal_location = 50)
plot_joyplot_frequencies(simulated_pop$frequencies, time_points = c(0, 50, 100))

## ----simulate selection-------------------------------------------------------
selection_matrix <- matrix(NA, nrow = 1, ncol = 5)
selection_matrix[1, ] <- c(50, 1, 1.1, 1.2, 1)

# we expect allele 'a' to do well at 0.5 Morgan:
selected_pop <- simulate_admixture(
                      module = sequence_module(molecular_data = fake_input_data1,
                                               morgan = 1,
                                               markers = chosen_markers),
                      pop_size = 1000,
                      total_runtime = 100,
                      select_matrix = selection_matrix)
plot_over_time(selected_pop$frequencies, focal_location = 50)

## ----migration----------------------------------------------------------------

migr_pop <-
  simulate_admixture(
        module = sequence_module(molecular_data = list(fake_input_data1,
                                                       fake_input_data2)),
        migration = migration_settings(migration_rate = 0.01,
                                       population_size = c(100, 100)),
        total_runtime = 100)

## ----from_simulation----------------------------------------------------------
simulated_pop <- simulate_admixture(
                      module = ancestry_module(), 
                      pop_size = 100,
                      total_runtime = 100)
prepared_pop  <-
  simulation_data_to_genomeadmixr_data(simulation_data = simulated_pop,
                                       markers = seq(0, 1, length.out = 100))

simulated_pop2 <- simulate_admixture(
                        module = sequence_module(molecular_data = prepared_pop),
                        pop_size = 100,
                        total_runtime = 100)

## ---- results = 'asis'--------------------------------------------------------
writeLines(print_mat(matrix(1 / 4, nrow = 4, ncol = 4)))

## ---- results = 'asis'--------------------------------------------------------
kappa <- 0.2
diag_entry <- 1 - (1 / 4 + 1 / 4 + 0.2 / 4)

writeLines(print_mat(matrix(c(diag_entry, kappa, 1 / 4, 1 / 4,
                              kappa, diag_entry, 1 / 4, 1 / 4,
                              1 / 4, 1 / 4, diag_entry, kappa,
                              1 / 4, 1 / 4, kappa, diag_entry),
                              nrow = 4, ncol = 4)))

## ----example mutation---------------------------------------------------------
  fake_input_data3 <-
  create_artificial_genomeadmixr_data(number_of_individuals = 100,
                                      marker_locations = 1:1000,
                                      used_nucleotides = 1)
  
  kappa <- 0.5
  rate_matrix <- matrix(c(0, kappa, 1, 1,
                          kappa, 0, 1, 1,
                          1, 1, 0, kappa,
                          1, 1, kappa, 0), nrow = 4, ncol = 4)
  
  mutation_rate <- 1e-5
  
  simulated_pop <- simulate_admixture(
                      module = sequence_module(molecular_data = fake_input_data3,
                                               morgan = 1,
                                               markers = chosen_markers,
                                               mutation_rate = mutation_rate,
                                               substitution_matrix = rate_matrix),
                      pop_size = 1000,
                      total_runtime = 100)
  
  freqs <- calculate_allele_frequencies(simulated_pop)
  
  require(magrittr, quietly = TRUE)
  freqs %>%
    dplyr::group_by(ancestor) %>%
    dplyr::summarise(mean(frequency))

