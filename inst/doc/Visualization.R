## ----setup, include=FALSE-----------------------------------------------------
library(GenomeAdmixR)
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width = 6)

## ----simulate-----------------------------------------------------------------
select_matrix <- matrix(nrow = 1, ncol = 5)
select_matrix[1, ] <- c(0.5, 1, 1 + 0.05, 1 + 0.1, 0)
population <- simulate_admixture(pop_size = 1000,
                                 number_of_founders = 2,
                                 total_runtime = 200,
                                 select_matrix = select_matrix,
                                 markers = c(0.5, seq(0, 1, length.out = 100)))

## ----plot over time-----------------------------------------------------------
plot_over_time(population$frequencies, focal_location = 0.500)

## ----plot frequencies---------------------------------------------------------
plot_frequencies(population, locations = seq(0, 1, length.out = 1000))

## ----plot difference frequencies----------------------------------------------
plot_difference_frequencies(population)

## ----plot start end-----------------------------------------------------------
plot_start_end(population)

## ----joyplot------------------------------------------------------------------
plot_joyplot_frequencies(population$frequencies,
                         time_points = c(0, 10, 25, 50, 100, 199))

## ----plot chromosome----------------------------------------------------------
plot(population$population[[1]])

plot_chromosome(population$population[[1]]$chromosome1, xmin = 0.45, xmax = 0.55)

