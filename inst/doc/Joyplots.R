## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width = 6)

## -----------------------------------------------------------------------------
library(GenomeAdmixR)
library(ggplot2)
packageVersion("GenomeAdmixR")

## ----setup select_matrix------------------------------------------------------
select_matrix <- matrix(ncol = 5, nrow = 2)
s <- 0.1
select_matrix[1, ] <- c(0.23, 1.0, 1 + 0.5 * s, 1 + s, 0)
select_matrix[2, ] <- c(0.27, 1.0, 1 + 0.5 * s, 1 + s, 1)

markers <- seq(from = 0.2, to = 0.3, length.out = 100)

## ----simulate-----------------------------------------------------------------
selected_pop <- simulate_admixture(pop_size = 1000,
                                   number_of_founders = 4,
                                   total_runtime = 1001,
                                   morgan = 1,
                                   select_matrix = select_matrix,
                                   markers = markers)

## ----joyplot all--------------------------------------------------------------
time_points <- seq(from = 0, to = 1000, by = 100)
plot_joyplot_frequencies(selected_pop$frequencies,
                         time_points,
                         picked_ancestor = "ALL")

## ----joyplot separate 0-------------------------------------------------------
p1 <- plot_joyplot_frequencies(selected_pop$frequencies,
                               time_points,
                               picked_ancestor = 0)
p1 + ggplot2::geom_vline(xintercept = 0.23, lty = 2)

## ----joyplot separate 1-------------------------------------------------------
p2 <- plot_joyplot_frequencies(selected_pop$frequencies,
                               time_points,
                               picked_ancestor = 1)
p2 + ggplot2::geom_vline(xintercept = 0.27, lty = 2)

