% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_sequence_migration.R
\name{simulate_sequence_migration}
\alias{simulate_sequence_migration}
\title{Individual based simulation of the breakdown of contiguous ancestry blocks in
two populations linked by migration}
\usage{
simulate_sequence_migration(
  input_data_population_1 = NA,
  input_data_population_2 = NA,
  pop_size = c(100, 100),
  total_runtime = 100,
  morgan = 1,
  recombination_rate = NA,
  num_threads = 1,
  select_matrix = NA,
  markers = NA,
  verbose = FALSE,
  multiplicative_selection = TRUE,
  migration_rate = 0,
  stop_at_critical_fst = FALSE,
  critical_fst = NA,
  generations_between_update = 100,
  sampled_individuals = 10,
  number_of_markers = 100,
  random_markers = TRUE,
  mutation_rate = 0,
  substitution_matrix = matrix(1/4, 4, 4)
)
}
\arguments{
\item{input_data_population_1}{Genomic data used as input, should be created
by the function \code{create_input_data} or by the function
\code{combine_input_data}}

\item{input_data_population_2}{Genomic data used as input, should be created
by thefunction \code{create_input_data} or by the function
\code{combine_input_data}}

\item{pop_size}{Vector containing the number of individuals in both
populations.}

\item{total_runtime}{Number of generations}

\item{morgan}{Length of the chromosome in Morgan (e.g. the number of
crossovers during meiosis)}

\item{recombination_rate}{rate in cM / Mbp, used to map recombination to the
markers. If the recombination_rate is not set, the value for morgan is used,
assuming that the markers included span an entire chromosome.}

\item{num_threads}{number of threads. Default is 1. Set to -1 to use all
available threads}

\item{select_matrix}{Selection matrix indicating the markers which are under
selection. If not provided by the user, the simulation proceeds neutrally. If
provided, each row in the matrix should contain five entries:
\itemize{
\item{location of the marker under selection (in Morgan) }
\item{fitness of wildtype (aa)}
\item{fitness of heterozygote (aA)}
\item{fitness of homozygote mutant (AA)}
\item{Ancestral type that represents the mutant allele A}
}}

\item{markers}{A vector of locations of markers (relative locations in
[0, 1]). If a vector is provided, ancestry at these marker positions is
tracked for every generation.}

\item{verbose}{Verbose output if TRUE. Default value is FALSE}

\item{multiplicative_selection}{Default: TRUE. If TRUE, fitness is calculated
for multiple markers by multiplying fitness values for each marker. If FALSE,
fitness is calculated by adding fitness values for each marker.}

\item{migration_rate}{Rate of migration between the two populations.
Migration is implemented such that with probability m (migration rate) one
of the two parents of a new offspring is from the other population, with
probability 1-m both parents are of the focal population.}

\item{stop_at_critical_fst}{option to stop at a critical FST value
, default is FALSE}

\item{critical_fst}{the critical fst value to stop, if
\code{stop_simulation_at_critical_fst} is TRUE}

\item{generations_between_update}{The number of generations after which the
simulation has to check again whether the critical Fst value is exceeded}

\item{sampled_individuals}{Number of individuals to be sampled at random from
the population to estimate Fst}

\item{number_of_markers}{Number of markers to be used to estimate Fst}

\item{random_markers}{Are the markers to estimate Fst randomly distributed,
or regularly distributed? Default is TRUE.}

\item{mutation_rate}{the per base probability of mutation. Default is 0.}

\item{substitution_matrix}{a 4x4 matrix representing the probability of
mutating to another base (where [1/2/3/4] = [a/c/t/g]), conditional on the
event of a mutation happening. Default is the JC69 matrix, with equal
probabilities for all transitions / transversions.}
}
\value{
A list with: \code{population_1}, \code{population_2} two population
objects, and three tibbles with allele frequencies (only contain values of a
vector was provided to the argument \code{markers}: \code{frequencies},
\code{initial_frequencies} and \code{final_frequencies}. Each tibble contains
five columns, \code{time}, \code{location}, \code{ancestor}, \code{frequency}
and \code{population}, which indicates the number of generations, the
location along the chromosome of the marker, the ancestral allele at that
location in that generation, the frequency of that allele and the population
in which it was recorded (1 or 2). If a critical fst value was used to
terminate the simulation, and object \code{FST} with the final FST estimate
is returned as well.
}
\description{
Individual based simulation of the breakdown of contiguous
ancestry blocks, with or without selection. Simulations can be started from
scratch, or from a predefined input population. Two populations are
simulated, connected by migration
}
