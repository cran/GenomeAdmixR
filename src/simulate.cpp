//
//  selection.cpp
//
//
//  Created by Thijs Janzen on 28/02/2018.
//
//
#include <vector>
#include <cstdlib>
#include <numeric>
#include <cmath>

#include <vector>
#include <algorithm>

#include <thread>

#include "Fish.h"
#include "random_functions.h"
#include "helper_functions.h"
#include "util.h"

#include <RcppParallel.h>

#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;

void update_pop(const std::vector<Fish>& Pop,
                std::vector<Fish>& new_generation,
                size_t pop_size,
                double morgan,
                const std::vector<double>& fitness,
                const double& maxFitness,
                bool use_selection,
                double multiplicative_selection,
                int num_threads) {

  if (Pop.size() != pop_size) {
    stop("wrong size pop");
  }
  if (new_generation.size() != pop_size) {
    stop("new_generation wrong size");
  }


  int num_seeds = num_threads * 10; // tbb might re-start threads due to the load-balancer
  if (num_threads == -1) {
    num_seeds = 200;
  }
  std::vector< int > seed_values(num_seeds);

  {
    rnd_t rndgen;
    for (int i = 0; i < num_seeds; ++i) {
      seed_values[i] = rndgen.random_number(INT_MAX); // large value
    }
  }

  if (num_threads == 1) {
    rnd_t rndgen2;
    for (unsigned i = 0; i < pop_size; ++i) {
      int index1 = 0;
      int index2 = 0;
      if (use_selection) {
        index1 = draw_prop_fitness(fitness, maxFitness, rndgen2);
        index2 = draw_prop_fitness(fitness, maxFitness, rndgen2);
        while(index2 == index1) index2 = draw_prop_fitness(fitness, maxFitness, rndgen2);
      } else {
        index1 = rndgen2.random_number( pop_size );
        index2 = rndgen2.random_number( pop_size );
        while(index2 == index1) index2 = rndgen2.random_number( pop_size );
      }

      new_generation[i] = mate(Pop[index1],
                               Pop[index2],
                                  morgan, rndgen2);
    }
  } else {

    set_num_threads();

    tbb::parallel_for(
      tbb::blocked_range<unsigned>(0, pop_size),
      [&](const tbb::blocked_range<unsigned>& r) {

      size_t local_seed = std::hash<std::thread::id>{}(std::this_thread::get_id());
      size_t local_time = static_cast<unsigned int>( time(NULL) );
      thread_local rnd_t rndgen2(local_seed + local_time);

        for (unsigned i = r.begin(); i < r.end(); ++i) {
          int index1 = 0;
          int index2 = 0;
          if (use_selection) {
            index1 = draw_prop_fitness(fitness, maxFitness, rndgen2);
            index2 = draw_prop_fitness(fitness, maxFitness, rndgen2);
            while(index2 == index1) index2 = draw_prop_fitness(fitness, maxFitness, rndgen2);
          } else {
            index1 = rndgen2.random_number( pop_size );
            index2 = rndgen2.random_number( pop_size );
            while(index2 == index1) index2 = rndgen2.random_number( pop_size );
          }

          new_generation[i] = mate(Pop[index1],
                                   Pop[index2],
                                      morgan, rndgen2);
        }
      });
  }
  return;
}

std::vector< Fish > simulate_Population(const std::vector< Fish>& sourcePop,
                                        const NumericMatrix& select_r,
                                        size_t pop_size,
                                        int total_runtime,
                                        double morgan,
                                        bool verbose,
                                        arma::mat& frequencies,
                                        bool track_frequency,
                                        const NumericVector& track_markers,
                                        bool track_junctions,
                                        std::vector<double>& junctions,
                                        bool multiplicative_selection,
                                        int num_alleles,
                                        const std::vector<int>& founder_labels,
                                        rnd_t& rndgen,
                                        int num_threads) {

  bool use_selection = false;
  if(select_r(0, 1) >= 0) use_selection = true;

  std::vector<Fish> Pop = sourcePop;
  std::vector<double> fitness;
  double maxFitness = -1;
  std::vector< std::array<double, 5>> select;
  if(use_selection) {
    for (size_t i = 0; i < select_r.nrow(); ++i) {
      std::array<double, 5> entry;
      for (size_t j = 0; j < 5; ++j) {
        entry[j] = select_r(i, j);
      }
      select.push_back(entry);
    }

    for(auto it = Pop.begin(); it != Pop.end(); ++it){
      double fit = calculate_fitness((*it), select, multiplicative_selection);
      if(fit > maxFitness) maxFitness = fit;

      fitness.push_back(fit);
    }
  }

  int updateFreq = total_runtime / 20;
  if(updateFreq < 1) updateFreq = 1;

  if(verbose) {
    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";
  }

  for(int t = 0; t < total_runtime; ++t) {
    if(track_junctions) junctions.push_back(calc_mean_junctions(Pop));

    if(track_frequency) {
      for(int i = 0; i < track_markers.size(); ++i) {
        if(track_markers[i] < 0) break;
        arma::mat local_mat = update_frequency_tibble(Pop,
                                                      track_markers[i],
                                                                   founder_labels,
                                                                   t,
                                                                   morgan);

        // now we have to find where to copy local_mat into frequencies
        int time_block = track_markers.size() * founder_labels.size(); // number of markers times number of alleles

        int start_add_time = t * time_block;
        int start_add_marker = i * founder_labels.size() + start_add_time;

        for(size_t j = 0; j < founder_labels.size(); ++j) {
          for(size_t k = 0; k < 4; ++k) {
            frequencies(start_add_marker + j, k)  = local_mat(j, k);
          }
        }
      }
    }

    std::vector<Fish> new_generation(pop_size);

    update_pop(Pop, new_generation, pop_size,
               morgan, fitness, maxFitness,  use_selection,
               multiplicative_selection, num_threads);

    if (t % updateFreq == 0 && verbose) {
      Rcout << "**";
    }

    if (t > 2 && is_fixed(Pop)) {
      if (verbose) Rcout << "\n After " << t << " generations, the population has become completely homozygous and fixed\n";
      R_FlushConsole();
      return(Pop);
    }

    Rcpp::checkUserInterrupt();
    Pop.swap(new_generation);
    if (use_selection) {
      for (size_t i = 0; i < pop_size; ++i) {
        fitness[i] = calculate_fitness(Pop[i], select, multiplicative_selection);
      }
      maxFitness = *std::max_element(fitness.begin(), fitness.end());
    }
  }
  if(verbose) Rcout << "\n";
  return(Pop);
}

// [[Rcpp::export]]
List simulate_cpp(Rcpp::NumericVector input_population,
                  NumericMatrix select,
                  size_t pop_size,
                  size_t number_of_founders,
                  Rcpp::NumericVector starting_proportions,
                  size_t total_runtime,
                  double morgan,
                  bool verbose,
                  bool track_frequency,
                  NumericVector track_markers,
                  bool track_junctions,
                  bool multiplicative_selection,
                  int num_threads) {
  try {
    rnd_t rndgen;

    std::vector< Fish > Pop;
    int number_of_alleles = number_of_founders;
    std::vector<int> founder_labels;

    track_markers = scale_markers(track_markers, morgan);

    if (input_population[0] > -1e4) {
      if (verbose) { Rcout << "found input population, converting\n"; force_output(); }

      Pop = convert_NumericVector_to_fishVector(input_population);

      number_of_founders = 0;
      for (auto it = Pop.begin(); it != Pop.end(); ++it) {
        update_founder_labels((*it).chromosome1, founder_labels);
        update_founder_labels((*it).chromosome2, founder_labels);
      }

      number_of_alleles = founder_labels.size();

      if (Pop.size() != pop_size) {
        // the new population has to be seeded from the input!
        std::vector< Fish > Pop_new;
        for (size_t j = 0; j < pop_size; ++j) {
          int index = rndgen.random_number(Pop.size());
          Pop_new.push_back(Pop[index]);
        }
        Pop = Pop_new;
      }
    } else {
      if (verbose) { Rcout << "starting generating random starting population\n"; force_output(); }
      for (size_t i = 0; i < pop_size; ++i) {
        int founder_1 = draw_random_founder(starting_proportions, rndgen);
        int founder_2 = draw_random_founder(starting_proportions, rndgen);

        Fish p1 = Fish( founder_1 );
        Fish p2 = Fish( founder_2 );

        Pop.push_back(mate(p1,p2, morgan, rndgen));
      }
      for (int i = 0; i < number_of_alleles; ++i) {
        founder_labels.push_back(i);
      }
    }

    arma::mat frequencies_table;

    if (track_frequency) {
      int number_of_markers = track_markers.size();
      arma::mat x(number_of_markers * number_of_alleles * total_runtime, 4); // 4 columns: time, loc, anc, type
      frequencies_table = x;
    }

    arma::mat initial_frequencies = update_all_frequencies_tibble(Pop,
                                                                  track_markers,
                                                                  founder_labels,
                                                                  0,
                                                                  morgan);

    std::vector<double> junctions;
    std::vector<Fish> outputPop = simulate_Population(Pop,
                                                      select,
                                                      pop_size,
                                                      total_runtime,
                                                      morgan,
                                                      verbose,
                                                      frequencies_table,
                                                      track_frequency,
                                                      track_markers,
                                                      track_junctions,
                                                      junctions,
                                                      multiplicative_selection,
                                                      number_of_alleles,
                                                      founder_labels,
                                                      rndgen,
                                                      num_threads);

    if (verbose) {
      Rcout << "done simulating\n";
    }
    arma::mat final_frequencies = update_all_frequencies_tibble(outputPop,
                                                                track_markers,
                                                                founder_labels,
                                                                total_runtime,
                                                                morgan);

    return List::create( Named("population") = convert_to_list(outputPop),
                         Named("frequencies") = frequencies_table,
                         Named("initial_frequencies") = initial_frequencies,
                         Named("final_frequencies") = final_frequencies,
                         Named("junctions") = junctions);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}
