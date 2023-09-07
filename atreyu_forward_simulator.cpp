#include "forward_simulator.h"

int main(int argc, char* argv[]) {

  int n_segments;
  unsigned long long popsize[100];
  float deathrate[100];
  float mutationrate[100];
  float fmutant = 1;
  float pfrac = 1;
  int nsim = 1;
  int seed;
  int verbose = 0;
  int constant_mutation_rate = 0;

  parse_arguments(argc, argv, &n_segments, popsize, deathrate, mutationrate, &fmutant, &pfrac, &nsim, &seed, &verbose, &constant_mutation_rate);
  
  // Random number generator: initialize
  std::mt19937 gen(seed);
  //std::minstd_rand gen(seed);
  std::uniform_real_distribution<> u(0., 1.); // uniform distribution between 0 and 1: will be used for Gillespie

  // Run the simulations
  unsigned long long target_size, mutants, sampled_mutants, wt, nevents, remaining_birth_wt;
  float death, mrate;

  unsigned long long n_birth_wt, n_death_wt;
  double propensity_birth_mutant, propensity_death_mutant, propensity_total, time_to_mutation, ellapsed_time;
  int sim, index, growth_direction;

  for (sim = 0; sim < nsim; sim++) {
    mutants = 0;
    nevents = 0;
    remaining_birth_wt = 0;
    wt = popsize[0];
    
    index = 0;
    while (index < n_segments) {
      target_size = popsize[index+1];
      death = deathrate[index];
      mrate = mutationrate[index];
       
      // Simulate until reaching target population size

      // Check consistency of user-provided death rates and population sizes
      if (death<1.){
        if (target_size < mutants + wt) errx(1,"in step %d of population growth, population size is decreasing while death rate is lower than 1",index);
        growth_direction = 1;
      }
      if (death>1.){
        errx(1,"not implemented yet");
        if (target_size > mutants + wt) errx(1,"in step %d of population growth, population size is increasing while death rate is higher than 1",index);
        growth_direction = -1;
      }
      if (death==1.) warnx("in step %d of population growth, death rate is 1, which will result in a quasi infinite loop",index);
      while (1) {
        
        /* Draw the number of wild-type divisions until the next mutation
        Ce sample from an exponential distribution with parameter mrate (waiting time to next event in a Poisson process)
        Contrary to the canonical Gillespie implementation, the waiting time is here not expressed as real time (a continuous variable) but as number of divisions of the wild-type (a discrete variable), and thus a geometric distribution would be more accurate
        */
        if (remaining_birth_wt==0) n_birth_wt = ceil(-log(u(gen))/mrate); else{
          // remaining WT divisions from last growth phase before next mutation (see below for more explanations)
          n_birth_wt = remaining_birth_wt;
          if (verbose) printf("  (remaining from last step)");
        }
        if (verbose) printf("  %llu WT birth until next mutation\n", n_birth_wt);
        
        // Compute the number of death of the WT during this time step
        #ifdef STOCHASTIC_WT_DEATH
        std::poisson_distribution<unsigned long long> dWT(death * n_birth_wt);
        n_death_wt = dWT(gen);
        #else
        n_death_wt = ceil(death * n_birth_wt);  // no need for poisson sampling, fully deterministic computations is sufficient for the WT
        #endif
        if (verbose) printf("  %llu WT death until next mutation\n", n_death_wt);

        /* Are we going to cross the next population size threshold? 
        If so, we need to stop *before* this mutation: only part of the WT divisions are performed, and the mutation is not performed
        But since the time of occurence of this mutation (number of remaining division events) has already been determined, we may want to keep it for next growth phase (if any).
        If we just discard it and draw again time to occurence during next growth phase, this could introduce a bias: long periods before occurence would have more chances of being discarded and drawn again than short ones. 
        However, if mutation rate changes between growth phases, then we have no other choice than discard and draw again. Comparaisons with Gillespie simulations show that the resulting potential bias is to small to be detected.
        */
        if ( growth_direction * (wt + mutants + n_birth_wt - n_death_wt) > growth_direction * target_size ) {
          unsigned long long new_n_birth_wt = ceil((target_size - wt - mutants) / (1-death));
          remaining_birth_wt = n_birth_wt - new_n_birth_wt;
          n_birth_wt = new_n_birth_wt;
          #ifdef STOCHASTIC_WT_DEATH
          std::poisson_distribution<unsigned long long> dWT(death * n_birth_wt);
          n_death_wt = dWT(gen);
          #else
          n_death_wt = ceil(death * n_birth_wt);
          #endif
          if (verbose) printf("  only performing %llu WT birth and %llu WT death as next target WT size is reached\n", n_birth_wt, n_death_wt);
        } else remaining_birth_wt = 0;
        
        // Compute the ellapsed time scaled by division rate of the WT, ie the number of synchronous generations for the WT
        // This implicitely implies a deterministic growth for the WT
        #ifdef STOCHASTIC_WT_DEATH
        #error "STOCHASTIC_WT_DEATH not fully implemented"
        #endif
        time_to_mutation = ( log(wt+n_birth_wt-n_death_wt)-log(wt) ) / ( 1 - death );
        if (verbose) printf("  ellasped time: %f\n", time_to_mutation);
        
        // Simulate population dynamics of the mutant with gillespie until we reach this time
        ellapsed_time = 0;
        while (1){
          if (mutants == 0) break; // Mutant population extinct or first mutation not yet arised
          // Compute the propensity of each reaction (and normalize by total reaction propensity)
          // All propensities are relative to growth rate of the wild-type
          propensity_birth_mutant = mutants * fmutant;
          propensity_death_mutant = mutants * fmutant * death;
          propensity_total = propensity_birth_mutant + propensity_death_mutant;
          propensity_birth_mutant /= propensity_total;
          propensity_death_mutant /= propensity_total;
          ellapsed_time += -log(u(gen)) / propensity_total; // draw waiting time to next event
          if(ellapsed_time < time_to_mutation){
            if (u(gen) < propensity_birth_mutant) mutants++; else mutants--;
          } else break; // no further event before the end of the inner simulation
        }
        
        // Actually update the variables 
        wt += n_birth_wt - n_death_wt;
        if (verbose) printf("  new number of mutants before the mutation: %llu\n", mutants);

        // Perform the mutation if we did not decide to stop before (when population threshold is reached)
        if (remaining_birth_wt != 0){
          assert(remaining_birth_wt > 0);
          if (!constant_mutation_rate) remaining_birth_wt = 0; // This may introduce a small bias, but otherwise we have a problem when mutation rate changes between two segments
          break;
        } else {
            mutants++;
            nevents++;
        }
      }
      index++;
    }

  // Perform binomial sampling if pfrac < 1
  sampled_mutants = mutants;
  if (pfrac < 1){
    std::binomial_distribution<int> bin(mutants, pfrac);
    sampled_mutants = bin(gen);
    if (verbose) printf("sampling %llu mutants of a %llu mutant population\n", sampled_mutants, mutants);
  } else if (pfrac > 1) errx(1, "pfrac is higher than 1");

  // Output the result of this simulation
  printf("%llu %llu %llu %llu\n",sampled_mutants,nevents,wt,wt+mutants);
  }

  return 0;
}

