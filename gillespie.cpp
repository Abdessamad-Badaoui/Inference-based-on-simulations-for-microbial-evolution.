#include "forward_simulator.h"
#define UPDATE_RATE_STEPS 1 // DO NOT PUT ANOTHER VALUE WITHOUT CHECKING EDGE CASES SUCH AS INITIAL LOW POPULATION SIZE AND **MUTANT EXTINCTION** WHICH HAS YET TO BE IMPLEMENTED

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
  std::uniform_real_distribution<> u(0., 1.);

  // Run the simulations
  unsigned long long target_size, mutants, sampled_mutants, wt, nevents;
  float death, mrate;
  int sim, index;
  unsigned long long step;
  
  double propensity_birth_wt, propensity_death_wt, propensity_birth_m, propensity_death_m, propensity_mut, propensity_total, rnext;
  
  for (sim = 0; sim < nsim; sim++) {
    mutants = 0;
    nevents = 0;
    wt = popsize[0];
    
    index = 0;
    while (index < n_segments) {
      target_size = popsize[index+1];
      death = deathrate[index];
      mrate = mutationrate[index];
       
      // Simulate until reaching target population size
      if (death<1.){
        if (target_size < mutants + wt) errx(1,"in step %d of population growth, population size is decreasing while death rate is lower than 1",index);
      }
      if (death>1.){
        if (target_size > mutants + wt) errx(1,"in step %d of population growth, population size is increasing while death rate is higher than 1",index);
      }
      if (death==1.) warnx("in step %d of population growth, death rate is 1, which will result in a quasi infinite loop",index);
      step=0;
      do {
        if ((wt==0)&&(mutants==0)){
            //errx(1,"All extinct\n");
            break;
        }
        if (step%UPDATE_RATE_STEPS==0){
            // Compute the relative rate of each process
            propensity_birth_wt = wt;
            propensity_death_wt = wt * death;
            propensity_birth_m = mutants * fmutant;
            propensity_death_m = mutants * fmutant * death;
            propensity_mut = wt * mrate;
            propensity_total = propensity_birth_wt + propensity_death_wt + propensity_birth_m + propensity_death_m + propensity_mut;
            assert(propensity_total > 0.);
            propensity_birth_wt /= propensity_total;
            propensity_death_wt /= propensity_total;
            propensity_birth_m /= propensity_total;
            propensity_death_m /= propensity_total;
            propensity_mut /= propensity_total;
        }
        step++;
        rnext = u(gen); // Which reaction occurs next?
        if (rnext < propensity_birth_wt) wt++;
        else if (rnext < propensity_birth_wt + propensity_death_wt) wt--;
        else if (rnext < propensity_birth_wt + propensity_death_wt + propensity_birth_m) mutants++;
        else if (rnext < propensity_birth_wt + propensity_death_wt + propensity_birth_m + propensity_death_m) mutants--;
        else {
            mutants++;
            nevents++;
        }
      } while (mutants + wt != target_size );
      index++;
    }
  
  // Perform binomial sampling if pfrac < 1
  sampled_mutants = mutants;
  if (pfrac < 1){
    std::binomial_distribution<int> bin(mutants, pfrac);
    sampled_mutants = bin(gen);
  } else if (pfrac > 1) errx(1, "pfrac is higher than 1");

  // Output the result of this simulation
  printf("%llu %llu %llu %llu\n",sampled_mutants,nevents,wt,wt+mutants);
  }

  return 0;
}

