#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <random>
#include <cassert>
#include <cstring>

void print_help(const char* progname){
    printf("Forward simulator for atreyu\n"
    "Simulates the growth of a population with birth, death, and accumulation of a mutation\n"
    "Usage:\n\t%s 'N0,N1,..,Nf' 'd0,d1,..,df-1' 'm0,m2,..,mf-1' [fmutant] [pfrac] [nrep] [seed] ['verbose']\n"
    "\t\tNi is population size at step i\n"
    "\t\tdi and mi are relative death rate and mutation rate between steps i and i+1\n"
    "\t\tfmutant is the fitness of the mutant relative to the wild-type (default = 1)\n"
    "\t\tpfrac is the plating (sampling) fraction (default = 1)\n"
    "\t\tnrep is the number of replicate simulations to run (default = 1)\n"
    "\t\tseed can be provided for perfect replication\n"
    "\t or %s help\n"
    "Output: one line per replicate simulation, containing\n"
    "\t 'number of mutants' 'number of mutational events' 'number of wt' 'total population size'\n", progname, progname);
}

void parse_arguments(int argc, char* argv[], int* n_segments, unsigned long long popsize[], float deathrate[], float mutationrate[], float* fmutant, float* pfrac, int* nsim, int* seed, int* verbose, int* constant_mutation_rate){
    // Check that required number of mandatory arguments is provided
    if (argc < 4) {
        print_help(argv[0]);
        exit(0);
    }

    // Parse the first argument into a list of population size
    char* a1 = argv[1];
    int index = 0;
    const char separator[] = ","; 
    char* comma = strtok(a1, separator);
    while (comma) {
        unsigned long long value = (long long) atof(comma);
        popsize[index++] = value;
        comma = strtok(NULL, separator);
    }
    int nsizes = index;


    // Parse the second argument into a list of death rates
    char* a2 = argv[2];
    index = 0;
    comma = strtok(a2, separator);
    while (comma) {
        double value = atof(comma);
        deathrate[index++] = value;
        comma = strtok(NULL, separator);
    }
    int ndrates = index;

    // Parse the third argument into a list of mutation rates
    char* a3 = (argv[3]);
    index = 0;
    comma = strtok(a3, separator);
    while (comma) {
        double value = atof(comma);
        mutationrate[index++] = value;
        comma = strtok(NULL, separator);
    }
    *n_segments = index;
    #ifdef CONSTANT_MUTATION_RATE
    double m = mutationrate[0];
    *constant_mutation_rate = 1;
    for (int s=1; s<*n_segments; s++){
        if(mutationrate[s]!=m)
        {
            *constant_mutation_rate = 0;
            fprintf(stderr,"Mutation rate changes between growth phases, so CONSTANT_MUTATION_RATE will be ignored\n");
            break;
        }
    }
    #endif

    // Check consistency between the first 3 arguments (population size, death rates and mutation rates)
    if (nsizes != ndrates + 1) errx(1,"the number of population sizes must be the number of death rates + 1");
    if (nsizes != (*n_segments) + 1) errx(1,"the number of population sizes must be the number of mutation rates + 1");
    if (nsizes <= 1) errx(1,"at least 2 times points (2 population sizes and 1 death rate) must be provided");

    // Parse the fourth argument (if existing) into relative fitness of the mutant
    if (argc >= 5) *fmutant = atof(argv[4]);

    // Parse the fifth argument (if existing) into plating fraction
    if (argc >= 6) *pfrac = atof(argv[5]);

    // Parse the sixth argument (if existing) into a number of simulations to run
    if (argc >= 7) *nsim = atoi(argv[6]);

    // Parse the seventh argument (if existing) into the random seed
    if (argc >= 8) *seed = atoi(argv[7]); else {
        std::random_device rd;
        *seed = rd();
    }

    // Parse "verbose" and extra arguments
    if (argc == 9) {
        if (memchr(argv[8], 'v', strlen(argv[8]))) *verbose = 1; else warnx("extra arguments ignored");
    }
    if (argc > 9) warnx("extra arguments ignored");

    // Display read parameters when 'verbose'
    if (*verbose){
        printf("Population sizes: ");
        for (int i = 0; i < nsizes; i++) printf("%llu ",popsize[i]);
        printf("\n");

        printf("Death rates: ");
        for (int i = 0; i < ndrates; i++) printf("%.3f ",deathrate[i]);
        printf("\n");

        printf("Mutation rate: ");
        for (int i = 0; i < *n_segments; i++) printf("%.2e ", mutationrate[i]);
        printf("\n");

        printf("Fitness of the mutant: %f\n", *fmutant);
        printf("Plating fraction: %f\n", *pfrac);
        printf("Number of simulations: %d\n", *nsim);
        printf("Random seed: %d\n", *seed);
    }
}