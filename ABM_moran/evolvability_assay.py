from landscape_evolution import *
import matplotlib as mpl
import cProfile
import pstats
import time
import random
import seaborn as sns
import pandas as pd
import statistics
import pickle

try:
    import itertools.izip as zip
except ImportError:
    import itertools


# This is taken from Trevor Bedford's github. It is a Wright-Fisher simulation which includes: Mutation, Selection and Genetic Drift

def main():

    # Length of bit-sequence or number of mutational sites.
    seq_length = 6

    # How many generations should the population evolve in a neutral landscape?
    generations = 100

    # How many times would you like to simulate evolution?
    repeats = 10

    # Possible mutations at any site.
    alphabet = ['0', '1']

    #Use bin/zfill to produce every possible binary haplotype.
    base_haplotypes = []
    for i in range(2**seq_length):
        base_haplotypes.append(bin(i)[2:].zfill(seq_length))

    #Use itertools to produce every possible binary haplotype.
    genotypes = [''.join(seq) for seq in itertools.product("01", repeat=seq_length)]

    # Values we will be looping over to study evolvability
    # Population size and mutation rate.
    pop_sizes = [1, 1e1, 1e2, 1e3, 1e4]#, 1e5]
    mutation_rates = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4]#, 1e-3]#, .01]

    # Landscape to evolve in with no selection -- generating mutants.
    neutral_ls = [1]*(2**seq_length)

    # Dictionaries to store populations and fitnesses tied to haplotypes.
    pop = {}
    fitness = {}
    fitness_drug = {}

    #Fill in the neutral landscape with equal fitnesses. There's a faster way to do this.
    for i in range(len(genotypes)):
        fitness[genotypes[i]] = neutral_ls[i]

    #Create array to store evolvability data.
    # mean = evolvability average over each starting haplotype.
    # std = evolvability std over each starting haplotype.
    evolvability_mean = np.zeros((len(pop_sizes), len(mutation_rates)))
    evolvability_std = np.zeros((len(pop_sizes), len(mutation_rates)))

    # Epistasis for drug-landscape.
    sigma = 1

    # Loop over all repeats identified at the top of the file.
    for i in range(repeats):
        print("Repeat: {}".format(i))

        # Generate a drug landscape to test in.
        # Rough mount fuji with epistasis = sigma.
        # Normalize landscape between 0 and 1, add small number to avoid 0 replication for global minimum.
        drug_ls = Landscape(seq_length, sigma)
        drug_ls.ls = ((drug_ls.ls-np.min(drug_ls.ls))/(np.max(drug_ls.ls)-np.min(drug_ls.ls))) + 1e-5

        # Clear fitness of previous loops data.
        fitness_drug.clear()

        # Assign each haplotype the corresponding fitness from the drug landscape.
        for i in range(len(genotypes)):
            fitness_drug[genotypes[i]] = drug_ls.ls[i]

        # Loop over all population sizes selected.
        pop_size_counter = 0
        for pop_size in pop_sizes:
            print("Population size: {}".format(pop_size))

            # Loop over all mutation rates selected.
            mut_rate_counter = 0
            for mutation_rate in mutation_rates:
                print("Mutation rate: {}".format(mutation_rate))

                # Loop over all haplotypes in the landscape.
                # Create an emply list to append all haplotype-specific data to be averaged.
                base_haplo_evo = []
                for base_haplotype in base_haplotypes:
                    print("Haplotype: {}".format(base_haplotype))

                    # Clear pop dictionary from previous loop.
                    # Assign population to begin in specified base_haplotype.
                    pop.clear()
                    pop[base_haplotype] = pop_size

                    # Clear history before simulation
                    # Run a Wright-Fisher simulation for the specified values.
                    history = []
                    simulate(pop, history, generations, mutation_rate, pop_size, seq_length, fitness, alphabet)

                    # Find the haplotype with the maximum fitness that is occupied by the simulation.
                    max_geno = max(pop.keys(), key=(lambda k : fitness_drug[k]))

                    # Find the associated fitness in the drug-landscape for that genotype and append to storage list.
                    base_haplo_evo.append(fitness_drug[max_geno])

                # Take the mean and standard deviation of the haplotype-specific data and store it for future.
                evolvability_mean[pop_size_counter, mut_rate_counter] += statistics.mean(base_haplo_evo)
                evolvability_std[pop_size_counter, mut_rate_counter] += statistics.pstdev(base_haplo_evo)

                mut_rate_counter += 1
            pop_size_counter +=1

    # Average evolvability data over the number of selected repeats.
    evolvability_mean = evolvability_mean/repeats
    evolvability_std = evolvability_std/repeats

    # Dump the datas into a pickle file for future analysis.
    #with open('mean_Evolvability.pkl', 'wb') as f:
    #    pickle.dump(evolvabilty_mean, f)
    #with open('std_Evolvability.pkl', 'wb') as f:
    #    pickle.dump(evolvability_std, f)

    # Plot basic values for the selected data.

    # Plots the mean haplotype data.
    plt.figure()
    plt.subplot(1,2,1)
    ax1 = sns.heatmap(evolvability_mean, xticklabels=mutation_rates, yticklabels=pop_sizes, cmap='YlGnBu', annot=True)
    ax1.set_xlabel('mutation rate')
    ax1.set_ylabel('pop_size')
    ax1.invert_yaxis()

    # Plots the std haplotype data.
    plt.subplot(1,2,2)
    ax2 = sns.heatmap(evolvability_std, xticklabels=mutation_rates, yticklabels=pop_sizes, cmap='YlGnBu', annot=True)
    ax2.set_xlabel('mutation rate')
    ax2.set_ylabel('pop_size')
    ax2.invert_yaxis()
    plt.show()




if __name__ == '__main__':
    main()
