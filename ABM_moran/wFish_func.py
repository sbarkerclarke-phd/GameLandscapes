import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import random
import cProfile
import pstats
from landscape_evolution import * 
import time
try:
    import itertools.izip as zip
except ImportError:
    import itertools



"""
Simulate Wright-Fisher evolution for some amount of generations
"""
 ## All games same type ##
def initialise():
    return(0)

def simulate_game(pop, game, generations, mutation_rate, pop_size, seq_length, A, alphabet):
    pop = {}
    fitness = {}
    #fitnessB = {}
    #base_haplotype = ''.join(["0" for i in range(seq_length)])
    base_haplotype = '01'+'0'*(seq_length-2)
    genotypes = [''.join(seq) for seq in itertools.product("01", repeat=seq_length)]

    pop[base_haplotype] = pop_size
    #pop2 = copy.deepcopy(pop)

    #Allocate landscape values to fitness dictionary
    for i in range(len(genotypes)):
        fitness[genotypes[i]] = A.ls[i]
        #fitnessB[genotypes[i]] = Bs[46].ls[i]
        
    #Generate random landscape
    #print(fitness)
    history = []
    simulate(pop, game, history, generations, mutation_rate, pop_size, seq_length, fitness, alphabet)
    print(fitness)
    return(history)
    

def modify_fitness(pop, fitness, game):

        #Define game type
        quadrant = game

        N_geno = len(fitness.keys())
        N_pairs = int(N_geno*(N_geno-1)/2)
        #N_T = sum(counts)
        N=0
        game_landscapes = {}

        res = list(itertools.combinations(fitness, 2))
        #print(res)

        for i in fitness.keys():
            game_landscapes[i]={}
            if pop.get(i)==None:
                pop[i]=0
        
        #print(res)
        
        for x,y in res:
            A,D = (fitness[x], fitness[y])
            if (quadrant==None):
                A,B,C,D = (0,0,0,0)
            elif(quadrant==1):
                A,B,C,D = (A,2*D,A/2,D)
            elif(quadrant==2):
                A,B,C,D = (A,D/2,A/2,D)
            elif(quadrant==3):
                A,B,C,D = (A,2*D,2*A,D)
            elif(quadrant==4):
                A,B,C,D = (A,D/2,2*A,D)

            #print(A,B,C,D)
            
            N_T = pop[x]+pop[y]
            if (N_T>0):
                game_landscapes[x][y] = A*pop[x]/N_T+ B*pop[y]/N_T
                game_landscapes[y][x] = C*pop[x]/N_T + D*pop[y]/N_T
                N=N+1
            else:
                game_landscapes[x][y]= game_landscapes[y][x]=0
                
        #print(game_landscapes['01'])

        fitness_update = {}
        
        for i in fitness.keys():
            fitness_update[i]= sum(game_landscapes[i].values())/ (sum(fitness.values())*2)

        if (quadrant==None):
            fitness_update = fitness
                
        return(fitness_update)
    

def simulate(pop, game, history ,generations, mutation_rate, pop_size, seq_length, fitness, alphabet):
    clone_pop = dict(pop) #Population distribution
    history.append(clone_pop) #Append old pop to history
    for i in range(generations):
        fitness = modify_fitness(pop,fitness, game)
        #print(pop)
        print(fitness)
        time_step(pop, mutation_rate, pop_size, seq_length, fitness, alphabet)
        clone_pop = dict(pop)
        history.append(clone_pop)
    print(pop)

def simulateSwitch(pop, history ,generations, mutation_rate, pop_size, seq_length, fitnessA, fitnessB, alphabet):
    clone_pop = dict(pop)
    history.append(clone_pop)

    count = 0
    for i in range(generations):
        if (count < 1):
            time_step(pop, mutation_rate, pop_size, seq_length, fitnessA, alphabet)
            clone_pop = dict(pop)
            history.append(clone_pop)
            count += 1
        else:
            time_step(pop, mutation_rate, pop_size, seq_length, fitnessB, alphabet)
            clone_pop = dict(pop)
            history.append(clone_pop)
            count += 1
            if (count > 2):
                count = 0

"""
Execuse one generation of the Wright-Fisher
"""
def time_step(pop, mutation_rate, pop_size, seq_length, fitness, alphabet):
    mutation_step(pop, mutation_rate, pop_size, seq_length, alphabet)
    offspring_step(pop, pop_size, fitness)


#################################################################################################################################
"""
Below is the code responsible for mutation in the Wright-Fisher model, in order of function calls.
"""

"""
First step of mutation -- get a count of the mutations that occur.
"""
def mutation_step(pop, mutation_rate, pop_size, seq_length, alphabet):
    mutation_count = get_mutation_count(mutation_rate, pop_size, seq_length)
    for i in range(mutation_count):
        mutation_event(pop, pop_size, seq_length, alphabet)

"""
Draw mutation count from a poisson distribution with mean equal to average number of expected mutations.
"""
def get_mutation_count(mutation_rate, pop_size, seq_length):
    mean = mutation_rate * pop_size * seq_length
    return np.random.poisson(mean)


"""
Function that find a random haplotype to mutate and adds that new mutant to the population. Reduces mutated population by 1.
"""
def mutation_event(pop, pop_size, seq_length, alphabet):
    haplotype = get_random_haplotype(pop, pop_size)
    if pop[haplotype] > 1:
        pop[haplotype] -= 1
        new_haplotype = get_mutant(haplotype, seq_length, alphabet)
        if new_haplotype in pop:
            pop[new_haplotype] += 1
        else:
            pop[new_haplotype] = 1

"""
Chooses a random haplotype in the population that will be returned.
"""
def get_random_haplotype(pop, pop_size):
    haplotypes = list(pop.keys())
    frequencies = [x/pop_size for x in pop.values()]
    total = sum(frequencies)
    frequencies = [x / total for x in frequencies]
    return fast_choice(haplotypes, frequencies)
    #return random.choices(haplotypes, weights=frequencies)[0]

"""
Receives the haplotype to be mutated and returns a new haplotype with a mutation with all neighbor mutations equally probable.
"""
def get_mutant(haplotype, seq_length, alphabet):
    site = int(random.random()*seq_length)
    possible_mutations = list(alphabet)
    possible_mutations.remove(haplotype[site])
    mutation = random.choice(possible_mutations)
    new_haplotype = haplotype[:site] + mutation + haplotype[site+1:]
    return new_haplotype


#################################################################################################################################
"""
Below is the code responsible for offspring in the Wright-Fisher model, in order of function calls.
"""


"""
Gets the number of counts after an offspring step and stores them in the haplotype. If a population is reduced to zero then delete it.
"""
def offspring_step(pop, pop_size, fitness):
    haplotypes = list(pop.keys())
    counts = get_offspring_counts(pop, pop_size, fitness)
    for (haplotype, count) in zip(haplotypes, counts):
        if (count > 0):
            pop[haplotype] = count
        else:
            del pop[haplotype]

"""
Returns the new population count for each haplotype given offspring counts weighted by fitness of haplotype
"""
def get_offspring_counts(pop, pop_size, fitness):
    haplotypes = list(pop.keys())
    frequencies = [pop[haplotype]/pop_size for haplotype in haplotypes]
    fitnesses = [fitness[haplotype] for haplotype in haplotypes]
    weights = [x * y for x,y in zip(frequencies, fitnesses)]
    total = sum(weights)
    weights = [x / total for x in weights]
    return list(np.random.multinomial(pop_size, weights))



"""
This is faster than numpy because numpy is dog water. Self-written code to choose an item from a list with some prob.
"""
def fast_choice(options, probs):
    x = random.random()
    cum = 0
    for i, p in enumerate(probs):
        cum += p
        if x < cum:
            return options[i]
    return options[-1]


#################################################################################################################################
"""
Below are helper functions for plotting//quantification
"""
def get_distance(seq_a, seq_b):
    diffs = 0
    length = len(seq_a)
    assert len(seq_a) == len(seq_b)
    for chr_a, chr_b in zip(seq_a, seq_b):
        if chr_a != chr_b:
            diffs += 1
    return diffs / float(length)

def get_diversity(population, pop_size):
    haplotypes = list(population.keys())
    haplotype_count = len(haplotypes)
    diversity = 0
    for i in range(haplotype_count):
        for j in range(haplotype_count):
            haplotype_a = haplotypes[i]
            haplotype_b = haplotypes[j]
            frequency_a = population[haplotype_a] / float(pop_size)
            frequency_b = population[haplotype_b] / float(pop_size)
            frequency_pair = frequency_a * frequency_b
            diversity += frequency_pair * get_distance(haplotype_a, haplotype_b)
    return diversity

def get_diversity_trajectory(history, pop_size):
    trajectory = [get_diversity(generation, pop_size) for generation in history]
    return trajectory

def diversity_plot(history, pop_size):
    mpl.rcParams['font.size']=14
    trajectory = get_diversity_trajectory(history, pop_size)
    plt.plot(trajectory, "#447CCD")
    plt.ylabel("diversity")
    plt.xlabel("generation")

def get_divergence(population, base_haplotype, pop_size):
    haplotypes = population.keys()
    divergence = 0
    for haplotype in haplotypes:
        frequency = population[haplotype] / float(pop_size)
        divergence += frequency * get_distance(base_haplotype, haplotype)
    return divergence

def get_divergence_trajectory(history, base_haplotype, pop_size):
    trajectory = [get_divergence(generation, base_haplotype, pop_size) for generation in history]
    return trajectory

def divergence_plot(history, base_haplotype, pop_size):
    mpl.rcParams['font.size']=14
    trajectory = get_divergence_trajectory(history, base_haplotype, pop_size)
    plt.plot(trajectory, "#447CCD")
    plt.ylabel("divergence")
    plt.xlabel("generation")

def get_frequency(haplotype, generation, history, pop_size):
    pop_at_generation = history[generation]
    if haplotype in pop_at_generation:
        return pop_at_generation[haplotype]/float(pop_size)
    else:
        return 0

def get_trajectory(haplotype, generations, history, pop_size):
    trajectory = [get_frequency(haplotype, gen, history, pop_size) for gen in range(generations)]
    return trajectory

def get_all_haplotypes(history):
    haplotypes = set()
    for generation in history:
        for haplotype in generation:
            haplotypes.add(haplotype)
    return haplotypes

def stacked_trajectory_plot(history, generations, pop_size, xlabel="generation"):
    colors_lighter = ["#A567AF", "#8F69C1", "#8474D1", "#7F85DB", "#7F97DF", "#82A8DD", "#88B5D5", "#8FC0C9", "#97C8BC", "#A1CDAD", "#ACD1A0", "#B9D395", "#C6D38C", "#D3D285", "#DECE81", "#E8C77D", "#EDBB7A", "#EEAB77", "#ED9773", "#EA816F", "#E76B6B"]
    mpl.rcParams['font.size']=18
    haplotypes = get_all_haplotypes(history)
    trajectories = [get_trajectory(haplotype, generations, history, pop_size) for haplotype in haplotypes]
    plt.stackplot(range(generations), trajectories, colors=colors_lighter)
    plt.ylim(0, 1)
    plt.ylabel("frequency")
    plt.xlabel(xlabel)

def get_snp_frequency(site, generation, history, pop_size):
    minor_allele_frequency = 0.0
    pop_at_generation = history[generation]
    for haplotype in pop_at_generation.keys():
        allele = haplotype[site]
        frequency = pop_at_generation[haplotype] / float(pop_size)
        if allele != "0":
            minor_allele_frequency += frequency
    return minor_allele_frequency

def get_snp_trajectory(site, generations, history, pop_size):
    trajectory = [get_snp_frequency(site, gen, history, pop_size) for gen in range(generations)]
    return trajectory

def get_all_snps(history, seq_length):
    snps = set()
    for generation in history:
        for haplotype in generation:
            for site in range(seq_length):
                if haplotype[site] != "0":
                    snps.add(site)
    return snps

def snp_trajectory_plot(history, seq_length, generations, pop_size, xlabel="generation"):
    colors = ["#781C86", "#571EA2", "#462EB9", "#3F47C9", "#3F63CF", "#447CCD", "#4C90C0", "#56A0AE", "#63AC9A", "#72B485", "#83BA70", "#96BD60", "#AABD52", "#BDBB48", "#CEB541", "#DCAB3C", "#E49938", "#E68133", "#E4632E", "#DF4327", "#DB2122"]
    mpl.rcParams['font.size']=18
    snps = get_all_snps(history, seq_length)
    trajectories = [get_snp_trajectory(snp, generations, history, pop_size) for snp in snps]
    data = []
    for trajectory, color in zip(trajectories, itertools.cycle(colors)):
        data.append(range(generations))
        data.append(trajectory)
        data.append(color)
    fig = plt.plot(*data)
    plt.ylim(0, 1)
    plt.ylabel("frequency")
    plt.xlabel(xlabel)


