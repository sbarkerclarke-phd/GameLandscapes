from landscape_evolution import *
from wFish_func import *
import matplotlib.pyplot as mpl
import matplotlib.pyplot as plt
import numpy as np
import cProfile
import pstats
import time
import random
import operator
try:
    import itertools.izip as zip
except ImportError:
    import itertools


# This is taken from Trevor Bedford's github. It is a Wright-Fisher simulation which includes: Mutation, Selection and Genetic Drift

def main_func():
    #Define things that remain same throughout simulation

    pop_size = 10000         # Size of Population
    seq_length = 4       # This is "N"
    generations = 200      # How long it runs
    mutation_rate = 0.0001     # per gen per individual per site
    repeats = 3            # Number of landscape replicates

    
    # Possible mutations at any site.
    alphabet = ['0', '1']

    #Population dictionary
    pop = {}

    #Landscape
    A = Landscape(seq_length, 2)
    print(A)
    A.ls = ((A.ls-np.min(A.ls))/(np.max(A.ls)-np.min(A.ls))) + 1e-5   # Normalize to fitness 0 --> 1 for fitness proportionality

    #genotypes = [''.join(seq) for seq in itertools.product("01", repeat=seq_length)]
    
    history_store = {}
    for i in range(repeats):
        history_store[i] = simulate_game(pop, None, generations, mutation_rate, pop_size, seq_length, A, alphabet)
        #print(max(history_store[i][99].items(), key=operator.itemgetter(1))[0])
        
    history_store1 = {}
    #for i in range(repeats):
        #history_store1[i] = simulate_game(pop, 1, generations, mutation_rate, pop_size, seq_length, A, alphabet)
        #print(max(history_store1[i][99].items(), key=operator.itemgetter(1))[0])

    #sum_list = [0]*seq_length**2   
    #for j in range(repeats):
    #    sum_list = [a + b for a, b in zip(sum_list, history_store[j])]

    #sum_list_f = [a/repeats for a in sum_list]
    #print(sum_list_f)

    history = history_store[0]
    history_game1=simulate_game(pop, 1, generations, mutation_rate, pop_size, seq_length, A, alphabet)
    history_game2=simulate_game(pop, 2, generations, mutation_rate, pop_size, seq_length, A, alphabet)
    history_game3=simulate_game(pop, 3, generations, mutation_rate, pop_size, seq_length, A, alphabet)
    history_game4=simulate_game(pop, 4, generations, mutation_rate, pop_size, seq_length, A, alphabet)

    #print(history_game1[100])
    #calculate fitness
    #avgFit = 0
    #for i in pop.keys():
        #avgFit += (pop[i]/pop_size)*fitness[i]


    #simulateSwitch(pop2, history ,generations, mutation_rate, pop_size, seq_length, fitness, fitnessB, alphabet)

    #Generate plot object
    fig = mpl.figure(num=None, figsize=(14, 14), dpi=80, facecolor='w', edgecolor='k')
    h = 5 
    w = 2
    #No game
    mpl.subplot2grid((h,w), (0,0))
    stacked_trajectory_plot(history, generations, pop_size)

    mpl.subplot2grid((h,w), (0,1))
    snp_trajectory_plot(history, seq_length, generations, pop_size)

    #Quadrant 1
    plt.subplot2grid((h,w), (1,0))
    stacked_trajectory_plot(history_game1, generations, pop_size)
    
    plt.subplot2grid((h,w), (1,1))
    snp_trajectory_plot(history_game1, seq_length, generations, pop_size)

    #Quadrant 2
    plt.subplot2grid((h,w), (2,0))
    stacked_trajectory_plot(history_game2, generations, pop_size)
    
    plt.subplot2grid((h,w), (2,1))
    snp_trajectory_plot(history_game2, seq_length, generations, pop_size)

    #Quadrant 3
    plt.subplot2grid((h,w), (3,0))
    stacked_trajectory_plot(history_game3, generations, pop_size)
    
    plt.subplot2grid((h,w), (3,1))
    snp_trajectory_plot(history_game3, seq_length, generations, pop_size)

    #Quadrant 4
    plt.subplot2grid((h,w), (4,0))
    stacked_trajectory_plot(history_game4, generations, pop_size)
    
    plt.subplot2grid((h,w), (4,1))
    snp_trajectory_plot(history_game4, seq_length, generations, pop_size)

    #ax1 = plt.subplot2grid((h,w), (4,0))
    #ax1.text(0.5, 0.5, "Quadrant 4", va="center", ha="center")

    #ax1.tick_params(labelbottom=False, labelleft=False)
    mpl.show()

main_func()
