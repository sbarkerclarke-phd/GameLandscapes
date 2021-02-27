from landscape_evolution import *
from wFish_func_fullmatrix import *
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

    pop_size = 1000         # Size of Population
    seq_length = 2       # This is "N"
    generations = 200      # How long it runs
    mutation_rate = 0.01     # per gen per individual per site
    repeats = 1            # Number of landscape replicates

    # Possible mutations at any site.
    alphabet = ['0', '1']

    #Population dictionary
    pop = {}
    
    #Genotypes
    genotypes = [''.join(seq) for seq in itertools.product("01", repeat=seq_length)]


    #Landscape (random or determined)
    A = Landscape(seq_length, 2, ls=np.array([0,0.25,0.25,0.5]))
    A.ls = ((A.ls-np.min(A.ls))/(np.max(A.ls)-np.min(A.ls))) + 1e-5   # Normalize to fitness 0 --> 1 for fitness proportionality
    
    #Fitness vector
    fitness0 = {}
    for i in range(len(genotypes)):
        fitness0[genotypes[i]] = A.ls[i]
    
    #Generate random payoff matrix
    payoff = get_payoff(fitness0, 1, genotypes)
    
    #Get pairwise game coordinates
    print(fitness0)
    game_coords = get_game_coords(payoff, genotypes, fitness0)
    
    #print(payoff)
    #print(game_coords)
    
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

    
    history=simulate_game(pop, None, generations, mutation_rate, pop_size, seq_length, A, alphabet)
    game1, history_game1=simulate_game(pop, 1, generations, mutation_rate, pop_size, seq_length, A, alphabet, True)
    game1pts = get_game_points(game1)
    print(game1)
    game2, history_game2=simulate_game(pop, 2, generations, mutation_rate, pop_size, seq_length, A, alphabet, True)
    game2pts = get_game_points(game2)
    game3, history_game3=simulate_game(pop, 3, generations, mutation_rate, pop_size, seq_length, A, alphabet, True)
    game3pts = get_game_points(game3)
    game4, history_game4=simulate_game(pop, 4, generations, mutation_rate, pop_size, seq_length, A, alphabet, True)
    game4pts = get_game_points(game4)
    
    #history_game3=simulate_game(pop, 3, generations, mutation_rate, pop_size, seq_length, A, alphabet)
    #history_game4=simulate_game(pop, 4, generations, mutation_rate, pop_size, seq_length, A, alphabet)

    #simulateSwitch(pop2, history ,generations, mutation_rate, pop_size, seq_length, fitness, fitnessB, alphabet)

    plt.rcParams['font.size'] = '10'

    #Generate plot object
    fig = mpl.figure(num=None, figsize=(14, 14), dpi=80, facecolor='w', edgecolor='k')
    h = 5 
    w = 3
    #No game
    mpl.subplot2grid((h,w), (0,0))
    stacked_trajectory_plot(history, generations, pop_size)
    plt.rcParams['font.size'] = '10'

    mpl.subplot2grid((h,w), (0,1))
    snp_trajectory_plot(history, seq_length, generations, pop_size)
    plt.rcParams['font.size'] = '10'
    
    #Fitness plot
    plt.subplot2grid((h,w),(0,2))
    plt.ylabel("Fitness")
    plt.rcParams['font.size'] = '10'
    #plt.xlim(0,1)
    
    xcoords = [0.2,0.5,0.5,0.8]
    ycoords = [0.5, 0.75, 0.25, 0.5]
    plt.xlim(-0.2,1.2)
    plt.ylim(-0.2,1.2)
    f = [i for i in fitness0.values()]
    size = [(i+0.1)*1000 for i in fitness0.values()]
     
    scatter = plt.scatter(xcoords, ycoords, c=size, s=size,  cmap="Blues", edgecolors="black")
    handles=scatter.legend_elements()[0]
    labels = ['00', "01,10",'11']
    plt.legend(handles,labels)
    plt.axis('off')        

    #Quadrant 1
    plt.subplot2grid((h,w), (1,0))
    stacked_trajectory_plot(history_game1, generations, pop_size)
    plt.rcParams['font.size'] = '10'

    plt.subplot2grid((h,w), (1,1))
    snp_trajectory_plot(history_game1, seq_length, generations, pop_size)
    #plt.legend()
    plt.rcParams['font.size'] = '10'
    
    plt.subplot2grid((h,w),(1,2))
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    plt.grid()

    for key in game1pts:
        a, b = game1pts[key]
        plt.plot(a,b, 'o', label=key)
    plt.legend(game1pts.keys(),bbox_to_anchor=(1, 1), loc='upper left', prop={'size': 12})
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, color='k')


    #Quadrant 2
    plt.subplot2grid((h,w), (2,0))
    stacked_trajectory_plot(history_game2, generations, pop_size)
    plt.rcParams['font.size'] = '10'
    
    plt.subplot2grid((h,w), (2,1))
    snp_trajectory_plot(history_game2, seq_length, generations, pop_size)
    plt.rcParams['font.size'] = '10'
    
    plt.subplot2grid((h,w),(2,2))
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    plt.grid()
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, color='k')
    for key in game2pts:
        x, y = game2pts[key]
        plt.plot(x,y, 'o')

    #Quadrant 3
    plt.subplot2grid((h,w), (3,0))
    stacked_trajectory_plot(history_game3, generations, pop_size)
    plt.rcParams['font.size'] = '10'

    plt.subplot2grid((h,w), (3,1))
    snp_trajectory_plot(history_game3, seq_length, generations, pop_size)
    plt.rcParams['font.size'] = '10'

    plt.subplot2grid((h,w),(3,2))
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    plt.grid()
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, color='k')
    for key in game3pts:
        x, y = game3pts[key]
        plt.plot(x,y, 'o')
       
    #Quadrant 4
    plt.subplot2grid((h,w), (4,0))
    stacked_trajectory_plot(history_game4, generations, pop_size)
    plt.rcParams['font.size'] = '10'
    
    plt.subplot2grid((h,w), (4,1))
    snp_trajectory_plot(history_game4, seq_length, generations, pop_size)
    plt.rcParams['font.size'] = '10'
    
    plt.subplot2grid((h,w),(4,2))
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    plt.grid()
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, color='k')
    for key in game4pts:
        x, y = game4pts[key]
        plt.plot(x,y, 'o')
    #ax1 = plt.subplot2grid((h,w), (4,0))
    #ax1.text(0.5, 0.5, "Quadrant 4", va="center", ha="center")

    #ax1.tick_params(labelbottom=False, labelleft=False)
    mpl.show()

main_func()
