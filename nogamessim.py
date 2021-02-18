import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
import seaborn as sns
import scipy as sp
import math
from population_basicabmonly import Population
from pymuller import muller



p = Population(fitness_data="manual",
               landscape_path="./test_fitness_data2.csv",
               n_sims=2,
               mut_rate=0.01,
               death_rate=0.05,
               init_counts=[0,0,200,0],
               n_gen=300)

print(p.random_mutations(p.n_allele)) #

#counts = p.run_abm_v2()

#print (counts[200])

#counts, n_survive = p.simulate()
#np.savetxt("nogame.csv",counts, delimiter=",")
print("NoGame")
p = Population(fitness_data="manual",
               landscape_path="./test_fitness_data2.csv",
               n_sims=2,
               mut_rate=0.01,
               death_rate=0.05,
               init_counts=[200,200,200,0],
               n_gen=300)

print(p.random_mutations(p.n_allele)) #
counts, n_survive = p.simulate()
#np.savetxt("nogame.csv",counts, delimiter=",")

print("Game1")

p = Population(fitness_data="game",
               landscape_path="./test_fitness_data2.csv",
               n_sims=2,
               mut_rate=0.01,
               death_rate=0.05,
               init_counts=[200,200,200,0],
               n_gen=300,
               game_quadrant=1)

print(p.random_mutations(p.n_allele)) #
counts, n_survive = p.simulate()
#np.savetxt("game1.csv",counts, delimiter=",")


print("Game2")

p = Population(fitness_data="game",
               landscape_path="./test_fitness_data2.csv",
               n_sims=2,
               mut_rate=0.01,
               death_rate=0.05,
               init_counts=[200,200,200,0],
               n_gen=300,
               game_quadrant=2)

print(p.random_mutations(p.n_allele)) #
counts, n_survive = p.simulate()
#np.savetxt("game2.csv",counts, delimiter=",")

print("Game3")

p = Population(fitness_data="game",
               landscape_path="./test_fitness_data2.csv",
               n_sims=2,
               mut_rate=0.01,
               death_rate=0.05,
               init_counts=[200,200,200,0],
               n_gen=300,
               game_quadrant=3)

print(p.random_mutations(p.n_allele)) #
counts, n_survive = p.simulate()
#np.savetxt("game3.csv",counts, delimiter=",")

print("Game4")

p = Population(fitness_data="game",
               landscape_path="./test_fitness_data2.csv",
               n_sims=2,
               mut_rate=0.01,
               death_rate=0.05,
               init_counts=[200,200,200,0],
               n_gen=300,
               game_quadrant=4)

print(p.random_mutations(p.n_allele)) #
counts, n_survive = p.simulate()
#np.savetxt("game4.csv",counts, delimiter=",")


