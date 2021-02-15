import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
import seaborn as sns
import scipy as sp
import math
from population_basicabmonly import Population

p = Population(fitness_data="manual",
               landscape_path="./test_fitness_data2.csv",
               n_sims=2,
               mut_rate=0.1,
               death_rate=0.2,
               n_gen=1000)

print(p.random_mutations(p.n_allele)) #

counts = p.run_abm_v2()
#counts, n_survive = p.simulate()

print(counts[0], counts[10], counts[20])

counts = p.run_abm_v2()
#counts, n_survive = p.simulate()
p = Population(fitness_data="manual",
               landscape_path="./test_fitness_data2.csv",
               n_sims=2,
               mut_rate=0.1,
               death_rate=0.3,
               n_gen=1000)

print(counts[0], counts[10], counts[20])
