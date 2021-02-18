import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
import seaborn as sns
import scipy as sp
import math
# import warnings

class Population:
###############################################################################    
    # Initializer
    def __init__(self,
                 carrying_cap = False,
                 curve_type='constant', # drug concentration curve
                 counts_log_scale = False, # plot counts on log scale
                 constant_pop = False, # normalize to a constant population size
                 drugless_path = None, # file path for the drugless growth rates
                 death_rate = 0.15, 
                 div_scale = 1, # scale the division rate to model different organisms
                 drug_log_scale = False, # plot the drug concentration curve on a log scale
                 drug_curve = None, # input a custom drug concentration curve
                 drug_regimen = None, # for modeling drug regimens
                 debug=False, # print the current time step
                 entropy_lim = None, # entropy plotting limits
                 fig_title = '',
                 fitness_data = 'generate', # 'generate' = generate fitness data using drugless growth rates, ic50, and drug concentration. 'manual' = input fitness landscape from csv
                 h_step = 500,
                 ic50_path = None, 
                 init_counts = None, # default is 10,000 wild type cells
                 k_elim = 0.001, # form modeling pharmacokinetics
                 k_abs = 0.07,
                 landscape_path = None, # path for custom fitness landscape
                 min_dose = 0,
                 mut_rate = 0.01, # mutation rate
                 max_cells = 10**6, # carrying capacity
                 max_dose = 1, 
                 n_gen=1000, # number of generations
                 n_impulse = 2, # for modeling dosage regimens
                 normalize = False, 
                 n_sims = 1, # number of simulations to average together
                 pad_right = False,
                 plot=True, # plot the result of simulate()
                 plot_entropy = False, # plot the entropy of the population over time underneath the timecourse
                 prob_drop=0, 
                 slope = None, 
                 thresh = 5000,
                 three_undefined = False, # this is a hack because allele #3 in the pyrimethamine data is not quantified
                 x_lim = None, # plotting
                 y_lim = None # plotting
                 ):
                
        # Evolutionary parameters
        
        # Number of generations (time steps)
        self.n_gen = n_gen
        self.max_cells = max_cells
        
        # ABM parameters
        self.mut_rate = mut_rate
        self.death_rate = death_rate
        
        # Timecourse (set after running self.simulate)
        self.counts = np.zeros([self.n_gen,16])
        self.counts_extinct = np.zeros([self.n_gen,16])
        self.counts_survive = np.zeros([self.n_gen,16])
                                            
        # Model parameters
        
        self.carrying_cap = True
        self.thresh = thresh # Threshold for hybrid model (not yet implemented)
        self.div_scale = div_scale # Scale the division rate to simulate different organisms
        self.n_sims = n_sims # number of simulations to average together in self.simulate
        self.constant_pop = constant_pop
        # self.v2 = v2
        self.debug = debug
        
        self.fitness_data = fitness_data
        
        # Generate fitness data from IC50 and drugless growth rate data
        if fitness_data == 'generate':
            # Data paths
            if drugless_path is None:
                # self.drugless_path = "C:\\Users\\Eshan\\Documents\\python scripts\\theory division\\abm_variable_fitness\\data\\ogbunugafor_drugless.csv"
                self.drugless_path = 'ogbunugafor_drugless.csv'
            else:
                self.drugless_path = drugless_path
                
            if ic50_path is None:
                # self.ic50_path = "C:\\Users\\Eshan\\Documents\\python scripts\\theory division\\abm_variable_fitness\\data\\pyrimethamine_ic50.csv"
                self.ic50_path = 'pyrimethamine_ic50.csv'
            else:
                self.ic50_path = ic50_path
            
            # load the data
            self.drugless_rates = self.load_fitness(self.drugless_path)
            self.ic50 = self.load_fitness(self.ic50_path)
            
            self.timestep_scale = max(self.drugless_rates)
            # determine number of alleles from data (not yet implemented)
            self.n_allele = self.drugless_rates.shape[0]
        
        # load fitness landscape from excel file
        elif fitness_data == 'manual':
            self.landscape_path = landscape_path
            self.landscape_data = self.load_fitness(self.landscape_path)
            
            self.timestep_scale = max(self.landscape_data)
            self.n_allele = self.landscape_data.shape[0]
            
        # Initial number of cells (default = 100 at 0000)
        if init_counts is None:
            self.init_counts = np.zeros(self.n_allele)
            self.init_counts[0] = 10**2
        else:
            self.init_counts = init_counts
        
        # Dose parameters
        self.curve_type = curve_type # linear, constant, heaviside, pharm, pulsed
        
        # Pharmacological paramters
        if k_abs < k_elim:
            raise Exception('Inappropriate pharmacokinetic values: k_abs < k_elim.')            
        
        self.k_elim = k_elim
        self.k_abs = k_abs 
        self.pad_right = pad_right
        self.max_dose = max_dose
        
        if slope is None:
            self.slope = self.max_dose/self.n_gen # Ramped parameter (what time step to reach maximum dose, determines slope)
        else:
            self.slope = slope
        

        self.n_impulse = n_impulse # number of impulses for a pulsed dose
        self.prob_drop = prob_drop # probability of dropping a dose
        self.h_step = h_step # when to turn on heaviside function
        self.min_dose = min_dose 
        
        # Generate drug dosage curves if one is not specified
        if drug_curve is None:
            self.drug_curve,u = self.gen_curves()
        else:
            self.drug_curve = drug_curve
            
        if drug_regimen is None:
            self.drug_regimen = u
        else:
            self.drug_regimen = drug_regimen
        
        # Visualization parameters
        self.plot = plot # boolean
        self.plot_entropy = plot_entropy
        self.drug_log_scale = drug_log_scale # plot drugs on log scale
        self.counts_log_scale = counts_log_scale # plot counts on log scale
        self.fig_title = fig_title
        self.normalize = normalize
        self.counts_log_scale = counts_log_scale
        if x_lim is None:
            self.x_lim = n_gen
        else:
           self.x_lim = x_lim 
        self.y_lim = y_lim
        self.entropy_lim = entropy_lim
###############################################################################       
    
    # Load data
    def load_fitness(self,data_path):
        # also use to load ic50 and drugless growth rate
        fitness = pd.read_csv(data_path)
        cols = list(fitness.columns)
        fit_array = np.array(cols)
        fit_array = fit_array.astype(np.float)
        return fit_array
        
###############################################################################
    # ABM helper methods
    
    # converts decimals to binary
    def int_to_binary(self,num):
        pad = int(math.log(self.n_allele,2))
        return bin(num)[2:].zfill(pad)
    
    # computes hamming distance between two genotypes
    def hammingDistance(self,s1,s2):
        assert len(s1) == len(s2)
        return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
    
    # converts an integer to a genotype and padding to the left by 0s
    def convertIntToGenotype(self,anInt,pad):
    	offset = 2**pad
    	return [int(x) for x in bin(offset+anInt)[3:]]
    
    def random_mutations(self,N):
        trans_mat = np.zeros([N,N])
        for mm in range(N):
            for nn in range(N):
                trans_mat[mm, nn] = self.hammingDistance( self.int_to_binary(mm) , self.int_to_binary(nn))
        # Don't include mutant no. 4
        
        trans_mat[trans_mat>1] = 0
        trans_mat = trans_mat/trans_mat.sum(axis=1)
        
    #    trans_mat[3,:] = 0
    #    trans_mat[:,3] = 0
    #    print(str(trans_mat))
        return trans_mat
    
    # compute fitness given a drug concentration
    def gen_fitness(self,allele,conc,drugless_rate,ic50):        
        c = -.6824968 # empirical curve fit

        # logistic equation from Ogbunugafor 2016
        conc = conc/10**6 # concentration in uM, convert to M
        
        # ic50 is already log-ed in the dataset
        log_eqn = lambda d,i: d/(1+np.exp((i-np.log10(conc))/c))
        if conc <= 0:
            fitness = drugless_rate[allele]
        else:
            fitness = log_eqn(drugless_rate[allele],ic50[allele])

        return fitness
    
    def gen_fit_land(self,conc):
        
        fit_land = np.zeros(self.n_allele)
        
        for allele in range(self.n_allele):
            fit_land[allele] = self.gen_fitness(allele,conc,self.drugless_rates,self.ic50)
        
        return fit_land
    
###############################################################################
    # Methods for generating drug curves
    
    # Equation for a simple 1 compartment pharmacokinetic model
    def pharm_eqn(self,t,k_elim=0.01,k_abs=0.1,max_dose=1):
        conc = np.exp(-k_elim*t)-np.exp(-k_abs*t)
        t_max = np.log(k_elim/k_abs)/(k_elim-k_abs)
        conc = conc/(np.exp(-k_elim*t_max)-np.exp(-k_abs*t_max))
        conc = conc*max_dose
        return conc
    
    # Old convolution algorithm (very slow)
    # Convolve the arbitrary curve u with the pharmacokinetic model
    # def convolve_pharm(self,u):
    #                    # k_elim=0.01,
    #                    # k_abs=0.1,
    #                    # max_dose=1):
    #     k_elim = self.k_elim
    #     k_abs = self.k_abs
    #     max_dose = self.max_dose
        
    #     # k_lim and k_abs are the absorption and elimination rate constants for the pharmacokinetic model
    #     # t_max is the max length of the output curve
    #     # algorithm is at best O(n^2)...
        
    #     conv = np.zeros(self.n_gen)
    #     for t in range(self.n_gen):
    #         for tau in range(self.n_gen):
    #             if t-tau >= 0 and t-tau<u.shape[0]:
    #                 conv[t] += u[t-tau]*self.pharm_eqn(tau,k_elim=k_elim,k_abs=k_abs,max_dose=max_dose)
    #     # conv = np.convolve()
    #     return conv

    # New convolution method (much faster with numpy)
    def convolve_pharm(self,u):
                       # k_elim=0.01,
                       # k_abs=0.1,
                       # max_dose=1):
        k_elim = self.k_elim
        k_abs = self.k_abs
        max_dose = self.max_dose
        
        # k_lim and k_abs are the absorption and elimination rate constants for the pharmacokinetic model
        # t_max is the max length of the output curve
        # algorithm is at best O(n^2)...
        
        # conv = np.zeros(self.n_gen)
        # for t in range(self.n_gen):
        #     for tau in range(self.n_gen):
        #         if t-tau >= 0 and t-tau<u.shape[0]:
        #             conv[t] += u[t-tau]*self.pharm_eqn(tau,k_elim=k_elim,k_abs=k_abs,max_dose=max_dose)
        
        pharm = np.zeros(self.n_gen)
        for i in range(self.n_gen):
            pharm[i] = self.pharm_eqn(i,k_elim=k_elim,k_abs=k_abs,max_dose=max_dose)
        
        # using FFT turns out to be much faster for convolutions!
        conv = np.convolve(u,pharm)
        conv = conv[0:self.n_gen]
        return conv
    
    # Generates an impulse train to input to convolve_pharm()
    def gen_impulses(self):
        gap = np.floor(self.n_gen/self.n_impulse)
        u = np.zeros(self.n_gen)
        if self.pad_right:
            impulse_indx = np.arange(self.n_impulse)*gap
        else:
            impulse_indx = np.arange(self.n_impulse+1)*gap-1
            
        # eliminate random doses
        keep_indx = np.random.rand(len(impulse_indx)) > self.prob_drop
        
        impulse_indx = impulse_indx[keep_indx]
        
        impulse_indx = impulse_indx.astype(int)
        u[impulse_indx]=1 
        return u
    
    # generates drug concentration curves
    def gen_curves(self):
        curve = np.zeros(self.n_gen)
        # print('hi')
        u = None
        if self.curve_type == 'linear': # aka ramp linearly till timestep defined by steepness
            # cur_dose = 0
            for i in range(self.n_gen):
                # if i <= self.steepness:
                #     slope = (self.max_dose-10**(-3))/self.steepness
                #     conc = slope*i+10**-3
                # else:
                #     # step = self.steepness
                #     slope = (self.max_dose-10**(-3))/self.steepness
                # if cur_dose < self.max_dose:
                conc = self.slope*i
                
                if conc > self.max_dose:
                    conc=self.max_dose
                # else:
                #     conc = self.max_dose
                # cur_dose = conc
                    # conc = slope*i+10**-3
                curve[i]=conc
                
        elif self.curve_type == 'constant':
            curve[:] = self.max_dose
            # print('here')

        elif self.curve_type == 'heaviside':
            for i in range(self.n_gen):
                if i <= self.h_step:
                    curve[i] = self.min_dose
                else:
                    curve[i] = self.max_dose 
        
        # Two compartment pharmacokinetic model
        elif self.curve_type == 'pharm':
            t_max = np.log(self.k_elim/self.k_abs)/(self.k_elim-self.k_abs)
            for i in range(self.n_gen):
                conc = np.exp(-self.k_elim*i)-np.exp(-self.k_abs*i)
                conc = conc/(np.exp(-self.k_elim*t_max)-np.exp(-self.k_abs*t_max))
                conc = conc*self.max_dose
                curve[i] = conc
        
        # Pulsed convolves an impulse train with the 1-compartment model (models patient taking a maintenence dose)
        elif self.curve_type == 'pulsed':
            u = self.gen_impulses()
            curve = self.convolve_pharm(u)
        return curve, u

###############################################################################
    # Run one abm simulation (ignores n_sim)
    # def run_abm(self):
        
    #     n_allele = len(self.drugless_rates)

    #     # Obtain transition matrix for mutations
    #     P = self.random_mutations( n_allele )

    #     # Keeps track of cell counts at each generation
    #     counts = np.zeros([self.n_gen, n_allele], dtype=int)
    #     # drug_curve = np.zeros(self.n_gen)

    #     counts[0,:] = self.init_counts
    
    #     for mm in range(self.n_gen-1):
    #         # Normalize to constant population
                            
    #         conc = self.drug_curve[mm]
            
    #         fit_land = np.zeros(self.n_allele)
            
    #         # fitness of allele 0010 is not quantified in the dataset - set to zero
    #         for kk in range(self.n_allele):
    #             if kk == 3:
    #                 fit_land[kk] = 0 # fitness of allele 0010 is not quantified in the dataset
    #             elif kk < 3:
    #                 fit_land[kk] = self.gen_fitness(kk,conc,self.drugless_rates,self.ic50)
    #             elif kk > 3:
    #                 fit_land[kk] = self.gen_fitness(kk,conc,self.drugless_rates,self.ic50)
            
    #         fit_land = fit_land*self.div_scale
    #         n_cells = np.sum( counts[mm] )
    
    #         # Scale division rates based on carrying capacity
    #         if self.carrying_cap:
    #             division_scale = 1 / (1+(2*np.sum(counts[mm])/self.max_cells)**4)
    #         else:
    #             division_scale = 1
    
    #         if counts[mm].sum()>self.max_cells:
    #             division_scale = 0
    
    #         div_rate = np.repeat( fit_land*division_scale, counts[mm] )
    #         cell_types = np.repeat( np.arange(n_allele) , counts[mm] )
    
    #         # Death of cells
    #         death_rates = np.random.rand(n_cells)
    #         surv_ind = death_rates > self.death_rate
    #         div_rate = div_rate[surv_ind]
    #         cell_types = cell_types[surv_ind]
    #         n_cells = len(cell_types)

    #         counts[mm+1] = np.bincount(cell_types, minlength=n_allele)
    
    #         #Divide and mutate cells
    #         div_ind = np.random.rand(n_cells) < div_rate
    
    #         # Mutate cells
    #         # initial state of allele types
    #         daughter_types = cell_types[div_ind].copy()
    
    #         # Generate random numbers to check for mutation
    #         daughter_counts = np.bincount( daughter_types , minlength=n_allele)
    
    #         # Mutate cells of each allele type
    #         for allele in np.random.permutation(np.arange(n_allele)):
    #             n_mut = np.sum( np.random.rand( daughter_counts[allele] ) < self.mut_rate )
    
    #             # note that columns in P are normalized to probability densities (columns sum to 1)
    #             mutations = np.random.choice(n_allele, size=n_mut, p=P[:,allele]).astype(np.uint8)
    
    #             #Add mutating cell to their final types
    #             counts[mm+1] +=np.bincount( mutations , minlength=n_allele)
    #             counts[:,3] =  0
                
    #             #Substract mutating cells from that allele
    #             daughter_counts[allele] -=n_mut
    
    #         counts[mm+1] += daughter_counts
            
    #         # Normalize to constant population            
    #         if self.constant_pop:
    #             cur_size = np.sum(counts[mm+1])
    #             counts[mm+1] = counts[mm+1]*self.init_counts[0]/cur_size
    #             counts[mm+1] = np.floor(counts[mm+1])

    #     return counts
    
    def run_abm_v2(self):
        
        n_allele = self.n_allele

        # Obtain transition matrix for mutations
        P = self.random_mutations( n_allele )

        # Keeps track of cell counts at each generation
        counts = np.zeros([self.n_gen, n_allele], dtype=int)
        # drug_curve = np.zeros(self.n_gen)

        counts[0,:] = self.init_counts
    
        for mm in range(self.n_gen-1):
            
            if self.debug:
                if np.mod(mm,10) == 0:
                    print(str(mm))
                            
            conc = self.drug_curve[mm]
            
            fit_land = np.zeros(self.n_allele)
            
            if self.fitness_data == 'generate':
            
                # fitness of allele 0010 is not quantified in the pyrimethamine dataset - set to zero
                if self.three_undefined == True:
                    for kk in range(self.n_allele):
                        if kk == 3:
                            fit_land[kk] = 0 # fitness of allele 0010 is not quantified in the dataset
                        elif kk < 3:
                            fit_land[kk] = self.gen_fitness(kk,conc,self.drugless_rates,self.ic50)
                        elif kk > 3:
                            fit_land[kk] = self.gen_fitness(kk,conc,self.drugless_rates,self.ic50)
                
                else:
                    for kk in range(self.n_allele):
                        if kk == 3:
                            fit_land[kk] = 0 # fitness of allele 0010 is not quantified in the dataset
                        elif kk < 3:
                            fit_land[kk] = self.gen_fitness(kk,conc,self.drugless_rates,self.ic50)
                        elif kk > 3:
                            fit_land[kk] = self.gen_fitness(kk,conc,self.drugless_rates,self.ic50)
                            
            elif self.fitness_data == 'manual':
                fit_land = self.landscape_data
                    
            fit_land = fit_land*self.div_scale
    
            # Scale division rates based on carrying capacity
            if self.carrying_cap:
                division_scale = 1 / (1+(2*np.sum(counts[mm])/self.max_cells)**4)
            else:
                division_scale = 1
    
            if counts[mm].sum()>self.max_cells:
                division_scale = 0
            
            fit_land = fit_land*division_scale
            
            # fit_land = fit_land/self.timestep_scale # ensures fitness is never greater than 1
            # death_rate = self.death_rate/self.timestep_scale # scales all other rates proportionally
            # mut_rate = self.mut_rate/self.timestep_scale
            
            death_rate = self.death_rate
            mut_rate = self.mut_rate
            
            counts[mm+1] = counts[mm]
    
            # Kill cells
            
            counts[mm+1] = counts[mm+1]-np.random.binomial(counts[mm],death_rate)
    
            # Divide cells
            
            divide = np.random.binomial((counts[mm+1]*self.timestep_scale).astype(int),fit_land/self.timestep_scale)
            
            # Mutate cells
            
            daughter_types = np.repeat( np.arange(n_allele) , divide )
            daughter_counts = np.bincount( daughter_types , minlength=n_allele )
    
            # Mutate cells of each allele type
            for allele in np.random.permutation(np.arange(n_allele)):
                n_mut = np.random.binomial(daughter_counts[allele],mut_rate)
    
                mutations = np.random.choice(n_allele, size=n_mut, p=P[:,allele]).astype(np.uint8)
    
                # Add mutating cell to their final types
                counts[mm+1] +=np.bincount( mutations , minlength=n_allele)
                counts[:,3] =  0
                # Substract mutating cells from that allele
                daughter_counts[allele] -=n_mut
    
            counts[mm+1] += daughter_counts
            
            # Normalize to constant population            
            if self.constant_pop:
                cur_size = np.sum(counts[mm+1])
                counts[mm+1] = counts[mm+1]*self.init_counts[0]/cur_size
                counts[mm+1] = np.floor(counts[mm+1])

        return counts

    # Runs abm simulation n_sim times and averages results. Then sets self.counts to final result. Also quantifies survival number
    def simulate(self):
        
        counts_t = np.zeros([self.n_gen,self.n_allele])
        counts = np.zeros([self.n_gen,self.n_allele])
        counts_survive = np.zeros([self.n_gen,self.n_allele])
        counts_extinct = np.zeros([self.n_gen,self.n_allele])
        
        n_survive = 0
        for i in range(self.n_sims):
            
            # if self.v2:
            #     counts_t = self.run_abm_v2()
            # else:
            #     counts_t = self.run_abm()
            
            counts_t = self.run_abm_v2()
                
            if any(counts_t[self.n_gen-1,:]>0.1*self.max_cells):
                n_survive+=1
                counts_survive += counts_t
                if self.plot is True:
                    title_t = 'Dose = ' + str(self.max_dose) + ' uM, survived'
                    self.plot_timecourse(counts_t = counts_t,
                                         title_t = title_t)
            else:
                counts_extinct += counts_t
                if self.plot is True:
                    title_t = 'Dose = ' + str(self.max_dose) + ' , extinct'
                    self.plot_timecourse(counts_t = counts_t,
                                         title_t = title_t)           
           
                
            counts+=counts_t
                                                                                                      
        counts = counts/self.n_sims
        counts_survive = counts_survive/n_survive
        if  (self.n_sims - n_survive) > 0:
            counts_extinct = counts_extinct/(self.n_sims-n_survive)
        
        self.counts = counts
        self.counts_survive = counts_survive
        self.counts_extinct = counts_extinct
        
        return counts, n_survive

    # def plot_timecourse(self,counts_t=None,title_t=None):
        
    #     if (self.counts == 0).all() and counts_t is None:
    #         print('No data to plot!')
    #         return
    #     elif counts_t is None:
    #         counts = self.counts
    #     else:
    #         counts = counts_t # an input other than self overrides self
    #     if title_t is not None:
    #         title = title_t
    #     else:
    #         title = self.fig_title    
            
    #     fig, ax = plt.subplots(figsize = (6,4))
    # #    plt.rcParams.update({'font.size': 22})
    #     counts_total = np.sum(counts,axis=0)
        
    #     sorted_index = counts_total.argsort()
    #     sorted_index_big = sorted_index[8:]
        
    #     colors = sns.color_palette('bright')
    #     colors = np.concatenate((colors[0:9],colors[0:7]),axis=0)
    #     # shuffle colors

    #     colors[[14,15]] = colors[[15,14]]
        
    #     cc = (cycler(color=colors) + 
    #           cycler(linestyle=['-', '-','-','-','-','-','-','-','-',
    #                             '--','--','--','--','--','--','--']))
        
    #     ax.set_prop_cycle(cc)
        
    #     ax2 = ax.twinx()

    #     color = [0.5,0.5,0.5]
    #     ax2.set_ylabel('Drug Concentration (uM)', color=color,fontsize=20) # we already handled the x-label with ax1
        
    #     if self.drug_log_scale:
    #         if all(self.drug_curve>0):
    #             drug_curve = np.log10(self.drug_curve)
    #         yticks = np.log10([10**-4,10**-3,10**-2,10**-1,10**0,10**1,10**2,10**3])    
    #         ax2.set_yticks(yticks)
    #         ax2.set_yticklabels(['0','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$',
    #                          '$10^1$','$10^2$','$10^3$'])
    #         ax2.set_ylim(-4,3)
    #     else:
    #         drug_curve = self.drug_curve
    #         ax2.set_ylim(0,1.1*max(drug_curve))
    
    # #    ax2.plot(drug_curve, color=color, linewidth=3.0, linestyle = 'dashed')
    #     ax2.plot(drug_curve, color=color, linewidth=2.0)
    #     ax2.tick_params(axis='y', labelcolor=color)
            
    #     ax2.legend(['Drug Conc.'],loc=(1.25,0.85),frameon=False,fontsize=15)
        
    #     ax2.tick_params(labelsize=18)
    # #    plt.yticks(fontsize=18)
    #     ax2.set_title(title,fontsize=20)
        
    #     if self.normalize:
    #         counts = counts/np.max(counts)
            
    #     for allele in range(counts.shape[1]):
    #         if allele in sorted_index_big:
    #             ax.plot(counts[:,allele],linewidth=3.0,label=str(self.int_to_binary(allele)))
    # #            print(str(allele))
    #         else:
    #             ax.plot(counts[:,allele],linewidth=3.0,label=None)
    #     ax.legend(loc=(1.25,.03),frameon=False,fontsize=15)
    # #    ax.legend(frameon=False,fontsize=15)
    # #        ax.legend([str(int_to_binary(allele))])
            
    #     ax.set_xlim(0,self.x_lim)
    #     ax.set_facecolor(color='w')
    #     ax.grid(False)
    
    #     ax.set_xlabel('Time',fontsize=20)
    #     ax.set_ylabel('Cells',fontsize=20)
    #     ax.tick_params(labelsize=20)
        
    #     if self.y_lim is not None:
    #         y_lim = self.y_lim
    #     else:
    #         y_lim = np.max(counts)
        
    #     if self.counts_log_scale:
    #         ax.set_yscale('log')
    #         ax.set_ylim(1,5*10**5)
    #     else:
    #         ax.set_ylim(0,y_lim)
    
    #     plt.show()
    #     return fig,ax
    
    def plot_timecourse(self,counts_t=None,title_t=None):
        
        if (self.counts == 0).all() and counts_t is None:
            print('No data to plot!')
            return
        elif counts_t is None:
            counts = self.counts
        else:
            counts = counts_t # an input other than self overrides self
        if title_t is not None:
            title = title_t
        else:
            title = self.fig_title    
            
        # fig, ax = plt.subplots(2,1,
        #                        figsize = (6,4),
        #                        sharex=True,)
    #    plt.rcParams.update({'font.size': 22})
        left = 0.1
        width = 0.8
        
        if self.plot_entropy == True:
            fig,(ax1,ax3) = plt.subplots(2,1,figsize=(6,4),sharex=True) 
            ax3.set_position([left, 0.2, width, 0.2]) # ax3 is entropy
        else:
            fig,ax1 = plt.subplots(1,1,figsize=(6,4),sharex=True)
        
        ax1.set_position([left, 0.5, width, 0.6]) # ax1 is the timecourse
                
        counts_total = np.sum(counts,axis=0)
        
        sorted_index = counts_total.argsort()
        sorted_index_big = sorted_index[-8:]
        
        colors = sns.color_palette('bright')
        colors = np.concatenate((colors[0:9],colors[0:7]),axis=0)
        # shuffle colors

        colors[[14,15]] = colors[[15,14]]
        
        cc = (cycler(color=colors) + 
              cycler(linestyle=['-', '-','-','-','-','-','-','-','-',
                                '--','--','--','--','--','--','--']))
        
        ax1.set_prop_cycle(cc)

        color = [0.5,0.5,0.5]
        
        if self.fitness_data == 'generate':
            ax2 = ax1.twinx() # ax2 is the drug timecourse
            ax2.set_position([left, 0.5, width, 0.6])
            ax2.set_ylabel('Drug Concentration (uM)', color=color,fontsize=20) # we already handled the x-label with ax1
            
            if self.drug_log_scale:
                if all(self.drug_curve>0):
                    drug_curve = np.log10(self.drug_curve)
                yticks = np.log10([10**-4,10**-3,10**-2,10**-1,10**0,10**1,10**2,10**3])    
                ax2.set_yticks(yticks)
                ax2.set_yticklabels(['0','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$',
                                 '$10^1$','$10^2$','$10^3$'])
                ax2.set_ylim(-4,3)
            else:
                drug_curve = self.drug_curve
                ax2.set_ylim(0,1.1*max(drug_curve))
        
        #    ax2.plot(drug_curve, color=color, linewidth=3.0, linestyle = 'dashed')
            ax2.plot(drug_curve, color=color, linewidth=2.0)
            ax2.tick_params(axis='y', labelcolor=color)
                
            ax2.legend(['Drug Conc.'],loc=(1.25,0.93),frameon=False,fontsize=15)
            
            ax2.tick_params(labelsize=15)
        #    plt.yticks(fontsize=18)
            ax2.set_title(title,fontsize=20)
        
        if self.normalize:
            counts = counts/np.max(counts)
            
        # print(str(sorted_index_big))
        for allele in range(counts.shape[1]):
            if allele in sorted_index_big:
                ax1.plot(counts[:,allele],linewidth=3.0,label=str(self.int_to_binary(allele)))
    #            print(str(allele))
            else:
                ax1.plot(counts[:,allele],linewidth=3.0,label=None)
                
        ax1.legend(loc=(1.25,-.12),frameon=False,fontsize=15)
    #    ax.legend(frameon=False,fontsize=15)
    #        ax.legend([str(int_to_binary(allele))])
            
        ax1.set_xlim(0,self.x_lim)
        ax1.set_facecolor(color='w')
        ax1.grid(False)
    
        # ax1.set_xlabel('Time',fontsize=20)
        ax1.set_ylabel('Cells',fontsize=20)
        ax1.tick_params(labelsize=15)
        
        if self.plot_entropy == True:
            e = self.entropy(counts)
            
            ax3.plot(e,color='black')
            ax3.set_xlabel('Time',fontsize=20)
            ax3.set_ylabel('Entropy',fontsize=20)
            if self.entropy_lim is not None:
                ax3.set_ylim(0,self.entropy_lim)
            ax3.tick_params(labelsize=15)
        
        if self.y_lim is not None:
            y_lim = self.y_lim
        else:
            y_lim = np.max(counts)
        
        if self.counts_log_scale:
            ax1.set_yscale('log')
            ax1.set_ylim(1,5*10**5)
        else:
            ax1.set_ylim(0,y_lim)
        
        plt.show()
        return fig
    
    # Calculate the shannon-gibbs entropy (normalized population size)
    def entropy(self,counts=None):
        if counts is None:
            counts = self.counts

        k = np.sum(counts,1)
        entropy = np.zeros(counts.shape[0])
        counts_t = np.zeros(counts.shape)
        for i in range(counts.shape[0]):
            counts_t[i,:] = np.divide(counts[i,:],k[i])
            # counts_log = np.log(counts_t[i,:]+1)
            # entropy[i] = np.dot(counts_t[i,:],counts_log)
            entropy[i] = sp.stats.entropy(counts_t[i,:])

        return entropy
    
    def plot_fitness_curves(self,fig_title=''):
    
        drugless_rates = self.drugless_rates
        ic50 = self.ic50
        
        fig, ax = plt.subplots(figsize = (12,6))
        
        powers = np.linspace(-3,5,20)
        conc = np.power(10*np.ones(powers.shape[0]),powers)
        fit = np.zeros(conc.shape[0])
        
        # colors = sns.color_palette('bright')
        # colors = np.concatenate((colors[0:9],colors[0:7]),axis=0)
        # colors[[14,15]] = colors[[15,14]]
        
        # cc = (cycler(color=colors) + 
              # cycler(linestyle=['-', '-','-','-','-','-','-','-','-',
                                # '--','--','--','--','--','--','--']))
        # ax.set_prop_cycle(cc) 
        
        for allele in range(16):
            if allele == 3:
                fit = np.zeros(conc.shape[0])
            if allele > 3:
                for j in range(conc.shape[0]):
                    fit[j] = self.gen_fitness(allele,conc[j],drugless_rates,ic50)
            else:
                for j in range(conc.shape[0]):
                    fit[j] = self.gen_fitness(allele,conc[j],drugless_rates,ic50)
            ax.plot(powers,fit,linewidth=3,label=str(self.int_to_binary(allele)))
    #    ind = np.arange(9)
        ax.legend(fontsize=15,frameon=False,loc=(1,-.10))
        ax.set_xticks([-3,-2,-1,0,1,2,3,4,5])
        ax.set_xticklabels(['$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$',
                             '$10^1$','$10^2$','$10^3$','$10^4$','$10^5$'])
        
        plt.title(fig_title,fontsize=20)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        
        plt.xlabel('Drug concentration ($\mathrm{\mu}$M)',fontsize=20)
        plt.ylabel('Growth Rate',fontsize=20)
        ax.set_frame_on(False)
        
        return ax
