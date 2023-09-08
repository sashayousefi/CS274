"""
This file contains the main fragment sampling class, which performs a Monte Carlo simulated annealing procedure to fold a protein.
"""

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.scoring import *
init(extra_options='-mute all  -constant_seed')

import numpy as np
import pandas as pd
import random
import utils
from Protein import Protein
from FragmentSet import FragmentSet

import os
import random


class MCMCSampler(object):
    def __init__(self, protein, fragset, nfrag, nmer, temp, temp_end, anneal_rate):
        """
        TO DO: initialize necessary variables
        The score function is given to you (Rosetta centroid score function)
        """
        self.scorefxn = create_score_function('score3')

        #parameters loaded in from main
        self.protein = protein
        self.length = protein.length
        self.nfrag = nfrag
        self.nmer = nmer
        self.temp = temp
        self.temp_end = temp_end
        self.anneal_rate = anneal_rate
        
        #store each group of lowRMS fragment sets for each position 
        self.frag_dict = {}
        for i in np.arange(1, self.length - self.nmer + 2):
            self.frag_dict[i] = fragset.get_lowRMS_fragments(i, self.nfrag)
        
        #values to write to log
        self.lowest_energy_conf = protein
        self.iteration = 0
        self.energy = self.compute_energy(protein)


        
    def compute_energy(self, protein):
        """
        TO DO
        Compute energy of protein.
        Hint: look at utils.py
        --------
        Params:
            - protein (Protein object): protein to score
        Return:
            - energy of conformation (float)
        """        
        return utils.score_pose(protein.pose, self.scorefxn)

    def perturb_fragment(self, protein, pos, samp): # you may want to add more arguments
        """
        TO DO
        Sample from possible fragments for a position, and replace torsion angles of that fragment in the protein.
        ---------
        Params:
            - protein (Protein object): protein to score
            - pos (int): fragment position in chain (1-indexed)
            - samp (int): sample fragment index 
        Returns:
            - pert_prot (Protein object)L perturbed protein with phi/psi angles changed to fragment phi/psi angles
        """
        pert_prot = Protein(pose=protein.pose) #copy the protein object to a duplicate object
        frag = self.frag_dict[pos][samp - 1]

        for i in np.arange(self.nmer):
            angles_pos = frag[i] #get angles for each position in the fragment
            pert_prot.set_torsion(int(pos + i), float(angles_pos[0]), float(angles_pos[1])) #replace torsion angles with sampled angles
        return pert_prot



    def metropolis_accept(self, protein, pert_protein, temp): # you may want to add more arguments
        """
        TO DO
        Calculate probability of accepting or rejecting move based on Metropolis criterion.
        --------
        Params:
            - protein (Protein object): protein with pre-pertubation conformation
            - pert_protein (Protein object): protein with pre-pertubation conformation
            - temp (float): temperature of the system
        Returns:
            - (float) the probability of accepting or rejecting a move based on the Metropolis criterion
        """
        delta_E = self.compute_energy(pert_protein) - self.compute_energy(protein)
        if delta_E <= 0:
            return 1
        else:
            return np.exp(-delta_E/temp) #probability of accepting an energy increasing move 



    def anneal_temp(self):
        """
        TO DO
        Anneal temperature using exponential annealing schedule. Consider kT to be a single variable (i.e. ignore Boltzmann constant)
        --------
        Params:
            - None: use class defined temperature and anneal rate variables
        Returns:
            - annealed temperature (float): if the move is accepted, the temperature is annealed by the anneal rate
        """
        return self.temp * self.anneal_rate

    def step(self):
        """
        TO DO
        Take a single MCMC step. Each step should do the following:
        1. sample position in chain
            - Note: think about positions you can sample a k-mer fragment from. 
        2. sample fragment at that position and replace torsions in a *copied version* of the protein
        3. measure energy after replacing fragment
        4. accept or reject based on Metropolis criterion
            - if accept: incorporate proposed insertion and anneal temperature
            - if reject: sample new fragment (go to step 2)
        """
        sample_pos = random.randint(1, self.length - self.nmer + 1)
        samp_list = np.random.choice(np.arange(1, self.nfrag + 1), self.nfrag, replace = False)
        for samp in samp_list:
            pert_protein = self.perturb_fragment(self.protein, sample_pos, samp)
            energy_post = self.compute_energy(pert_protein)
            #store the lowest energy structure
            if energy_post < self.compute_energy(self.protein):
                self.lowest_energy_conf = pert_protein
            #run metropolis
            prob = self.metropolis_accept(self.protein, pert_protein, self.temp)
            unif_val = np.random.uniform(0, 1)
            if prob >= unif_val: #if the probability of accepting is greater than a uniform random value
                self.protein = pert_protein
                self.temp = self.anneal_temp()
                self.iteration += 1
                self.energy = energy_post
                break #break once we accept a move and anneal the temperature
                
            
    def simulate(self):
        """
        TO DO
        Run full MCMC simulation from start_temp to end_temp. 
        Be sure to save the best (lowest-energy) structure, so you can access it after.
        It is also a good idea to track certain variables during the simulation (temp, energy, and more).
        -------- 
        Params:
            - None: use the step function and the class defined temperature, energy, and iteration number to run/log the simulation
        Returns:
            - lowest energy conformation (Protein): the lowest energy conformation of the protein determined by the simulation 
            - temperature (list): list of temperature values at each iteration 
            - iteration (list): list of iteration numbers
            - energy (list): list of energy values at each iteration  
        """
        iteration, temperature, energy = [], [], [] #lists to hold iteration, temp, and energy values for the log files
        while self.temp >= self.temp_end: #while temperature is greater than the stop temp value
            self.step() #take an MCMC step
            iteration.append(self.iteration)
            temperature.append(self.temp) 
            energy.append(self.energy)
        return self.lowest_energy_conf, iteration, temperature, energy


    
