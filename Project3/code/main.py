"""
This is the master file, which you should use to set up and run the simulations.
You may define functions or classes as necessary

For an input sequence, do the following (see project page for details):
	1. load sequence from fasta file and initialize protein into extended configuration
	2. Run a series of simulations (n=5 or 10):
		- Perform MCMC sampling with 9-mer fragments from kT=100 to kT=1 (assembly stage)
		- Perform MCMC sampling with 3-mer fragments from kT=1 to kT=0.1 (refinement stage)
		- Take best (lowest-energy) structure after refinement and perform energy minimization (see utils.relax)
		- Log energy and RMSD to native structure after minimization
	3. Visualize lowest-RMSD structure in PyMol

"""
from pyrosetta import *
from pyrosetta.rosetta import *
init(extra_options='-mute all -constant_seed')

from Protein import Protein
from FragmentSampler import MCMCSampler
from FragmentSet import FragmentSet

import utils
import argparse
import numpy as np
import pandas as pd
import time
import os

def get_args():
	"""	
	Function to get the arguments from the command line input. Takes in a fasta file, log directory (default "."),
	number of simulations (default 10), number of fragments (default 3), and anneal_rate (default 0.999).
	"""
	parser = argparse.ArgumentParser(description='Split FASTA files')
	parser.add_argument('--fasta', help='.fasta file containing sequence', metavar='FILE')
	parser.add_argument('--logdir', help='directory to save all log files', metavar='DIR', default = '.')
	parser.add_argument('--nsims', help='number of simulations', metavar='NUM', default = 10)
	parser.add_argument('--nfrags', help='number of fragments to sample from at each iteration', metavar='NUM', default = 3)
	parser.add_argument('--anneal_rate', help='temperature and annealing parameter', metavar='NUM', default = 0.999)
	return parser.parse_args()

def main():
	"""	
	Main function to protein folding run simulation.
	"""
	#get arguments and create log directory
	args = get_args()
	fasta, logdir, nsims, nfrags, anneal_rate = args.fasta, args.logdir, args.nsims, args.nfrags, args.anneal_rate
	if not os.path.isdir(logdir):
		os.mkdir(logdir)
		
	#read in protein sequence from fasta file
	with open(fasta, 'r') as f:
		protein_name = f.readline().strip()[1:].split('_')[0]
		protein_seq = f.readline().strip()
	protein = Protein(sequence = protein_seq) 

	#create a simulation log dataframe
	fasta_path = os.path.dirname(str(fasta))
	simulation_log = pd.DataFrame(columns = ['sim_number', 'energy', 'rmsd'])

	for sim in np.arange(1, int(args.nsims) + 1):
		sim_path = logdir + "/sim_" + str(sim)
		if not os.path.isdir(sim_path):
			os.mkdir(sim_path)
		#create fragsets for each kmer
		frag_path_3 = os.path.join(fasta_path, protein_name + '_3mers.frag')
		rmsd_path_3 = os.path.join(fasta_path, protein_name + '_3mers.rmsd')
		frag_path_9 = os.path.join(fasta_path, protein_name + '_9mers.frag')
		rmsd_path_9 = os.path.join(fasta_path, protein_name + '_9mers.rmsd')		
		fragset_3mer = FragmentSet(frag_path_3, rmsd_path_3)
		fragset_9mer = FragmentSet(frag_path_9, rmsd_path_9)

		#stage 1 = run simulation with the 9mer using the input sequence
		samp_9mer = MCMCSampler(protein, fragset_9mer, nfrag = int(nfrags), nmer = 9, temp = 100, temp_end = 1, anneal_rate = float(args.anneal_rate))
		lowest_energy_conf_9mer, iteration_9mer, temperature_9mer, energy_9mer = samp_9mer.simulate()

		#stage_2 = run simulation with the 3mer using the lowest energy conformation from stage 1 folding
		samp_3mer = MCMCSampler(lowest_energy_conf_9mer, fragset_3mer, nfrag = int(nfrags), nmer = 3, temp = 1, temp_end = 0.1, anneal_rate = float(args.anneal_rate))
		lowest_energy_conf_3mer, iteration_3mer, temperature_3mer, energy_3mer = samp_3mer.simulate()
		
		#paths to save initial and best structures
		initial_path = os.path.join(sim_path, 'initial.pdb')
		best_path = os.path.join(sim_path, 'best.pdb')

		#logging values per simulation 
		df_log = pd.DataFrame(data = {'iteration': iteration_9mer + iteration_3mer, \
			'temperature': temperature_9mer + temperature_3mer, 'energy': energy_9mer + energy_3mer})
		df_log.to_csv(os.path.join(sim_path,'df_log.txt'), sep="\t", index=False)

		#saving proteins to simulation directory
		lowest_energy_conf_3mer.save_pdb(best_path)
		protein.save_pdb(initial_path)

		#performing energy minimization on the best structure from stage 2
		relaxed_prot, rmsd, score = utils.relax(best_path, os.path.join(fasta_path, protein_name + '.pdb'))
		simulation_log.loc[len(simulation_log.index)] = [sim, score, rmsd]

	#log values for each simulation (simulation number, the energy after relaxation, and the RMSD to target after relaxation)
	simulation_log.to_csv(os.path.join(args.logdir, "simulation_summary.txt"), sep="\t", index=False)	



if __name__=='__main__':
	main()

