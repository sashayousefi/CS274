import utils
import numpy as np
import pandas as pd
import os
import re

class FragmentSet(object):
	def __init__(self, fragfile, rmsdfile):
		"""
		This class contains the fragment library for the input protein. It must do the following:
		- Read in fragment file and parse fragments at each position. Fragment files are of the form <protein>_<frag_len>mers.frag
		- Read in RMSD file containing pre-calculated RMSD to native structure for each fragment at each position.
		- Based on fragments and their corresponding RMSDs, rank fragments at each position by RMSD
		
		"""
		self.pos_dict = {} #dictionary to store fragments at each position
		with open(fragfile, 'r') as f:
			curr_position = 0
			for line in f:
				line = line.strip().split()
				if line == []: #if we reach the end of the file, end parsing
					continue;
				elif str(line[0]) == 'position:': #read in the position in the protein
					curr_position = int(line[1])
					self.pos_dict[curr_position] = {}
				else:
					frag_idx = int(re.sub("[^0-9]", "", line[-1])) #get the fragment index
					if frag_idx not in self.pos_dict[curr_position].keys():
						self.pos_dict[curr_position][frag_idx] = []
					phi = float(line[5]) #store the fragment phi angle
					psi = float(line[6]) #store the fragment psi angle. 
					self.pos_dict[curr_position][frag_idx].append((phi, psi)) 
		self.rmsd = pd.read_csv(rmsdfile, sep = '\t',  names=["position", "fragment", "rmsd"]) #read in fragment rmsd values
		self.rmsd = self.rmsd.sort_values(['position', 'rmsd']) #sort rmsd values from smalest to largest for each position



	def get_lowRMS_fragments(self, pos, N):
		"""
		Returns the top-ranked fragments by RMSD at a defined position in the chain
		--------
		Params
			- pos (int): fragment position in chain (1-indexed)
			- N (int): number of fragments to return
		Returns
			- lowRMS_fragments (list): top N fragments at pos by RMSD. This should be a list of lists of (phi, psi) tuples. 
			  For example, a 3-mer fragment could be represented as the following: [(-60.892, 142.456), (-72.281, 128.933), (-132.337, -175.477)]
		"""
		lowRMS_fragments = []
		cand_frag = self.rmsd[self.rmsd['position'] == pos]['fragment'].values
		for i in np.arange(int(N)):
			lowRMS_fragments.append(self.pos_dict[pos][cand_frag[i]])
		return lowRMS_fragments #returns the top N candidate fragments with the lowest rmsd
			






