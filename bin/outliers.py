#!/usr/local/bin/python

from optparse import OptionParser
import os

import h5py
import numpy as np
import pandas as pd

from operator import itemgetter


def main():
	usage = 'usage: <scores_file>'
	parser = OptionParser(usage)
	parser.add_option('-t', dest = "tissues",
		default = None, type = int,
		help = "Number of target indexes")
	parser.add_option('-g', dest = "genes",
		default = None, type = int,
		help = "Number of target genes")
	parser.add_option('-o', dest = "out_dir",
		default = "outliers_results")
	parser.add_option('-n', dest = "name",
		default = "outliers.tsv")
	parser.add_option('--sd', dest = 'sd',
		default = 6, type = int,
		help = "Standard deviation coeficient")
	(options, args) = parser.parse_args()

	h5_file = h5py.File(args[0], 'r')
	
	if not os.path.isdir(options.out_dir):
		os.mkdir(options.out_dir)

	if options.genes : 
		genes = range(options.genes)
	else : 
		genes = range(h5_file['seqs'].shape[0])

	if options.tissues : 
		tissues = range(options.tissues)
	else : 
		tissues = range(h5_file['sum'].shape[3])

	### Genome nt sequence from 1hot
	genomes_1hot = h5_file['seqs'][min(genes):max(genes)+1]
	genomes = []
	for genome_1hot in genomes_1hot:
		genome = []
		for nt in genome_1hot:
			if nt[0] :
				genome.append("A")
			elif nt[1] :
				genome.append("C")
			elif nt[2] :
				genome.append("G")
			else :
				genome.append("T")

		genomes.append(genome)


	### Un-normalize scores
	scores = h5_file['sum'][min(genes):max(genes)+1,:,:,min(tissues):max(tissues)+1]
		
	norm_scores = []	

	for tissue in range(scores.shape[3]):
		ref = scores[:,:,:,tissue] * genomes_1hot
		ref = np.sum(ref, axis = 2)
		tmp = scores[:,:,:,tissue]

		for i in range(4):
			tmp[:,:,i] = tmp[:,:,i] - ref

		norm_scores.append(tmp)

	norm_scores = np.array(norm_scores)


	### Outliers
	mean_tg = np.mean(norm_scores, axis = (2,3)) ### mean by gene and by tissue
	sd_tg = np.std(norm_scores, axis = (2,3))  ### standard deviation by gene and by tissue

	outliers = []
	positions = []
	mutations = []
	genes = []
	tissues = []
	reference = []
	distance = []
	chrom = []
	start = []
	stop = []
		
	for gene in range(mean_tg.shape[1]) :

		for tissue in range(mean_tg.shape[0]):		

			### for sup outliers :
			sup = (norm_scores[tissue, gene, :, :] > 
			(mean_tg[tissue, gene] + options.sd*sd_tg[tissue, gene]))
			
			pos_sup = np.where(sup)[0]
			mut_sup = np.where(sup)[1]
			out_sup = norm_scores[tissue, gene,:,:][sup]			
			dist_sup = abs(out_sup-mean_tg[tissue, gene])/sd_tg[tissue, gene]			

			### for inf outliers
			inf = (norm_scores[tissue, gene, :, :] < 
			(mean_tg[tissue, gene] - options.sd*sd_tg[tissue, gene]))

			pos_inf = np.where(inf)[0]
			mut_inf = np.where(inf)[1]
			out_inf = norm_scores[tissue, gene,:,:][inf]			
			dist_inf = abs(out_inf-mean_tg[tissue, gene])/sd_tg[tissue, gene]

			### concat things
			positions = [*positions, *pos_sup, *pos_inf]
			mutations = [*mutations, *mut_sup, *mut_inf]
			outliers = [*outliers, *out_sup, *out_inf]
			distance = [*distance, *dist_sup, *dist_inf]			
												
			if len(pos_sup) != 0 :				
				reference.extend(itemgetter(*pos_sup)(genomes[gene]))
			else : 
				reference = reference
			
			if len(pos_inf) != 0:				
				reference.extend(itemgetter(*pos_inf)(genomes[gene]))
			else:
				reference = reference

			tissues.extend([tissue]*(len(pos_inf) + len(pos_sup)))
			genes.extend([gene]*(len(pos_inf) + len(pos_sup)))
			chrom.extend([h5_file['chrom'][gene]]*(len(pos_inf) + len(pos_sup)))
			start.extend([h5_file['start'][gene]]*(len(pos_inf) + len(pos_sup)))
			stop.extend([h5_file['end'][gene]]*(len(pos_inf) + len(pos_sup)))
							

	### Output
	for mut in range(len(mutations)):
		if mutations[mut] == 0 :
			mutations[mut] = "A"
		elif mutations[mut] == 1 :
			mutations[mut] = "C"
		elif mutations[mut] == 2 :
			mutations[mut] = "G"
		else :
			mutations[mut] = "T"

	df = pd.DataFrame({
	"position":positions,
	"mutation":mutations, 
	"outlier":outliers,
	"tissue":tissues,
	"gene":genes,
	"reference":reference,
	"chrom":chrom,
	"start":start,
	"stop":stop,
	"distance":distance
	})

	df.to_csv("%s/%s" % (options.out_dir, options.name),
		sep = "\t",
		index = False)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
  main()
