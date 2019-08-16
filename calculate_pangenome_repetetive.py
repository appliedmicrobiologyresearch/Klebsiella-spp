#!/usr/bin/env python


import sys
import numpy as np
import random


#usage: python calculate_pangenome_repetetive.py gene_presence_absence.Rtab

def main():



	inputOptions = sys.argv[1:]
	
	
	strain2genes={}
	input_file = [n for n in open(inputOptions[0],'r').read().replace("\r","").split("\n") if len(n)>0]

	for line in input_file:
		if line[0:5]=="Gene	":
			strains=line.replace("Gene	","").split("\t")
			for strain in strains:
				strain2genes[strain]=list()
		else:
			gene=line.split("\t")[0]
			
			for strain,gene_present in zip(strains,line.split("\t")[1:]):
				if gene_present != "0":
					strain2genes[strain].append(gene)
		
	strains=strain2genes.keys()

	repetitions=100
		
	pan_genomes={}
	core_genomes={}
	for i in range(1,len(strains)+1):
		pan_genomes[i]=list()
		core_genomes[i]=list()
		
	for i in range(0,repetitions):
		random.shuffle(strains)
		
		pangenome=list(strain2genes[strains[0]])
		coregenome=list(strain2genes[strains[0]])
		counter=0
			
		for strain in strains:	
				
			counter+=1
		
			coregenome=set(coregenome).intersection(set(strain2genes[strain]))
			pangenome=set(pangenome).union(set(strain2genes[strain]))
		
		
			pan_genomes[counter].append(len(pangenome))
			core_genomes[counter].append(len(coregenome))
		
	outfile = open('result.txt', 'w')		
	outfile.write("replicate\tnumber_of_strains\tnumber_of_orthologous_cluster\tsd\n")
	for i in range(1,len(strains)+1):
			
		outfile.write("result_Core"+"\t"+str(i)+"\t"+str(np.mean(core_genomes[i]))+"\t"+str(np.std(core_genomes[i]))+"\n")
		outfile.write("result_Pan"+"\t"+str(i)+"\t"+str(np.mean(pan_genomes[i]))+"\t"+str(np.std(pan_genomes[i]))+"\n")
		



main()
