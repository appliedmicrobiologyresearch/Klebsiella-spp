# Daniel Wüthrich

import sys
import os
import networkx as nx
import matplotlib.pyplot as plt
import math
import random, string
import operator


#usage: ref_protein1.faa protein2.faa ref_annotation1.gff annotation2.gff strain1.strain2

inputOptions = sys.argv[1:]

def main():

	ANIdata = [n for n in open("ANIm_percentage_identity.tab",'r').read().replace("\r","").split("\n") if len(n)>0]

	edges=[]

	strains=ANIdata[0].split("\t")[1:]
	for line in ANIdata[1:]:
		strain1=line.split("\t")[0]

		for strain2, value in zip(strains, line.split("\t")[1:]):
			if float(value) >= 0.999:
				edges.append([strain1,strain2])
			
	graph=nx.Graph()
	graph.add_edges_from(edges)
	
	counter=0

	for cluster in nx.connected_components(graph):
		print list(cluster)[0]
		counter+=1
		


main()