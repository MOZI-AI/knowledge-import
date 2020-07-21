__author__ = "Hedra"
__email__ = "hedra@singularitynet.io"


# The following script imports the Physical Entity (PE) Identifier mapping files from https://reactome.org/download-data

# Requires: NCBI2Reactome_PE_Pathway.txt
#	    UniProt2Reactome_PE_Pathway.txt
# 	    ChEBI2Reactome_PE_Pathway.txt

# from https://reactome.org/download/current/

import pandas as pd
import wget
import os
import sys
import metadata
from datetime import date
from atomwrappers import *
import argparse

# Get each of the files first

# URL's

ncbi = "https://reactome.org/download/current/NCBI2Reactome_PE_Pathway.txt"
uniprot = "https://reactome.org/download/current/UniProt2Reactome_PE_Pathway.txt"
chebi = "https://reactome.org/download/current/ChEBI2Reactome_PE_Pathway.txt"
script = "https://github.com/MOZI-AI/knowledge-import/PE_Identifier_mapping.py"

# If you have the files downloaded, make sure the file names are the same 
# Or modify the file names in this code to match yours.

def get_data(name):
	for data in name:
		if(not os.path.isfile('raw_data/{}'.format(data.split('/')[-1]))):
			print("Downloading the datasets, It might take a while") 
			wget.download(data, "raw_data/")
			print("Done")

# The column 'R_PE_name' contains the Gene Symbol and its location information, so we need to split it
# Example: A1BG [extracellular region]
# A1BG is the Gene symbol and 'extracellular region' is the gene location
# some has extra symbols which needs preprocessing e.g. CCL5(24-91) [extracellular region], p-S472-AKT3 [plasma membrane]

def find_location(PEname, filter=False):
	if "[" in PEname and "]" in PEname:
		loc = PEname[PEname.find("[")+1:PEname.find("]")]
		gene = PEname.split("[" +loc +"]")[0]
	else:
		loc = ""
		gene = PEname
	gene = gene.replace(gene[gene.find("("):PEname.find(")")+1], "").replace(")", "").replace("(","")
	if "-" in gene:
		gene = [i for i in gene.split("-") if not i.strip().isdigit()][-1]
	gene = gene.strip()
	if filter:
		return gene
	return gene,loc

def import_dataset(dataset, delim, without_location=False):
	print("Started importing " + dataset)
	if "UniProt" in dataset or "ChEBI" in dataset:
		data = pd.read_csv(dataset, low_memory=False, delimiter=delim, names=["db_id", "R_PE_id", "R_PE_name","pathway","url","event_name", "evidence_code", "species","un1","un2","un3","un4","un5","un6"])

	else:	
		data = pd.read_csv(dataset, low_memory=False, delimiter=delim, names=["db_id", "R_PE_id", "R_PE_name","pathway","url","event_name", "evidence_code", "species"])
	mapping_entrez = pd.read_csv("raw_data/entrez.txt", low_memory=False, sep="\t")
	# Take only symbols of Human species
	data_human = data[data['species'] == 'Homo sapiens'][['db_id','R_PE_name','pathway']]

	if without_location:
		if not os.path.exists(os.path.join(os.getcwd(), 'gene-level-without-location')):
			os.makedirs('gene-level-without-location')  
		file_name = open("gene-level-without-location/"+dataset.split("/")[-1]+"_without_location_{}.scm".format(str(date.today())), "w")

	if not os.path.exists(os.path.join(os.getcwd(), 'dataset')):
		os.makedirs('dataset')
	with open("dataset/"+dataset.split("/")[-1]+"_{}.scm".format(str(date.today())), 'w') as f:
		if "NCBI" in dataset:
			genes = []
			pathways = []
			infered = {}
			gene_symbols = mapping_entrez["Approved symbol"].values
			for i in range(len(data_human)):
				gene_sym, location = find_location(data_human.iloc[i]['R_PE_name'])
				pathway = data_human.iloc[i]['pathway']
				db_id = data_human.iloc[i]['db_id']
				try:
					gene = mapping_entrez[mapping_entrez["NCBI Gene ID"] == int(db_id)]["Approved symbol"].values[0]
				except:
					if len(gene_sym.split(" ")) > 1:
						if str(db_id) in infered.keys():
							gene = infered[str(db_id)]
						else:
							# non_exist.append(gene_sym + '\t' +str(db_id))
							continue
					else:
						if gene_sym in gene_symbols:
							gene = gene_sym
							infered[str(db_id)] = gene
						else:
							continue
				if not gene.isdigit() and not len(gene) == 1 and not gene in ["", " "]:
					gene = gene.strip()
					member = CMemberLink(CGeneNode(gene),ReactomeNode(pathway))
					eva = CEvaluationLink(CPredicateNode("has_location"), CListLink(CGeneNode(gene), CConceptNode(location)))
					cont = CContextLink(member, eva)
					f.write(cont.recursive_print())
					if without_location:
						file_name.write(member.recursive_print())
					if not gene in genes:
						genes.append(gene)
					if not pathway in pathways:
						pathways.append(pathway) 
			version = "NCBI2reactome_pathway_mapping:latest"
			num_pathways = {"Reactome Pathway": len(pathways)}
			metadata.update_meta(version,ncbi,script,genes=len(genes),pathways=num_pathways)
		elif "UniProt" in dataset:
			molecules = []
			pathways = []
			for i in range(len(data_human)):
				prot = str(data_human.iloc[i]['R_PE_name'])
				loc = prot[prot.find("[")+1:prot.find("]")]
				prot_name = prot.split("[" +loc +"]")[0]
				pathway = data_human.iloc[i]['pathway']
				protein = [i for i in str(data_human.iloc[i]['db_id']).split("-") if not i.strip().isdigit()][-1]
				protein = protein.strip()
				member = CMemberLink(ProteinNode(protein), ReactomeNode(pathway))
				eva_loc = CEvaluationLink(CPredicateNode("has_location"), CListLink(ProteinNode(protein), CConceptNode(loc)))
				eva_name = CEvaluationLink(CPredicateNode("has_name"), CListLink(ProteinNode(protein), CConceptNode(prot_name)))
				cont = CContextLink(member, eva_loc)
				f.write(cont.recursive_print())
				if without_location:
					file_name.write(member.recursive_print())
				if not protein in molecules:
					molecules.append(protein)
					f.write(eva_name.recursive_print())
				if not pathway in pathways:
					pathways.append(pathway)
			version = "Uniprot2reactome_pathway_mapping:latest"
			num_pathways = {"Reactome Pathway": len(pathways)}
			metadata.update_meta(version,ncbi,script,prot=len(molecules),pathways=num_pathways)
		elif "ChEBI" in dataset:
			molecules = []
			pathways = []
			for i in range(len(data_human)):
				chebi = str(data_human.iloc[i]['R_PE_name'])
				loc = chebi[chebi.find("[")+1:chebi.find("]")]
				chebi_name = chebi.split("[" +loc +"]")[0].replace('"',"")
				chebi_id = str(data_human.iloc[i]['db_id'])
				if not chebi_id is "nan":
					chebi_id = chebi_id.strip()
					pathway = data_human.iloc[i]['pathway']
					member = CMemberLink(ChebiNode(chebi_id), ReactomeNode(pathway))
					eva_loc = CEvaluationLink(CPredicateNode("has_location"), CListLink(ChebiNode(chebi_id), CConceptNode(loc)))
					eva_name = CEvaluationLink(CPredicateNode("has_name"), CListLink(ChebiNode(chebi_id), CConceptNode(chebi_name)))
					cont = CContextLink(member, eva_loc)
					f.write(cont.recursive_print())
					if without_location:
						file_name.write(member.recursive_print())
					if not chebi_id in molecules:
						molecules.append(chebi_id)
						f.write(eva_name.recursive_print())
					if not pathway in pathways:
						pathways.append(pathway)
			version = "Chebi2reactome_pathway_mapping:latest"
			num_pathways = {"Reactome Pathway": len(pathways)}
			metadata.update_meta(version,ncbi,script,chebi=len(molecules),pathways=num_pathways)
	print("Done")

def parse_arg():
	parser = argparse.ArgumentParser(description='Import Physical entities to pathway mapping from https://reactome.org')
	parser.add_argument('--option', type=str, default='all',
                        help='which dataset to import: prot, chebi, ncbi')
	return parser.parse_args()

if __name__ == "__main__":

	option = parse_arg().option
	if option == "ncbi":
		get_data([ncbi])
		import_dataset('raw_data/NCBI2Reactome_PE_Pathway.txt', '\t', without_location=True)

	elif option == "prot":
		get_data([uniprot])
		import_dataset('raw_data/UniProt2Reactome_PE_Pathway.txt', '\t')

	elif option == "chebi":
		get_data([chebi])
		import_dataset('raw_data/ChEBI2Reactome_PE_Pathway.txt', '\t', without_location=True)

	elif option == "all":
		get_data([chebi, uniprot, ncbi])
		import_dataset('raw_data/NCBI2Reactome_PE_Pathway.txt', '\t', without_location=True)

		import_dataset('raw_data/UniProt2Reactome_PE_Pathway.txt', '\t')

		import_dataset('raw_data/ChEBI2Reactome_PE_Pathway.txt', '\t', without_location=True)
	else:
	    print("Specify the correct option, Try again")
