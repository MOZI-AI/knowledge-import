__author__ = "Hedra"
__email__ = "hedra@singularitynet.io"

# The following script imports the following files from Small molecule database at http://smpdb.ca/

	#1 Metabolite names linked to SMPDB pathways CSV (includes KEGG and ChEBI IDs)	
	#2 Protein names linked to SMPDB pathways CSV (includes UniProt IDs)

# Requires: smpdb_metabolites.csv.zip
#	    smpdb_proteins.csv.zip

# from 	http://smpdb.ca/downloads/smpdb_metabolites.csv.zip
# 	http://smpdb.ca/downloads/smpdb_proteins.csv.zip


import pandas as pd
import os
from urllib.request import urlopen
from zipfile import ZipFile
from io import BytesIO
import wget
import metadata
from datetime import date
from atomwrappers import *
import argparse

script = "https://github.com/MOZI-AI/knowledge-import/SMPDB_pathway.py"

def import_metabolites(gene_level=False):
    pathways = []
    chebis = []
    source = "http://smpdb.ca/downloads/smpdb_metabolites.csv.zip"

    if not "smpdb_metabolites.csv.zip" in os.listdir("raw_data/"):

        print("Started downloading smpdb_metabolites.csv, it will take some time to download")
        wget.download(source, "raw_data/")
        
    ZipFile("raw_data/smpdb_metabolites.csv.zip").extractall("raw_data/smpdb_chebi")
    pathway_chebi = os.listdir("raw_data/smpdb_chebi")
     
    print("Started importing {} files of smpdb_metabolites".format(len(pathway_chebi)))

    # For a gene level dataset, excelude the name
    if gene_level:
        if not os.path.exists(os.path.join(os.getcwd(), 'gene-level')):
            os.makedirs('gene-level')  
        g = open("gene-level/smpdb_chebi_{}.scm".format(str(date.today())), "w")

    with open("dataset/smpdb_chebi_{}.scm".format(str(date.today())), 'w') as f:
        for filename in pathway_chebi:
            data = pd.read_csv("raw_data/smpdb_chebi/"+filename, low_memory=False)

            for r,c in data.iterrows():
                chebi_id = filter_nan(str(data.iloc[r]['ChEBI ID']).split(".")[0].strip())
                smpdb_id = filter_nan(str(data.iloc[r]['SMPDB ID']).strip())
                chebi_name = filter_nan(str(data.iloc[r]['IUPAC']).strip())
                try:
                    if chebi_id: 
                        chebi_id= "ChEBI:" + chebi_id 
                    member = CMemberLink(ChebiNode(chebi_id), SMPNode(smpdb_id))
                    f.write(member.recursive_print() + "\n")
                    if gene_level:
                        g.write(member.recursive_print() + "\n")
                    if not chebi_id in chebis:
                        ch_name = CEvaluationLink(CPredicateNode("has_name"), CListLink(ChebiNode(chebi_id), CConceptNode(chebi_name)))
                        f.write(ch_name.recursive_print() + "\n")
                        chebis.append(chebi_id)
                    if not smpdb_id in pathways:
                        pathways.append(smpdb_id)
                except AttributeError:
                    print("Null value detected")
                    continue

    num_pathways = {"SMPDB Pathway": len(pathways)} 
    metadata.update_meta("smpdb_metabolites: Latest",source, script,chebi=len(chebis), pathways=num_pathways)        
    print("Done. Check dataset/smpdb_chebi.scm")

def import_proteins(gene_level=False):
    pathways = []
    proteins = []
    genes = []
    source = "http://smpdb.ca/downloads/smpdb_proteins.csv.zip"

    if not "smpdb_proteins.csv.zip" in os.listdir("raw_data/"):
        print("Started downloading smpdb_proteins.csv, It will take some time to download \n")
        wget.download(source, "raw_data")

    ZipFile("raw_data/smpdb_proteins.csv.zip").extractall("raw_data/smpdb_prot")
    pathway_prot = os.listdir("raw_data/smpdb_prot")

    print("Started importing {} files of smpdb_proteins".format(len(pathway_prot)))
    
    if gene_level:
        g = open("gene-level/smpdb_gene_{}.scm".format(str(date.today())), "w")

    with open("dataset/smpdb_protein_{}.scm".format(str(date.today())), 'w') as f:
        for filename in pathway_prot:
            data = pd.read_csv("raw_data/smpdb_prot/"+filename, low_memory=False)
            for r,c in data.iterrows():
                protein = filter_nan(str(data.iloc[r]['Uniprot ID']).split(".")[0].strip())
                protein_name = filter_nan(str(data.iloc[r]['Protein Name']).strip()) 
                gene = filter_nan(str(data.iloc[r]['Gene Name']).upper().strip())
                smpdb_id = filter_nan(str(data.iloc[r]['SMPDB ID']).strip())
                smpdb_name = filter_nan(str(data.iloc[r]['Pathway Name']).strip())
                try:
                    member = CMemberLink(CGeneNode(gene), SMPNode(smpdb_id))
                    f.write(member.recursive_print() + "\n")
                    expression = CEvaluationLink(CPredicateNode("expresses"), CListLink(CGeneNode(gene), ProteinNode(protein)))
                    f.write(expression.recursive_print() + "\n")                    
                    if gene_level:
                        g.write(member.recursive_print() + "\n")
                    if not smpdb_id in pathways:
                        smp_name = CEvaluationLink(CPredicateNode("has_name"), CListLink(SMPNode(smpdb_id), CConceptNode(smpdb_name)))
                        f.write(smp_name.recursive_print() + "\n")
                        pathways.append(smpdb_id)
                    if not protein in proteins:
                        prot_name = CEvaluationLink(CPredicateNode("has_name"), CListLink(ProteinNode(protein), CConceptNode(protein_name)))
                        f.write(prot_name.recursive_print() + "\n")
                        proteins.append(protein)
                    if not gene in genes:
                        genes.append(gene)
                except AttributeError:
                    print("Null value detected")
                    continue
    
    num_pathways = {"SMPDB Pathway": len(pathways)} 
    metadata.update_meta("smpdb_proteins: Latest",source, script,genes=len(genes), prot=len(proteins),pathways=num_pathways)
    print("Done. Check dataset/smpdb_protein.scm and gene-level/smpdb_gene.scm")

def filter_nan(value):
    if str(value).lower() == "nan":
        return False
    else:
        return str(value)

def parse_arg():
	parser = argparse.ArgumentParser(description='Imports metabolite and protein sets of SMPDB pathway from http://smpdb.ca/downloads')
	parser.add_argument('--option', type=str, default='all',
                        help='which dataset to import: P for proteins, M for metabolites')
	return parser.parse_args()

if __name__ == "__main__":

	option = parse_arg().option
	if option == "P" or option == "p":
		import_proteins(gene_level=True)
	elif option == "M" or option == "m":
		import_metabolites(gene_level=True)
	elif option == "B" or option == "all":
		import_proteins(gene_level=True)
		import_metabolites(gene_level=True)
	else:
	    print("Incorect option, Try again")
