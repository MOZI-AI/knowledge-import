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

script = "https://github.com/MOZI-AI/knowledge-import/SMPDB_pathway.py"

def atomese(node1, node1_type, node2, node2_type, node1_prefix="", node2_prefix="", predicate=False):
    if node1.lower() != "nan" and node2.lower() != "nan":
        atom1 = '({} "{}{}")'.format(node1_type, node1_prefix, node1)
        atom2 = '({} "{}{}")'.format(node2_type, node2_prefix, node2)
        if predicate:
            return '(EvaluationLink \n'+'\t(PredicateNode "'+ predicate +'")\n'+'\t(ListLink \n\t\t'+ atom1 +'\n\t\t'+ atom2 +'))\n'
        else:
            return '(MemberLink \n'+'\t'+ atom1 +'\n\t'+ atom2 +')\n'
    else:
        return ""

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
                chebi_id = str(data.iloc[r]['ChEBI ID']).split(".")[0].strip()
                smpdb_id = str(data.iloc[r]['SMPDB ID']).strip()
                chebi_name = str(data.iloc[r]['IUPAC']).strip()

                if not chebi_id in chebis:
                    chebis.append(chebi_id)
                if not smpdb_id in pathways:
                    pathways.append(smpdb_id) 

                f.write(atomese(chebi_id, 'MoleculeNode', smpdb_id, 'ConceptNode', node1_prefix='ChEBI:') )
                g.write(atomese(chebi_id, 'MoleculeNode', smpdb_id, 'ConceptNode', node1_prefix='ChEBI:') )
                f.write(atomese(chebi_id, 'MoleculeNode', chebi_name, 'ConceptNode', node1_prefix='ChEBI:', predicate='has_name') )

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
                protein = str(data.iloc[r]['Uniprot ID']).split(".")[0].strip()
                protein_name = str(data.iloc[r]['Protein Name']).strip() 
                gene = str(data.iloc[r]['Gene Name']).upper().strip()
                smpdb_id = str(data.iloc[r]['SMPDB ID']).strip()
                smpdb_name = str(data.iloc[r]['Pathway Name']).strip()

                if not protein in proteins:
                   proteins.append(protein)
                if not gene in genes:
                   genes.append(gene)
                if not smpdb_id in pathways:
                   pathways.append(smpdb_id)   

                f.write(atomese(gene, 'GeneNode', protein, 'MoleculeNode',node2_prefix='Uniprot:', predicate='expresses') )
                f.write(atomese(gene, 'GeneNode', smpdb_id, 'ConceptNode'))
                if gene_level:
                   g.write(atomese(gene, 'GeneNode', smpdb_id, 'ConceptNode'))
                f.write(atomese(protein, 'MoleculeNode', smpdb_id, 'ConceptNode', node1_prefix='Uniprot:'))
                f.write(atomese(smpdb_id, 'ConceptNode', smpdb_name, 'ConceptNode', predicate='has_name'))
                f.write(atomese(protein, 'MoleculeNode', protein_name, 'ConceptNode',node1_prefix='Uniprot:', predicate='has_name'))

            # print("Imported "+filename)
    
    num_pathways = {"SMPDB Pathway": len(pathways)} 
    metadata.update_meta("smpdb_proteins: Latest",source, script,genes=len(genes), prot=len(proteins),pathways=num_pathways)
    print("Done. Check dataset/smpdb_protein.scm and gene-level/smpdb_gene.scm")

## Import them
if __name__ == "__main__":
	print("Import the following files from Small molecule database \n" +
	      "Press M to import Metabolite names linked to SMPDB pathways \n"+
	      "Press P to import Protein names linked to SMPDB pathways \n"+
	      "Press B for both\n")
	option = input()
	if option == "P" or option == "p":
		import_proteins(gene_level=True)
	elif option == "M" or option == "m":
		import_metabolites(gene_level=True)
	elif option == "B" or option == "b":
		import_proteins(gene_level=True)
		import_metabolites(gene_level=True)
	else:
	    print("Incorect option, Try again")