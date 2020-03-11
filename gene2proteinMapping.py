__author__ = "Hedra"
__email__ = "hedra@singularitynet.io"


# The following script imports the id table gene2preteinIDs.csv  
# this file is extracted using uniprotIDmap.R  

# Run uniprotIDmap.R to get entrez2uniprot.csv or Download it here
# https://gitlab.com/opencog-bio/ocImport/raw/master/data/entrez2uniprot.csv.xz

# Requires: entrez2uniprot.csv

import pandas as pd
import os
import math
import metadata

# Define helper functions 

def expres(predicate, node1, node2):
    return ""+'\n(EvaluationLink \n'+'\t(PredicateNode "'+ predicate +'")\n'+'\t\t(ListLink \n'+ "\t\t"+ node1 + "\t\t"+ node2 +'))\n'

script = "https://github.com/MOZI-AI/knowledge-import/gene2proteinMapping.py"

if not "entrez2uniprot.csv" in os.listdir("raw_data/"):
	print("Generate the entres to protein ID mapping table first \n" )

else:
    data = pd.read_csv("raw_data/entrez2uniprot.csv", dtype={'uniprot': str, 'entrez': float, 'symbol': str})
    print("Started importing")
    prot = []
    genes = []
    if not os.path.exists(os.getcwd()+'/dataset'):
        os.makedirs('dataset')
    with open("dataset/entrez_to_protein.scm", 'w') as f:
        for i in range(len(data)):
            try:
                g = data.iloc[i]['symbol']
                p = data.iloc[i]['uniprot'].strip()
                genes.append(g)
                prot.append(p)
                f.write(expres("expresses", '(GeneNode '+ '"' + g +'")\n', '(MoleculeNode "'+'Uniprot:'+ p+'")\n'))
            except:
                continue
            if not math.isnan(data.iloc[i]['entrez']):
                f.write(expres("has_entrez_id", '(GeneNode '+ '"' + g +'")\n', '(ConceptNode "'+'entrez:'+ str(int(data.iloc[i]['entrez']))+'")\n'))
        metadata.update_meta("gene2proteinMapping:latest", 
        "entrez2uniprot.csv",script,genes=len(set(genes)),prot=len(set(prot)))

        print("Done")

