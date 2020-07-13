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
from datetime import date
from atomwrappers import *

# Define helper functions 

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
    with open("dataset/entrez_to_protein_{}.scm".format(str(date.today())), 'w') as f:
        for i in range(len(data)):
            try:
                g = data.iloc[i]['symbol']
                p = data.iloc[i]['uniprot'].strip()
                genes.append(g)
                prot.append(p)
                expresion = CEvaluationLink(CPredicateNode("expresses"), CListLink(CGeneNode(g),ProteinNode(p)))
                f.write(expresion.recursive_print())
            except:
                continue
            if not math.isnan(data.iloc[i]['entrez']):
                entrez_id = str(int(data.iloc[i]['entrez']))
                has_entrez = CEvaluationLink(CPredicateNode("has_entrez_id"), CListLink(CGeneNode(g),CConceptNode("entrez:"+entrez_id)))
                f.write(has_entrez.recursive_print())

        metadata.update_meta("gene2proteinMapping:latest", 
        "entrez2uniprot.csv",script,genes=len(set(genes)),prot=len(set(prot)))

        print("Done")

