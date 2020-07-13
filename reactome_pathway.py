__author__ = "Hedra"
__email__ = "hedra@singularitynet.io"

# The following script imports Reactome pathways and their relationship from https://reactome.org

# Requires: 
	# https://reactome.org/download/current/ReactomePathwaysRelation.txt
	# https://reactome.org/download/current/ReactomePathways.txt

import pandas as pd
from urllib.request import urlopen
import os
import metadata
from datetime import date
from atomwrappers import *

# URL's
pathway_rln = "https://reactome.org/download/current/ReactomePathwaysRelation.txt"
pathway = "https://reactome.org/download/current/ReactomePathways.txt"

def homosapien(pathway_id):
	if pathway_id.startswith("R-HSA"):
		return True
	else:
		return False

if not os.path.isfile('ReactomePathwaysRelation.txt'):
	print("Downloading ReactomePathwaysRelation.txt")
	pathway_relation = pd.read_csv(urlopen(pathway_rln), low_memory=False, delimiter='\t', names=["parent", "child"])
	print("Done")
else:
	pathway_relation = pd.read_csv('ReactomePathwaysRelation.txt', low_memory=False, delimiter='\t', names=["parent", "child"])  

if not os.path.isfile('ReactomePathways.txt'):
	print("Downloading ReactomePathways.txt")
	pathway_list = pd.read_csv(urlopen(pathway), low_memory=False, delimiter='\t', names=["ID", "name", "Species"])
	print("Done")
else:
	pathway_list = pd.read_csv('ReactomePathways.txt', low_memory=False, delimiter='\t', names=["ID", "name", "Species"]) 

pathway_list = pathway_list[pathway_list['Species']=='Homo sapiens'] 
max_len = max(len(pathway_list), len(pathway_relation)) 

print("Started importing")

script = "https://github.com/MOZI-AI/knowledge-import/reactome_pathway.py"
pathways = pathway_relation['parent'].values + pathway_relation['child'].values
if not os.path.exists(os.path.join(os.getcwd(), 'dataset')):
    os.makedirs('dataset')
output = "dataset/reactome_{}.scm".format(str(date.today()))
with open(output, 'w') as f:
    for i in range(len(pathway_list)):
        pw_name = pathway_list.iloc[i]['name']
        pw_id = pathway_list.iloc[i]['ID']
        eva_name = CEvaluationLink(CPredicateNode("has_name"), CListLink(ReactomeNode(pw_id), CConceptNode(pw_name)))
        f.write(eva_name.recursive_print() + "\n")

    for i in range(len(pathway_relation)):
        pw_parent = pathway_relation.iloc[i]['parent']
        pw_child = pathway_relation.iloc[i]['child']										
        if homosapien(pw_child) and homosapien(pw_parent): 
            inherit = CInheritanceLink(ReactomeNode(pw_child), ReactomeNode(pw_parent))
            f.write(inherit.recursive_print() + "\n")						

num_pathways = {"Reactome Pathway": len(set(pathways))}
metadata.update_meta("Reactome Pathways relationship:latest", 
        pathway_rln+" "+pathway,script,pathways=num_pathways)

print("Done, check {}".format(output))
