__author__ = "Hedra"
__email__ = "hedra@singularitynet.io"


# The following script maps Uniprot to GO 
# Requires:  goa_human_isoform_valid.gaf
# source: http://current.geneontology.org/annotations/goa_human_isoform.gaf.gz

import os
import wget
import gzip
import metadata
from datetime import date
from atomwrappers import *
import json
import find_gons

dataset_url = "http://current.geneontology.org/annotations/goa_human_isoform.gaf.gz"
lines = []
prot = []
go = []
if not os.path.isfile('raw_data/goa_human_isoform_valid.gaf'):
    print("Downloading dataset")
    lines = gzip.open(wget.download(dataset_url, "raw_data/")).readlines()
    lines = [l.decode("utf-8") for l in lines]
else:
    lines = open('raw_data/goa_human_isoform_valid.gaf').readlines()

with open("raw_data/go-namespace.json", "r") as ns:
    go_namespace = json.load(ns)
init_namespace = len(go_namespace)

with open("dataset/uniprot2GO_{}.scm".format(str(date.today())), 'w') as f:
    print("\nStarted importing")
    for i in lines:
        if 'UniProtKB' in i:
            go_namespace, go_term = find_gons.find_type(i.split('\t')[4], go_namespace)
            protein = ProteinNode(i.split('\t')[1])
            if go_term:
                f.write(CMemberLink(protein,go_term).recursive_print() + "\n")
            prot.append(i.split('\t')[1])
            go.append(go_term)

if len(go_namespace) > init_namespace:
    with open("raw_data/go-namespace.json", "w") as ns:
        json.dump(go_namespace, ns, indent=2)

script = "https://github.com/MOZI-AI/knowledge-import/uniprot2GO.py"
metadata.update_meta("Uniprot-GO:latest", dataset_url,script,prot=len(set(prot)), goterms={"go-terms":len(set(go))})
print("Done, check dataset/uniprot2GO.scm")