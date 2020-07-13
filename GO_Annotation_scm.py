# imports Go Annotation to atomese
# Requires: file gene_association.goa_ref_human.gz from http://geneontology.org/gene-associations/gene_association.goa_ref_human.gz
import wget
import gzip
import os
import metadata
from datetime import date
from atomwrappers import *
import pandas as pd
import json
import find_gons

source = "http://current.geneontology.org/annotations/goa_human.gaf.gz"
if not os.path.exists("raw_data/goa_human.gaf.gz"):
    dataset = wget.download(source, "raw_data")
with gzip.open("raw_data/goa_human.gaf.gz", "rb") as f:
    lines = f.readlines()
lines = [l.decode("utf-8") for l in lines]
line_no = []

with open("namespace_R.json", "r") as ns:
    go_namespace = json.load(ns)
init_namespace = len(go_namespace)

for num, line in enumerate(lines , 1):
  if "UniProtKB" in line :
      line_no.append(num)

lines_annotate = lines[line_no[0] -1 : len(lines)]

#open file to write
if not os.path.exists(os.path.join(os.getcwd(), 'gene-level')):
      os.makedirs('gene-level')
if not os.path.exists(os.path.join(os.getcwd(), 'dataset')):
      os.makedirs('dataset') 
scm_output = open('dataset/GO_annotation_{}.scm'.format(str(date.today())), 'w')
scm_gene_level = open('gene-level/GO_annotation_gene-level_{}.scm'.format(str(date.today())), 'w')

#add GOC Validation Date
scm_output.write(";"+((lines[0]).split('!')[1]).split('$')[0]+ "\n")
scm_output.write(";"+((lines[1]).split('!')[1]).split('$')[0]+ "\n\n")

genes = []
go = []
#loop through lines
for l in lines_annotate:
    gene_symbol =l.split('\t')[2]
    go_id = (l.split('\t')[4])
    gene_name = l.split('\t')[9]
    go_namespace, go_ns = find_gons.find_type(go_id, go_namespace)
    if go_ns:
        member = CMemberLink(CGeneNode(gene_symbol.upper()), go_ns)
        scm_output.write(member.recursive_print() + "\n")
        scm_gene_level.write(member.recursive_print() + "\n")
        go.append(go_id)
    else: 
        print("Unknown namespace: {}".format(go_id))
    if not gene_symbol in genes:
        genes.append(gene_symbol)
        eval_name = CEvaluationLink(CPredicateNode("has_name"), CListLink(CGeneNode(gene_symbol), CConceptNode(gene_name)))
        scm_output.write(eval_name.recursive_print() + "\n")
scm_output.close()
scm_gene_level.close()

if len(go_namespace) > init_namespace:
    with open("raw_data/go-namespace.json", "w") as ns:
        json.dump(go_namespace, ns, indent=2)

script = "https://github.com/MOZI-AI/knowledge-import/GO_Annotation_scm.py"
metadata.update_meta("GO_Annotation:latest", source,script,genes=len(genes), goterms={"go-terms":len(set(go))})
print("Done, check dataset/GO_annotation.scm and gene-level/GO_annotation.scm")