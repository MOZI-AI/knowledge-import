#!/usr/bin/env python3
# 2017-12-03
# Edited March 2020 
# Script to convert go.obo to atomspace representation in scheme
# Requires: file go.obo from http://www.berkeleybop.org/ontologies/go.obo or http://snapshot.geneontology.org/ontology/go.obo

import re
import wget
import metadata
import os
from datetime import date
import json
from atomwrappers import *
import find_gons

source = "http://snapshot.geneontology.org/ontology/go.obo"
output_data = 'dataset/GO_{}.scm'.format(str(date.today()))

if not os.path.exists('raw_data/go.obo'):
    wget.download(source,"raw_data/")
dataset_file = open('raw_data/go.obo')
lines = dataset_file.readlines()

# store line of number --- "[Terms]" and [Typedef]
line_no = []
print("\nStarted importing\n")
for num, line in enumerate(lines, 1):
    if "[Term]" in line or "[Typedef]" in line:
        line_no.append(num)

line_no.sort()

with open("raw_data/go-namespace.json", "r") as ns:
    go_namespace = json.load(ns)
init_namespace = len(go_namespace)

# open file to write
scm_output = open(output_data, 'w')
goterm = {"biological_process":[],"molecular_function":[],"cellular_component":[]}
i = 0
# partition each line and call functions
while i < len(line_no):
    if i + 1 == len(line_no):
        part = lines[line_no[i] : len(lines)]
    else:
        part = lines[line_no[i] : line_no[i+1] - 1]
    test = [l.partition(':') for l in part]
    k = 0
    is_a = []
    idd =""
    name= ""
    namespace=""
    obsolete =""
    while k < len(test):
        if (test[k][0] == 'is_obsolete'):
            obsolete = (test[k][2].partition('\n')[0]).partition(' ')[2].replace('\\', '\\\\')
        elif (test[k][0] == 'id'):
            idd = (test[k][2].partition('\n')[0]).partition(' ')[2].replace('\\', '\\\\')
        elif (test[k][0] == 'name'):
            name = (test[k][2].partition('\n')[0]).partition(' ')[2].replace('\\', '\\\\')
        elif (test[k][0] == 'namespace'):
            namespace = (test[k][2].partition('\n')[0]).partition(' ')[2].replace('\\', '\\\\')
        elif (test[k][0] == 'is_a'):
            is_a.append(((test[k][2].partition('\n')[0]).partition('!')[0]).partition(' ')[2].replace('\\', '\\\\').strip())
        k = k +1

    if obsolete != 'true' and "GO:" in idd:
        go_namespace, go_term = find_gons.find_type(idd, go_namespace, go_ns=namespace)
        go_name = CEvaluationLink(CPredicateNode("has_name"),CListLink(go_term, CConceptNode(name)))
        scm_output.write(go_name.recursive_print() + "\n")
        if namespace in goterm.keys():
            goterm[namespace].append(idd)
        if len(is_a) != 0:
            isa_len = 0
            while isa_len < len(is_a):
                go_namespace, parent_term = find_gons.find_type(is_a[isa_len], go_namespace)
                if parent_term:
                    inherit = CInheritanceLink(go_term, parent_term) 
                    scm_output.write(inherit.recursive_print() + "\n")
                else:
                    print("Unknown namespace: {}".format(is_a[isa_len]))
                isa_len = isa_len + 1
    i= i + 1

if len(go_namespace) > init_namespace:
    with open("raw_data/go-namespace.json", "w") as ns:
        json.dump(go_namespace, ns, indent=2)
ns = {}
for k in goterm.keys():
    ns[k] = len(set(goterm[k]))
script = "https://github.com/MOZI-AI/knowledge-import/GO_scm.py"
metadata.update_meta("GO Obo:latest", source,script,goterms=ns)
print("Done, check dataset/GO.scm")