# Script to convert go-plus csv to atomese representation in scheme
# Requires: file go-plus from https://bioportal.bioontology.org/ontologies/GO-PLUS
# You can also get the file from https://gitlab.com/opencog-bio/pln_mozi/blob/master/raw_data/GO-PLUS.csv.gz

import re
import wget
import metadata
import os
import pandas as pd
from datetime import date
import json
from atomwrappers import *
import find_gons 

source = "https://bioportal.bioontology.org/ontologies/GO-PLUS"
source_csv = "https://gitlab.com/opencog-bio/pln_mozi/blob/master/raw_data/GO-PLUS.csv.gz"
source_csv_latest = "http://data.bioontology.org/ontologies/GO-PLUS/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb&download_format=csv"

if not os.path.exists("raw_data/GO-PLUS.csv.gz"):
    dataset = wget.download(source_csv_latest, "raw_data")
df = pd.read_csv("raw_data/GO-PLUS.csv.gz", dtype=str)

with open("raw_data/go-namespace.json", "r") as ns:
    go_namespace = json.load(ns)
init_namespace = len(go_namespace)

def get_term(class_id):
    if str(class_id) == "nan":
        term = False
    else:
        term = class_id.split("/")[-1]
        if term.startswith("GO") or term.startswith("CL") or term.startswith("UBERON") or term.startswith("CHEBI"):
            term = term.replace("_",":")
            term = term.replace("CHEBI", "ChEBI")
    return term

def get_type(term, parent=[]):
    if "UBERON" in term:
        term = UberonNode(term)
    elif "CL" in term:
        term = CelltypeNode(term)
    elif "GO:" in term:
        term = find_namespace(term)
    elif "ChEBI" in term:
        if term in parent:
            term = ChebiOntology(term)
        else:
            term = ChebiNode(term)

    else:
        term = CConceptNode(term)
    return term

def find_namespace(go_term):
    global go_namespace
    go_namespace, go_term = find_gons.find_type(go_term, go_namespace)
    return go_term

# Parent CHEBI's should be a ConceptNode, not a MoleculeNode
parents = df["Parents"]
parent_chebis = []
for i in [i.split("|") for i in parents if str(i) != "nan"]:
    for c in i:
        term = get_term(c)
        if "ChEBI" in term:
            parent_chebis.append(term)
all_col = df.columns            
go_columns = ["negatively regulated by","negatively regulates", "positively regulated by", "positively regulates", "regulated by", "regulates", "has part", "part of"]
uberon_columns = open("raw_data/uberon_columns.txt", "r").read().splitlines()
uberon_columns = [c for c in uberon_columns if c in all_col]
cl_columns = ['has part', 'Parents', 'has role']

if not os.path.exists("dataset/go-plus"):
    os.mkdir("dataset/go-plus/")

go = open("dataset/go-plus/Go-Plus-GO_{}.scm".format(str(date.today())),"w")
uberon = open("dataset/go-plus/Go-Plus-UBERON_{}.scm".format(str(date.today())),"w")
cl = open("dataset/go-plus/Go-Plus-CL_{}.scm".format(str(date.today())),"w")
chebi = open("dataset/go-plus/Go-Plus-CHEBI_{}.scm".format(str(date.today())),"w")

go_with_def = open("dataset/go-plus/Go-Plus-GO_with_definition_{}.scm".format(str(date.today())),"w")
uberon_with_def = open("dataset/go-plus/Go-Plus-UBERON_with_definition_{}.scm".format(str(date.today())),"w")
cl_with_def = open("dataset/go-plus/Go-Plus-CL_with_definition_{}.scm".format(str(date.today())),"w")
chebi_with_def = open("dataset/go-plus/Go-Plus-CHEBI_with_definition_{}.scm".format(str(date.today())),"w")

meta = {}

print("Started importing")
for i in range(len(df)):
    try:
        term = get_term(df.iloc[i]["Class ID"])
        obsolete = df.iloc[i]["Obsolete"]
        definition = CConceptNode(str(df.iloc[i]["definition"]))
        if term and obsolete != "true" and "GO:" in term:
            term = find_namespace(term) 
            if term:
                name = CConceptNode(get_term(df.iloc[i]["Preferred Label"]))
                eva_name = CEvaluationLink(CPredicateNode("has_name"), CListLink(term,name))
                eva_defn = CEvaluationLink(CPredicateNode("has_definition"), CListLink(term,definition))
                go.write(eva_name.recursive_print() + "\n")
                go_with_def.write(eva_defn.recursive_print() + "\n" + eva_name.recursive_print() + "\n")
                for col in go_columns:
                    """
                    positive/negatively regulated by is inverse of positive/negatively regulates
                    has part is inverse of part of, keep the predicate the same with reverse order
                    """
                    if col.endswith("regulated by"):
                        term2 = find_namespace(get_term(df.iloc[i][col]))
                        if term2:
                            col_pred = col.replace("regulated by", "regulates")
                            col_pred = "GO_{}".format(col_pred.replace(" ", "_"))
                            eva = CEvaluationLink(CPredicateNode(col_pred), CListLink(term, term2))
                            go.write(eva.recursive_print() + "\n")
                            go_with_def.write(eva.recursive_print() + "\n")
                    elif col == "part of":
                        term2 = find_namespace(get_term(df.iloc[i][col]))
                        if term2:
                            col_pred = "GO_has_part"
                            eva = CEvaluationLink(CPredicateNode(col_pred), CListLink(term, term2)) 
                            go.write(eva.recursive_print() + "\n")
                            go_with_def.write(eva.recursive_print() + "\n")
                    else:
                        term2 = find_namespace(get_term(df.iloc[i][col]))
                        if term2:
                            col_pred = "GO_{}".format(col.replace(" ", "_"))
                            eva = CEvaluationLink(CPredicateNode(col_pred), CListLink(term, term2))                 
                            go.write(eva.recursive_print() + "\n")
                            go_with_def.write(eva.recursive_print() + "\n")

        elif term and obsolete != "true" and "UBERON" in term:
            term = UberonNode(term)
            name = get_term(df.iloc[i]["Preferred Label"])
            eva_name = CEvaluationLink(CPredicateNode("has_name"), CListLink(term, CConceptNode(name)))
            eva_defn = CEvaluationLink(CPredicateNode("has_definition"), CListLink(term, definition))
            uberon.write(eva_name.recursive_print() + "\n")
            uberon_with_def.write(eva_name.recursive_print() + "\n"+ eva_defn.recursive_print() + "\n")
            for col in uberon_columns:
                term_2 = get_term(df.iloc[i][col])
                if term_2:
                    term_2 = get_type(term_2)
                    col_pred = "UBERON_{}".format(col.replace(" ", "_"))
                    eva = CEvaluationLink(CPredicateNode(col_pred), CListLink(term, term_2))
                    uberon.write(eva.recursive_print() + "\n")

        elif term and obsolete != "true" and "CL" in term or "ChEBI" in term:
            if "CL" in term:
                file_name = cl
                file_name_with_def = cl_with_def
            else:
                file_name = chebi
                file_name_with_def = chebi_with_def

            term = get_type(term, parent=parent_chebis)
            name = get_term(df.iloc[i]["Preferred Label"])
            eva_name = CEvaluationLink(CPredicateNode("has_name"), CListLink(term, CConceptNode(name)))
            eva_defn = CEvaluationLink(CPredicateNode("has_definition"), CListLink(term, definition))
            file_name.write(eva_name.recursive_print() + "\n")
            file_name_with_def.write(eva_name.recursive_print() + "\n" + eva_defn.recursive_print() + "\n")
            for col in cl_columns:
                if col == "Parents":
                    parents = df.iloc[i][col]
                    if str(parents) != "nan":
                        for p in parents.split("|"):
                            term2 = get_type(get_term(p), parent=parent_chebis)
                            inherit = CInheritanceLink(term, term2) 
                            file_name.write(inherit.recursive_print() + "\n")
                            file_name_with_def.write(inherit.recursive_print() + "\n")
                else:
                    term2 = get_term(df.iloc[i][col])
                    if term2: 
                        term2 = get_type(term2, parent=parent_chebis)
                        eva_link = CEvaluationLink(CPredicateNode(col.replace(" ", "_")), CListLink(term, term2))
                        file_name.write(eva_link.recursive_print() + "\n")
                        file_name_with_def.write(eva_link.recursive_print() + "\n")
    except Exception as e:
        print("Exception {} at row {} ".format(e, i))
        continue
print("Done")
            