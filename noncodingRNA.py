# imports mapping from Gene to RNA transcribes into atomese
# Requires: file GCF_000001405.25_GRCh37.p13_feature_table.txt from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/
import wget
import gzip
import os
import metadata
import pandas as pd
from datetime import date

source = "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_feature_table.txt.gz"

def expres(predicate, node1, node2):
    return ""+'\n(EvaluationLink \n'+'\t(PredicateNode "'+ predicate +'")\n'+'\t(ListLink \n\t\t'+ node1 +'\n\t\t'+ node2 +'))\n'

dataset = "GCF_000001405.25_GRCh37.p13_feature_table.txt.gz"
if not dataset in os.listdir("raw_data/"):
	wget.download(source, "raw_data")

data = pd.read_csv("raw_data/"+dataset, sep="\t",dtype=str)
col = ["product_accession","name","symbol"]
data = data[col].dropna()

# RefSeq accession numbers and molecule types
# https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/
# non_coding_refseq starts with NR_ or XR_

print("Started importing")
rnas = []
genes = []
with open("dataset/noncodingRNA_{}.scm".format(str(date.today())), 'w') as f:
    for i in range(len(data)):
        rna = data.iloc[i]["product_accession"].split(".")[0]
        if rna.split("_")[0] in ["NR", "XR"]:
            gene = data.iloc[i]["symbol"].strip()
            name = data.iloc[i]["name"]
            rnas.append(rna)
            genes.append(gene)
            f.write(expres("transcribed_to", '(GeneNode "{}")'.format(gene), '(MoleculeNode "{}")'.format(rna)))
            f.write(expres("has_name", '(MoleculeNode "{}")'.format(rna), '(ConceptNode "{}")'.format(name)))

version = dataset.split(".")[1]
script = "https://github.com/MOZI-AI/knowledge-import/noncodingRNA.py"

metadata.update_meta("noncodingRNA:{}".format(version),dataset,script,genes=len(set(genes)),ncrna=len(set(rnas)))

print("Done")