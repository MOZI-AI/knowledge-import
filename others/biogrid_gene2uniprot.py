__author__ = "Hedra"
__email__ = "hedra@singularitynet.io"


# The following script imports the biogrid_genes mapped to their coding uniprots though biogrid_id  

# Requires: uniprot to biogrid_id mapping file uniprot2biogrid.csv and
# Biogrid gene to biogrid_id mapping file gene2biogrid.csv (run biogrid_genes.py to get this file)

import os
from datetime import date
import pandas as pd
import metadata
from atomwrappers import *

script = "https://github.com/MOZI-AI/knowledge-import/biogrid_gene2uniprot.py"

def to_atomese(data):
    print("importing the data")
    df = data.dropna()
    genes = []
    proteins = []
    if not os.path.exists(os.path.join(os.getcwd(), 'dataset')):
        os.makedirs('dataset')
    output_file = "dataset/biogridgene2uniprot_{}.scm".format(str(date.today()))   
    with open(output_file, 'w') as f:
        for i in range(df.shape[0]):
            gene = df.iloc[i]['gene_symbol'].upper().strip()
            biogrid_id = str(df.iloc[i]['biogrid_id'])
            prot = df.iloc[i]['uniprot'].strip()
            if gene and biogrid_id and prot:
                if gene not in genes:
                    genes.append(gene)
                if prot not in proteins:
                    proteins.append(prot)
                expresion = CEvaluationLink(CPredicateNode("expresses"), CListLink(CGeneNode(gene),ProteinNode(prot)))
                bio_prot = CEvaluationLink(CPredicateNode("has_biogridID"), CListLink(ProteinNode(prot), CConceptNode("Bio:"+biogrid_id)))
                bio_gene = CEvaluationLink(CPredicateNode("has_biogridID"), CListLink(CGeneNode(gene), CConceptNode("Bio:"+biogrid_id)))
                f.write("\n".join([x.recursive_print() for x in [expresion, bio_prot, bio_gene]]))
                
    metadata.update_meta("Biogrid-Gene2uniprot:latest", 
        "uniprot2biogrid.csv, gene2biogrid.csv",script,genes=str(len(genes)),prot=len(proteins))
    print("Done, check {}".format(output_file))
 
if __name__ == "__main__":
    '''
        Requires: uniprot to biogrid_id mapping file uniprot2biogrid.csv and
        Biogrid gene symbold to biogrid_id mapping file gene2biogrid.csv 
        (run biogrid_genes2id.py to get gene2biogrid.csv)
    '''
    print("imports the biogrid_genes mapped to their coding uniprots though biogrid_id\n")
    try:
        bio = pd.read_csv("raw_data/gene2biogrid.csv", sep="\t")
        uniprot = pd.read_csv("raw_data/uniprot2biogrid.csv", sep=",")
    except Exception as e:
        print(e)
    for i in range(uniprot.shape[0]):
        biogrid_id = uniprot.iloc[i]['biogrid']
        prot = uniprot.iloc[i]['uniprot']
        # some uniprots has morethan one biogrid ID separated by comma
        for b in biogrid_id.split(","):
                bio.loc[bio['biogrid_id']==int(b), 'uniprot'] = prot
    to_atomese(bio)
