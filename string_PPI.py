# Uniprot to string mapping https://string-db.org/mapping_files/uniprot/human.uniprot_2_string.2018.tsv.gz
# String PPI dataset https://stringdb-static.org/download/protein.actions.v11.0/9606.protein.actions.v11.0.txt.gz
# Columns definition http://www.string-db.org/help/faq/#what-does-the-columns-in-proteinsactions-file-mean

import pandas as pd
import wget
import os
import sys
import metadata
import datetime
from atomwrappers import *

source = "https://stringdb-static.org/download/protein.actions.v11.0/9606.protein.actions.v11.0.txt.gz"
mapping = "https://string-db.org/mapping_files/uniprot/human.uniprot_2_string.2018.tsv.gz"

def import_string():
    print("started at " + str(datetime.datetime.now()))
    if not os.path.exists('raw_data/9606.protein.actions.v11.0.txt.gz'):
        wget.download(source,"raw_data/")
    if not os.path.exists('raw_data/human.uniprot_2_string.2018.tsv.gz'):
        wget.download(mapping,"raw_data/")
    
    df_data = pd.read_csv("raw_data/9606.protein.actions.v11.0.txt.gz", dtype=str, sep="\t")
    df_data_symmetric = df_data[df_data['is_directional'] == "f"]
    df_data_asymmetric = df_data[df_data['is_directional'] == "t"]
    df_mapping = pd.read_csv("raw_data/human.uniprot_2_string.2018.tsv.gz", dtype=str, sep="\t", names=["code", "uniprot", "ensembl","num1","num2"])
   
    # create a mapping dictionary
    mapping_dict = {} 
    for e in df_mapping["ensembl"]:
        if not e in mapping_dict.keys():
            mapping_dict[e] = df_mapping[df_mapping["ensembl"] == e]["uniprot"].values[0]
    print("Done with the Dict, importing into atomese")
    print(len(df_data))
    notmapped = []

    """
        If the directionality of the interaction is true and a is acting, use ListLink and keep the order. Otherwise use SetLink
        * is_directional - describes if the diractionality of the particular interaction is known.
        * a_is_acting - the directionality of the action if applicable ('t' gives that item_id_a is acting upon item_id_b)
        Example:
        item_id_a   item_id_b   mode    is_directional  a_is_acting
        ENSP00000000233	ENSP00000216366 reaction    f   f   
        <=> EvaluationLink 
                PredicateNode "reaction"
                SetLink ENSP00000000233 ENSP00000216366

        ENSP00000000233	ENSP00000216366 reaction    t   f
        <=> EvaluationLink 
                PredicateNode "reaction"
                ListLink ENSP00000216366 ENSP00000000233
        ENSP00000000233	ENSP00000216366	reaction    t   t
        <=> EvaluationLink 
                PredicateNode "reaction"
                ListLink ENSP00000000233 ENSP00000216366
    
        Keep symmetric relations and ignore if the same relation happens to be asymmetric
    """
    symmetric = {}
    if not os.path.exists(os.path.join(os.getcwd(), 'string_dataset')):
        os.makedirs('string_dataset')
    with open("string_dataset/string_ppi_{}.scm".format(str(datetime.date.today())), "w") as f, open('string_dataset/string_ggi_{}.scm'.format(str(datetime.date.today())), 'w') as g:
        for i in range(len(df_data_symmetric)):
            try:
                prot1 = df_data_symmetric.iloc[i]['item_id_a']
                prot2 = df_data_symmetric.iloc[i]['item_id_b']
                mode = df_data_symmetric.iloc[i]['mode']
                score = int(df_data_symmetric.iloc[i]['score'])               

                if prot1 in mapping_dict.keys() and prot2 in mapping_dict.keys():
                    prot1 = mapping_dict[prot1]
                    prot2 = mapping_dict[prot2]
                else:
                    if not prot1 in mapping_dict.keys():
                        notmapped.append(prot1)
                    else:
                        notmapped.append(prot2)
                    continue
                
                protein1 = ProteinNode(prot1.split("|")[0])
                gene1 = CGeneNode(prot1.split("|")[1].split("_")[0].upper()) 
                protein2 = ProteinNode(prot2.split("|")[0])
                gene2 = CGeneNode(prot2.split("|")[1].split("_")[0].upper())
                stv = "(stv {} {})".format(1.0, score/1000)

                f.write(str(CEvaluationLink(CPredicateNode(mode), CSetLink(protein1, protein2), stv=stv)) + "\n")
                f.write(str(CEvaluationLink(CPredicateNode(mode), CSetLink(gene1, gene2), stv=stv)) + "\n")
                symmetric[gene1 + gene2] = mode

            except Exception as e:
                print(e)
        for i in range(len(df_data_asymmetric)):
            try:
                prot1 = df_data_asymmetric.iloc[i]['item_id_a']
                prot2 = df_data_asymmetric.iloc[i]['item_id_b']
                mode = df_data_asymmetric.iloc[i]['mode']
                a_is_acting = df_data_asymmetric.iloc[i]['a_is_acting'] 
                score = int(df_data_asymmetric.iloc[i]['score'])               

                if prot1 in mapping_dict.keys() and prot2 in mapping_dict.keys():
                    prot1 = mapping_dict[prot1]
                    prot2 = mapping_dict[prot2]
                else:
                    if not prot1 in mapping_dict.keys():
                        notmapped.append(prot1)
                    else:
                        notmapped.append(prot2)
                    continue
                
                protein1 = ProteinNode(prot1.split("|")[0])
                gene1 = CGeneNode(prot1.split("|")[1].split("_")[0].upper()) 
                protein2 = ProteinNode(prot2.split("|")[0])
                gene2 = CGeneNode(prot2.split("|")[1].split("_")[0].upper())
                stv = "(stv {} {})".format(1.0, score/1000)

                if not (gene1+gene2 in symmetric.keys() and symmetric[gene1+gene2] == mode):
                    if a_is_acting is "t": 
                        f.write(str(CEvaluationLink(CPredicateNode(mode), CListLink(protein1, protein2), stv=stv)) + "\n")
                        g.write(str(CEvaluationLink(CPredicateNode(mode), CListLink(gene1, gene2), stv=stv)) + "\n")
                    else:
                        f.write(str(CEvaluationLink(CPredicateNode(mode), CListLink(protein2, protein1), stv=stv)) + "\n")
                        g.write(str(CEvaluationLink(CPredicateNode(mode), CListLink(gene2, gene1), stv=stv)) + "\n")

            except Exception as e:
                print(e)
        
        print("Done " + str(datetime.datetime.now()))
        with open("string_dataset/notmapped_ensembles.txt", "w") as n:
            n.write("\n".join(set(notmapped)))

if __name__ == "__main__":
    import_string()