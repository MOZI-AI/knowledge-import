import argparse
import wget
import pandas as pd
import datetime
from atomwrappers import *

# source: ftp://ftp.bgee.org/current/download/calls/expr_calls/Homo_sapiens_expr_simple_development.tsv.gz

def parse_args():
    parser = argparse.ArgumentParser(description='Converts human Genes to Anatomy mapping data into atomese')
    parser.add_argument('--dbsource', type=str, default='',
                        help='the source data file')
    parser.add_argument('--output_file', type=str, default='',
                        help='output file')
    return parser.parse_args()

def gene2anatomy(dbsource, output_file):
    try:
        df = pd.read_csv(dbsource, sep="\t", dtype=str)
    except:
        raise("\nError parsing the data source\n") 

    df = df[["Gene name", "Anatomical entity ID","Anatomical entity name"]]
    df = df[df["Anatomical entity ID"].str.startswith("CL:")]
    df = df.drop_duplicates()
    print(len(df))
    genes = df["Gene name"]
    genes = [g.split(".")[0].upper() for g in genes]
    df["Gene name"] = genes
    cell_types = []
    for i in range(len(df)):
        gene = df.iloc[i]["Gene name"]
        anatomy = df.iloc[i]["Anatomical entity ID"]
        name = df.iloc[i]["Anatomical entity name"]
        output_file.write(CMemberLink(CGeneNode(gene), CelltypeNode(anatomy)).recursive_print() + "\n")
        if not anatomy in cell_types:
            cell_types.append(anatomy)
            output_file.write(CEvaluationLink(CPredicateNode("has_name"), CListLink(CelltypeNode(anatomy), CConceptNode(name))).recursive_print() + "\n")

def main():
    args = parse_args()
    if(not(args.dbsource)):
        dbsource = wget.download("ftp://ftp.bgee.org/current/download/calls/expr_calls/Homo_sapiens_expr_simple_development.tsv.gz", "raw_data/")
    else:
        dbsource = args.dbsource
    if(not(args.output_file)):
        output_file = open("dataset/gene2anatomy_{}.scm".format(str(datetime.date.today())), "w")
    else:
        output_file = open(args.output_file, 'w')
    
    gene2anatomy(dbsource, output_file)

if __name__ == "__main__":
    main()