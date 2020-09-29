import argparse
import wget
import pandas as pd
import datetime
from atomwrappers import *

# source: https://www.proteinatlas.org/download/normal_tissue.tsv.zip

def parse_args():
    parser = argparse.ArgumentParser(description='Converts proteinatlas normal tissue data into atomese')
    parser.add_argument('--dbsource', type=str, default='',
                        help='the source data file')
    parser.add_argument('--clmapping', type=str, default='',
                        help='the source data file')
    parser.add_argument('--output_file', type=str, default='',
                        help='output file')
    return parser.parse_args()

def preprocess(mapping_file):
    mapping_df = pd.read_csv(mapping_file, sep=",", dtype=str)
    mapping_dict = {}
    for cl in mapping_df["cl_id"]:
        dic = {}
        sub_df = mapping_df[mapping_df["cl_id"] == cl]
        dic['tissue'] = sub_df["Tissue"].values[0]
        dic['celltype'] = sub_df["Cell type"].values[0]
        dic['clname'] = sub_df["cl_name"].values[0]
        mapping_dict[cl] = dic
    return mapping_dict

def find_type(anatomy):
    if "UBERON" in anatomy:
        return UberonNode(anatomy)
    else:
        return CelltypeNode(anatomy)

def gene2anatomy(dbsource, mapping_file, output_file):
    try:
        df = pd.read_csv(dbsource, sep="\t", dtype=str)
        # df.fillna("N/A")
        # mapping_dict = preprocess(mapping_file)
        mapping_df = pd.read_csv(mapping_file, sep=",", dtype=str)
        # mapping_df.fillna("N/A")
    except:
        raise("\nError parsing the data source\n") 

    df = df[["Gene name", "Tissue","Cell type"]]
    df = df.drop_duplicates()
    mapping_df = mapping_df.drop_duplicates()
    merged_df = pd.merge(df, mapping_df, how="left", on=["Tissue","Cell type"])
    merged_df = merged_df.fillna("N/A")
    print("Input data length: {}, Merged Data length {}".format(len(df),len(merged_df)))
    genes = df["Gene name"]
    genes = [g.split(".")[0].upper() for g in genes]
    df["Gene name"] = genes
    cell_types = []

    for r,v in merged_df.iterrows():
        gene = v["Gene name"]
        cl_id = v["cl_id"]
        cl_name = v["cl_name"]
        if not cl_id == "N/A":
            member = CMemberLink(CGeneNode(gene), find_type(cl_id))
            output_file.write(member.recursive_print() + "\n") 
            if not cl_id in cell_types:
                cell_types.append(cl_id)           
                name = CEvaluationLink(CPredicateNode("has_name"), CListLink(find_type(cl_id), CConceptNode(cl_name)))
                output_file.write(name.recursive_print() + "\n")

def main():
    args = parse_args()
    if(not(args.dbsource)):
        dbsource = wget.download("https://www.proteinatlas.org/download/normal_tissue.tsv.zip", "raw_data/")
    else:
        dbsource = args.dbsource
    if(not(args.clmapping)):
        mapping = "raw_data/PA_normal_tissue2cl.csv"
    else:
        mapping = args.dbsource
    if(not(args.output_file)):
        output_file = open("dataset/PA_gene2celltype_{}.scm".format(str(datetime.date.today())), "w")
    else:
        output_file = open(args.output_file, 'w')
    
    gene2anatomy(dbsource, mapping, output_file)

if __name__ == "__main__":
    main()