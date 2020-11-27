__author__ = 'Abdulrahman Semrie<hsamireh@gmail.com>'

import pandas as pd
import argparse
import os
from atomwrappers import *
from zipfile import ZipFile
from io import BytesIO
from pharmagkb import build_request

def create_gene_expr_ln(patient, gene, value, fp, overexpr):
    if not pd.isna(value):
        if overexpr:
            gene_exec_ln = CLazyExecutionOutputLink(CSchemaNode("make-overexpression-schema-for-gene"), CGeneNode(gene))
            exec_ln = CExecutionLink(gene_exec_ln, patient,
                                     CNumberNode(str(value)))
        else:
            gene_exec_ln = CLazyExecutionOutputLink(CSchemaNode("make-underexpression-schema-for-gene"), CGeneNode(gene))
            exec_ln = CExecutionLink(gene_exec_ln, patient,
                                     CNumberNode(str(value)))

        fp.write(exec_ln.recursive_print() + "\n")


def expr_to_int(gene, val, expr_dict, overexpr):
    if overexpr:
        if val <= expr_dict[gene]: #For the overexpression,  we consider underexpression to the value 0
            return 0
        else:
            return val - expr_dict[gene]

    else:
        if val >= expr_dict[gene]: #For underexpression, we consider overexpression to have value 0
            return 0
        else:
            return expr_dict[gene] - val

def create_quantitative_predicate_ln(gene, overexpr, fp):
    if overexpr:
        quant_ln = CQuantitativePredicateLink(
                CLazyExecutionOutputLink(CSchemaNode("make-overexpression-schema-for-gene"), CGeneNode(gene)),
                CLazyExecutionOutputLink(CSchemaNode("make-overexpression-predicate-for-gene"), CGeneNode(gene)))
    else:
        quant_ln = CQuantitativePredicateLink(
            CLazyExecutionOutputLink(CSchemaNode("make-underexpression-schema-for-gene"), CGeneNode(gene)),
            CLazyExecutionOutputLink(CSchemaNode("make-underexpression-predicate-for-gene"), CGeneNode(gene)))

    fp.write(quant_ln.recursive_print() + "\n")

def create_member_ln(gene, fp):
    member_ln = CMemberLink(CGeneNode(gene), CConceptNode("profiled-genes"))
    fp.write(member_ln.recursive_print() + "\n")

def import_gene_expr(patient_df, overexpr_fp, underexpr_fp):
    df = patient_df.dropna(axis=1, how='all')
    mean_dict = {}
    for col in df.columns.values:
        if col == "patient_ID": continue
        mean_dict[str(col)] = df[col].median()

    for i in range(df.shape[0]):
        patient = CPatientNode(str(int(df.iloc[i]["patient_ID"])))
        for k in mean_dict:
            create_member_ln(k, overexpr_fp)
            create_member_ln(k, underexpr_fp)
            create_quantitative_predicate_ln(k, True, overexpr_fp)
            create_quantitative_predicate_ln(k, False, underexpr_fp)
            create_gene_expr_ln(patient, k, expr_to_int(k, df.iloc[i][k], mean_dict, True), overexpr_fp, True)
            create_gene_expr_ln(patient, k, expr_to_int(k, df.iloc[i][k], mean_dict, False), underexpr_fp, False)

def parse_args():
    parser = argparse.ArgumentParser(description="convert clinical trial patient data to atomese")
    parser.add_argument("--table", type=str, default='',
                        help="Path to clinical trial info table in csv format")
    parser.add_argument("--path", type=str, default='',
                        help="Path to save the output atomese file")
    return parser.parse_args()

def download_data():
    merged_file_ln = "https://snet-bio-data.s3-us-west-2.amazonaws.com/example15bmc/merged-combat15.csv.xz"
    return build_request(merged_file_ln).read()

def main():
    print("Importing data")
    args = parse_args()
    if args.table:
        patient_df = pd.read_csv(args.table)
    else:
        merged_file = download_data()
        patient_df = pd.read_csv(merged_file, compression="xz")

    if args.path:
        save_path = os.path.abspath(args.path)
    else:
        save_path = os.getcwd()

    overexpr_file = os.path.join(save_path, "patient_gene_over_expr.scm")
    underexpr_file = os.path.join(save_path, "patient_gene_under_expr.scm")

    with open(overexpr_file, "w") as f1:
        with open(underexpr_file, "w") as f2:
            import_gene_expr(patient_df, f1, f2)

if __name__ == "__main__":
    main()
    print("Done")