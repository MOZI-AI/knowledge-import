__author__ = 'Abdulrahman Semrie<hsamireh@gmail.com>'

import pandas as pd
import argparse
import os
from atomwrappers import *
from zipfile import ZipFile
from io import BytesIO
from pharmagkb import build_request

gene_dict = {"p53_mutation": "TP53",  "PTEN_pos": "PTEN", "cytokeratin5_pos": "CK5"}
ihc_dict = {"ER_preTrt": "ER_status", "Erbeta_preTrt": "Erbeta_status", "PR_preTrt": "PR_status", "HER2_preTrt": "HER2_status"}
misc_dict = {"preTrt_lymph_node_status": "pre_node_status", "postTrt_lymph_node_status": "post_node_status", "tumor_stage_preTrt": "pre_tumor_stage", "tumor_stage_postTrt": "post_tumor_stage",
             "clinical_AJCC_stage": "cancer_stage"}

def create_ihc_marker_ln(patient, marker, value, fp):
    """
    Creates & writes EvaluationLinks for patient ihc biomarker status
       eg.  (EvaluationLink (stv 1 1) (Predicate "ER_status")
                    (ListLink (Concept "1234") (Concept "positive")))
    """
    if not pd.isna(value):
        marker_name = ihc_dict[marker]
        if value == 0:
            eval_ln = CEvaluationLink(CPredicateNode(marker_name), CListLink(patient, CConceptNode("negative")), stv=CStv(1.0, 1.0))
        else:
            eval_ln = CEvaluationLink(CPredicateNode(marker_name), CListLink(patient, CConceptNode("positive")), stv=CStv(1.0, 1.0))

        fp.write(eval_ln.recursive_print() + "\n")


def create_mutation_ln(patient, mut, value, fp):
    """
    Creates & writes EvaluationLinks for mutation of patient genes
       eg.  (EvaluationLink (stv 1 1) (Predicate "has_mutation")
                    (ListLink (Gene "TP53") (Concept "1234")))
    """
    if not pd.isna(value) and value != 0:
        gene_name = gene_dict[mut]
        eval_ln = CEvaluationLink(CPredicateNode("has_mutation"), CListLink(CGeneNode(gene_name), patient), stv=CStv(value, 1.0))
        fp.write(eval_ln.recursive_print() + "\n")


def create_outcome_ln(patient, outcome, value, fp):
    """
    Creates & writes EvaluationLinks for post treatment outcome of patients
       eg.  (EvaluationLink (stv 1 1) (Predicate "DFS_outcome")
                    (ListLink (Concept "1234") (Concept "positive") ))
    """
    if not pd.isna(value):
        if value == 0:
            eval_ln = CEvaluationLink(CPredicateNode(outcome), CListLink(patient, CConceptNode("negative")), stv=CStv(1.0, 1.0))
        else:
            eval_ln = CEvaluationLink(CPredicateNode(outcome), CListLink(patient, CConceptNode("positive")), stv=CStv(1.0, 1.0))

        fp.write(eval_ln.recursive_print() + "\n")


def create_treatment_ln(patient, treatment, value, drugs_df, fp):
    """
    Creates & writes EvaluationLinks for treatments received by patients
       eg.  (EvaluationLink (stv 1 1) (Predicate "treatment")
                    (ListLink (Concept "1234") (Concept "anthracycline") ))
    """
    if not pd.isna(value) and value != 0:
        if treatment == "surgery_type":
            eval_ln = CEvaluationLink(CPredicateNode("surgery"), CListLink(patient, CConceptNode(treatment)), stv=CStv(1.0, 1.0))
            surgery_type_ln = CEvaluationLink(CPredicateNode("surgery_type"), CListLink(patient, CConceptNode(value)), stv=CStv(1.0, 1.0))
            fp.write(surgery_type_ln.recursive_print() + "\n")
        else:
            try:
                drug_series = drugs_df[drugs_df.Name == treatment].iloc[0]
                if drug_series["Type"] == "Drug":
                    cross_ref = drug_series["Cross-references"]
                    cross_ref = cross_ref.split(",")
                    chebi = [x for x in cross_ref if x.startswith("ChEBI")]
                    if len(chebi) > 0: #look for a chebi id first
                        chebi_id = (chebi[0]).split(":")[-1]
                        eval_ln = CEvaluationLink(CPredicateNode("has_treatment"), CListLink(patient, CChebiNode("ChEBI:" + chebi_id)),
                                              stv=CStv(1.0, 1.0))
                        fp.write(eval_ln.recursive_print() + "\n")
                    else: #look for pubchem
                        pubchem = [x for x in cross_ref if x.startswith("PubChem")]
                        if len(pubchem) > 0:
                            pubchem_id = (pubchem[0]).split(":")[-1]
                            eval_ln = CEvaluationLink(CPredicateNode("has_treatment"),
                                      CListLink(patient, CPubchemNode("PubChem:" + pubchem_id)),
                                      stv=CStv(1.0, 1.0))
                            fp.write(eval_ln.recursive_print() + "\n")

            except IndexError:
                pass

def create_misc_ln(patient, name, value, fp):
    if not pd.isna(value):
        eval_ln = CEvaluationLink(CPredicateNode(misc_dict[name]), CListLink(patient, CConceptNode(value)), stv=CStv(1.0, 1.0))
        fp.write(eval_ln.recursive_print() + "\n")


def import_info(df, drugs_df, fp):
    ihc_marker_cols = ['ER_preTrt','Erbeta_preTrt',
                       'PR_preTrt', 'HER2_preTrt']
    mutation_cols = ['p53_mutation', 'PTEN_pos', 'cytokeratin5_pos']
    outcome_cols = ['OS', 'DFS', 'RFS', 'metastasis', 'pCR']
    treatment_cols = df.loc[:, 'radiotherapyClass':'neoadjuvant_or_adjuvant'].columns.tolist()
    treatment_cols.append('surgery_type')
    misc_status_col = ["preTrt_lymph_node_status", "postTrt_lymph_node_status", "clinical_AJCC_stage", "tumor_stage_preTrt", "tumor_stage_postTrt"]

    for i in range(df.shape[0]):
        patient = CConceptNode(str(df.iloc[i]["patient_ID"]))
        for marker in ihc_marker_cols:
            create_ihc_marker_ln(patient, marker, df.iloc[i][marker], fp)

        for mutation in mutation_cols:
            create_mutation_ln(patient, mutation, df.iloc[i][mutation], fp)

        for outcome in outcome_cols:
            create_outcome_ln(patient, outcome, df.iloc[i][outcome], fp)

        for treatment in treatment_cols:
            create_treatment_ln(patient, treatment, df.iloc[i][treatment], drugs_df, fp)

        for ns in misc_status_col:
            create_misc_ln(patient, ns, df.iloc[i][ns], fp)

def parse_args():
    parser = argparse.ArgumentParser(description="convert clinical trial patient data to atomese")
    parser.add_argument("--table", type=str, default='',
                        help="Path to clinical trial info table in csv format")
    parser.add_argument("--drugs", type=str, default='',
                        help="Path to PharmaGKB drugs zip file")
    parser.add_argument("--output", type=str, default='',
                        help="Path to save the output atomese file")
    return parser.parse_args()

def download_data():
    drugs = "https://s3.pgkb.org/data/drugs.zip"
    clinical_table = "https://raw.githubusercontent.com/singnet/cancer/master/data/curatedBreastData/bcClinicalTable.csv"

    drugs_zip = ZipFile(BytesIO(build_request(drugs).read()))
    clinical_file = BytesIO(build_request(clinical_table).read())
    return clinical_file, drugs_zip

def import_data():
    print("Importing data")
    args = parse_args()
    if args.drugs and args.table:
        clinical_df = pd.read_csv(args.table)
        drug_df = pd.read_csv(ZipFile(args.drugs).open("drugs.tsv"), sep="\t")
    else:
        clinical_file, drug_file = download_data()
        clinical_df = pd.read_csv(clinical_file)
        drug_df = pd.read_csv(drug_file.open("drugs.tsv"), sep="\t")

    if args.output:
        outfile = args.output
    else:
        outfile = "patient_data.scm"

    clinical_df = clinical_df.dropna(axis=1, how='all') #Remove columns that have null for all values

    with open(outfile, "w") as f:
        import_info(clinical_df, drug_df, f)


if __name__ == "__main__":
    import_data()
    print("Done!")