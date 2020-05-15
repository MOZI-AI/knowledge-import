__author__ = "Hedra"
__email__ = "hedra@singularitynet.io"

# The following script imports SARS-CoV-2 (COVID-19) and Coronavirus-Related Interactions from thebiogrid.com

# Requires: BIOGRID-CORONAVIRUS-3.5.183.tab3.zip

# from https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-3.5.183/BIOGRID-CORONAVIRUS-3.5.183.tab3.zip

# The version 183 is first release (March 25 2020), It can also be any of the latest versions with the same format

import argparse
import pandas as pd
import wget
import os
import sys
import metadata
from datetime import date
import zipfile
from atomwrappers import *


def checkdisc(diction, key, value):
    try:
        diction.setdefault(key, []).append(value)
    except KeyError:
        return "key error"


def evaLink(node1, node1_type, node2, node2_type, predicate, prefix1="", prefix2="", symmetric=False, stv=""):
    if not (str(node1) in ["-", "nan"] or str(node2) in ["-", "nan"]):
        if symmetric:
            list_type = "SetLink"
        else:
            list_type = "ListLink"
        return ("(EvaluationLink {}\n".format(stv) +
                "\t (PredicateNode \"" + predicate + "\")\n" +
                "\t ({} \n".format(list_type) +
                "\t\t ({}".format(node1_type) + " \"" + prefix1 + str(node1) + "\")\n" +
                "\t\t ({}".format(node2_type) + " \"" + prefix2 + str(node2) + "\")))\n")
    else:
        return ""


def member(node1, node1_type, node2, node2_type, prefix1="", prefix2=""):
    if not (str(node1) in ["-", "nan"] or str(node2) in ["-", "nan"]):
        return ('(MemberLink\n' +
                '\t({} "'.format(node1_type) + prefix1 + str(node1) + '")\n' +
                '\t({} "'.format(node2_type) + prefix2 + str(node2) + '"))\n')
    else:
        return ""


def add_taxonomy(taxonomy_id, node, fp):
    member_ln = CMemberLink(node, CConceptNode("ncbi:" + str(taxonomy_id)))
    fp.write(member_ln.recursive_print() + "\n")


def add_protein_interaction(proteins_lst, prot_node_1, gene_node_1, prot_node_2, gene_node_2, bio_id_1, bio_id_2, fp):
    if not prot_node_1.name in proteins_lst:
        express_ln = CEvaluationLink(CPredicateNode("expresses"), CListLink(gene_node_1, prot_node_1))
        gene_bio = CEvaluationLink(CPredicateNode("has_biogridID"),
                                   CListLink(gene_node_1, CConceptNode("Bio:" + bio_id_1)))
        prot_bio = CEvaluationLink(CPredicateNode("has_biogridID"),
                                   CListLink(prot_node_1, CConceptNode("Bio:" + bio_id_1)))

        fp.write(express_ln.recursive_print() + "\n")
        fp.write(gene_bio.recursive_print() + "\n")
        fp.write(prot_bio.recursive_print() + "\n")

        proteins_lst.append(prot_node_1.name)

    if not prot_node_2.name in proteins_lst:
        express_ln = CEvaluationLink(CPredicateNode("expresses"), CListLink(gene_node_2, prot_node_2))
        gene_bio = CEvaluationLink(CPredicateNode("has_biogridID"),
                                   CListLink(gene_node_2, CConceptNode("Bio:" + bio_id_2)))
        prot_bio = CEvaluationLink(CPredicateNode("has_biogridID"),
                                   CListLink(prot_node_2, CConceptNode("Bio:" + bio_id_2)))

        fp.write(express_ln.recursive_print() + "\n")
        fp.write(gene_bio.recursive_print() + "\n")
        fp.write(prot_bio.recursive_print() + "\n")

        proteins_lst.append(prot_node_2.name)


def process_data(version, file_path):
    if file_path:
        try:
            data = pd.read_csv(file_path, low_memory=False, delimiter='\t')
            version = file_path.split('-')[-1].replace(".tab3.txt", "")
            import_data(data, file_path, version, gene_level=True)
        except Exception as e:
            print(e)
    else:
        if version:
            source = 'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-' + version + '/BIOGRID-CORONAVIRUS-' + version + '.tab3.zip'
        else:
            source = 'https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-CORONAVIRUS-LATEST.tab3.zip'
        try:
            dataset = wget.download(source, "raw_data")
            version = zipfile.ZipFile(dataset).namelist()[0].split('-')[-1].replace(".tab3.txt", "")
            print(version)
            data = pd.read_csv(dataset, low_memory=False, delimiter='\t')
        except:
            print("Error processing biogrid version {0}".format(version))
            raise

        import_data(data, source, version, gene_level=True)


def import_data(data, source, version, gene_level=False, form='tab2'):
    # Set the gene_level to True to get only the GGI without extra entrez and pubmedID info
    print("started importing")
    if not os.path.exists(os.path.join(os.getcwd(), 'dataset')):
        os.makedirs('dataset')

    if gene_level:
        if not os.path.exists(os.path.join(os.getcwd(), 'gene-level')):
            os.makedirs('gene-level')
        g = open('gene-level/COVID-19-biogrid_' + version + "_gene-level_" + str(date.today()) + '.scm', 'w')

    with open('dataset/COVID-19-biogrid_' + version + "_" + str(date.today()) + '.scm', 'w') as f:
        gene_pairs = []
        protein_pairs = []
        entrez = []
        covid_genes = []
        proteins = []
        for i in range(len(data)):
            if not (pd.isnull(data.iloc[i]['Official Symbol Interactor A']) or pd.isnull(
                    data.iloc[i]['Official Symbol Interactor B'])):
                gene1 = str(data.iloc[i]['Official Symbol Interactor A']).upper().strip()
                gene2 = str(data.iloc[i]['Official Symbol Interactor B']).upper().strip()
                prot1 = str(data.iloc[i]['SWISS-PROT Accessions Interactor A']).strip()
                prot2 = str(data.iloc[i]['SWISS-PROT Accessions Interactor B']).strip()
                score = data.iloc[i]['Score']
                entrez1 = str(data.iloc[i]['Entrez Gene Interactor A']).strip()
                entrez2 = str(data.iloc[i]['Entrez Gene Interactor B']).strip()
                taxonomy_id_1 = int(data.iloc[i]['Organism ID Interactor A'])
                taxonomy_id_2 = int(data.iloc[i]['Organism ID Interactor B'])

                gene_node_1 = CGeneNode(gene1)
                gene_node_2 = CGeneNode(gene2)

                prot_node_1 = CMoleculeNode("Uniprot:" + prot1)
                prot_node_2 = CMoleculeNode("Uniprot:" + prot2)

                stv_node = None
                if not str(score) in ["-", "nan"]:
                    stv_node = CStv(1.0, round(float(score), 3))

                if (gene1, gene2) not in gene_pairs:

                    if not gene1 in entrez:
                        entrez_ln_1 = CEvaluationLink(CPredicateNode("has_entrez_id"),
                                                      CListLink(gene_node_1, CConceptNode("entrez:" + entrez1)))
                        f.write(entrez_ln_1.recursive_print() + "\n")
                        entrez.append(gene1)

                    if not gene2 in entrez:
                        eval_ln_2 = CEvaluationLink(CPredicateNode("has_entrez_id"),
                                                    CListLink(gene_node_2, CConceptNode("entrez:" + entrez2)))
                        f.write(eval_ln_2.recursive_print() + "\n")
                        entrez.append(gene2)

                    interacts_ln = CEvaluationLink(CPredicateNode("interacts_with"),
                                                   CSetLink(gene_node_1, gene_node_2), stv=stv_node)
                    f.write(interacts_ln.recursive_print() + "\n")

                    if gene_level:
                        g.write(interacts_ln.recursive_print() + "\n")

                    if taxonomy_id_1 == 2697049:
                        covid_genes.append(gene1)
                        add_taxonomy(taxonomy_id_1, gene_node_1, f)
                        add_taxonomy(taxonomy_id_1, prot_node_1, f)
                        if gene_level:
                            add_taxonomy(taxonomy_id_1, gene_node_1, g)
                    if taxonomy_id_2 == 2697049:
                        covid_genes.append(gene2)

                        add_taxonomy(taxonomy_id_2, gene_node_2, f)
                        add_taxonomy(taxonomy_id_2, prot_node_2, f)

                        if gene_level:
                            add_taxonomy(taxonomy_id_2, gene_node_2, g)

                    gene_pairs.append((gene1, gene2))

                if (prot1, prot2) not in protein_pairs:
                    interacts_ln = CEvaluationLink(CPredicateNode("interacts_with"),
                                                   CSetLink(prot_node_1, prot_node_2), stv=stv_node)

                    f.write(interacts_ln.recursive_print() + "\n")

                    bio_1 = str(data.iloc[i]['BioGRID ID Interactor A']).strip()
                    bio_2 = str(data.iloc[i]['BioGRID ID Interactor B']).strip()
                    add_protein_interaction(proteins, prot_node_1, gene_node_1, prot_node_2, gene_node_2, bio_1, bio_2,
                                            f)

                    protein_pairs.append((prot1, prot2))

        f.write(evaLink("2697049", "ConceptNode", "SARS-CoV-2", "ConceptNode", "has_name", prefix1="ncbi:"))
        g.write(evaLink("2697049", "ConceptNode", "SARS-CoV-2", "ConceptNode", "has_name", prefix1="ncbi:"))
    gene_pairs = set((a, b) if a <= b else (b, a) for a, b in gene_pairs)
    number_of_interactions = len(gene_pairs)
    script = "https://github.com/MOZI-AI/knowledge-import/coronavirus_biogrid.py"
    metadata.update_meta("Coronavirus Biogrid:" + version, source, script, genes=str(len(set(entrez))),
                         prot=len(set(proteins)), interactions=str(number_of_interactions))
    print("Done, check " + 'dataset/COVID-19-biogrid_' + version + "_" + str(date.today()) + '.scm')
    with open("Covid19-genes", "w") as co:
        co.write("\n".join(list(set(covid_genes))))


def parse_args():
    parser = argparse.ArgumentParser(description='convert biogrid db to atomese')
    parser.add_argument('--path', type=str, default='',
                        help='process local file in biogrid format')
    parser.add_argument('--download', action='store_true', default=True,
                        help='download and process db from biogrid')
    parser.add_argument('--version', type=str, default='',
                        help='version to download(by default lastest is used)')
    return parser.parse_args()


if __name__ == "__main__":
    """
  usage:
  run the script with the path to the source data (if downloaded)
        python coronavirus_biogrid.py --path=path/to/the/source_data 
  Or run the script and specify a version number you wanted or just hit enter (to get the latest)
  """
    arguments = parse_args()
    process_data(arguments.version, arguments.path)
