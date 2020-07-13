__author__ = "Hedra"
__email__ = "hedra@singularitynet.io"

# The following script imports Biogrid interaction of Human genes from thebiogrid.com

# Requires: BIOGRID-ORGANISM-Homo_sapiens-3.5.169.tab2.txt

# from https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.169/BIOGRID-ORGANISM-3.5.169.tab2.zip

# Or any of the latest version with the same type

import argparse
import pandas as pd
from urllib.request import urlopen
from zipfile import ZipFile
from io import BytesIO
import os
import collections
import metadata
from datetime import date
from atomwrappers import *



def import_data_from_web(version, form='tab2'):
    if form not in ('tab2', 'tab3'):
        raise RuntimeError("format {0} is not supported".format(form))
    if version:
        source = 'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-' + version + '/BIOGRID-ORGANISM-' + version + '.{0}.zip'.format(
            form)
    else:
        source = 'https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.{0}.zip'.format(
            form)
    print("Started downloading the source:\n {0}".format(source))
    try:
        thefile = urlopen(source)
        extracted_files = ZipFile(BytesIO(thefile.read()))
        if version:
            dataset = 'BIOGRID-ORGANISM-Homo_sapiens-' + version + '.{0}.txt'.format(form)
        else:
            dataset = [i for i in extracted_files.namelist() if "BIOGRID-ORGANISM-Homo_sapiens" in i][0]
            version = dataset.split('-')[-1].replace(".{0}.txt".format(form), "")
        data = pd.read_csv(extracted_files.open(dataset), low_memory=False, delimiter='\t')
    except:
        print("Error processing biogrid version {0}".format(version))
        raise
    import_data(data, source, version, gene_level=True, form=form)
    data.to_csv("raw_data/" + dataset, sep='\t', index=False)


def import_local_data(file, form='tab2'):
    path = os.path.abspath(file)
    if os.path.isfile(path):
        try:
            data = pd.read_csv(path, low_memory=False, delimiter='\t')
            version = file.split('-')[-1].replace(".{0}.txt".format(form), "")
            import_data(data, path, version, gene_level=False, form=form)
        except Exception as e:
            print(e)


def import_data(data, source, version, gene_level=False, form='tab2'):
    # Set the gene_level to True to get only the GGI without extra entrez and pubmedID info
    if form == 'tab2':
        pubsource = 'Pubmed ID'
    elif form == 'tab3':
        pubsource = 'Publication Source'
    else:
        raise RuntimeError("format {0} is not yet supported".format(form))
    data = data[['Entrez Gene Interactor A', 'Entrez Gene Interactor B', 'Official Symbol Interactor A',
                 'Official Symbol Interactor B', pubsource]]
    print("started importing")
    dataset_path = os.path.join(os.getcwd(), 'dataset')
    if not os.path.exists(dataset_path):
        os.makedirs('dataset')

    if gene_level:
        if not os.path.exists(os.path.join(os.getcwd(), 'gene-level')):
            os.makedirs('gene-level')
        g = open('gene-level/biogrid_gene_gene_' + version + "_gene-level_" + str(date.today()) + '.scm', 'w')

    biogrid_path = os.path.join(dataset_path, 'biogrid_gene_gene_' + version + '_' + str(date.today()) + '.scm')
    with open(biogrid_path, 'w') as f:
        pairs = collections.defaultdict(list)
        entrez = []
        for i in range(len(data)):
            if not (pd.isnull(data.iloc[i]['Official Symbol Interactor A']) or pd.isnull(
                    data.iloc[i]['Official Symbol Interactor B'])):
                node1 = str(data.iloc[i]['Official Symbol Interactor A']).upper().strip()
                node2 = str(data.iloc[i]['Official Symbol Interactor B']).upper().strip()
                gene_1 = CGeneNode(node1)
                gene_2 = CGeneNode(node2)
                pubmed = data.iloc[i][pubsource]
                if node1 > node2:
                    interactors = node1 + ':' + node2
                else:
                    interactors = node2 + ':' + node1

                pubmed_link = 'https://www.ncbi.nlm.nih.gov/pubmed/?term=' + str(pubmed)
                pairs[interactors].append(CConceptNode(pubmed_link))

                if not node1 in entrez:
                    entrez_1 = CConceptNode("entrez:" + str(data.iloc[i]['Entrez Gene Interactor A']))
                    eval_ln = CEvaluationLink(CPredicateNode("has_entrez_id"), CListLink(gene_1, entrez_1))
                    f.write(eval_ln.recursive_print() + "\n")
                    entrez.append(node1)

                if not node2 in entrez:
                    entrez_2 = CConceptNode("entrez:" + str(data.iloc[i]['Entrez Gene Interactor B']))
                    eval_ln = CEvaluationLink(CPredicateNode("has_entrez_id"), CListLink(gene_2, entrez_2))
                    f.write(eval_ln.recursive_print() + "\n")
                    entrez.append(node2)

        number_of_genes = []
        for p in pairs.keys():
            gene_1 = CGeneNode(p.split(':')[0])
            gene_2 = CGeneNode(p.split(':')[1])
            interacts_ln = CEvaluationLink(CPredicateNode("interacts_with"), CSetLink(gene_1, gene_2))
            eval_ln = CEvaluationLink(CPredicateNode("has_pubmedID"), CListLink(
                interacts_ln,
                CListLink(*pairs[p])
            ))
            f.write(eval_ln.recursive_print() + "\n")
            if gene_level:
                g.write(eval_ln.recursive_print() + "\n")

            number_of_genes.append(gene_1.name)
            number_of_genes.append(gene_2.name)

    number_of_interactions = len(set(pairs.keys()))
    script = "https://github.com/MOZI-AI/knowledge-import/biogrid.py"
    metadata.update_meta("Biogrid:" + version, source, script, genes=str(len(set(number_of_genes))),
                         interactions=str(number_of_interactions))
    print("Done, check " + 'dataset/biogrid_gene_gene_' + version + "_" + str(date.today()) + '.scm')


def parse_args():
    parser = argparse.ArgumentParser(description='convert biogrid db to atomese')
    parser.add_argument('--format', type=str, default='tab2',
                        help='biogrid format of the data: tab2 or tab3')
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
        python biogrid.py --path /path/to/the/source_data
  Or run the script and specify a version number you wanted or just hit enter (to get the latest)
  """
    arguments = parse_args()
    form = arguments.format
    version = arguments.version
    if arguments.path:
        dataset_path = arguments.path
        if form:
            import_local_data(dataset_path, form=form)
        else:
            import_local_data(dataset_path)
    else:
        print("Imports interaction between genes (Homo_sapiens) from thebiogrid.com")
        import_data_from_web(version, form=form)
