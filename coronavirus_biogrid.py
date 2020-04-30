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

def checkdisc(diction, key, value):
  try:
    diction.setdefault(key,[]).append(value)
  except KeyError:
    return "key error"

def evaLink(node1, node1_type, node2, node2_type, predicate, prefix1="", prefix2="",symmetric=False, stv=""):
    if not (str(node1) in ["-", "nan"] or str(node2) in ["-", "nan"]):
        if symmetric:
            list_type = "SetLink"
        else:
            list_type = "ListLink"
        return ("(EvaluationLink {}\n".format(stv) +
            "\t (PredicateNode \""+ predicate + "\")\n" +
            "\t ({} \n".format(list_type) +
            "\t\t ({}".format(node1_type) + " \"" + prefix1 + str(node1) + "\")\n" +
            "\t\t ({}".format(node2_type) + " \"" + prefix2 + str(node2) + "\")))\n" )
    else:
        return ""

def member(node1, node1_type, node2, node2_type, prefix1="", prefix2=""):
    if not (str(node1) in ["-", "nan"] or str(node2) in ["-", "nan"]):
        return ('(MemberLink\n' + 
                '\t({} "'.format(node1_type) + prefix1 + str(node1) + '")\n'+
                '\t({} "'.format(node2_type) + prefix2 + str(node2) + '"))\n')
    else:
        return ""

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
            source = 'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-'+ version +'/BIOGRID-CORONAVIRUS-'+ version +'.tab3.zip'
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
        g = open('gene-level/COVID-19-biogrid_'+version+"_gene-level_"+str(date.today())+'.scm','w')

    with open('dataset/COVID-19-biogrid_'+version+"_"+str(date.today())+'.scm','w') as f:
        pairs = []
        entrez = []
        genes = []
        covid_genes = []
        proteins = []
        for i in range(len(data)):
            if not (pd.isnull(data.iloc[i]['Official Symbol Interactor A']) or pd.isnull(data.iloc[i]['Official Symbol Interactor B'])):
                gene1 = str(data.iloc[i]['Official Symbol Interactor A']).upper().strip()
                gene2 = str(data.iloc[i]['Official Symbol Interactor B']).upper().strip()
                prot1 = str(data.iloc[i]['SWISS-PROT Accessions Interactor A']).strip()
                prot2 = str(data.iloc[i]['SWISS-PROT Accessions Interactor B']).strip()
                score = data.iloc[i]['Score']
                entrez1 = str(data.iloc[i]['Entrez Gene Interactor A']).strip()
                entrez2 = str(data.iloc[i]['Entrez Gene Interactor B']).strip()
                stv = ""
                if not str(score) in ["-", "nan"]:
                    stv = '(stv 1.0 {})'.format(round(float(score),3))
                taxonomy_id_1 = int(data.iloc[i]['Organism Interactor A'])
                taxonomy_id_2 = int(data.iloc[i]['Organism Interactor B'])

                if (gene1, gene2) not in pairs or (gene2, gene1) not in genes:
                    if not gene1 in entrez:
                        f.write(evaLink(gene1, "GeneNode", entrez1,"ConceptNode", "has_entrez_id",prefix2="entrez:"))
                        entrez.append(gene1)

                    if not gene2 in entrez:
                        f.write(evaLink(gene2, "GeneNode", entrez2,"ConceptNode", "has_entrez_id",prefix2="entrez:"))
                        entrez.append(gene2)

                    f.write(evaLink(gene1, "GeneNode",gene2,"GeneNode", "interacts_with", symmetric=True, stv=stv))
                    if gene_level:
                        g.write(evaLink(gene1, "GeneNode",gene2,"GeneNode", "interacts_with", symmetric=True, stv=stv))

                    if not prot1 in proteins:
                        bio = str(data.iloc[i]['BioGRID ID Interactor A']).strip()
                        f.write(evaLink(gene1, "GeneNode", prot1,"MoleculeNode", "expresses", prefix2="Uniprot:"))
                        f.write(evaLink(gene1, "GeneNode", bio,"ConceptNode", "has_biogridID", prefix2="Bio:"))
                        f.write(evaLink(prot1, "MoleculeNode", bio,"ConceptNode", "has_biogridID", prefix1="Uniprot:",prefix2="Bio:"))

                    if not prot2 in proteins:
                        bio = str(data.iloc[i]['BioGRID ID Interactor B']).strip()
                        f.write(evaLink(gene2, "GeneNode", prot2,"MoleculeNode", "expresses", prefix2="Uniprot:"))
                        f.write(evaLink(gene2, "GeneNode", bio,"ConceptNode", "has_biogridID", prefix2="Bio:"))
                        f.write(evaLink(prot2, "MoleculeNode", bio,"ConceptNode", "has_biogridID", prefix1="Uniprot:",prefix2="Bio:"))

                    if taxonomy_id_1 == 2697049: 
                        covid_genes.append(gene1)
                        f.write(evaLink(gene1, "GeneNode", taxonomy_id_1 , "ConceptNode","from_organism", prefix2="TaxonomyID:"))
                        f.write(evaLink(prot1, "MoleculeNode", taxonomy_id_1, "ConceptNode","from_organism", prefix1="Uniprot:",prefix2="TaxonomyID:"))
                        if gene_level:
                            g.write(evaLink(gene1, "GeneNode", taxonomy_id_1 , "ConceptNode","from_organism",prefix2="TaxonomyID:"))
                            g.write(evaLink(prot1, "MoleculeNode", taxonomy_id_1, "ConceptNode","from_organism", prefix1="Uniprot:", prefix2="TaxonomyID:"))
                    if taxonomy_id_2 == 2697049:
                        covid_genes.append(gene2)
                        f.write(evaLink(gene2, "GeneNode", taxonomy_id_2 , "ConceptNode","from_organism",prefix2="TaxonomyID:"))
                        f.write(evaLink(prot2, "MoleculeNode", taxonomy_id_2, "ConceptNode","from_organism", prefix1="Uniprot:",prefix2="TaxonomyID:"))
                        if gene_level:
                            g.write(evaLink(gene2, "GeneNode", taxonomy_id_2 , "ConceptNode", "from_organism", prefix2="TaxonomyID:"))
                            g.write(evaLink(prot2, "MoleculeNode", taxonomy_id_2, "ConceptNode", "from_organism", prefix1="Uniprot:", prefix2="TaxonomyID:"))                            
                    
                    genes.append(gene1)
                    genes.append(gene2)
                    proteins.append(prot1)
                    proteins.append(prot2)
                    pairs.append((gene1,gene2))
        f.write(evaLink("2697049", "ConceptNode", "SARS-CoV-2", "ConceptNode","has_name",prefix1="TaxonomyID:"))
        g.write(evaLink("2697049", "ConceptNode", "SARS-CoV-2", "ConceptNode","has_name",prefix1="TaxonomyID:"))
    pairs = set((a,b) if a<=b else (b,a) for a,b in pairs)
    number_of_interactions = len(pairs)
    script = "https://github.com/MOZI-AI/knowledge-import/coronavirus_biogrid.py"
    metadata.update_meta("Coronavirus Biogrid:"+version, source,script,genes=str(len(set(genes))),prot=len(set(proteins)), interactions=str(number_of_interactions))
    print("Done, check "+'dataset/COVID-19-biogrid'+version+"_"+str(date.today())+'.scm')
    with open("Covid19-genes","w") as co:
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
  