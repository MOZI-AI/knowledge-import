"""
converter from pharmagkb to atomese
http://pharmgkb.org
"""

__author__ = "Anatoly Belikov"
__email__ = "abelikov@singularitynet.io"

import re
import urllib.request
import argparse
import xml.etree.ElementTree as ET
from zipfile import ZipFile
from io import BytesIO
import pandas
from atomwrappers import *


chebi_re = re.compile(".*ChEBI:CHEBI:(\d+).*")
pubchem_re = re.compile(".*PubChem Compound:(\d+).*")
pubchem_re_sub = re.compile(".*PubChem Substance:(\d+).*")
drugbank_re = re.compile('.*DrugBank:D?B?(\d+).*')
re_dict = dict()
re_dict['ChEBI'] = [chebi_re]
re_dict['PubChem'] = [pubchem_re, pubchem_re_sub]
re_dict['DrugBank'] = [drugbank_re]


def pharma_to_id(chem_table, name):
    """
    extract references to the substance from chem_table

    Parameters:
    -----------
    chem_table: pandas.DataFrame
        pharagkb chemicals
    name: str
        pharmagkb id for the substance
    
    Returns
    -------
    dict
        database name: id pairs
    """
    chem = chem_table[chem_table['PharmGKB Accession Id'] == name]
    if not len(chem):
        print("Not found chemical row for {0}".format(name))
        return dict() 
    assert len(chem) == 1
    references = chem['Cross-references'].tolist()
    chemical_id = dict()

    for r in references:
        if isinstance(r, str):
            for id_type, reglist in re_dict.items():
                for regex in reglist:
                    match = regex.match(r)
                    if match:
                        for num in match.groups():
                            chemical_id[id_type] = num
    # todo: convert PubChem to ChEBI if possible
    if not chemical_id:
        print("Not found pubchem or chebi id for {0}".format(name))
    return chemical_id


def gen_chemical_members(mol_id_map, pathway_id):
    """
    Generate member links for the substance pointed mol_id_map
    
    Parameters
    ----------
    mol_id_map: dict
        Database name, id pairs
    pathway_id: str
        pharma gkb id for pathway
        
    Returns
    -------
    tuple[list, list]
        member links connecting substance to pathway
        MoleculeNodes for the substance
    """
    tmp = list()
    molecules = list()
    for (id_type, mol_id) in mol_id_map.items():
        molecule = CMoleculeNode('{0}:{1}'.format(id_type, mol_id))
        member = CMemberLink(molecule,
                             CConceptNode(pathway_id))
        tmp.append(member)
        molecules.append(molecule)
    return tmp, molecules


gene_re = re.compile('([A-Z0-9-]*).*')
def gen_gene_member(gene, pathway_id, organism=None):
    match = gene_re.match(gene)
    assert match is not None
    gene = match.group(1)
    result = []
    member = CMemberLink(CGeneNode(gene), CConceptNode(pathway_id))
    result.append(member)
    if organism is not None:
        org_link = CMemberLink(CGeneNode(gene), CConceptNode("oranism:NCBI{0}".format(organism)))
        result.append(org_link)
    return result


def process_genes(genes_str, pathway_id, organism=None):
    """
    Add member link for comma-separated string of genes
    Parameters
    ----------
    genes_str: str
        genes
    pathway_id: str
        pharmagkb id
    organism: str
        optional parameter, organism genes belong to
        default is that genes belong to human sapiens - ncbi 9606
    """
    tmp = []
    if isinstance(genes_str, str):
        for gene in genes_str.split(','):
            tmp += gen_gene_member(gene.strip(), pathway_id, organism=organism)
    return tmp

uniprot_re = re.compile('.*UniProtKB:([A-Za-z0-9-]+).*')
def process_proteins(pharma_id, pathway_id, genes_data):
    """
    Generate member links for the protein pointed by pharmagkb id
    
    Parameters
    ----------
    pharma_id: str
        pharma gkb id for protein
    pathway_id: str
        pharma gkb id for pathway
    genes_data: pandas.DataFrame
        table of genes data
        
    Returns
    -------
    tuple[list, list]
        member links connecting protein to pathway
        MoleculeNodes for proteins
    """
    tmp = []
    proteins = []
    gene = genes_data[genes_data['PharmGKB Accession Id'] == pharma_id]
    if len(gene) == 0:
        print("no gene for {0} in genes table".format(pharma_id))
        return tmp, proteins
    assert len(gene) == 1
    for ref in gene['Cross-references'].tolist()[0].split(','):
        match = uniprot_re.match(ref)
        if match:
            for prot_id in match.groups():
                members, molecules = gen_chemical_members({'Uniprot': prot_id}, pathway_id)
                tmp += members
                proteins += molecules
        else:
            assert 'uniprot' not in ref.lower()
    return tmp, proteins


def gen_location(protein_node, location_node, pathway_id):
    return CContextLink(CConceptNode(pathway_id),
               CEvaluationLink(CPredicateNode("has_location"),
                               CListLink(protein_node,
                                         location_node)))
                               

def generate_locations(elem, ns, chemical_nodes, pathway_id):
    result = []
    for location_elem in elem.findall('./bp:cellularLocation', ns):
           for location in location_elem.attrib.values():
               match = go_location_re.match(location) 
               if match:
                   go_location = match.group(1)
                   location_node = CConceptNode(go_location)
                   for chemical_node in chemical_nodes:
                       context = gen_location(chemical_node, location_node, pathway_id)
                       result.append(context)
               else:
                   assert 'go:' not in location.lower()
    return result


protein_ref_re = re.compile('pgkb.[a-z0-9]+.[A-Za-z0-9]+.*(PA\d+).*')
go_location_re = re.compile('.*(GO:\d+).*')
def convert_pathway(pathway, chem_data, genes_data, pathway_id, pathway_name, ns):
    ev_name = CEvaluationLink(
                 CPredicateNode("has_name"),
                 CListLink(CConceptNode(pathway_id),
                           CConceptNode(pathway_name)))
    tmp = [ev_name]
    # properties often don't have valid attributes 
    for protein in pathway.findall('./bp:Protein', ns):
        name = protein.find('./bp:standardName', ns).text
        organism = None
        if name.strip().startswith('HIV'):
            organism = '12721'
        tmp += process_genes(name, pathway_id, organism) 
        protein_ref_id = None
        for key, value in protein.attrib.items():
            match = protein_ref_re.match(value)
            if match:
                protein_ref_id = match.group(1)
        if protein_ref_id is None:
            for ent_ref in protein.findall('./bp:entityReference', ns):
                for value in ent_ref.attrib.values():
                    match = protein_ref_re.match(value)
                    if match:
                        protein_ref_id = match.group(1)
                    else:
                        assert 'PA' not in value
        if protein_ref_id is None:
           # give up
           name = protein.find('./bp:standardName', ns).text
           print("can't map protein to uniprot for {0}".format(name))
           continue
        members, protein_nodes = process_proteins(protein_ref_id, pathway_id, genes_data)
        tmp += generate_locations(protein, ns, protein_nodes, pathway_id)
        tmp += members
    for smallmolecule in pathway.findall('./bp:SmallMolecule', ns):
        reference = smallmolecule.find('bp:entityReference', ns)
        assert reference is not None
        for value in reference.attrib.values():
            pharma_pkg_id = re.match('.*\.(PA\d+)\.?.*', value)
            if pharma_pkg_id is None:
                # try by standard name
                name = smallmolecule.find('./bp:standardName', ns).text.lower()
                row = chem_data[chem_data.Name == name]
                if len(row):
                    assert len(row) == 1
                    pharma_pkg_id = row.iloc[0]['PharmGKB Accession Id']
                else:
                    print("no pharmapkg id for {0}".format(value))
                    continue
            else:
                pharma_pkg_id = pharma_pkg_id.group(1)
            molecule_drug = pharma_to_id(chem_data, pharma_pkg_id)
            members, chemical_nodes = gen_chemical_members(molecule_drug, pathway_id)
            tmp += generate_locations(smallmolecule, ns, chemical_nodes, pathway_id)
            tmp += members
    return '\n'.join([x.recursive_print() for x in tmp])


# https://effbot.org/zone/element-namespaces.htm
def parse_map(source_file):
    """
    Extract namespaces from xml file
    """

    events = "start", "start-ns", "end-ns"

    root = None
    ns_map = []
    result = dict()
    for event, elem in ET.iterparse(source_file, events):
        if event == "start-ns":
            ns_map.append(elem)
        elif event == "end-ns":
            ns_map.pop()
        elif event == "start":
            if root is None:
                root = elem
            result.update(dict(ns_map))

    return result


PATHWAY_RE = re.compile('(PA\d+)-(\w+).owl')
def get_pathway_id_name(root, ns):
    name = root.findall('./bp:Pathway/bp:displayName', ns)[0].text
    tmp = root.findall('./bp:Pathway[@rdf:about]', ns)
    res = []
    for x in tmp:
        if x.findall('./bp:pathwayComponent', ns):
            name = x.find('./bp:displayName', ns).text
            res.append(x)
            continue
    assert len(res) == 1
    res = res[0].find('./bp:xref[@rdf:resource]', ns)
    pathway_id = None
    if res is not None:
        for k,v in res.attrib.items():
            if k.endswith('resource'):
                pathway_id = v.split('/')[-1].split('.')[-1]
    return pathway_id, name


def build_request(url):
    headers = dict()
    # seems to work with wget headers 
    headers['User-Agent'] = "Wget/1.19.5 (linux-gnu)"
    headers['Accept'] = '*/*'
    headers['Accept-Encoding'] = 'identity'
    headers['Connection'] = 'Keep-Alive'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    response = urllib.request.urlopen(req)
    return response


def download():
    pathway = 'https://s3.pgkb.org/data/pathways-biopax.zip' 
    genes = 'https://s3.pgkb.org/data/genes.zip'
    chemicals = 'https://s3.pgkb.org/data/chemicals.zip'
    pathway_zip = ZipFile(BytesIO(build_request(pathway).read()))
    genes_zip = ZipFile(BytesIO(build_request(genes).read()))
    chem_zip = ZipFile(BytesIO(build_request(chemicals).read()))
    return pathway_zip, genes_zip, chem_zip


def parse_args():
    parser = argparse.ArgumentParser(description='convert biogrid db to atomese')
    parser.add_argument('--pathways', type=str, default='',
                        help='zip archive with pathways in owl format')
    parser.add_argument('--chemicals', type=str, default='',
                        help='zip archive with chemicals in tsv format')
    parser.add_argument('--genes', type=str, default='',
                        help='zip archive with genes data in tsv format')
    parser.add_argument('--output', type=str, default='/tmp/pharmagkb.scm',
                        help='path to output file')
    return parser.parse_args()


def main():
    args = parse_args()
    if (args.pathways and args.genes and args):
        pathway_file = ZipFile(BytesIO(open(args.pathways, 'rb').read()))
        chemicals_file = ZipFile(BytesIO(open(args.chemicals, 'rb').read()))
        genes_file = ZipFile(BytesIO(open(args.genes, 'rb').read()))
    else:
        pathway_file, genes_file, chemicals_file = download()
    chem_tsv = chemicals_file.open('chemicals.tsv')
    genes_tsv = genes_file.open('genes.tsv')

    genes_data = pandas.read_csv(genes_tsv, sep="\t")
    chem_data = pandas.read_csv(chem_tsv, sep="\t")

    pathway_files = [x for x in pathway_file.namelist() if x.endswith('.owl')]

    out_path = args.output
    output = open(out_path, 'wt')
    for filename in pathway_files:
        extracted_file = pathway_file.open(filename).read()
        tree = ET.fromstring(extracted_file)
        ns = parse_map(pathway_file.open(filename))
        pathway_id, pathway_name = get_pathway_id_name(tree, ns)
        if pathway_id is None:
            pathway_id = filename.split('-')[0]
        res = convert_pathway(tree, chem_data, genes_data, pathway_id, pathway_name, ns)
        output.write(res)
        output.write('\n' * 3)
    output.close()


if __name__ == '__main__':
    main()
