import re
import subprocess
import urllib.request
import argparse
import xml.etree.ElementTree as ET
from zipfile import ZipFile
from io import BytesIO
import pandas


class CAtom:
    pass


class CNode(CAtom):
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return '({0} "{1}")'.format(self.atom_type, self.name) 

    def recursive_print(self, result='', indent=''):
        return result + indent + str(self)


class CLink(CAtom):
    def __init__(self, *atoms):
        self.outgoing = atoms

    def __str__(self):
        outgoing = '\n'.join([str(x) for x in self.outgoing])
        return '({0} {1})'.format(self.atom_type, outgoing)

    def recursive_print(self, result='', indent=''):
        result += '({0}'.format(self.atom_type)
        indent = indent + '    '
        for x in self.outgoing:
            result = x.recursive_print(result + '\n', indent)
        result += ')'
        return result


class CEvaluationLink(CLink):
    atom_type = 'EvaluationLink'

class CPredicateNode(CNode):
    atom_type = 'PredicateNode'

class CConceptNode(CNode):
    atom_type = 'ConceptNode'

class CMoleculeNode(CNode):
    atom_type = 'MoleculeNode'

class CMemberLink(CLink):
    atom_type = 'MemberLink'

class CListLink(CLink):
    atom_type = 'ListLink'

class CGeneNode(CNode):
    atom_type = 'GeneNode'

reaction_names = dict()
reaction_names['Biochemical Reaction'] = 'biochemical_reaction'
reaction_names['Activation'] = 'activation_of'
reaction_names['Transport'] = 'transport_of'


pharma2chebi_map = dict()
pharma2chebi_map['serotonin'] = 28790
pharma2chebi_map['tropisetron'] = 32269

chebi_re = re.compile(".*ChEBI:CHEBI:(\d+).*")
pubchem_re = re.compile(".*PubChem Compound:(\d+).*")
pubchem_re_sub = re.compile(".*PubChem Substance:(\d+).*")
drugbank_re = re.compile('.*DrugBank:D?B?(\d+).*')
re_dict = dict()
re_dict['ChEBI'] = [chebi_re]
re_dict['PubChem'] = [pubchem_re, pubchem_re_sub]
re_dict['DrugBank'] = [drugbank_re]

def pharma_to_id(chem_table, name):
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
    tmp = []
    for (id_type, mol_id) in mol_id_map.items():
        member = CMemberLink(CMoleculeNode('{0}:{1}'.format(id_type, mol_id)),
                             CConceptNode(pathway_id))
        tmp.append(member)
    return tmp


gene_re = re.compile('^([A-Z0-9]*)$')
def gen_gene_member(gene, pathway_id):
    member = CMemberLink(CGeneNode(gene), CConceptNode(pathway_id))
    return [member]


def process_genes(genes_str, pathway_id):
    tmp = []
    if isinstance(genes_str, str):
        for gene in genes_str.split(','):
            tmp += gen_gene_member(gene.strip(), pathway_id)
    return tmp

uniprot_re = re.compile('.*UniProtKB:([A-Za-z0-9-]+).*')
def process_proteins(pharma_id, pathway_id, genes_data):
    tmp = []
    gene = genes_data[genes_data['PharmGKB Accession Id'] == pharma_id]
    if len(gene) == 0:
        print("no gene for {0} in genes table".format(pharma_id))
        return tmp
    assert len(gene) == 1
    for ref in gene['Cross-references'].tolist()[0].split(','):
        match = uniprot_re.match(ref)
        if match:
            for prot_id in match.groups():
                tmp += gen_chemical_members({'Uniprot': prot_id}, pathway_id)
        else:
            assert 'uniprot' not in ref.lower()
    return tmp


protein_ref_re = re.compile('pgkb.[a-z0-9]+.[A-Za-z0-9]+.*(PA\d+).*')
def convert_pathway(pathway, chem_data, genes_data, pathway_id, pathway_name, ns):
    ev_name = CEvaluationLink(
                 CPredicateNode("has_name"),
                 CListLink(CConceptNode(pathway_id),
                           CConceptNode(pathway_name)))
    tmp = [ev_name]
    # properties ofter don't have valid attributes 
    for protein in pathway.findall('./bp:Protein', ns):
        tmp += process_genes(protein.find('./bp:standardName', ns).text, pathway_id) 
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
        tmp += process_proteins(protein_ref_id, pathway_id, genes_data)
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
            tmp += gen_chemical_members(molecule_drug, pathway_id)
    return '\n'.join([x.recursive_print() for x in tmp])
        

def parse_map(tree):

    events = "start", "start-ns", "end-ns"

    root = None
    ns_map = []
    result = dict()
    for event, elem in ET.iterparse(tree, events):
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
        output.write('\n')
    output.close()


if __name__ == '__main__':
    main()


