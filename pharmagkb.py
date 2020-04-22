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
from gzip import GzipFile
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



def gen_gene_member(gene, pathway_id, organism=None):
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
def gen_proteins(pharma_id, name, pathway_id, pharma2uniprot):
    """
    Generate member links for the protein pointed by pharmagkb id
    
    Parameters
    ----------
    pharma_id: List[str]
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
    for prot in pharma_id:
        entry = pharma2uniprot[pharma2uniprot.pharma_id == prot + ';']
        if not len(entry):
            print("not found uniprot id for {0}".format(name))
            continue
        for i in range(len(entry)):
            prot_id = entry.iloc[i].Entry
            molecule = CMoleculeNode('Uniprot:{0}'.format(prot_id))
            member = CMemberLink(molecule,
                                CConceptNode(pathway_id))
            tmp.append(member)
            proteins.append(molecule)
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


def parse_molecule(smallmolecule, ns, chem_data):
    """
    Parse SmallMolecule element
    
    Parameters
    ----------
    smallmolecule: xml.etree.ElementTree.Element
        SmallMolecule
    ns: dict
        namespaces from the owl file
    chem_data: pandas.DataFrame
        table chemicals from pharagkb
    
    Returns
    -------
    dict, str
        external db name: id pairs
        human readable name
    """
    reference = smallmolecule.find('bp:entityReference', ns)
    molecule_drug = dict()
    assert reference is not None
    name = smallmolecule.find('./bp:standardName', ns).text.lower()
    for value in reference.attrib.values():
        pharma_pkg_id = re.match('.*\.(PA\d+)\.?.*', value)
        if pharma_pkg_id is None:
            # try by standard name

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
    return molecule_drug, name


def process_small_molecules(pathway, ns, pathway_id, chem_data, id_map):
    tmp = list()
    for smallmolecule in pathway.findall('./bp:SmallMolecule', ns):
        molecule_drug, name = parse_molecule(smallmolecule, ns, chem_data)
        for db_name in ('ChEBI', 'PubChem', 'DrugBank'):
            if db_name in molecule_drug:
                molecule = CMoleculeNode("{0}:{1}".format(db_name, molecule_drug[db_name]))
                member = CMemberLink(molecule,
                                    CConceptNode(pathway_id))
                break
        if not molecule_drug:
            # something unknown
            molecule = CConceptNode(name)
            member = CMemberLink(molecule,
                                CConceptNode(pathway_id))
        tmp.append(member)
        id_map[about(smallmolecule, ns)] = [molecule]
        tmp += generate_locations(smallmolecule, ns, [molecule], pathway_id)
            
    return tmp
 

def match_protein_id(value):
    match = protein_ref_re.match(value)
    if match:
        return match.group(2)
    else:
        assert 'PA' not in value
    return None


def parse_protein(protein, pathway, ns, pathway_id, pharma2uniprot, elem_chemical_map, tmp=None):
    name = protein.find('./bp:standardName', ns).text
    protein_elem_id = about(protein, ns)
    if protein_elem_id in elem_chemical_map:
        return elem_chemical_map[protein_elem_id]

    organism = None
    if name.strip().startswith('HIV'):
        organism = '12721'
    reference_elem = protein.find('./bp:entityReference', ns)
    if reference_elem is None:
        print("failed to find reference for {0}".format(protein_elem_id))
        elem_chemical_map[protein_elem_id] = []
        return []
    ent_ref_id = resource(reference_elem, ns)
    ent_ref = pathway.find('./*[@rdf:about="{0}"]'.format(ent_ref_id), ns)
    protein_ref_id = []
    if ent_ref.tag.endswith('ProteinReference'):
        # it is ether protein group or a protein
        xref = ent_ref.find('./bp:xref', ns)
        if xref is not None:
            value = resource(xref, ns)
            prot_id = match_protein_id(value)
            if prot_id:
                protein_ref_id.append(prot_id)
        else:
            # protein group
            for ent_mem in ent_ref.findall('./bp:memberEntityReference', ns):
                value = resource(ent_mem, ns)
                prot_id = match_protein_id(value)
                if prot_id:
                    protein_ref_id.append(prot_id)
    else:
        import pdb;pdb.set_trace()
    if protein_ref_id is not None:
        members, protein_nodes = gen_proteins(protein_ref_id, name, pathway_id, pharma2uniprot)
        elem_chemical_map[protein_elem_id] = protein_nodes
        if tmp is not None:
            tmp += generate_locations(protein, ns, protein_nodes, pathway_id)
            tmp += members
    else:
        # give up
        name = protein.find('./bp:standardName', ns).text
        print("can't map protein to uniprot for {0}".format(name))
        return None
    return elem_chemical_map.get(protein_elem_id, None)


def process_proteins(pathway, ns, pathway_id, genes_data, pharma2uniprot, elem_chemical_map):
    tmp = list()
    # properties often don't have valid attributes 
    for protein in pathway.findall('./bp:Protein', ns):
        parse_protein(protein, pathway, ns, pathway_id, pharma2uniprot, elem_chemical_map, tmp=tmp)
    return tmp



protein_ref_re = re.compile('.*(\.|/)(PA\d+)')
go_location_re = re.compile('.*(GO:\d+).*')
drug_re = re.compile('pgkb.drug.*(PA\d+)')



def wrap_set(molecules):
    if len(molecules) == 1:
        result = molecules[0]
    else:
        result = CSetLink(*molecules)
    return result


def gen_interaction(interaction, pathway, pathway_id, ns, id_map, interaction_name):
    result = list()
    left_elem = interaction.find('./bp:left', ns)
    right_elem = interaction.find('./bp:right', ns)
    if left_elem is None or right_elem is None:
        print("din't find \"from\" participant in interaction {0}".format(about(interaction, ns)))
        id_map[about(interaction, ns)] = result
        return result
    left = resource(left_elem, ns)
    right = resource(right_elem, ns)
    parse_elem(pathway.find('./*[@rdf:about="{0}"]'.format(left), ns), pathway, ns, pathway_id, id_map)
    parse_elem(pathway.find('./*[@rdf:about="{0}"]'.format(right), ns), pathway, ns, pathway_id, id_map)      
    left_mol = id_map.get(left, ())
    right_mol = id_map.get(right, ())
    if left_mol and right_mol:
        ev = CEvaluationLink(
            CPredicateNode(interaction_name),
            CListLink(wrap_set(left_mol),
                        wrap_set(right_mol)))
        result.append(ev)
    if not result:
        print("failed to parse {0}".format(about(interaction, ns)))
    id_map[about(interaction, ns)] = result
    return result


def parse_subelements(interaction, xpath, pathway, pathway_id, ns, id_map):
    left_items = list()
    for left in interaction.findall(xpath, ns):
        left_id = resource(left, ns)
        parse_elem(pathway.find('./*[@rdf:about="{0}"]'.format(left_id), ns), pathway, ns, pathway_id, id_map)
        left_items += id_map.get(left_id, [])
    return left_items


def gen_conversion(interaction, pathway, pathway_id, ns, id_map):
    result = list()
    left = parse_subelements(interaction, './bp:left', pathway, pathway_id, ns, id_map)
    right = parse_subelements(interaction, './bp:right', pathway, pathway_id, ns, id_map)
    if left and right:
        ev = CEvaluationLink(
                CPredicateNode('conversion_of'),
                CListLink(wrap_set(left),
                        wrap_set(right)))
        result.append(ev)
        id_map[about(interaction, ns)] = [ev]
    else:
        id_map[about(interaction, ns)] = []
    return result


control_name_map = dict()
control_name_map['ACTIVATION'] = 'activation_of'
control_name_map['INHIBITION'] = 'inhibition_of'
control_name_map['leads_to'] = 'leads_to'


def parse_control(control, pathway, ns, pathway_id, id_map):
    result = list()
    controller_el = control.find('./bp:controller', ns)
    if controller_el is None:
        print("no controller in {0}".format(about(control, ns)))
        id_map[about(control, ns)] = result
        return result
    controller_id = resource(controller_el, ns)
    controller = process_component(pathway.find('./*[@rdf:about="{0}"]'.format(controller_id), ns), 
                                   pathway, ns, pathway_id, id_map, [])
    # controlled is a some interaction or control
    controlled_id = control.find('./bp:controlled[@rdf:resource]', ns).attrib['{{{0}}}resource'.format(ns['rdf'])]

    controlled = process_component(pathway.find('./*[@rdf:about="{0}"]'.format(controlled_id), ns), 
                                   pathway, ns, pathway_id, id_map, [])
    control_type = None
    if about(control, ns).startswith('pgkb.leadsTo'):
        print("leadsTo control is not implemented")
    elif about(control, ns).startswith('pgkb.control.transport'):
        print("transport control is not implemented")
        id_map[about(control, ns)] = result
        return result
    else:
        control_type_elem = control.find('./bp:controlType', ns)
        if control_type_elem:
            control_type = control_type_elem.text
        else:
            print("no control type in {0}".format(about(control, ns)))
    if controlled and controller and control_type:
        res = CEvaluationLink(
                    CPredicateNode(control_name_map[control_type]),
                    CListLink(
                        wrap_set(controller),
                        wrap_set(controlled)))
        result.append(res)
    id_map[about(control, ns)] = result
    return result
    

def parse_catalysis(element, pathway, ns, pathway_id, id_map):
    result = list()
    controller_el = element.find('./bp:controller', ns)
    if controller_el is None:
        print("failed to parse - no controller in Catalysis {0}".format(about(element, ns)))
        id_map[about(element, ns)] = []
        return result
    controller_id = resource(controller_el, ns)
    controller = process_component(pathway.find('./*[@rdf:about="{0}"]'.format(controller_id), ns), 
                                   pathway, ns, pathway_id, id_map, result=result)
    # controled is a some interaction
    controlled_id = resource(element.find('./bp:controlled[@rdf:resource]', ns), ns)
    controlled_elem = find_about_element(pathway, ns, controlled_id)
    controlled = process_component(controlled_elem, pathway, ns, pathway_id, id_map, result=result)
    if controlled and controller:
        res = CEvaluationLink(
                CPredicateNode("catalysys_of"),
                CListLink(
                    wrap_set(controller),
                    wrap_set(controlled)))
        result.append(res)
        id_map[about(element, ns)] = [res]
    else:
        id_map[about(element, ns)] = []
    return result


def find_about_element(pathway, ns, elem_id):
    return pathway.find('./*[@rdf:about="{0}"]'.format(elem_id), ns)


def parse_elem(elem, pathway, ns, pathway_id, id_map):
    result = list()
    if about(elem, ns) in id_map:
        return id_map[about(elem, ns)]
    if elem.tag.endswith('Complex'):
        name = elem.find('./bp:standardName', ns).text
        node = CConceptNode(name)
        member = CMemberLink(node,
                        CConceptNode(pathway_id))
        id_map[about(elem, ns)] = [node]
        result.append(node)
        result.append(member)
    elif elem.tag.endswith('PhysicalEntity'):
        print('PhysicalEntity as part of interaction is not supported {0}'.format(elem.find('./bp:standardName', ns).text))
        id_map[about(elem, ns)] = []
    elif elem.tag.endswith('Pathway'):
        name = elem.find('./bp:standardName', ns).text
        print("Pathway as component of interaction is not supported: {0}".format(name))
        id_map[about(elem, ns)] = []
        return result
    elif elem.tag.endswith('Dna') or elem.tag.endswith('Rna'):
        print("Dna and Rna as component of interaction is not supported: {0}".format(about(elem, ns)))
        id_map[about(elem, ns)] = []
        return result
    else:
        import pdb;pdb.set_trace()
    return result
        
        
        
def about(elem, ns):
    ab = '{{{0}}}about'.format(ns['rdf'])
    return elem.attrib[ab]

def resource(elem, ns):
    res = '{{{0}}}resource'.format(ns['rdf'])
    return elem.attrib[res]


def parse_interaction(interaction, pathway, ns, pathway_id, id_map):
    result = []
    perticipant_el = interaction.findall('./bp:participant', ns)
    if not perticipant_el:
        print("no participant in interaction: {0}".format(about(interaction, ns)))
        id_map[about(interaction, ns)] = []
    else:
        participant_id = [resource(p, ns) for p in perticipant_el]
        for par_id in participant_id:
            participant = find_about_element(pathway, ns, par_id)
            parse_elem(participant, pathway, ns, pathway_id, id_map)
            result += id_map.get(par_id, [])
        # it is subiteraction, replace it with it's participant
        id_map[about(interaction, ns)] = result
    return result


def process_component(interaction, pathway, ns, pathway_id, id_map, result):
    interaction_name = interaction.tag.split('}')[-1]
    if about(interaction, ns) in id_map:
        return id_map[about(interaction, ns)]
    if interaction_name == 'BiochemicalReaction':
        result += gen_interaction(interaction, pathway, pathway_id, ns, id_map, 'reaction')
    elif interaction_name == 'Transport':
        result += gen_interaction(interaction, pathway, pathway_id, ns, id_map, 'transport_of')
    elif interaction_name == 'Catalysis':
        result += parse_catalysis(interaction, pathway, ns, pathway_id, id_map)
    elif interaction_name == 'Interaction':
        result += parse_interaction(interaction, pathway, ns, pathway_id, id_map)
    elif interaction_name == 'Control':
        result += parse_control(interaction, pathway, ns, pathway_id, id_map)
    elif interaction_name == 'Conversion':
        result += gen_conversion(interaction, pathway, pathway_id, ns, id_map)
    elif interaction_name == 'Pathway':
        result += parse_elem(interaction, pathway, ns, pathway_id, id_map)
    elif interaction_name == 'ComplexAssembly':
        print("ComplexAssembly parsing is not yet implemented")
        id_map[about(interaction, ns)] = []
    elif interaction_name == 'Complex':
        print("Complex parsing is not yet implemented")
        id_map[about(interaction, ns)] = []
    elif interaction_name == 'Degradation':
        print("Degradation parsing is not yet implemented")
        id_map[about(interaction, ns)] = []
    elif interaction_name == 'TemplateReactionRegulation':
        print("TemplateReactionRegulation parsing is not yet implemented")
        id_map[about(interaction, ns)] = []
    else:
        import pdb;pdb.set_trace()
    return id_map[about(interaction, ns)]


def process_components(pathway, ns, pathway_id, id_map):
    result = list()
    for component in pathway.findall('bp:Pathway/bp:pathwayComponent', ns):
        for comp in component.attrib.values():
            interaction = pathway.find('./*[@rdf:about="{0}"]'.format(comp), ns)
            process_component(interaction, pathway, ns, pathway_id, id_map, result=result)
    return result
        
        
def convert_pathway(pathway, chem_data, genes_data, pharma2uniprot, pathway_id, pathway_name, ns):
    print("processing pathway {0} {1}".format(pathway_id, pathway_name))
    ev_name = CEvaluationLink(
                 CPredicateNode("has_name"),
                 CListLink(CConceptNode(pathway_id),
                           CConceptNode(pathway_name)))
    tmp = [ev_name]
    id_map = dict()
    tmp += process_proteins(pathway, ns, pathway_id, genes_data, pharma2uniprot, id_map)
    tmp += process_small_molecules(pathway, ns, pathway_id, chem_data, id_map)
    tmp += process_components(pathway, ns, pathway_id, id_map)
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


pharma2uniprot_url = 'https://github.com/noskill/knowledge-import/raw/master/uniprot2pharmagkb.tab.gz'
def download():
    pathway = 'https://s3.pgkb.org/data/pathways-biopax.zip' 
    genes = 'https://s3.pgkb.org/data/genes.zip'
    chemicals = 'https://s3.pgkb.org/data/chemicals.zip'

    pathway_zip = ZipFile(BytesIO(build_request(pathway).read()))
    genes_zip = ZipFile(BytesIO(build_request(genes).read()))
    chem_zip = ZipFile(BytesIO(build_request(chemicals).read()))
    pharma2uniprot = GzipFile(fileobj=BytesIO(build_request(pharma2uniprot_url).read()))
    return pathway_zip, genes_zip, chem_zip, pharma2uniprot


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
    parser.add_argument('--pharma2uniprot', type=str, default='',
                        help='path to pharma2uniprot file')
    return parser.parse_args()


def remove_duplicates(pharma2uniprot):
    loc_entry = ([(i,pharma2uniprot.iloc[i]) for (i,x) in enumerate(pharma2uniprot['Cross-reference (PharmGKB)'].tolist()) if x.count('PA') > 1])
    table = pharma2uniprot.iloc[0:0]
    prev_iloc = 0
    for (iloc, entry) in loc_entry:
        table = table.append(pharma2uniprot.iloc[prev_iloc: iloc])
        for pharma_gkb_id in entry['Cross-reference (PharmGKB)'].split(';'):
            if pharma_gkb_id:
                new_entry = entry.copy()
                new_entry['Cross-reference (PharmGKB)'] = pharma_gkb_id + ';'
                table = table.append(new_entry)
        prev_iloc = iloc + 1
    table = table.append(pharma2uniprot.iloc[prev_iloc:])
    table = table.rename(columns={'Cross-reference (PharmGKB)': 'pharma_id'})
    return table


def main():
    args = parse_args()
    if (args.pathways and args.genes and args):
        pathway_file = ZipFile(BytesIO(open(args.pathways, 'rb').read()))
        chemicals_file = ZipFile(BytesIO(open(args.chemicals, 'rb').read()))
        genes_file = ZipFile(BytesIO(open(args.genes, 'rb').read()))
        # file is small, handle it separately
        if args.pharma2uniprot:
            pharma2uniprot_file = GzipFile(args.pharma2uniprot)
        else:
            pharma2uniprot_file = GzipFile(fileobj=BytesIO(urllib.request.urlopen(pharma2uniprot_url).read()))
    else:
        pathway_file, genes_file, chemicals_file, pharma2uniprot_file = download()

    chem_tsv = chemicals_file.open('chemicals.tsv')
    genes_tsv = genes_file.open('genes.tsv')

    genes_data = pandas.read_csv(genes_tsv, sep="\t")
    chem_data = pandas.read_csv(chem_tsv, sep="\t")
    pharma2uniprot = remove_duplicates(pandas.read_csv(pharma2uniprot_file, sep='\t'))
    
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
        res = convert_pathway(tree, chem_data, genes_data, pharma2uniprot, pathway_id, pathway_name, ns)
        output.write(res)
        output.write('\n' * 3)
    output.close()


if __name__ == '__main__':
    main()
