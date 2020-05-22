"""
Import CAUSAL ACTIVITY MODELS (CAMS) from https://geneontology.cloud/home
"""

import os
import re
import urllib.request
import argparse

from zipfile import ZipFile

from io import BytesIO
import pandas
from atomwrappers import *
from sif import SIF

tab_re = re.compile('.*(\t).*')


def download():
    url = 'https://s3.amazonaws.com/geneontology-public/gocam/GO-CAMs.sif.zip'
    return BytesIO(urllib.request.urlopen(url).read())


class RGDError(RuntimeError):
    pass

non_human_db = ('RGD', 'ZFIN', 'MGI', 'FB', 'WB', 'Xenbase', 'EMAPA', 'SGD')

def format_id(node_id, protomapping):
    if node_id.count(':') == 2:
        node_id = ':'.join(node_id.split(':')[1:])
    if node_id.count(':'):
        db_name, num = node_id.split(':')
        if db_name in ('GO', 'CL', 'PR'):
            return node_id
        if db_name == 'UBERON':
            return node_id
        elif db_name == 'UniProtKB':
            return 'Uniprot:' + num
        elif db_name in non_human_db:
            # it's about rats, zebrafish, mice, fly, worm, frog, yeast
            raise RGDError('RGD')
        elif db_name == 'CHEBI':
            return 'ChEBI:' + num
        elif db_name == 'PR':
            data = protomapping[protomapping.PR == node_id]
            if len(data) == 1:
                ref = data.iloc[0].external
                uniprot = [x for x in ref.split(',') if x.strip().startswith('UniProtKB')]
                if uniprot:
                    return 'Uniprot:' + uniprot[0].strip().split(':')[1]
                import pdb;pdb.set_trace()
            else:
                import pdb;pdb.set_trace()
        else:
            import pdb;pdb.set_trace()
    return node_id



def generate_has_name(concept_id, concept_name):
    if invalid_name.match(concept_name.name) is not None:
        raise RGDError('invalid')
    ev = CEvaluationLink(
            CPredicateNode("has_name"),
            CListLink(concept_id, concept_name))
    return [ev]


invalid_name = re.compile('^({0}):.*'.format('|'.join(non_human_db)))


def process_sif_file(sif_human, sif_ref, protomapping):
    result = []
    for line_h, line_ref in zip(sif_human.lines, sif_ref.lines):
        tmp = []
        try:
            node_a = CConceptNode(format_id(line_ref.node_a, protomapping))
            relation = CPredicateNode(line_ref.relation)
            tmp += generate_has_name(relation, CConceptNode(line_h.relation))
            tmp += generate_has_name(node_a, CConceptNode(line_h.node_a))
            for id_b, name_b_h in zip(line_ref.nodes_b, line_h.nodes_b):
                node_b = CConceptNode(format_id(id_b, protomapping))
                ev = CEvaluationLink(
                        relation,
                        CListLink(node_a,
                                  node_b))
                tmp.append(ev)
                tmp += generate_has_name(node_b, CConceptNode(name_b_h))
        except RGDError:
            continue
        result += tmp
    return result
 

def main():
    args = parse_args()

    protomapping = pandas.read_csv(args.protomapping, sep='\t',
                                   names=['PR', 'external', 'relation'])
    result = []
    sif_human = SIF(args.readable_sif)
    sif_ref = SIF(args.db_ref_sif)
    result += process_sif_file(sif_human, sif_ref, protomapping)

    with open(args.output, 'wt') as f:
        for item in result:
            f.write(item.recursive_print() + '\n')


def parse_args():
    usage = """
    this script will generate one scm file from two sif files.
    First run ttl2sif.py from gocam-sif-pyexport to get two sif files from one ttl:
    one with human readable names and one with ids
    """
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('--readable-sif', type=str, required=True,
                        help='path to sif file with human-readable names')
    parser.add_argument('--output', type=str, required=True,
                        help='path to output file')
    parser.add_argument('--db-ref-sif', type=str, required=True,
                        help='sif with references to databases')
    parser.add_argument('--protomapping', type=str, required=True,
                        help='path to promapping.txt from protoconsortium.org')
    return parser.parse_args()


main()
