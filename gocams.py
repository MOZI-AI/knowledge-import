import pandas
import os
import re
import urllib.request
import argparse

from zipfile import ZipFile

from io import BytesIO
from atomwrappers import *
from sif import SIF

tab_re = re.compile('.*(\t).*')


def download():
    url = 'https://s3.amazonaws.com/geneontology-public/gocam/GO-CAMs.sif.zip'
    return BytesIO(urllib.request.urlopen(url).read())


def main():
    args = parse_args()
    if not args.output:
        print('missing output path')
        return 1
    if args.input:
        fileobj = BytesIO(open(args.input, 'rb').read())
    else:
        fileobj = download()
    archive = ZipFile(fileobj)

    
    tmp = []
    for fpath in archive.filelist:
        filename = fpath.orig_filename.split('/')[1]
        pathway_name = filename.split('.')[0]
        tmp.append(CInheritanceLink(
                        CConceptNode(pathway_name),
                        CConceptNode('pathway')))
        sif_file = SIF(archive.open(fpath))
        for line in sif_file.lines:
            
            node_a = CConceptNode(line.node_a)
            relation = CPredicateNode(line.relation)
            for node_b in line.nodes_b:
                ev = CEvaluationLink(
                        relation,
                        CListLink(node_a,
                                  CConceptNode(node_b)))
                ctx = CContextLink(CConceptNode(pathway_name),
                                   ev)
                tmp.append(ctx)

    with open(args.output, 'wt') as f:
        for item in tmp:
            f.write(item.recursive_print() + '\n')


def parse_args():
    parser = argparse.ArgumentParser(description='convert biogrid db to atomese')
    parser.add_argument('--input', type=str, default='',
                        help='path to zip archive with pathways in sif format')
    parser.add_argument('--output', type=str, default='',
                        help='path to output file')
    return parser.parse_args()


main()
