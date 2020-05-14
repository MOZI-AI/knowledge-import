import argparse
import xml.etree.ElementTree as ET
from atomwrappers import *
import requests

def parse_args():
    parser = argparse.ArgumentParser(description='convert obofoundry db to atomese')
    parser.add_argument('--dbowl', type=str, default='',
                        help='owl file with ontology')
    parser.add_argument('--iao', type=str, default='',
                        help='owl file with iao definitions (bfo_classes_only.owl for example)')
    parser.add_argument('--rofile', type=str, default='',
                        help='owl file with ro definitions (ro.owl for example)')
    parser.add_argument('--output', type=str, default='',
                        help='output schema file with atomspace (result.scm for example)')
    parser.add_argument('--promap', type=str, default='',
                        help='file which contains mapping of PR (PRotein ontology) data to other databases')
    return parser.parse_args()

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

def parseIAO(filename):
    with open(filename, 'rb') as file:
        iao_tree = ET.fromstring(file.read())
    iao_ns = parse_map(filename)
    iao_dict = dict()
    properties = iao_tree.findall("./owl:AnnotationProperty", iao_ns)
    for element in properties:
        if(element.text == None):
            continue
        label = element.findall("./rdfs:label", iao_ns)[0]
        iao_dict[next(iter(element.attrib.values())).split("/")[-1]] = label.text
    return iao_dict

def parseRO(filename):
    with open(filename, 'rb') as file:
        ro_tree = ET.fromstring(file.read())
    ro_ns = parse_map(filename)
    ro_dict = dict()
    properties = ro_tree.findall("./owl:AnnotationProperty", ro_ns)
    for element in properties:
        try:
            label = element.findall("./rdfs:label", ro_ns)[0]
        except Exception:
            continue
        ro_dict[next(iter(element.attrib.values())).split("/")[-1]] = label.text
    object_properties = ro_tree.findall("./owl:ObjectProperty", ro_ns)
    for element in object_properties:
        try:
            label = element.findall("./rdfs:label", ro_ns)[0]
        except Exception:
            continue
        ro_dict[next(iter(element.attrib.values())).split("/")[-1]] = label.text
    return ro_dict

def parseMapping(filename):
    is_a_dict = {}
    exact_dict = {}
    error_list = []
    with open(filename) as file:
        for line in file:
            splitted = line.split("\t")
            inheritance_type = splitted[-1].split("\n")[0]
            if (inheritance_type == 'is_a'):
                is_a_dict[splitted[0]] = splitted[1]
            elif(inheritance_type == 'exact'):
                exact_dict[splitted[0]] = splitted[1]
            else:
                error_list.append(inheritance_type)
    return is_a_dict, exact_dict

def ParseOntology(filename, exact_dict, ro_dict, iao_dict):
    with open(filename, 'rb') as file:
        ontology_tree = ET.fromstring(file.read())
    ontology_ns = parse_map(filename)
    classes = ontology_tree.findall("./owl:Class", ontology_ns)
    tmp = []
    for cls in classes:
        id = cls.findall("./oboInOwl:id", ontology_ns)
        if(not id):
            if (cls.attrib.values()):
                class_name = next(iter(cls.attrib.values())).split("/")[-1].split("_")
                if (len(class_name) == 1):
                    main_node_id = class_name[0]
                else:
                    main_node_id = class_name[0] + ":" + class_name[1]
            else:
                continue
        else:
            main_node_id = id[0].text
        if (main_node_id in exact_dict):
            main_node_id = exact_dict[main_node_id]
        parent_classes = cls.findall("./rdfs:subClassOf", ontology_ns)
        label = cls.findall("./rdfs:label", ontology_ns)
        if(label):
            tmp.append(CEvaluationLink(CPredicateNode("has_name"), CListLink(CConceptNode(main_node_id), CConceptNode(label[0].text))))
        if (parent_classes):
            for parent in parent_classes:
                parent_with_restriction = parent.findall("./owl:Restriction", ontology_ns)
                if(parent_with_restriction):
                    continue
                collection = parent.findall("./owl:Class", ontology_ns)
                if(collection):
                    collection_descr = collection[0].findall("./owl:intersectionOf", ontology_ns)[0].findall("./rdf:Description", ontology_ns)
                    if not collection_descr:
                        continue
                    parent_id = next(iter(collection_descr[0].attrib.values())).split("/")[-1].split("_")
                    if (len(parent_id) == 1):
                        parentid = parent_id[0]
                    else:
                        parentid = parent_id[0] + ":" + parent_id[1]
                    if (parentid in exact_dict):
                        parentid = exact_dict[parentid]
                    restriction = collection[0].findall("./owl:intersectionOf", ontology_ns)[0].findall("./owl:Restriction", ontology_ns)
                    if(restriction):
                        restriction_property = next(iter(restriction[0].findall("./owl:onProperty", ontology_ns)[0].attrib.values())).split("/")[-1]
                        if (restriction_property in ro_dict):
                            restriction_label = ro_dict[restriction_property]
                        elif(restriction_property in iao_dict):
                            restriction_label = iao_dict[restriction_property]
                        else:
                            restriction_label = restriction_property
                        restriction_svf = next(iter(restriction[0].findall("./owl:someValuesFrom", ontology_ns)[0].attrib.values())).split("/")[-1].split("_")
                        if (len(restriction_svf) == 1):
                            restrictionsvf = restriction_svf[1]
                        else:
                            restrictionsvf = restriction_svf[0] + ":" + restriction_svf[1]
                        tmp_eval_link = CEvaluationLink(CPredicateNode(restriction_label), CConceptNode(restrictionsvf))
                        tmp_simil_link = CSimilarityLink(CConceptNode(main_node_id), CConceptNode(parentid))
                        tmp.append(CContextLink(tmp_eval_link, tmp_simil_link))
                else:
                    parent_id = next(iter(parent.attrib.values())).split("/")[-1].split("_")
                    if(len(parent_id) == 1):
                        parentid = parent_id[0]
                    else:
                        parentid = parent_id[0] + ":" + parent_id[1]
                    tmp.append(CInheritanceLink(CConceptNode(main_node_id), CConceptNode(parentid)))
    return tmp

def main():
    args = parse_args()
    if(not(args.dbowl)):
        resp = requests.get("https://raw.githubusercontent.com/obophenotype/cell-ontology/master/cl.owl")
        dbowl_filename = "/tmp/cl.owl"
        open(dbowl_filename, 'wb').write(resp.content)
    else:
        dbowl_filename = args.dbowl
    if(not(args.iao)):
        resp = requests.get("https://raw.githubusercontent.com/BFO-ontology/BFO/v2019-08-26/bfo_classes_only.owl")
        iao_filename = "/tmp/bfo_classes_only.owl"
        open(iao_filename, 'wb').write(resp.content)
    else:
        iao_filename = args.iao
    if (not (args.rofile)):
        resp = requests.get("https://raw.githubusercontent.com/oborel/obo-relations/master/ro.owl")
        ro_filename = "/tmp/ro.owl"
        open(ro_filename, 'wb').write(resp.content)
    else:
        ro_filename = args.rofile
    if(not(args.promap)):
        resp = requests.get("https://proconsortium.org/download/current/promapping.txt")
        promap_filename = "/tmp/promap.txt"
        open(promap_filename, 'wb').write(resp.content)
    else:
        promap_filename = args.promap

    iao_dict = parseIAO(iao_filename)
    ro_dict = parseRO(ro_filename)
    _, exact_dict = parseMapping(promap_filename)
    classes_dict = ParseOntology(dbowl_filename, exact_dict, ro_dict, iao_dict)

    if(not(args.output)):
        output = open("result.scm", 'wt')
    else:
        output = open(args.output, 'wt')

    classes_to_print =  '\n'.join([x.recursive_print() for x in classes_dict])
    output.write(classes_to_print)
    output.close()


if __name__ == '__main__':
    main()
