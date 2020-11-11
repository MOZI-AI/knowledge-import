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

variables_list = ["$A", "$B", "$C", "$D", "$E", "$F", "$G", "$H", "$I", "$J", "$K", "$L", "$M", "$N", "$O", "$P",
                  "$Q", "$R", "$S", "$T", "$V", "$U", "$W", "$X", "$Y", "$Z"]

def transform_todots(id):
    try:
        _id = id.split("_")
        _id = _id[0] + ":" + _id[1]
    except Exception:
        return id
    return _id

def checkIfInProps(id, props1, props2):
    if id in props1:
        id = props1[id]
    elif id in props2:
        id = props2[id]
    else:
        id = transform_todots(id)
    return id

def checkIfInExact(id, exact_dict):
    id = transform_todots(id)
    if id in exact_dict:
        return exact_dict[id]
    else:
        return id

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

def makeSatisfLink(variable, restr, list_member):
    tmp_eval_link = CEvaluationLink(CPredicateNode(restr), CListLink(CVariable(variable), list_member))
    return CSatisfyingSetScopeLink(CVariable(variable), tmp_eval_link)

def parseInceptionClass(cls, ontology_ns, ro_dict, iao_dict, exact_dict):
    result_links = []
    for owl_cls in cls:
        owl_intersect = owl_cls.findall("./owl:intersectionOf", ontology_ns)
        owl_union = owl_cls.findall("./owl:unionOf", ontology_ns)
        if (owl_intersect):
            for intersect in owl_intersect:
                restrictions = intersect.findall("./owl:Restriction", ontology_ns)
                links_to_intersect = []
                if (restrictions):
                    parsed_restrictions = parseRestrictions(restrictions, ontology_ns, ro_dict, iao_dict, exact_dict)
                    if(not (parsed_restrictions == -1)):
                        for x in parsed_restrictions:
                            links_to_intersect.append(x)
                descriptions = intersect.findall("./rdf:Description", ontology_ns)
                if (descriptions):
                    for descr in descriptions:
                        links_to_intersect.append(parseDescription(descr, exact_dict))
            else:
                result_links.append(CAndLink(*links_to_intersect))
        elif (owl_union):
            for union in owl_union:
                restrictions = union.findall("./owl:Restriction", ontology_ns)
                links_to_union = []
                if (restrictions):
                    parsed_restrictions = parseRestrictions(restrictions, ontology_ns, ro_dict, iao_dict, exact_dict)
                    if(not (parsed_restrictions == -1)):
                        for x in parsed_restrictions:
                            links_to_union.append(x)
                descriptions = union.findall("./rdf:Description", ontology_ns)
                if (descriptions):
                    for descr in descriptions:
                        links_to_union.append(parseDescription(descr, exact_dict))
                result_links.append(COrLink(*links_to_union))
    return result_links

def parseRestrictions(restrictions, ontology_ns, ro_dict, iao_dict, exact_dict):
    var_count = 0
    satisf_links = []
    for rest in restrictions:
        prop = rest.findall("./owl:onProperty", ontology_ns)
        prop_value = checkIfInProps(next(iter(prop[0].attrib.values())).split("/")[-1], ro_dict, iao_dict)
        svf = rest.findall("./owl:someValuesFrom", ontology_ns)
        avf = rest.findall("./owl:allValuesFrom", ontology_ns)
        min_cardinality = rest.findall("./owl:minQualifiedCardinality", ontology_ns) #currently don't know how to deal with those 3
        max_cardinality = rest.findall("./owl:maxQualifiedCardinality", ontology_ns)
        has_self = rest.findall("./owl:hasSelf", ontology_ns)
        if(max_cardinality):
            continue
        if(has_self):
            continue
        if(min_cardinality):
            continue
        if(svf):
            if(len(svf) == 1):
                inception_class = svf[0].findall("./owl:Class", ontology_ns)
                inception_restriction = svf[0].findall("./owl:Restriction", ontology_ns)
                if(inception_class):
                    parsed_inception_class = parseInceptionClass(inception_class, ontology_ns, ro_dict, iao_dict, exact_dict)
                    if(parsed_inception_class == -1):
                        return -1
                    else:
                        svf_value = parsed_inception_class[0]
                elif (inception_restriction):
                    parsed_inception_restrictions = parseRestrictions(inception_restriction, ontology_ns, ro_dict, iao_dict, exact_dict)
                    if (parsed_inception_restrictions == -1):
                        return -1
                    else:
                        svf_value = parsed_inception_restrictions[0]
                else:
                    svf_value = CConceptNode(checkIfInExact(next(iter(svf[0].attrib.values())).split("/")[-1], exact_dict))
        elif(avf):
            svf_value = CConceptNode(checkIfInExact(next(iter(avf[0].attrib.values())).split("/")[-1], exact_dict))
        satisf_links.append(makeSatisfLink(variables_list[var_count], prop_value, svf_value))
        var_count+=1
    if(satisf_links):
        return satisf_links
    else:
        return -1

def parseDescription(description, exact_dict):
    return CConceptNode(checkIfInExact(next(iter(description.attrib.values())).split("/")[-1], exact_dict))

def parseEquivalentClass(equiv_classes, ontology_ns, class_id, exact_dict, ro_dict, iao_dict):
    for eq_cls in equiv_classes:
        owl_class = eq_cls.findall("./owl:Class", ontology_ns)
        if(owl_class):
            return parseInceptionClass(owl_class, ontology_ns, ro_dict, iao_dict, exact_dict)
        else:
            equiv_id = checkIfInExact(next(iter(eq_cls.attrib.values())).split("/")[-1], exact_dict)
            if(equiv_id):
                return CSimilarityLink(CConceptNode(class_id), CConceptNode(equiv_id))
    return -1

def parseSubclass(subclasses, ontology_ns, class_id, exact_dict, ro_dict, iao_dict):
    subclass_links = []
    for parent in subclasses:
        parent_with_restriction = parent.findall("./owl:Restriction", ontology_ns)
        collection = parent.findall("./owl:Class", ontology_ns)
        if (parent_with_restriction):
            restriction_link = parseRestrictions(parent_with_restriction, ontology_ns, ro_dict, iao_dict, exact_dict)
            if(not (restriction_link == -1)):
                for x in restriction_link:
                    subclass_links.append(CSubsetLink(CConceptNode(class_id), x))
        elif(collection):
            collection_link = parseInceptionClass(collection, ontology_ns, ro_dict, iao_dict, exact_dict)
            if (not (collection_link == -1)):
                for x in collection_link:
                    subclass_links.append(CSubsetLink(CConceptNode(class_id), x))
        else:
            parent_id = checkIfInExact(next(iter(parent.attrib.values())).split("/")[-1], exact_dict)
            subclass_links.append(CSubsetLink(CConceptNode(class_id), CConceptNode(parent_id)))
    return subclass_links

def ParseOntology(filename, exact_dict, ro_dict, iao_dict):
    with open(filename, 'rb') as file:
        ontology_tree = ET.fromstring(file.read())
    ontology_ns = parse_map(filename)
    classes = ontology_tree.findall("./owl:Class", ontology_ns)
    tmp = []
    for cls in classes:
        id = cls.findall("./oboInOwl:id", ontology_ns)
        if (not id):
            if (cls.attrib.values()):
                main_node_id = checkIfInExact(next(iter(cls.attrib.values())).split("/")[-1], exact_dict)
            else:
                continue
        else:
            main_node_id = checkIfInExact(id[0].text, exact_dict)
        equiv_classes = cls.findall("./owl:equivalentClass", ontology_ns)
        if (main_node_id in exact_dict):
            main_node_id = exact_dict[main_node_id]
        if (equiv_classes):
            parsed_equiv = parseEquivalentClass(equiv_classes, ontology_ns, main_node_id, exact_dict, ro_dict, iao_dict)
            if(not (parsed_equiv == -1)):
                if(isinstance(parsed_equiv, list)):
                    tmp.extend(parsed_equiv)
                else:
                    tmp.append(parsed_equiv)
        parent_classes = cls.findall("./rdfs:subClassOf", ontology_ns)
        label = cls.findall("./rdfs:label", ontology_ns)
        if(label):
            tmp.append(CEvaluationLink(CPredicateNode("has_name"), CListLink(CConceptNode(main_node_id), CConceptNode(label[0].text))))
        if (parent_classes):
            parsed_subclasses = parseSubclass(parent_classes, ontology_ns, main_node_id, exact_dict, ro_dict, iao_dict)
            if (not (parsed_subclasses == -1)):
                if (isinstance(parsed_subclasses, list)):
                    tmp.extend(parsed_subclasses)
                else:
                    tmp.append(parsed_subclasses)

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
