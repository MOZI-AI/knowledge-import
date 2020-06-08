import argparse
import xml.etree.ElementTree as ET
from atomwrappers import *
import requests

def parse_args():
    parser = argparse.ArgumentParser(description='convert human view db to atomese')
    parser.add_argument('--dbowl', type=str, default='',
                        help='owl file with ontology')
    parser.add_argument('--output', type=str, default='',
                        help='output schema file with atomspace (result.scm for example)')
    return parser.parse_args()

def parseProperties(filename, property):
    props = {}
    with open(filename, 'rb') as file:
        for line in file:
            if line[:17].decode("utf-8") == property:
                line = line.decode("utf-8")
                splitted = line.split(":")[2].split("(")
                props[splitted[0].split(" ")[0]] = splitted[1].split(")")[0]
    return props

def transform_todots(id):
    id = id.split("_")
    id = id[0] + ":" + id[1]
    return id

def checkIfInProps(id, props1, props2):
    if id in props1:
        id = props1[id]
    elif id in props2:
        id = props2[id]
    else:
        id = transform_todots(id)
    return id

def makeSatisfLink(variable, restr, list_member):
    tmp_eval_link = CEvaluationLink(CPredicateNode(restr), CListLink(CVariable(variable), list_member))
    return CSatisfyingSetScopeLink(CVariable(variable), tmp_eval_link)

def makeRestrictionParentLink(id, parent, restriction_label, restriction_svf):
    tmp_eval_link = CEvaluationLink(CPredicateNode(restriction_label), CConceptNode(restriction_svf))
    tmp_simil_link = CSimilarityLink(CConceptNode(id), CConceptNode(parent))
    return CContextLink(tmp_eval_link, tmp_simil_link)

def makeRestrictionParentLink2(id, parent, restriction_svf):
    tmp_satisf_link = makeSatisfLink("$X", restriction_svf, CConceptNode(parent))
    return CSubsetLink(CConceptNode(id), tmp_satisf_link)

def makeTwoRestrictionParentLink(id, parent, restriction_svf, restriction_svf2):
    tmp_satisf_link = makeSatisfLink("$X", restriction_svf2, CConceptNode(parent))
    tmp_satisf_link2 = makeSatisfLink("$Y", restriction_svf, tmp_satisf_link)
    return CSubsetLink(CConceptNode(id), tmp_satisf_link2)

def makeDifficultCase(class_id_1, class_id_2, class_id_3, svf_prop_1, svf_prop_2):
    tmp_satisf_link = makeSatisfLink("$X", svf_prop_1, CConceptNode(class_id_2))
    part_1 = CAndLink(CConceptNode(class_id_1), tmp_satisf_link)
    part_2 = makeSatisfLink("$Y", svf_prop_2, CConceptNode(class_id_3))
    return CSubsetLink(part_1, part_2)

def makeRareCase(class_id_1, class_id_2, class_id_3, restriction_svf):
    tmp_satisf_link = makeSatisfLink("$X", restriction_svf, CConceptNode(class_id_3))
    tmp_and_link = CAndLink(CConceptNode(class_id_2), tmp_satisf_link)
    return CSubsetLink(CConceptNode(class_id_1), tmp_and_link)

def makeRareCase2(classes, restrictions):
    tmp_satisf_link = makeSatisfLink("$X", restrictions[-1], CConceptNode(classes[-1]))
    tmp_and_link = CAndLink(CConceptNode(classes[-2]), tmp_satisf_link)
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrictions[-2], tmp_and_link)
    tmp_and_link_2 = CAndLink(CConceptNode(classes[-3]), tmp_satisf_link_2)
    tmp_satisf_link_3 = makeSatisfLink("$Z", restrictions[-3], tmp_and_link_2)
    return CSubsetLink(CConceptNode(classes[-4]), tmp_satisf_link_3)

def makeDifficultCase2(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    part2 = makeSatisfLink("$Y", restrs[-2], tmp_satisf_link)
    tmp_satisf_link_2 = makeSatisfLink("$Z", restrs[-3], CConceptNode(classes[-2]))
    part1 = CAndLink(CConceptNode(classes[0]), tmp_satisf_link_2)
    return CSubsetLink(part1, part2)

def makeDifficultCaseCore(classes, restrs):
    part2 = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link = makeSatisfLink("$Y", restrs[-2], CConceptNode(classes[-2]))
    part1 = CAndLink(CConceptNode(classes[0]), tmp_satisf_link)
    return CSubsetLink(part1, part2)

def findAllInLine(line, toFind):
    founded_ids = []
    res = 0
    prev_res = 0
    while not (res == -1):
        res = line[res:].find(toFind)
        if (not (res == -1)):
            res = res + len(toFind) + prev_res
            prev_res = res
            founded_ids.append(res)
    return founded_ids



def parseSubclassLine(line, props1, props2):
    obo_ids = findAllInLine(line, "obo:")
    core_ids = findAllInLine(line, "core:")
    svf_ids = findAllInLine(line, "ObjectSomeValuesFrom(")
    intersection_ids = findAllInLine(line, "ObjectIntersectionOf(")
    min_cardinality_id = findAllInLine(line, "ObjectMinCardinality(")
    max_cardinality_id = findAllInLine(line, "ObjectMaxCardinality(")
    classes = []
    restrs = []
    if(len(min_cardinality_id) > 0) or (len(max_cardinality_id) > 0):
        return None #currently dont know how to deal with cardinalities
    if ((len(obo_ids) == 2) & (len(core_ids) == 0) & (len(svf_ids) == 0)):
        class_id = transform_todots(line[obo_ids[0]:obo_ids[1]].split(" ")[0])
        parent_id = transform_todots(line[obo_ids[1]:].split(")")[0])
        return CSubsetLink(CConceptNode(class_id), CConceptNode(parent_id))
    elif ((len(obo_ids) == 3) & (len(core_ids) == 0) & (len(svf_ids) == 1)):
        class_id = transform_todots(line[obo_ids[0]:obo_ids[1]].split(" ")[0])
        svf_prop = checkIfInProps(line[obo_ids[1]:obo_ids[2]].split(" ")[0], props1, props2)
        parent_id = transform_todots(line[obo_ids[2]:].split(")")[0])
        return makeRestrictionParentLink2(class_id, parent_id, svf_prop)
    elif ((len(obo_ids) == 2) & (len(core_ids) == 1) & (len(svf_ids) == 1)):
        class_id = transform_todots(line[obo_ids[0]:core_ids[0]].split(" ")[0])
        core = line[core_ids[0]:obo_ids[1]].split(" ")[0]
        parent_id = transform_todots(line[obo_ids[1]:].split(")")[0])
        return makeRestrictionParentLink2(class_id, parent_id, core)
    elif ((len(obo_ids) == 4) & (len(core_ids) == 0) & (len(svf_ids) == 2) & (len(intersection_ids) == 0)):
        if((svf_ids[0] > obo_ids[0]) & (svf_ids[0] < obo_ids[1]) & (svf_ids[1] > obo_ids[1]) & (svf_ids[1] < obo_ids[2]) & (svf_ids[1] < obo_ids[3])):
            class_id = transform_todots(line[obo_ids[0]:obo_ids[1]].split(" ")[0])
            svf_prop_1 = checkIfInProps(line[obo_ids[1]:obo_ids[2]].split(" ")[0], props1, props2)
            svf_prop_2 = checkIfInProps(line[obo_ids[2]:obo_ids[3]].split(" ")[0], props1, props2)
            parent_id = transform_todots(line[obo_ids[3]:].split(")")[0])
            return makeTwoRestrictionParentLink(class_id, parent_id, svf_prop_1, svf_prop_2)
    elif((len(obo_ids) == 5) & (len(svf_ids) == 2) & (len(intersection_ids) == 1)):
         if ((obo_ids[0] > intersection_ids[0]) & (svf_ids[0] > obo_ids[0]) & (obo_ids[1] > svf_ids[0]) &
                 (obo_ids[2] < svf_ids[1])):
            class_id_1 = transform_todots(line[obo_ids[0]:obo_ids[1]].split(" ")[0])
            svf_prop_1 = checkIfInProps(line[obo_ids[1]:obo_ids[2]].split(" ")[0], props1, props2)
            class_id_2 = transform_todots(line[obo_ids[2]:obo_ids[3]].split(")")[0])
            svf_prop_2 = checkIfInProps(line[obo_ids[3]:obo_ids[4]].split(" ")[0], props1, props2)
            class_id_3 = transform_todots(line[obo_ids[4]:].split(")")[0])
            return makeDifficultCase(class_id_1, class_id_2, class_id_3, svf_prop_1, svf_prop_2)
    elif((len(obo_ids) == 4) & (len(intersection_ids) == 1) & (len(svf_ids) == 1)):
        if((obo_ids[0] < intersection_ids[0]) & (obo_ids[1] > intersection_ids[0]) & (obo_ids[2] > svf_ids[0])):
            class_id_1 = transform_todots(line[obo_ids[0]:obo_ids[1]].split(" ")[0])
            class_id_2 = transform_todots(line[obo_ids[1]:obo_ids[2]].split(" ")[0])
            class_id_3 = transform_todots(line[obo_ids[3]:].split(")")[0])
            restr_svf = checkIfInProps(line[obo_ids[2]:obo_ids[3]].split(" ")[0], props1, props2)
            return makeRareCase(class_id_1, class_id_2, class_id_3, restr_svf)
    elif((len(obo_ids) == 7) & (len(intersection_ids) == 2) & (len(svf_ids) == 3)):
        classes.append(transform_todots(line[obo_ids[0]:obo_ids[1]].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[1]:obo_ids[2]].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[2]:obo_ids[3]].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[3]:obo_ids[4]].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[4]:obo_ids[5]].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[5]:obo_ids[6]].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[6]:].split(")")[0]))
        return makeRareCase2(classes, restrs)
    elif((len(intersection_ids) == 1) & (len(svf_ids) == 3) & (len(obo_ids) == 6)):
        classes.append(transform_todots(line[obo_ids[0]:obo_ids[1]].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[1]:obo_ids[2]].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[2]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[3]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[5]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[4]:].split(" ")[0], props1, props2))
        return makeDifficultCase2(classes, restrs)
    elif ((len(intersection_ids) == 1) & (len(svf_ids) == 2) & (len(obo_ids) == 4) & (len(core_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:obo_ids[1]].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[1]:obo_ids[2]].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[2]:].split(")")[0]))
        restrs.append(line[core_ids[0]:].split(" ")[0])
        classes.append(line[obo_ids[-1]:].split(")")[0])
        return makeDifficultCaseCore(classes, restrs)
    return None

def makeEquivClasses(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_and_link = CAndLink(CConceptNode(classes[-2]), tmp_satisf_link)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link)

def makeEquivClassesUnion(classes):
    return CSimilarityLink(CConceptNode(classes[0]), COrLink(CConceptNode(classes[1]), CConceptNode(classes[2])))

def makeEquivClasses2svfs1inter(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], CConceptNode(classes[-2]))
    part2 = CAndLink(CConceptNode(classes[-3]), tmp_satisf_link, tmp_satisf_link_2)
    return CSimilarityLink(CConceptNode(classes[0]), part2)

def makeEquivRareCase(classes):
    return CSimilarityLink(CConceptNode(classes[0]), CAndLink(CConceptNode(classes[1]), CConceptNode(classes[2])))

def makeEquivRareCase2(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_and_link = CAndLink(CConceptNode(classes[-2]), tmp_satisf_link)
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], tmp_and_link)
    tmp_and_link_2 = CAndLink(CConceptNode(classes[-3]), tmp_satisf_link_2)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link_2)

def makeEquivSVFCore(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[0], CConceptNode(classes[-1]))
    tmp_and_link = CAndLink(CConceptNode(classes[-2]), tmp_satisf_link)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link)

def makeEquivRareCase3(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[0], CConceptNode(classes[-1]))
    tmp_and_link = CAndLink(CConceptNode(classes[-2]), CConceptNode(classes[-3]), tmp_satisf_link)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link)

def makeEquiv2Core2Svf1Inters(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], CConceptNode(classes[-2]))
    tmp_and_link = CAndLink(CConceptNode(classes[-3]), tmp_satisf_link_2, tmp_satisf_link)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link)

def makeEquiv3Svf1Intersect(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], CConceptNode(classes[-2]))
    tmp_satisf_link_3 = makeSatisfLink("$Z", restrs[-3], CConceptNode(classes[-3]))
    tmp_and_link = CAndLink(CConceptNode(classes[-4]), tmp_satisf_link_3, tmp_satisf_link_2, tmp_satisf_link)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link)

def makeEquiv1Intersect2Svf1Core(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], CConceptNode(classes[-2]))
    tmp_and_link = CAndLink(CConceptNode(classes[-3]), tmp_satisf_link_2, tmp_satisf_link)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link)

def makeEquiv3Union(classes):
    tmp_or_link = COrLink(CConceptNode(classes[-1]), CConceptNode(classes[-2]), CConceptNode(classes[-3]))
    return CSimilarityLink(CConceptNode(classes[0]), tmp_or_link)

def makeEquiv4Svf1Inters(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], CConceptNode(classes[-2]))
    tmp_satisf_link_3 = makeSatisfLink("$Z", restrs[-3], CConceptNode(classes[-3]))
    tmp_satisf_link_4 = makeSatisfLink("$U", restrs[-4], CConceptNode(classes[-4]))
    tmp_and_link = CAndLink(CConceptNode(classes[-5]), tmp_satisf_link_4, tmp_satisf_link_3, tmp_satisf_link_2, tmp_satisf_link)
    return CSimilarityLink(CConceptNode(classes[-6]), tmp_and_link)

def makeEquiv4Svf1Inters2(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    and_part_3 = makeSatisfLink("$Y", restrs[-2], tmp_satisf_link)
    tmp_satisf_link_2 = makeSatisfLink("$Z", restrs[-3], CConceptNode(classes[-2]))
    and_part_2 = makeSatisfLink("$U", restrs[-4], tmp_satisf_link_2)
    tmp_and_link = CAndLink(CConceptNode(classes[-3]), and_part_2, and_part_3)
    return CSimilarityLink(CConceptNode(classes[-4]), tmp_and_link)

def makeEquiv2Svf1Intersect(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], tmp_satisf_link)
    tmp_and_link = CAndLink(CConceptNode(classes[-2]), tmp_satisf_link_2)
    return CSimilarityLink(CConceptNode(classes[-3]), tmp_and_link)

def makeEquiv3Svf1Union(classes, restrs):
    part_2 = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link = makeSatisfLink("$Y",  restrs[-2], CConceptNode(classes[-2]))
    tmp_satisf_link_2 = makeSatisfLink("$Z", restrs[-3], CConceptNode(classes[-3]))
    part_1 = CAndLink(tmp_satisf_link_2, tmp_satisf_link)
    return CSimilarityLink(part_1, part_2)

def makeEquivRareCase4(classes, restrs):
    part_2 = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link = makeSatisfLink("$Y", restrs[-2], CConceptNode(classes[-2]))
    tmp_satisf_link_2 = makeSatisfLink("$Z", restrs[-3], CConceptNode(classes[-3]))
    tmp_satisf_link_3 = makeSatisfLink("$U", restrs[-4], CConceptNode(classes[-4]))
    part_1 = COrLink(tmp_satisf_link_3, tmp_satisf_link_2, tmp_satisf_link)
    return CSimilarityLink(part_1, part_2)

def makeEquiv4Svf1Inter1Core(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], CConceptNode(classes[-2]))
    tmp_satisf_link_3 = makeSatisfLink("$Z", restrs[-3], CConceptNode(classes[-3]))
    tmp_satisf_link_4 = makeSatisfLink("$U", restrs[-4], CConceptNode(classes[-4]))
    tmp_and_link = CAndLink(CConceptNode(classes[-5]), tmp_satisf_link_4, tmp_satisf_link_3, tmp_satisf_link_2, tmp_satisf_link)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link)

def makeEquivRareCase5(classes):
    tmp_or_link = COrLink(CConceptNode(classes[-1]), CConceptNode(classes[-2]), CConceptNode(classes[-3]), CConceptNode(classes[-4]))
    return CSimilarityLink(CConceptNode(classes[0]), tmp_or_link)

def makeEquiv3Svf1Intersect1Core(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], CConceptNode(classes[-2]))
    tmp_satisf_link_3 = makeSatisfLink("$Z", restrs[-3], CConceptNode(classes[-3]))
    tmp_and_link = CAndLink(CConceptNode(classes[-4]), tmp_satisf_link, tmp_satisf_link_2, tmp_satisf_link_3)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link)

def makeEquiv3svf1Intersect2(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], tmp_satisf_link)
    tmp_satisf_link_3 = makeSatisfLink("$Z", restrs[-3], CConceptNode(classes[-2]))
    tmp_and_link = CAndLink(CConceptNode(classes[-3]), tmp_satisf_link_3, tmp_satisf_link_2)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link)

def makeEquivRareCase6(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], CConceptNode(classes[-2]))
    tmp_satisf_link_3 = makeSatisfLink("$Z", restrs[-3], CConceptNode(classes[-3]))
    tmp_and_link = CAndLink(CConceptNode(classes[-4]), tmp_satisf_link_3, tmp_satisf_link_2, tmp_satisf_link)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link)

def makeEquiv5Svf1Intersect(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], CConceptNode(classes[-2]))
    tmp_satisf_link_3 = makeSatisfLink("$Z", restrs[-3], CConceptNode(classes[-3]))
    tmp_satisf_link_4 = makeSatisfLink("$U", restrs[-4], CConceptNode(classes[-4]))
    tmp_satisf_link_5 = makeSatisfLink("$V", restrs[-5], CConceptNode(classes[-5]))
    tmp_and_link = CAndLink(CConceptNode(classes[-6]), tmp_satisf_link_5, tmp_satisf_link_4, tmp_satisf_link_3, tmp_satisf_link_2, tmp_satisf_link)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link)

def makeEquivRareCase7(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], tmp_satisf_link)
    tmp_satisf_link_3 = makeSatisfLink("$Z", restrs[-3], CConceptNode(classes[-2]))
    tmp_satisf_link_4 = makeSatisfLink("$U", restrs[-4], CConceptNode(classes[-3]))
    tmp_and_link = CAndLink(CConceptNode(classes[-4]), tmp_satisf_link_4, tmp_satisf_link_3, tmp_satisf_link_2)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link)

def makeEquivRareCase8(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], CConceptNode(classes[-2]))
    tmp_satisf_link_3 = makeSatisfLink("$Z", restrs[-3], CConceptNode(classes[-3]))
    tmp_satisf_link_4 = makeSatisfLink("$U", restrs[-4], CConceptNode(classes[-4]))
    tmp_satisf_link_5 = makeSatisfLink("$U", restrs[-5], CConceptNode(classes[-5]))
    tmp_and_link = CAndLink(CConceptNode(classes[-6]), tmp_satisf_link_5, tmp_satisf_link_4, tmp_satisf_link_3, tmp_satisf_link_2, tmp_satisf_link)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link)

def makeEquivRareCase9(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], CConceptNode(classes[-2]))
    tmp_satisf_link_3 = makeSatisfLink("$Z", restrs[-3], CConceptNode(classes[-3]))
    tmp_satisf_link_4 = makeSatisfLink("$U", restrs[-4], CConceptNode(classes[-4]))
    tmp_and_link = CAndLink(CConceptNode(classes[-5]), tmp_satisf_link_4, tmp_satisf_link_3, tmp_satisf_link_2, tmp_satisf_link)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link)

def makeEquivRareCase10(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], CConceptNode(classes[-2]))
    tmp_satisf_link_3 = makeSatisfLink("$Z", restrs[-3], CConceptNode(classes[-3]))
    tmp_and_link = CAndLink(CConceptNode(classes[-4]), tmp_satisf_link_3, tmp_satisf_link_2, tmp_satisf_link)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link)

def makeEquivRareCase11(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], CConceptNode(classes[-2]))
    tmp_satisf_link_3 = makeSatisfLink("$Z", restrs[-3], tmp_satisf_link_2)
    tmp_satisf_link_4 = makeSatisfLink("$U", restrs[-4], CConceptNode(classes[-3]))
    tmp_satisf_link_5 = makeSatisfLink("$V", restrs[-5], tmp_satisf_link_4)
    tmp_satisf_link_6 = makeSatisfLink("$W", restrs[-6], CConceptNode(classes[-4]))
    tmp_and_link = CAndLink(CConceptNode(classes[-5]), tmp_satisf_link_6, tmp_satisf_link_5, tmp_satisf_link_3, tmp_satisf_link)
    return CSimilarityLink(CConceptNode(classes[0]), tmp_and_link)

def makeEquivRareCase12(classes, restrs):
    tmp_satisf_link = makeSatisfLink("$X", restrs[-1], CConceptNode(classes[-1]))
    tmp_satisf_link_2 = makeSatisfLink("$Y", restrs[-2], CConceptNode(classes[-2]))
    tmp_satisf_link_3 = makeSatisfLink("$Z", restrs[-3], CConceptNode(classes[-3]))
    tmp_satisf_link_4 = makeSatisfLink("$U", restrs[-4], CConceptNode(classes[-4]))
    tmp_satisf_link_5 = makeSatisfLink("$V", restrs[-5], CConceptNode(classes[-5]))
    tmp_or_link = COrLink(tmp_satisf_link_5, tmp_satisf_link_4, tmp_satisf_link_3, tmp_satisf_link_2)
    return CSimilarityLink(tmp_or_link, tmp_satisf_link)

def parseEquvalentLine(line, props1, props2):
    obo_ids = findAllInLine(line, "obo:")
    core_ids = findAllInLine(line, "core:")
    svf_ids = findAllInLine(line, "ObjectSomeValuesFrom(")
    intersection_ids = findAllInLine(line, "ObjectIntersectionOf(")
    union_ids = findAllInLine(line, "ObjectUnionOf(")
    min_cardinality_id = findAllInLine(line, "ObjectMinCardinality(")
    max_cardinality_id = findAllInLine(line, "ObjectMaxCardinality(")
    classes = []
    restrs = []
    if ((len(min_cardinality_id) > 0) or (len(max_cardinality_id) > 0)):
        return None  # currently dont know how to deal with cardinalities
    elif((len(obo_ids) == 4) & (len(svf_ids) == 1) & (len(intersection_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        return makeEquivClasses(classes, restrs)
    elif ((len(obo_ids) == 3) & (len(union_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[2]:].split(")")[0]))
        return makeEquivClassesUnion(classes)
    elif ((len(obo_ids) == 6) & (len(intersection_ids) == 1) & (len(svf_ids) == 2)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[4]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[5]:].split(")")[0]))
        return makeEquivClasses2svfs1inter(classes, restrs)
    elif((len(obo_ids) == 3) & (len(intersection_ids) == 1) & (len(svf_ids) == 0)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[2]:].split(")")[0]))
        return makeEquivRareCase(classes)
    elif((len(obo_ids) == 6) & (len(svf_ids) == 2) & (len(intersection_ids) == 2)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[4]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[5]:].split(")")[0]))
        return makeEquivRareCase2(classes, restrs)
    elif((len(obo_ids) == 3) & (len(core_ids) == 1) & (len(intersection_ids) == 1) & (len(svf_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(line[core_ids[0]:].split(" ")[0])
        classes.append(transform_todots(line[obo_ids[2]:].split(")")[0]))
        return makeEquivSVFCore(classes, restrs)
    elif((len(obo_ids) == 5) & (len(svf_ids) == 1) & (len(intersection_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[2]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[3]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[4]:].split(")")[0]))
        return makeEquivRareCase3(classes, restrs)
    elif((len(obo_ids) == 4) & (len(core_ids) == 2) & (len(svf_ids) == 2) & (len(intersection_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(line[core_ids[0]:].split(" ")[0])
        classes.append(transform_todots(line[obo_ids[2]:].split(")")[0]))
        restrs.append(line[core_ids[1]:].split(" ")[0])
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        return makeEquiv2Core2Svf1Inters(classes, restrs)
    elif ((len(obo_ids) == 8) & (len(intersection_ids) == 1) & (len(svf_ids) == 3)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[4]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[5]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[6]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[7]:].split(")")[0]))
        return makeEquiv3Svf1Intersect(classes, restrs)
    elif((len(obo_ids) == 5) & (len(svf_ids) == 2) & (len(intersection_ids) == 1) & (len(core_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(line[core_ids[0]:].split(" ")[0])
        classes.append(transform_todots(line[obo_ids[4]:].split(")")[0]))
        return makeEquiv1Intersect2Svf1Core(classes, restrs)
    elif((len(obo_ids) == 4) & (len(union_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[2]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        return makeEquiv3Union(classes)
    elif((len(svf_ids) == 4) & (len(obo_ids) == 10) & (len(intersection_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[4]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[5]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[6]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[7]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[8]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[9]:].split(")")[0]))
        return makeEquiv4Svf1Inters(classes, restrs)
    elif((len(obo_ids) == 8) & (len(svf_ids) == 4) & (len(intersection_ids) == 1) & (len(core_ids) == 0)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        restrs.append(checkIfInProps(line[obo_ids[3]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[4]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[5]:].split(" ")[0], props1, props2))
        restrs.append(checkIfInProps(line[obo_ids[6]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[7]:].split(")")[0]))
        return makeEquiv4Svf1Inters2(classes, restrs)
    elif((len(obo_ids) == 5) & (len(svf_ids) == 2) & (len(intersection_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        restrs.append(checkIfInProps(line[obo_ids[3]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[4]:].split(")")[0]))
        return makeEquiv2Svf1Intersect(classes, restrs)
    elif((len(svf_ids) == 3) & (len(obo_ids) == 6) & (len(union_ids) == 1)):
        restrs.append(checkIfInProps(line[obo_ids[0]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[1]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[4]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[5]:].split(")")[0]))
        return makeEquiv3Svf1Union(classes, restrs)
    elif((len(svf_ids) == 4) & (len(union_ids) == 1) & (len(obo_ids) == 8)):
        restrs.append(checkIfInProps(line[obo_ids[0]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[1]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[4]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[5]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[6]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[7]:].split(")")[0]))
        return makeEquivRareCase4(classes, restrs)
    elif((len(obo_ids) == 9) & (len(core_ids) == 1) & (len(svf_ids) == 4) & (len(intersection_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[4]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[5]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[6]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[7]:].split(")")[0]))
        restrs.append(line[core_ids[0]:].split(" ")[0])
        classes.append(transform_todots(line[obo_ids[8]:].split(")")[0]))
        return makeEquiv4Svf1Inter1Core(classes, restrs)
    elif((len(obo_ids) == 5) & (len(union_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[2]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[3]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[4]:].split(")")[0]))
        return makeEquivRareCase5(classes)
    elif((len(obo_ids) == 7) & (len(core_ids) == 1) & (len(svf_ids) == 3) & (len(intersection_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[4]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[5]:].split(")")[0]))
        restrs.append(line[core_ids[0]:].split(" ")[0])
        classes.append(transform_todots(line[obo_ids[6]:].split(")")[0]))
        return makeEquiv3Svf1Intersect1Core(classes, restrs)
    elif((len(obo_ids) == 7) & (len(svf_ids) == 3) & (len(intersection_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[4]:].split(" ")[0], props1, props2))
        restrs.append(checkIfInProps(line[obo_ids[5]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[6]:].split(")")[0]))
        return makeEquiv3svf1Intersect2(classes, restrs)
    elif((len(obo_ids) == 6) & (len(core_ids) == 2) & (len(svf_ids) == 3) & (len(intersection_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(line[core_ids[0]:].split(" ")[0])
        classes.append(transform_todots(line[obo_ids[4]:].split(")")[0]))
        restrs.append(line[core_ids[1]:].split(" ")[0])
        classes.append(transform_todots(line[obo_ids[5]:].split(")")[0]))
        return makeEquivRareCase6(classes, restrs)
    elif((len(obo_ids) == 12) & (len(svf_ids) == 5) & (len(intersection_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[4]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[5]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[6]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[7]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[8]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[9]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[10]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[11]:].split(")")[0]))
        return makeEquiv5Svf1Intersect(classes, restrs)
    elif((len(obo_ids) == 9) & (len(svf_ids) == 4) & (len(intersection_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[4]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[5]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[6]:].split(" ")[0], props1, props2))
        restrs.append(checkIfInProps(line[obo_ids[7]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[8]:].split(")")[0]))
        return makeEquivRareCase7(classes, restrs)
    elif((len(obo_ids) == 10) & (len(core_ids) == 2) & (len(svf_ids) == 5) & (len(intersection_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[4]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[5]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[6]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[7]:].split(")")[0]))
        restrs.append(line[core_ids[0]:].split(" ")[0])
        classes.append(transform_todots(line[obo_ids[8]:].split(")")[0]))
        restrs.append(line[core_ids[1]:].split(" ")[0])
        classes.append(transform_todots(line[obo_ids[9]:].split(")")[0]))
        return makeEquivRareCase8(classes, restrs)
    elif((len(obo_ids) == 8) & (len(core_ids) == 2) & (len(svf_ids) == 4) & (len(intersection_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[4]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[5]:].split(")")[0]))
        restrs.append(line[core_ids[0]:].split(" ")[0])
        classes.append(transform_todots(line[obo_ids[6]:].split(")")[0]))
        restrs.append(line[core_ids[1]:].split(" ")[0])
        classes.append(transform_todots(line[obo_ids[7]:].split(")")[0]))
        return makeEquivRareCase9(classes, restrs)
    elif ((len(core_ids) == 3) & (len(obo_ids) == 5) & (len(svf_ids) == 3) & (len(intersection_ids) == 1)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(line[core_ids[0]:].split(" ")[0])
        classes.append(transform_todots(line[obo_ids[2]:].split(")")[0]))
        restrs.append(line[core_ids[1]:].split(" ")[0])
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(line[core_ids[2]:].split(" ")[0])
        classes.append(transform_todots(line[obo_ids[4]:].split(")")[0]))
        return makeEquivRareCase10(classes, restrs)
    elif((len(obo_ids) == 12) & (len(intersection_ids) == 1) & (len(svf_ids) == 6)):
        classes.append(transform_todots(line[obo_ids[0]:].split(" ")[0]))
        classes.append(transform_todots(line[obo_ids[1]:].split(" ")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[4]:].split(" ")[0], props1, props2))
        restrs.append(checkIfInProps(line[obo_ids[5]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[6]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[7]:].split(" ")[0], props1, props2))
        restrs.append(checkIfInProps(line[obo_ids[8]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[9]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[10]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[11]:].split(")")[0]))
        return makeEquivRareCase11(classes, restrs)
    elif((len(obo_ids) == 10) & (len(union_ids) == 1) & (len(svf_ids) == 5)):
        restrs.append(checkIfInProps(line[obo_ids[0]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[1]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[2]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[3]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[4]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[5]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[6]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[7]:].split(")")[0]))
        restrs.append(checkIfInProps(line[obo_ids[8]:].split(" ")[0], props1, props2))
        classes.append(transform_todots(line[obo_ids[9]:].split(")")[0]))
        return makeEquivRareCase12(classes, restrs)
    return None

def ParseOntology(filename):
    tmp = []
    annotation_props = parseProperties(filename, "# Annotation Property")
    object_props = parseProperties(filename, "# Object Property")
    with open (filename, 'rb') as file:
        for line in file:
            line = line.decode("utf-8")
            if line[:7] == "# Class":
                splitted = line.split(":")[2].split("(")
                try:
                    class_id = transform_todots(splitted[0]).split(" ")[0]
                    class_name = splitted[1].split(")")[0]
                    tmp.append(CEvaluationLink(CPredicateNode("has_name"),
                                               CListLink(CConceptNode(class_id), CConceptNode(class_name))))
                except Exception:
                    continue
            elif line[:10] == "SubClassOf":
                parsed_subclass_line = parseSubclassLine(line, annotation_props, object_props)
                if(parsed_subclass_line):
                    tmp.append(parsed_subclass_line)
            elif line[:17] == "EquivalentClasses":
                parsed_equivalent_line = parseEquvalentLine(line, annotation_props, object_props)
                if (parsed_equivalent_line):
                    tmp.append(parsed_equivalent_line)
    return tmp


def main():
    args = parse_args()
    if (not (args.dbowl)):
        resp = requests.get("http://ontologies.berkeleybop.org/uberon/subsets/human-view.owl")
        human_view_file = "/home/daddywesker/human-view.owl"
        open(human_view_file, 'wb').write(resp.content)
    else:
        human_view_file = args.dbowl

    db_dict = ParseOntology(human_view_file)

    if (not (args.output)):
        output = open("result.scm", 'wt')
    else:
        output = open(args.output, 'wt')

    classes_to_print = '\n'.join([x.recursive_print() for x in db_dict])
    output.write(classes_to_print)
    output.close()


if __name__ == '__main__':
    main()

