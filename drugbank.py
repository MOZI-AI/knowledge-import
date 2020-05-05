# To run:
# PYTHONIOENCODING=UTF-8 python3 drugbank.py

import os
import requests
import xml.etree.ElementTree as ET
from datetime import date

xml_file = "raw_data/drugbank/full database.xml"
tag_prefix = "{http://www.drugbank.ca}"
output_file = "dataset/drugbank_{}.scm".format(str(date.today()))

xml_root = ET.parse(xml_file).getroot()

if os.path.exists(os.path.join(os.getcwd(), output_file)):
  os.remove(output_file)
out_fp = open(output_file, "a", encoding = "utf8")

def find_tag(obj, tag):
  return obj.find(tag_prefix + tag)

def findall_tag(obj, tag):
  return obj.findall(tag_prefix + tag)

def get_child_tag_text(obj, tag):
  return find_tag(obj, tag).text

def evalink(pred, node_type1, node_type2, node1, node2):
  print("--- Creating EvaluationLink with:\npredicate = {}\nnode1 = {}\nnode2 = {}\n".format(pred, node1, node2))
  out_fp.write("(EvaluationLink\n")
  out_fp.write("\t(PredicateNode \"" + pred + "\")\n")
  out_fp.write("\t(ListLink\n")
  out_fp.write("\t\t(" + node_type1 + " \"" + node1 + "\")\n")
  out_fp.write("\t\t(" + node_type2 + " \"" + node2 + "\")\n")
  out_fp.write("\t)\n")
  out_fp.write(")\n")

def memblink(node_type1, node_type2, node1, node2):
  print("--- Creating MemberLink with:\nnode1 = {}\nnode2 = {}\n".format(node1, node2))
  out_fp.write("(MemberLink\n")
  out_fp.write("\t(" + node_type1 + " \"" + node1 + "\")\n")
  out_fp.write("\t(" + node_type2 + " \"" + node2 + "\")\n")
  out_fp.write(")\n")

def inhlink(node_type1, node_type2, node1, node2):
  print("--- Creating InheritanceLink with:\nnode1 = {}\nnode2 = {}\n".format(node1, node2))
  out_fp.write("(InheritanceLink\n")
  out_fp.write("\t(" + node_type1 + " \"" + node1 + "\")\n")
  out_fp.write("\t(" + node_type2 + " \"" + node2 + "\")\n")
  out_fp.write(")\n")

def get_pubchem_cid(sid):
  print("--- Getting PubChem CID for SID:{}\n".format(sid))
  response = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/" + sid + "/cids/txt")
  if response.status_code != 200:
    print("=== Failed to get PubChem CID for SID:{}\n".format(sid))
    return None
  else:
    return response.text.strip()

# Go through the whole file once, to get the external IDs
id_dict = {}
for drug in xml_root:
  drugbank_id = get_child_tag_text(drug, "drugbank-id")
  chebi = None
  pubchem_cid = None
  pubchem_sid = None
  for external_id in findall_tag(find_tag(drug, "external-identifiers"), "external-identifier"):
    resource = get_child_tag_text(external_id, "resource")
    identifier = get_child_tag_text(external_id, "identifier")
    if resource == "ChEBI":
      chebi = "ChEBI:" + identifier
    elif resource == "PubChem Compound":
      pubchem_cid = "PubChem:" + identifier
    elif resource == "PubChem Substance":
      cid = get_pubchem_cid(identifier)
      if cid == None:
        pubchem_sid = "PubChemSID:" + identifier
      else:
        pubchem_cid = "PubChem:" + cid
  if chebi != None:
    id_dict[drugbank_id] = chebi
  elif pubchem_cid != None:
    id_dict[drugbank_id] = pubchem_cid
  elif pubchem_sid != None:
    id_dict[drugbank_id] = pubchem_sid
  else:
    id_dict[drugbank_id] = "DrugBank:" + drugbank_id

for drug in xml_root:
  drugbank_id = get_child_tag_text(drug, "drugbank-id")
  standard_id = id_dict.get(drugbank_id)
  name = get_child_tag_text(drug, "name").lower()
  description = get_child_tag_text(drug, "description")

  evalink("has_name", "MoleculeNode", "ConceptNode", standard_id, name)

  if description != None:
    evalink("has_description", "MoleculeNode", "ConceptNode", standard_id, description)

  for group in findall_tag(find_tag(drug, "groups"), "group"):
    drug_group = group.text
    inhlink("MoleculeNode", "ConceptNode", standard_id, drug_group + " drug")

  for article in findall_tag(find_tag(find_tag(drug, "general-references"), "articles"), "article"):
    pubmed_id = get_child_tag_text(article, "pubmed-id")
    if pubmed_id != None:
      evalink("has_pubmedID", "MoleculeNode", "ConceptNode", standard_id, "https://www.ncbi.nlm.nih.gov/pubmed/?term=" + pubmed_id)

  for other_drug in findall_tag(find_tag(drug, "drug-interactions"), "drug-interaction"):
    other_drug_drugbank_id = get_child_tag_text(other_drug, "drugbank-id")
    other_drug_standard_id = id_dict.get(other_drug_drugbank_id)
    evalink("interacts_with", "MoleculeNode", "MoleculeNode", standard_id, other_drug_standard_id)

  for pathway in findall_tag(find_tag(drug, "pathways"), "pathway"):
    smpdb_id = get_child_tag_text(pathway, "smpdb-id")
    for involved_drug in findall_tag(find_tag(pathway, "drugs"), "drug"):
      involved_drug_drugbank_id = get_child_tag_text(involved_drug, "drugbank-id")
      involved_drug_standard_id = id_dict.get(involved_drug_drugbank_id)
      memblink("MoleculeNode", "ConceptNode", involved_drug_standard_id, smpdb_id)
    for uniprot_id in findall_tag(find_tag(pathway, "enzymes"), "uniprot-id"):
      uniprot_id = uniprot_id.text
      evalink("catalyzed_by", "ConceptNode", "MoleculeNode", smpdb_id, "Uniprot:" + uniprot_id)
