# To run:
# PYTHONIOENCODING=UTF-8 python3 drugbank.py

import os
import re
import requests
import wget
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
  try:
    response = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/" + sid + "/cids/txt", timeout=20)
  except:
    print("=== Connection error")
    return None

  if response.status_code != 200:
    print("=== Failed to find a PubChem CID for SID:{}\n".format(sid))
    return None
  else:
    return response.text.strip()

# Get ChEBI IDs for reference later
chebi_obo = "raw_data/chebi.obo"
chebi_url = "ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo"

if os.path.exists(chebi_obo):
  print("Removing file: {}".format(chebi_obo))
  os.remove(chebi_obo)

chebi_file = wget.download(chebi_url, "raw_data")
print("\nFile downloaded: {}".format(chebi_file))

chebi_fp = open(chebi_file, "r", errors="ignore")
chebi_dict = {}
chebi_name = []
chebi_id = None

for line in chebi_fp:
  line = line.replace("\n", "")
  if line == "[Term]":
    if len(chebi_name) > 0 and chebi_id != None:
      for name in chebi_name:
        chebi_dict[name.lower()] = chebi_id
      # print("ChEBI ID: {}\nName: {}\n".format(chebi_id, chebi_name))
    chebi_name = []
    chebi_id = None
  elif line.startswith("id: "):
    chebi_id = line.replace("id: CHEBI:", "")
  elif line.startswith("name: "):
    chebi_name.append(line.replace("name: ", ""))
  elif line.startswith("synonym: ") and "EXACT" in line:
    name = re.match(".+\"(.+)\".+", line).group(1)
    if name not in chebi_name:
      chebi_name.append(name)

chebi_fp.close()

# Then go through the whole file once, to get the external IDs
id_dict = {}
for drug in xml_root:
  drugbank_id = get_child_tag_text(drug, "drugbank-id")
  name = get_child_tag_text(drug, "name").lower()

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
      # Prefix will be added later
      pubchem_sid = identifier

  # Try to get the ChEBI ID from the official database if it's not found in DrugBank
  if chebi == None:
    chebi = chebi_dict.get(name)

  # Try to get the PubChem CID from the official database if it's not found in DrugBank
  if pubchem_cid == None and pubchem_sid != None:
    pubchem_cid = get_pubchem_cid(pubchem_sid)

  if chebi != None:
    id_dict[drugbank_id] = chebi
  elif pubchem_cid != None:
    id_dict[drugbank_id] = pubchem_cid
  elif pubchem_sid != None:
    id_dict[drugbank_id] = "PubChemSID:" + pubchem_sid
  else:
    # If no desired external IDs is found, use the DrugBank ID
    id_dict[drugbank_id] = "DrugBank:" + drugbank_id

# Finally do the conversion for each of the drugs
drug_groups = []
for drug in xml_root:
  drugbank_id = get_child_tag_text(drug, "drugbank-id")
  standard_id = id_dict.get(drugbank_id)
  name = get_child_tag_text(drug, "name").lower()
  description = get_child_tag_text(drug, "description")

  evalink("has_name", "MoleculeNode", "ConceptNode", standard_id, name)

  if description != None:
    evalink("has_description", "MoleculeNode", "ConceptNode", standard_id, description.replace("\"", "\\\"").strip())

  for group in findall_tag(find_tag(drug, "groups"), "group"):
    drug_group = group.text + " drug"
    inhlink("MoleculeNode", "ConceptNode", standard_id, drug_group)
    if drug_group not in drug_groups:
      inhlink("ConceptNode", "ConceptNode", drug_group, "drug")
      drug_groups.append(drug_group)

  for article in findall_tag(find_tag(find_tag(drug, "general-references"), "articles"), "article"):
    pubmed_id = get_child_tag_text(article, "pubmed-id")
    if pubmed_id != None:
      evalink("has_pubmedID", "MoleculeNode", "ConceptNode", standard_id, "https://www.ncbi.nlm.nih.gov/pubmed/?term=" + pubmed_id)

  for other_drug in findall_tag(find_tag(drug, "drug-interactions"), "drug-interaction"):
    other_drug_drugbank_id = get_child_tag_text(other_drug, "drugbank-id")
    other_drug_standard_id = id_dict.get(other_drug_drugbank_id)
    # For some reason a few of them are not in the 'full database' file?
    if other_drug_standard_id == None:
      other_drug_standard_id = other_drug_drugbank_id
    evalink("interacts_with", "MoleculeNode", "MoleculeNode", standard_id, other_drug_standard_id)

  for pathway in findall_tag(find_tag(drug, "pathways"), "pathway"):
    smpdb_id = get_child_tag_text(pathway, "smpdb-id")
    for involved_drug in findall_tag(find_tag(pathway, "drugs"), "drug"):
      involved_drug_drugbank_id = get_child_tag_text(involved_drug, "drugbank-id")
      involved_drug_standard_id = id_dict.get(involved_drug_drugbank_id)
      # For some reason a few of them are not in the 'full database' file?
      if involved_drug_standard_id == None:
        involved_drug_standard_id = involved_drug_drugbank_id
      memblink("MoleculeNode", "ConceptNode", involved_drug_standard_id, smpdb_id)
    for uniprot_id in findall_tag(find_tag(pathway, "enzymes"), "uniprot-id"):
      uniprot_id = uniprot_id.text
      evalink("catalyzed_by", "ConceptNode", "MoleculeNode", smpdb_id, "Uniprot:" + uniprot_id)
