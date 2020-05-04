# To run:
# PYTHONIOENCODING=UTF-8 python3 drugbank.py

import xml.etree.ElementTree as ET

xml_file = "raw_data/drugbank/full_database.xml"
tag_prefix = "{http://www.drugbank.ca}"

def find_tag(obj, tag):
  return obj.find(tag_prefix + tag)

def findall_tag(obj, tag):
  return obj.findall(tag_prefix + tag)

def get_child_tag_text(obj, tag):
  return find_tag(obj, tag).text

for drug in ET.parse(xml_file).getroot():
  drugbank_id = get_child_tag_text(drug, "drugbank-id")
  name = get_child_tag_text(drug, "name")
  description = get_child_tag_text(drug, "description")
  for group in findall_tag(find_tag(drug, "groups"), "group"):
    drug_group = group.text
  for article in findall_tag(find_tag(find_tag(drug, "general-references"), "articles"), "article"):
    pubmed_id = get_child_tag_text(article, "pubmed-id")
  for other_drug in findall_tag(find_tag(drug, "drug-interactions"), "drug-interaction"):
    # TODO: Need to get ChEBI ID for other_drug
    other_drug = get_child_tag_text(other_drug, "drugbank-id")
  for pathway in findall_tag(find_tag(drug, "pathways"), "pathway"):
    smpdb_id = get_child_tag_text(pathway, "smpdb-id")
    for involved_drug in findall_tag(find_tag(pathway, "drugs"), "drug"):
      # TODO: Need to get ChEBI ID for involved_drug
      involved_drug = get_child_tag_text(involved_drug, "drugbank-id")
    for uniprot_id in findall_tag(find_tag(pathway, "enzymes"), "uniprot-id"):
      uniprot_id = uniprot_id.text
