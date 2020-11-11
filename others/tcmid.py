# A script to convert TCMID to AtomSpace representation in Scheme
#
# Some additional dependencies may be required:
#   apt install unrar
#   python3 -m pip install rarfile
#
# To run the script:
#   PYTHONIOENCODING=utf-8 python3 tcmid.py

import os
import rarfile
import re
import wget
from datetime import date

rarfile.UNRAR_TOOL = "unrar"

output_file = "dataset/tcmid_{}.scm".format(str(date.today()))
tcmid_prescription = "prescription-TCMID.v2.01.rar"
tcmid_herb = "herb-TCMID.v2.01.rar"
tcmid_network = "ingredient_targets_disease_drug-TCMID.v2.03.rar"
tcmid_gnsp = "Ingredient_MS-TCMID.v2.01.rar"
tcmid_spectrum = "Herb_MS-TCMID.v2.01.rar"
tcmid_source_rars = [
  # Need to know the properties of the herb first, so tcmid_herb
  # should be processed before tcmid_prescription
  tcmid_herb,
  tcmid_prescription,
  tcmid_gnsp,
  tcmid_spectrum,
  tcmid_network
]
tcmid_base_url = "http://119.3.41.228:8000/static/download/"
chebi_obo_file = "raw_data/chebi.obo"
chebi_url = "ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo"

if os.path.exists(os.path.join(os.getcwd(), output_file)):
  os.remove(output_file)

def evalink(pred, node_type1, node_type2, node1, node2):
  print("--- Creating EvaluationLink with:\npredicate = {}\nnode1 = {}\nnode2 = {}\n".format(pred, node1, node2))
  out_fp.write("(EvaluationLink\n")
  out_fp.write("\t(PredicateNode \"" + pred + "\")\n")
  out_fp.write("\t(ListLink\n")
  out_fp.write("\t\t(" + node_type1 + " \"" + node1 + "\")\n")
  out_fp.write("\t\t(" + node_type2 + " \"" + node2 + "\")\n")
  out_fp.write("\t)\n")
  out_fp.write(")\n")

def memblink(node1, node2):
  print("--- Creating MemberLink with:\nnode1 = {}\nnode2 = {}\n".format(node1, node2))
  out_fp.write("(MemberLink\n")
  out_fp.write("\t(ConceptNode \"" + node1 + "\")\n")
  out_fp.write("\t(ConceptNode \"" + node2 + "\")\n")
  out_fp.write(")\n")

def is_available(entry):
  return False if entry == None or entry.strip() == "" or entry.strip().lower() == "na" or entry.strip().lower() == "n/a" else True

# ----------
# Keep a record of which part of a herb would be used in a formula
herb_part_dict = {}

out_fp = open(output_file, "a", encoding='utf8')

for rar_name in tcmid_source_rars:
  rar_path = "raw_data/{}".format(rar_name)

  if os.path.exists(rar_path):
    print("Removing file: {}".format(rar_path))
    os.remove(rar_path)

  rar_file = wget.download(tcmid_base_url + rar_name, "raw_data")
  print("\nFile downloaded: {}".format(rar_file))

  with rarfile.RarFile(rar_file) as rf:
    # There should only be one file per RAR file
    # Decode using UTF-8 for the Chinese characters
    lines = rf.read(rf.infolist()[0]).decode("utf-8", "ignore").split("\n")

    if rar_file.endswith(tcmid_herb):
      # Skip the first line (columns) in this file
      for line in lines[1:]:
        print("--- Reading line: " + line)
        if is_available(line):
          contents = line.split("\t")
          pinyin_name = contents[0]
          english_name = contents[2]
          properties = [x.lower().strip() for x in contents[4].split(",")]
          meridians = [x.lower().strip() for x in contents[5].split(",")]
          use_part = contents[6]
          if is_available(pinyin_name) and is_available(english_name):
            evalink("has_name", "ConceptNode", "ConceptNode", pinyin_name, english_name)
          if is_available(pinyin_name) and is_available(use_part):
            use_part_full_name = "{} {}".format(pinyin_name, use_part)
            herb_part_dict[pinyin_name] = use_part_full_name
            evalink("has_part", "ConceptNode", "ConceptNode", pinyin_name, use_part_full_name)
          for prop in properties:
            if is_available(pinyin_name) and is_available(prop):
              evalink("has_property", "ConceptNode", "ConceptNode", pinyin_name, "TCM:" + prop)
          for meri in meridians:
            if is_available(pinyin_name) and is_available(meri):
              evalink("meridian_affinity", "ConceptNode", "ConceptNode", pinyin_name, "TCM:" + meri)

    elif rar_file.endswith(tcmid_prescription):
      # Skip the first line (columns) in this file
      for line in lines[1:]:
        print("--- Reading line: " + line)
        if is_available(line):
          contents = line.split("\t")
          prescription = contents[0]
          composition = contents[3].split(",")
          for compo in composition:
            if is_available(compo) and is_available(prescription):
              compo_part = herb_part_dict[compo] if compo in herb_part_dict else compo
              evalink("composition", "ConceptNode", "ConceptNode", compo_part, prescription)
              memblink(compo, "herb")
          if is_available(prescription):
            memblink(prescription, "prescription")

    elif rar_file.endswith(tcmid_spectrum):
      # Skip the first line (columns) in this file
      for line in lines[1:]:
        print("--- Reading line: " + line)
        if is_available(line):
          contents = line.split("\t")
          pinyin_name = contents[1]
          spectrum_description = [x.lower().strip() for x in contents[6].split(";")]
          for sd in spectrum_description:
            if is_available(sd) and is_available(pinyin_name):
              evalink("has_hplc_description", "ConceptNode", "ConceptNode", pinyin_name, sd)

    elif rar_file.endswith(tcmid_gnsp):
      # Skip the first line (columns) in this file
      for line in lines[1:]:
        print("--- Reading line: " + line)
        if is_available(line):
          contents = line.split("\t")
          ingredient = contents[0].replace("\"", "").lower().strip()
          gnsp_id = contents[1].replace("\"", "").strip()
          if is_available(ingredient) and is_available(gnsp_id):
            evalink("has_gnsp_id", "MoleculeNode", "ConceptNode", ingredient, gnsp_id)

    elif rar_file.endswith(tcmid_network):
      ##### Get ChEBI info #####
      if os.path.exists(chebi_obo_file):
        print("Removing file: {}".format(chebi_obo_file))
        os.remove(chebi_obo_file)

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

      for line in lines:
        print("--- Reading line: " + line)
        if is_available(line):
          contents = line.split("\t")
          ingredient = contents[0].lower().strip()
          chebi_id = chebi_dict.get(ingredient)
          uniprot_id = contents[2]
          gene = contents[3]
          omim_ids = contents[4].split(";")
          drug_ids = contents[5].split(";")

          full_name = "ChEBI:" + chebi_id if is_available(chebi_id) else "TCM:" + ingredient

          if is_available(uniprot_id):
            evalink("interacts_with", "MoleculeNode", "MoleculeNode", full_name, "Uniprot:" + uniprot_id)
          if is_available(gene):
            evalink("interacts_with", "MoleculeNode", "GeneNode", full_name, gene)
          for omim in omim_ids:
            if is_available(omim):
              evalink("treats", "MoleculeNode", "ConceptNode", full_name, "OMIM:" + omim)
          for drug in drug_ids:
            if is_available(drug) and is_available(uniprot_id):
              evalink("targets", "MoleculeNode", "MoleculeNode", "DrugBank:" + drug, "Uniprot:" + uniprot_id)

out_fp.close()
