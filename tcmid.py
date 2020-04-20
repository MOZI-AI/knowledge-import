# python3 -m pip install rarfile
# apt install unrar
# PYTHONIOENCODING=utf-8 python3 tcmid.py

import os
import rarfile
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
#  tcmid_network,
  tcmid_gnsp,
  tcmid_spectrum
]
tcmid_base_url = "http://119.3.41.228:8000/static/download/"

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
  return False if entry.strip() == "" or entry.strip().lower() == "na" else True

# ----------
# Keep a record of which part of a herb would be used in a formula
herb_part_dict = {}

out_fp = open(output_file, "a", encoding='utf8')

for rar_name in tcmid_source_rars:
  rar_path = "raw_data/{}".format(rar_name)

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
              evalink("targets_meridian", "ConceptNode", "ConceptNode", pinyin_name, "TCM:" + meri)

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
          ingredient = contents[0].lower().strip()
          gnsp_id = contents[1].replace("\"", "").strip()
          if is_available(ingredient) and is_available(gnsp_id):
            evalink("has_gnsp_id", "MoleculeNode", "ConceptNode", ingredient, gnsp_id)

out_fp.close()
