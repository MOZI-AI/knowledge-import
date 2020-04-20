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
  tcmid_prescription,
  tcmid_herb,
  tcmid_network,
  tcmid_gnsp,
  tcmid_spectrum
]
tcmid_base_url = "http://119.3.41.228:8000/static/download/"

if os.path.exists(os.path.join(os.getcwd(), output_file)):
  os.remove(output_file)

def evalink(pred, node_type1, node_type2, node1, node2):
  out_fp.write("(EvaluationLink\n")
  out_fp.write("\t(PredicateNode \"" + pred + "\")\n")
  out_fp.write("\t(ListLink\n")
  out_fp.write("\t\t(" + node_type1 + " \"" + node1 + "\")\n")
  out_fp.write("\t\t(" + node_type2 + " \"" + node2 + "\")\n")
  out_fp.write("\t)\n")
  out_fp.write(")\n")

def memblink(node1, node2):
  out_fp.write("(MemberLink\n")
  out_fp.write("\t(ConceptNode \"" + node1 + "\")\n")
  out_fp.write("\t(ConceptNode \"" + node2 + "\")\n")
  out_fp.write(")\n")

out_fp = open(output_file, "a", encoding='utf8')

for rar_name in tcmid_source_rars:
  rar_path = "raw_data/{}".format(rar_name)

  if os.path.exists(rar_path):
    print("Removing file: {}".format(rar_path))
    os.remove(rar_path)

  rar_file = wget.download(tcmid_base_url + rar_name, "raw_data")
  print("\nFile downloaded: {}".format(rar_file))

out_fp.close()
