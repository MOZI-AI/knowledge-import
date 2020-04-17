# python3 -m pip install rarfile
# apt install unrar
# PYTHONIOENCODING=utf-8 python3 tcmid.py

import os
import rarfile
import wget
from datetime import date

rarfile.UNRAR_TOOL = "unrar"

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
#  tcmid_spectrum
]
tcmid_base_url = "http://119.3.41.228:8000/static/download/"

output_file = "dataset/tcmid_{}.scm".format(str(date.today()))
out_fp = open(output_file, "a")

if os.path.exists(os.path.join(os.getcwd(), output_file)):
  os.remove(output_file)
