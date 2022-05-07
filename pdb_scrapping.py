import json
from pickletools import float8
from numpy import float128, float64
import wget
import os
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

url_for_pdbs = "https://files.rcsb.org/download/"

# check if data folder exists
if not os.path.isdir('data'):
  os.makedirs('data')

# save the path of data folder
dest_folder = os.path.realpath("./data")

# Opening JSON file
f = open('primary_structures.json')

# returns JSON object as
# a dictionary
data = json.load(f)

pdbs_not_exist = ["6z66", "6z8n", "6z4s", "6zin", "6z4q", "6za8", "6yvr", "6z4v", "6orv", "7coy", "6nwa", "6kig", "6kif", "5zf0", "6kmx", "6kmw", "6k33", "6l4u", "6tcl", "6dhe", "4yuu", "6jlu", "6j3y", "6j40", "7b0n", "6x89", "5gup", "5xth", "5xti", "4wz7", "6h8k", "6qc6", "6q9b", "6qa9", "6qc4", "6q9e", "6tmj", "6yny", "6tmk", "6n2y", "6n2z", "6n30", "6n2d", "6j5k"]
# Iterating through the json and download pdb files
for i in data['objects']:
  if not os.path.exists(dest_folder + "/" + i['pdbid'] + ".pdb"):
    res = ''.join(c for c in i['resolution'] if (c.isdigit() or c =='.'))
    if res != "":
      if i['pdbid'] not in pdbs_not_exist and float(res) < 2.5:
        wget.download(url=url_for_pdbs + i['pdbid'] + ".pdb", out=dest_folder)

# Closing file
f.close()