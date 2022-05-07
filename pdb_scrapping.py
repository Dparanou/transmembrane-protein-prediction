import json
import wget
import os
import os.path
from os import path
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

url_for_pdbs = "https://files.rcsb.org/download/"

# check if data folder exists
if not path.isdir('data'):
  os.makedirs('data')

# save the path of data folder
dest_folder = os.path.realpath("./data")

# Opening JSON file
f = open('primary_structures.json')

# returns JSON object as
# a dictionary
data = json.load(f)

# Iterating through the json and download pdb files
for i in data['objects']:
  wget.download(url=url_for_pdbs + i['pdbid'] + ".pdb", out=dest_folder)

# Closing file
f.close()