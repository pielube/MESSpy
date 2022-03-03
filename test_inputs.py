import os
import json

"""
Reading input files
"""

path = r'./inputs'

file = 'structure.json'
filepath = os.path.join(path,file)

with open(filepath,'r') as f:
    structure = json.load(f)

file = 'economics.json'
filepath = os.path.join(path,file)

with open(filepath,'r') as f:
    economics = json.load(f)

file = 'general.json'
filepath = os.path.join(path,file)

with open(filepath,'r') as f:
    general = json.load(f)

file = 'refcase.json'
filepath = os.path.join(path,file)

with open(filepath,'r') as f:
    refcase = json.load(f)