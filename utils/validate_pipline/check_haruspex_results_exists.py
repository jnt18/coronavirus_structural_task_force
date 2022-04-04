#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script can be modified to detect any interesting information from a 
README.md in the pdb folders.
So far it only checks wether such a README.md exists and if so, it checks
for the method == Cryo EM and if Haruspex results are available.

Created on Mon Jun 28 17:04:23 2021

@author: Maximilian Edich
"""

look_method = "**Method**: Cryo-EM"
look_info = "## HARUSPEX results\n\nNot available"



import os

# check if script is located in utils
if os.path.basename(os.getcwd()) != "utils":
    exit("ERROR: place and launch the script from utils folder!")


pdb_dir = os.path.abspath(os.path.join(os.getcwd(), "..", "pdb"))

# init list gathering useful information from readmes
no_read_me_prots = []
hits = []

# iterate down the folder tree until correct README.md is found
for protein in os.listdir(pdb_dir):
    protein_dir = os.path.join(pdb_dir, protein)
    if not os.path.isdir(protein_dir):
        continue
    for taxonomy in os.listdir(protein_dir):
        taxonomy_dir = os.path.join(protein_dir, taxonomy)
        if not os.path.isdir(taxonomy_dir):
            continue
        for pdb_id in os.listdir(taxonomy_dir):
            if not os.path.isdir(os.path.join(taxonomy_dir, pdb_id)):
                continue
            # look into README and check if it is from cryo EM
            read_me_path = os.path.join(taxonomy_dir, pdb_id, "README.md")
            try:
                file = open(read_me_path, 'r')
                content = file.read()
                file.close()
            except FileNotFoundError:
                no_read_me_prots.append(os.path.join(protein,
                                                     taxonomy, pdb_id))
                continue
            
            # check if method is the one you are looking for
            if content.find(look_method) > -1:
                if content.find(look_info) > -1:
                    # open pdb file and look for additional info
                    """pdb_path = os.path.join(taxonomy_dir, pdb_id, pdb_id + ".pdb")
                    file = open (pdb_path, 'r')
                    pdb_content = file.read()
                    file.close()
                    emd_number_start = pdb_content.find("EMD-")
                    emd_number = pdb_content[emd_number_start : (emd_number_start + 10)].strip()
                    """
                    #hits.append(os.path.join(protein, taxonomy, pdb_id) + " | " + str(emd_number))
                    hits.append(os.path.join(protein, taxonomy, pdb_id))

print("Missing README.md files in:")
print("length: " + str(len(no_read_me_prots)))
for miss in no_read_me_prots:
    print(miss)

print()
print("Structures obtained from '"
      + look_method + "' and containing \n'" + look_info + "':")
print("length: " + str(len(hits)))
for hit in hits:
    print(hit)
                
            
    