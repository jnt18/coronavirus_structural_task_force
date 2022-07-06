#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 14:07:29 2021

This contains functions used to analyze the data frame or to fix bugs.

@author: Maximilian Edich
"""

import os
import pandas as pd

db_string_1 = "main_repo_database_SARS-CoV.pkl"
df_SC1 = pd.read_pickle(db_string_1)
db_string_2 = "main_repo_database_SARS-CoV-2.pkl"
df_SC2 = pd.read_pickle(db_string_2)



# # # Remove duplicates
def scan_for_duplicates(taxonomy):
    """
    Look for duplicates but do not modify the database.
    Returns: Number of relevant cases.
    """
    # specify used data frame
    if taxonomy == "SARS-CoV-2":
        df = df_SC2
    else:
        df = df_SC1
    
    duplicates = []
    for entry in df['pdb_id']:
        if len(df.loc[df['pdb_id'] == entry]) != 1:
            if not (entry in duplicates):
                duplicates.append(entry)
    print("List of pdb_ids with duplicates in data base:")
    print(duplicates)
    print("Number of duplicates: " + str(len(duplicates)))
    print("Total number of database entries: " + str(len(df)))
    
    return len(duplicates)
    

def remove_duplicates(taxonomy):
    """
    Look for duplicates, drop them from the database and apply the changes.
    Returns: None
    """
    # specify used data frame
    if taxonomy == "SARS-CoV-2":
        df = df_SC2
        db_string = db_string_2
    else:
        df = df_SC1
        db_string = db_string_1
    
    # locate duplicates
    duplicates = []
    for entry in df['pdb_id']:
        if len(df.loc[df['pdb_id'] == entry]) != 1:
            if not (entry in duplicates):
                duplicates.append(entry)
    print(duplicates)
    print(len(duplicates))
    
    print(len(df))
    
    # drop duplicates
    nf = df.drop_duplicates()
    
    for dupe in duplicates:
        prot = nf.loc[nf['pdb_id'] == dupe]
        print(prot)
    
    print(len(nf))
    nf.to_pickle(db_string)


# # # Change path of not assigned
def scan_for_not_assigned(taxonomy):
    """
    Look for entries containing 'not_assigned' in their path.
    Does not apply any changes to the database.
    Returns: Number of relevant cases.
    """
    # specify used data frame
    if taxonomy == "SARS-CoV-2":
        df = df_SC2
    else:
        df = df_SC1
        
    # find entries with 'not_assigned' in their path
    print("Entries containing 'not_assigned':")
    not_assigned = []
    for entry in df['path_in_repo']:
        if str(entry).find('not_assigned') >= 0:
            prot = df.loc[df['path_in_repo'] == entry, 'protein'].iloc[0]
            pdb_id = df.loc[df['path_in_repo'] == entry, 'pdb_id'].iloc[0]
            print(entry)
            not_assigned.append(entry)
    print("Number of such cases: " + str(len(not_assigned)))
    
    return len(not_assigned)


def update_path_of_not_assigned(taxonomy):
    """
    Look for entries containing 'not_assigned' in their path.
    Fix their path, if they are in fact assigned and apply changes to
    database.
    Returns: None
    """
    # specify used data frame
    if taxonomy == "SARS-CoV-2":
        df = df_SC2
        db_string = db_string_2
    else:
        df = df_SC1
        db_string = db_string_1
        
    # find entries with 'not_assigned' in their path
    not_assigned = []
    for entry in df['path_in_repo']:
        if str(entry).find('not_assigned') >= 0:
            prot = df.loc[df['path_in_repo'] == entry, 'protein'].iloc[0]
            pdb_id = df.loc[df['path_in_repo'] == entry, 'pdb_id'].iloc[0]
            print("entry:" + entry)
            print("prot:" + prot)
            
            if prot == "not_assigned":
                # pdb has to be assigned manually
                not_assigned.append(pdb_id)
                print(prot)
            else:
                # update path
                new_path = os.path.join("pdb", prot, taxonomy, pdb_id)
                print(new_path)
                df.loc[df["pdb_id"] == pdb_id,"path_in_repo"] = new_path
                print(df.loc[df["pdb_id"] == pdb_id,"path_in_repo"])
            print()
            
    print(len(not_assigned))
    print(not_assigned)
    print("these pdbs have to be assigned manually")
    df.to_pickle(db_string)
    

def assign_protein_to_data_frame_entry(pdb_id, protein, taxonomy):
    """
    Assign a protein to a certain database entry and save changes.
    Used for not_assigned proteins which have the files in correct place.
    Returns: None
    """
    # specify used data frame
    if taxonomy == "SARS-CoV-2":
        df = df_SC2
        db_string = db_string_2
    else:
        df = df_SC1
        db_string = db_string_1
        
    df.loc[df["pdb_id"] == pdb_id, "protein"] = protein
    df.to_pickle(db_string)

def delete_entry_from_database(pdb_id, taxonomy):
    """
    Use this function to delete entries which do not correspong to SARS-CoV
    or SARS-CoV-2.
    One example is 7f8l, which was somehow classified as SARS-CoV
    """
    # specify used data frame
    if taxonomy == "SARS-CoV-2":
        df = df_SC2
        db_string = db_string_2
    else:
        df = df_SC1
        db_string = db_string_1
        
    print("Total number of database entries: " + str(len(df)))
    # load old entry and drop it from old dataframe
    entry = df.loc[df['pdb_id'] == pdb_id]
    df_new = df.drop(entry.index)
    print("Total number of database entries: " + str(len(df_new)))
    df_new.to_pickle(db_string)
    
def print_entry(pdb_id, taxonomy):
    """
    Prints the entry with the given id.
    """
    # specify used data frame
    if taxonomy == "SARS-CoV-2":
        df = df_SC2
    else:
        df = df_SC1
        
    entry = df.loc[df['pdb_id'] == pdb_id]
    print(entry['protein'])


def run(taxonomy):
    print_entry('3i6k', 'SARS-CoV')
    assign_protein_to_data_frame_entry('3i6k', 'membrane_glycoprotein', taxonomy)
    print("Perform scans for errors in database.")
    errors = scan_for_duplicates(taxonomy)
    errors += scan_for_not_assigned(taxonomy)
    if errors == 0:
        print("Scan was succesful, no errors found!\n")
    else:
        print("Scan finished. Scan revealed " + str(errors) + " errors!\n")
        print("look into code of 'analyze_and_fix_dataframe.py' to perform manual "
              + "fixing with the respective functions!")
    # for fixing the errors, call the respective functions from above

run("SARS-CoV")