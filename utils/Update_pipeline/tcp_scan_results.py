#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 14:36:05 2022
This module iterates over all new structures and checks their assignments
on errors by checking also the title, taxonomy, and so on.

@author: localadmin
"""


def main(taxnonomy, c_new_pdb_lst):
    """
    Parameters
    ----------
    taxnonomy : String
        Value is either 'SARS-CoV' or 'SARS-CoV-2'.
    c_new_pdb_lst : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    print("scan!")
    print(c_new_pdb_lst)
    print("scan off.")
    
    # iterate over each entry and keept track of warnings
    warnings = []
    for entry in c_new_pdb_lst:
        # TODO get titel from PDB and look for taxonomy keyword
        # TODO if title taxonomy and real taxonomy do not agree, store warning
        pass
    
    # print out warnings
    print("Warnings in total: " + str(len(warnings)))
    for warn in warnings:
        print(warn)
    
    return
