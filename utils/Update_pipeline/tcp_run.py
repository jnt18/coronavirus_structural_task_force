import tcp_main
import analyze_and_fix_dataframe
import os
import RMSD
import mk_Alignment_strc_vs_seq as align
import argparse
import protein_and_domain_classifier as domain_classifier

"""
executable script for weekly update
"""

# handle input arguments
parser = argparse.ArgumentParser(description='Run weekly update.')
parser.add_argument('-t', '--taxonomy', type=str, required=True, help="Give the taxonomy, either 'SARS-CoV' or 'SARS-CoV-2'")
args = parser.parse_args()

taxo_sars_cov = "SARS-CoV"
taxo_sars_cov_2 = "SARS-CoV-2"

#These are used to only return pdb entries of taxonomy_id and exlcude everything in negate_taxonomy_id
#sars_cov_1_2_id is the id for both CoV-1 and CoV-2 to get only SARS-CoV-1 structures negate sars_cov_2_id
#sars_cov_2_id is the id only for CoV-2 so dont negate sars_cov_1_2_id but use the var: dont_negate (this is the id of an octopus)
sars_cov_1_2_id = "694009" #SARS-CoV-1 & SARS-CoV-2
sars_cov_2_id = "2697049" #SARS-CoV-2
dont_negate = "6645"


# set correct configurations based on input argument
if args.taxonomy == taxo_sars_cov:
    # SARS-CoV
    # taxonomy name used in PDB search
    taxonomy = "Severe acute respiratory syndrome coronavirus"
    taxonomy_id = sars_cov_1_2_id
    negate_taxonomy_id = sars_cov_2_id
elif args.taxonomy == taxo_sars_cov_2:
    # SARS-CoV-2
    # taxonomy name used in PDB search
    taxonomy = "Severe acute respiratory syndrome coronavirus 2"
    taxonomy_id = sars_cov_2_id
    negate_taxonomy_id = dont_negate
else:
    exit("ERROR: illegal taxonomy entered! Maybe a typo?")

# taxonomy name used to name files
taxo = args.taxonomy

# identify repo path
repo_path = os.path.abspath(os.path.join(__file__ ,"..", "..", "..", "pdb"))
reports_path = os.path.abspath(os.path.join(__file__ , "..",
    "weekly_reports", tcp_main.get_time() + "_update_report_" + taxo + ".txt"))
print(reports_path)


print("Searching for new and changed structures")
c_new_pdb_lst, changed_prot_list = tcp_main.main(taxonomy_id=taxonomy_id, negate_taxonomy_id=negate_taxonomy_id, taxo=taxo)
print("Doing sequence aligntment")
align.main(changed_prot_list, c_new_pdb_lst, repo_path, taxo)
print("Calculating RMSD")
RMSD.main(changed_prot_list, repo_path)

# check for common errors in the data base. These are currently:
# duplicate entries
# entries with wrong path containing 'not_assigned' after assignment
# entries not assigned to a protein at all
analyze_and_fix_dataframe.run()

# check if new nsp3 structures are available
try:
    # open file
    file = open(reports_path, 'r')
    content = file.read()
    file.close()
    # check for new nsp3
    if content.find('nsp3') > -1:
        print("New Nsp3 structure: perform domain classification...")
        domain_classifier.main(taxo, 'nsp3')
    
except FileNotFoundError:
    print("weekly reports file not found. No domain classification was made.")
except PermissionError():
    print("Permission denied for opening weekly reports file. No domain classification was made.")
