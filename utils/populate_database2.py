import csv
import os
import sqlite3
from gemmi import cif
import json

# build connection to SQL data base
conn = sqlite3.connect("stats.db")
cursor = conn.cursor()

# code to initialize database
'''
cursor.execute("""CREATE TABLE stats (pdbid text,
                                    datapath text,
                                    github text,
                                    protein text,
                                    virus text,
                                    method text,
                                    hasRerefinement int,
                                    resolution real,
                                    rmsd real,
                                    rwork real,
                                    rfree real)""")
'''

conn.commit()
url = 'https://github.com/thorn-lab/coronavirus_structural_task_force/tree/master'


def fill_database(workdir):
    """
    Collect relevant data of each pdb entry from local working directory and write this information into
    an SQL database.
    A fixed folder structure is required to find the  data where it is expected.

    :param workdir: The folder named 'pdb', containing sub folders of following structure:
    protein > virues > pdb_identifiers > pdb and cif files + validation folder
    :return: None.
    """

    # iterate through all folders
    for path, dirs, files in os.walk(workdir):
        # detect SARS-CoV and SARS-CoV-2 folders, since these contain sub-folders for each pdb id
        if path.endswith("SARS-CoV") or path.endswith("SARS-CoV-2"):

            # iterate through sub-folders named by pdb id
            for pdb_id in dirs:
                # build path to sub-folders
                folder = os.path.join(path, pdb_id)

                # extract virus name and protein name for this entry
                local_path_split = os.path.split(path)
                virus_name = local_path_split[1]
                protein_name = os.path.split(local_path_split[0])[1]

                # check if entry contains a folder, which indicates a re-refinement
                if (os.path.exists(os.path.join(folder, "isolde"))
                        or os.path.exists(os.path.join(folder, "refmac"))
                        or os.path.exists(os.path.join(folder, "ccpem"))
                        or os.path.exists(os.path.join(folder, "rerefinements"))
                        or os.path.exists(os.path.join(folder, "reintegration"))):
                    has_re_refinement = 1
                else:
                    has_re_refinement = 0

                # read cif file values
                try:
                    # extract block values from cif file and extract used method from those
                    cif_file = cif.read(os.path.join(folder, pdb_id + '.cif'))
                    block = cif_file.sole_block()
                    block_values = block.find_values('_exptl.method')
                    method = block_values.str(0)
                except Exception as e:
                    # no cif file found or readable
                    print("Error during reading out values from cif in " + folder)
                    print(e)
                    block = None
                    method = None

                # open readme file
                description = ""
                if os.path.exists(os.path.join(folder, 'README.md')):
                    try:
                        readme = open(os.path.join(folder, 'README.md')).read()
                        # extract description from readme file
                        description_lines = False
                        for line in readme.splitlines():
                            if line.startswith("## Basefolder"):
                                break
                            if description_lines:
                                if not line == "":
                                    description = '' + line
                            if line == "## Description":
                                description_lines = True
                    except UnicodeDecodeError as e:
                        print("Error during README reading in " + folder + ":")
                        print(e)

                # if molprobity data exists, read this data
                molprobity_score = None
                if os.path.exists(os.path.join(folder, 'validation', 'molprobity', 'molprobity.out')):
                    molprobity = open(os.path.join(folder, 'validation', 'molprobity', 'molprobity.out')).read()
                    molprobity = molprobity.splitlines()
                    molprobity_score_line = molprobity[-5]
                    molprobity_score = ''.join((ch if ch in '0123456789.' else ' ') for ch in molprobity_score_line)
                    molprobity_score = float(molprobity_score)

                # check used methods and set values according to method
                rfree = None
                rwork = None
                rmsd = None
                resolution = None
                if block is None:
                    # if it was not possible to read out block values from cif, then it was not possible to get
                    # the used method as well
                    pass
                elif method == 'X-RAY DIFFRACTION':
                    if block.find_value('_refine.ls_R_factor_R_free') is not None:
                        rfree = float(block.find_value('_refine.ls_R_factor_R_free'))
                    if block.find_value('_refine.ls_R_factor_R_work') is not None:
                        rwork = float(block.find_value('_refine.ls_R_factor_R_work'))
                    rmsd = None
                    # fsc = mmcif_dict['_entry_for_fsc']
                    try:
                        resolution = float(block.find_value('_reflns.d_resolution_high'))
                    except Exception as e:
                        resolution = None
                elif method == 'ELECTRON MICROSCOPY':
                    rfree = None
                    rwork = None
                    rmsd = None
                    # fsc = mmcif_dict['_entry_for_fsc']
                    resolution = block.find_value('_em_3d_reconstruction.resolution')
                elif method == 'SOLUTION NMR':
                    rfree = None
                    rwork = None
                    rmsd = block.find_value('pdbx_nmr_ensemble_rms')
                    # fsc = mmcif_dict['_entry_for_fsc']
                    resolution = None
                elif method == 'SOLID-STATE NMR':
                    pass
                else:
                    print("Error during value read out based on method in " + pdb_id)
                    print("method: " + method)
                
                # build relative path starting at 'pdb/'
                norm_path = os.path.normpath(folder)
                rel_path_list = norm_path.split(os.sep)
                # find index of 'pdb'
                pdb_index = 0
                while pdb_index < len(rel_path_list):
                    if rel_path_list[pdb_index] == 'pdb':
                        break
                    pdb_index += 1
                # build path from list
                rel_path = rel_path_list[pdb_index]
                for path_part in rel_path_list[(pdb_index + 1):]:
                    rel_path = os.path.join(rel_path, path_part)
                print(rel_path)
                
                # build github link
                github_link = url + "/" + rel_path
                github_link = github_link.replace('\\', '/')

                # parse all collected values
                parsed_values = (pdb_id, rel_path, github_link, protein_name, virus_name, method, has_re_refinement,
                                 resolution, rmsd, rwork, rfree)
                EM_values = (pdb_id, resolution)
                general_values = (pdb_id, description, method, rel_path, github_link, protein_name, virus_name,
                                  has_re_refinement, molprobity_score)
                MX_values = (pdb_id, rfree, resolution, rwork)
                NMR = (pdb_id, rmsd)

                # ensure safe working with file stream in case of exception
                with conn:
                    # execute SQL commands to fill in data into database
                    cursor.execute('SELECT ID FROM General WHERE ID = ?', [pdb_id])
                    data = cursor.fetchall()
                    if len(data) == 0:
                        cursor.execute('INSERT INTO stats VALUES(?,?,?,?,?,?,?,?,?,?,?)', parsed_values)
                        cursor.execute('INSERT INTO EM VALUES(?,?)', EM_values)
                        cursor.execute('INSERT INTO General VALUES(?,?,?,?,?,?,?,?,?)', general_values)
                        cursor.execute('INSERT INTO MX VALUES(?,?,?,?)', MX_values)
                        cursor.execute('INSERT INTO NMR VALUES(?,?)', NMR)
                    else:
                        update_id = (pdb_id,)
                        cursor.execute(
                            'UPDATE stats Set pdbid = ?, datapath = ?, github = ?, protein = ?, virus = ?, method = ?, hasRerefinement = ?, resolution = ?, rmsd = ?, rwork = ?, rfree = ? WHERE pdbid = ?',
                            (parsed_values + update_id))
                        cursor.execute('UPDATE EM Set ID = ?, Resolution = ? WHERE ID = ?', (EM_values + update_id))
                        cursor.execute(
                            'UPDATE General Set ID = ?, Description = ?, Method = ?, datapath = ?, github = ?, protein = ?, virus = ?, hasRerefinement = ?, Molprobity_score = ? WHERE ID = ?',
                            (general_values + update_id))
                        cursor.execute('UPDATE MX Set ID = ?, rfree = ?, Resolution = ?, Rwork = ? WHERE ID = ?',
                                       (MX_values + update_id))
                        if NMR[1] is not None:
                            cursor.execute('UPDATE NMR Set ID = ?, RMSD = ? WHERE ID = ?', (NMR + update_id))


# execute main function
pwd = os.path.join(os.getcwd(), os.pardir, "pdb")
pwd = os.path.abspath(pwd)
print(pwd)
fill_database(pwd)

# create excel textfiles
outlist = open('mxList.txt', 'w')
print('path, Rfree', file=outlist)
with conn:
    for row in cursor.execute('SELECT * FROM stats WHERE method=? ORDER BY protein, rfree', ('X-RAY DIFFRACTION',)):
        print(row[1], row[-1], file=outlist)
print('Done!')

emlist = open('emList.txt', 'w')
print('path, resolution', file=emlist)
with conn:
    for row in cursor.execute('SELECT * FROM stats WHERE method=? ORDER BY protein, resolution',
                              ('ELECTRON MICROSCOPY',)):
        print(row[1], row[-4], file=emlist)
print('Done!')

# create csv files
with open('stats.csv', 'w', newline='') as csvfile:
    full_db = csv.writer(csvfile, dialect='excel')
    full_db.writerow(['pdb', 'path', 'url', 'protein', 'virus', 'method', 'resolution', 'rmsd', 'rwork', 'rfree'])
    reader = csv.DictReader(csvfile)
    with conn:
        for row in cursor.execute('SELECT * FROM stats ORDER BY protein'):
            full_db.writerow(row)
print('Done!')

with open('EM.csv', 'w', newline='') as csvfile:
    full_db = csv.writer(csvfile, dialect='excel')
    full_db.writerow(['pdb', 'resolution'])
    with conn:
        for row in cursor.execute('SELECT * FROM EM ORDER BY resolution'):
            full_db.writerow(row)
print('Done!')

with open('General.csv', 'w', newline='') as csvfile:
    full_db = csv.writer(csvfile, dialect='excel')
    full_db.writerow(['pdb', 'description', 'method', 'datapath', 'github', 'protein', 'virus', 'hasRerefinement',
                      'Molprobity_score'])
    with conn:
        for row in cursor.execute('SELECT * FROM General ORDER BY method'):
            full_db.writerow(row)
print('Done!')

with open('MX.csv', 'w', newline='') as csvfile:
    full_db = csv.writer(csvfile, dialect='excel')
    full_db.writerow(['pdb', 'rfree', 'rwork', 'resolution'])
    with conn:
        for row in cursor.execute('SELECT * FROM MX ORDER BY resolution'):
            full_db.writerow(row)
print('Done!')

with open('NMR.csv', 'w', newline='') as csvfile:
    full_db = csv.writer(csvfile, dialect='excel')
    full_db.writerow(['pdb', 'rmsd'])
    with conn:
        for row in cursor.execute('SELECT * FROM NMR ORDER BY rmsd'):
            full_db.writerow(row)
print('Done!')
