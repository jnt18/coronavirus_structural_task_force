import seaborn as sb
from matplotlib import pyplot as plt
import numpy as np
import gemmi as gm
from pathlib import Path
import pandas as pd
import re

from lib.update_package import utils


def calculate_rmsd(df: pd.DataFrame, start: str, end: str, repo_path: str | Path):

    new_df = utils.get_current_df(df, start, end)
    old_df = df.loc[list(set(df.index) - set(new_df.index)), :]
    for protein, new_df_protein in new_df.groupby("protein"):
        results_df_protein = pd.DataFrame(
            columns=[
                "PDB-1",
                "PDB-2",
                "Chain-1",
                "Chain-2",
                "RMSD",
                "Aligned_Atoms",
            ]
        )
        for id_1, row_1 in new_df_protein.iterrows():
            polymers_1 = get_polymers(id_1, row_1, repo_path)
            for id_2, row_2 in old_df[old_df.protein == protein].iterrows():
                polymers_2 = get_polymers(id_2, row_2, repo_path)
                max_rmsd = 0.0
                for chain_1, polymer_1 in polymers_1.items():
                    for chain_2, polymer_2 in polymers_2.items():
                        ptype = polymer_1.check_polymer_type()
                        sup = gm.calculate_superposition(
                            polymer_1, polymer_2, ptype, gm.SupSelect.CaP
                        )
                        rmsd = sup.rmsd
                        if rmsd >= max_rmsd:
                            max_rmsd, aligned_atoms = rmsd, sup.count
                            best_pair = (chain_1, chain_2)

                results_df_protein.loc[len(results_df_protein), :] = (
                    id_1,
                    id_2,
                    best_pair[0],
                    best_pair[1],
                    round(max_rmsd, 3),
                    aligned_atoms,
                )

            old_df.loc[id_1, :] = row_1

        folder_path = Path(repo_path, new_df_protein.loc[id_1, "path_in_repo"]).parents[
            1
        ]
        save_file(results_df_protein, Path(folder_path, "RMSD.csv"))
        heat_map(folder_path, protein)


def get_polymers(id, row, repo_path):
    structure = gm.read_structure(f"{repo_path}/{row.path_in_repo}/{id}.pdb")[0]
    try:
        chains = [ord(chain) - 65 for chain in row.relevant_chains.split(", ")]
    except AttributeError:
        chains = range(len(structure))
    structure_by_chain = {}
    for chain in chains:
        try:
            structure_by_chain.update({chain: structure[chain].get_polymer()})
        except IndexError:
            print(f"error for {id}, no chain {chain}")
    return structure_by_chain


def save_file(results_df_protein, rmsd_path):
    if not rmsd_path.exists():
        header = (
            f"This document contains the RMSD of the best combination of chains of entries in this folder \n\n"
            + ",".join(results_df_protein.columns)
            + "\n"
        )
        with open(rmsd_path, "w") as doc:
            doc.write(header)

    results_df_protein.to_csv(rmsd_path, header=None, mode="a", index=False)


def heat_map(folder_path, protein):
    rmsd_path = Path(folder_path, "RMSD.csv")
    rmsd_df = pd.read_csv(rmsd_path, header=1)[["PDB-1", "PDB-2", "RMSD"]]
    if len(rmsd_df) > 1:
        hmap = get_distance_matrix(rmsd_df)
        fig, ax = plt.subplots()
        sb.heatmap(
            hmap,
            cmap="viridis",
            cbar=True,
            cbar_kws={"label": "[Å]", "orientation": "vertical"},
            xticklabels=True,
            yticklabels=True,
        )
        ax.set_title(protein + " best RMSD")
        plt.savefig(Path(folder_path, "heatmap.png"), dpi=800)
        plt.close(fig)


def get_distance_matrix(df):
    inv_df = df.copy().rename(columns={"PDB-1": "PDB-2", "PDB-2": "PDB-1"})
    full_rmsd_df = pd.concat([df, inv_df]).drop_duplicates()
    return full_rmsd_df.pivot(index="PDB-1", columns="PDB-2", values="RMSD")
