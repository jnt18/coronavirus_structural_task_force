from pathlib import Path
import pandas as pd
import gemmi

from lib.update_package import utils


def sequence_alignment(
    df: pd.DataFrame,
    start: str,
    end: str,
    repo_path: str | Path,
    fasta_path: str | Path,
):

    new_df = utils.get_current_df(df, start, end)
    fasta_sequences = utils.get_protein_seq_from_fasta(fasta_path)

    for _, group in new_df.groupby(["protein", "taxonomy"]):

        align_path = (
            Path(repo_path, group.path_in_repo.iloc[0]).parent
            / "structure_sequence_alignment.txt"
        )

        if not align_path.exists():
            header = (
                "This is the alignment of sequence in the pdb_file and the reference genome.\n"
                "For the alignment the python tool gemmi "
                "[https://gemmi.readthedocs.io/en/latest/index.html] was used.\n\n"
            )
            align_path.write_text(header)

        for id, row in group.iterrows():

            with open(align_path, "a") as doc:
                doc.write(f">>>{id}:\n")

                id_folder = Path(repo_path, row.path_in_repo, id)
                pdb_path = f"{id_folder}.pdb"
                cif_path = f"{id_folder}.cif"
                pdb = gemmi.read_pdb(pdb_path).entities
                cif = gemmi.read_structure(cif_path).entities

                if len(pdb) == 0:
                    print(f"No sequences found for {id}")
                    continue

                try:
                    chains = [
                        ord(chain) - 65 for chain in row.relevant_chains.split(", ")
                    ]
                except AttributeError:
                    chains = range(len(pdb))

                chain_sequences = extract_sequences(pdb, chains)

                protein_names = row.protein.split("-")

                remove_rna_sequences(chains, chain_sequences, protein_names)

                if len(protein_names) == len(chains):
                    candidates_per_chain = [[name] for name in protein_names]
                else:
                    # First and last chains are anchored to their respective names;
                    # middle chains are ambiguous and must be matched against all names.
                    candidates_per_chain = (
                        [[protein_names[0]]]
                        + [protein_names] * (len(chain_sequences) - 2)
                        + [[protein_names[-1]]]
                    )

                for chain, chain_seq, candidates in zip(
                    chains, chain_sequences, candidates_per_chain
                ):
                    result, fasta, name = best_alignment(
                        chain_seq, candidates, fasta_sequences
                    )
                    if not result:
                        print(f"No reference sequences for {candidates}")
                        continue
                    doc.write(
                        f">Reference genome ({name}) against deposited genome (chain {chr(chain+65)}):\n"
                    )
                    write_seq_report(doc, result, fasta, chain_seq)

                # Align deposited genome to structure sequence
                pdb_seq = "".join(extract_sequences(pdb, range(len(pdb))))
                cif_seq = "".join(extract_sequences(cif, range(len(cif))))

                result = gemmi.align_string_sequences(list(cif_seq), list(pdb_seq), [])

                doc.write(">Deposited genome against structure sequence:\n")

                write_seq_report(doc, result, cif_seq, pdb_seq)

                mismatches, gaps = analyze_mismatches(result)
                doc.write(f"Mismatch in {mismatches}\n")
                doc.write(f"Unmodelled in {gaps}\n\n")


def extract_sequences(entities, chains):
    """Return list of sequences for specified chains."""

    def get_string_seq(entities, chain_idx: int) -> str:
        """Return the one-letter sequence for a chain, replacing unknowns with X."""
        sequence = entities[chain_idx].full_sequence
        letters = (
            gemmi.find_tabulated_residue(res).one_letter_code for res in sequence
        )
        return "".join(code if code.isupper() else "X" for code in letters)

    return [get_string_seq(entities, int(chain)) for chain in chains]


def remove_rna_sequences(chains, chain_sequences, protein_names):
    try:  # No sequence alignment for rna sequences
        protein_names.remove("rna")
    except:
        ValueError
    for chain, chain_seq in zip(chains, chain_sequences):
        if set(chain_seq).issubset({"A", "C", "G", "U", "I", "N"}):
            chains.remove(chain)
            chain_sequences.remove(chain_seq)


def write_seq_report(doc, result, seq1: str, seq2: str):
    """Write alignment statistics and formatted alignment."""
    doc.write(f"Identity: {result.calculate_identity():.2f}\n")
    doc.write(f"Score: {result.score:.2f}\n")
    doc.write(f"Alignment:\n{result.formatted(seq1, seq2)}\n")


def best_alignment(chain_seq, protein_names, fasta_sequences):
    """Find the best alignment of a chain against reference FASTA sequences."""
    best_score = -float("inf")
    best_result = best_fasta = best_name = None

    for name in protein_names:
        fastas = [v for k, v in fasta_sequences if k == name]

        for fasta in fastas:
            result = gemmi.align_string_sequences(list(fasta), list(chain_seq), [])

            if result.score > best_score:
                best_score = result.score
                best_result = result
                best_fasta = fasta
                best_name = name

    return best_result, best_fasta, best_name


def analyze_mismatches(result):
    """Return mismatch and gap positions from match string."""
    mismatches = []
    gaps = []

    for i, letter in enumerate(result.match_string):
        pos = i + 1
        if letter == "|":
            continue
        elif letter == ".":
            mismatches.append(pos)
        else:
            gaps.append(pos)

    return mismatches, gaps
