from concurrent.futures import ProcessPoolExecutor
from functools import partial
from Bio import pairwise2
from Bio.Align import substitution_matrices
import pandas as pd
import re
from tqdm import tqdm
import swifter 
# swifter.config.set_defaults(allow_dask_on_strings=True, progress_bar=True)


def sanitize_sequence_advanced(sequence: str) -> str:
    if not isinstance(sequence, str):
        return "" 

    ptm_replacements = {
        "(MSE)": "M",  # Selenomethionine -> Methionine
        "(SEP)": "S",  # Phosphoserine -> Serine
        "(TPO)": "T",  # Phosphothreonine -> Threonine
        "(PTR)": "Y",  # Phosphotyrosine -> Tyrosine
        "(NEP)": "K",  # N-Epsilon-Phospholysine -> Lysine
        "(MLY)": "K",  # Monomethyllysine -> Lysine
        "(M2L)": "K",  # Dimethyllysine -> Lysine
        "(M3L)": "K",  # Trimethyllysine -> Lysine
        "(ALY)": "K",  # Acetyllysine -> Lysine
        "(HLY)": "K",  # Hydroxylysine -> Lysine
        "(M1G)": "R",  # Monomethylarginine -> Arginine
        "(M2G)": "R",  # Dimethylarginine -> Arginine
        "(CIR)": "R",  # Citrulline -> Arginine
        "(HYP)": "P",  # Hydroxyproline -> Proline
        "(CGU)": "E",  # Gamma-carboxyglutamate -> Glutamate
        "(NH2)": "",   # C-Terminal Amidation -> Remove
        "(ACE)": "",   # N-Acetyl Group -> Remove
    }

    processed_seq = sequence
    for mod_code, standard_aa in ptm_replacements.items():
        processed_seq = processed_seq.replace(mod_code, standard_aa)
    valid_chars = "ACDEFGHIKLMNPQRSTVWY"
    sanitized_seq = re.sub(f"[^{valid_chars}]", "X", processed_seq.upper())
    return sanitized_seq

def create_multi_class_mask(uniprot_sequence: str, pdb_construct_sequence: str, blosum62) -> list[int] | None:
    alignments = pairwise2.align.globalds(
        uniprot_sequence,
        pdb_construct_sequence,
        blosum62,
        -10,
        -0.5,
        penalize_end_gaps=False  # FIX: Added required argument
    )

    if not alignments:
        return None

    aligned_uniprot, aligned_pdb, score, begin, end = alignments[0]
    modification_mask = []

    for uniprot_char, pdb_char in zip(aligned_uniprot, aligned_pdb):
        if uniprot_char == '-':
            continue
        if pdb_char == '-':
            modification_mask.append(1)
        elif uniprot_char == pdb_char:
            modification_mask.append(0)
        else:
            modification_mask.append(2)

    if len(modification_mask) != len(uniprot_sequence):
        return None

    return modification_mask


def process_group(group_df, blosum62):
    ret = {}
    for i, row in group_df.iterrows():
        ret[i] = create_multi_class_mask(row['uniprot_seq'], row['pdb_sequence_sanitized'], blosum62)
    return ret

def process_groups(groups, blosum62):
    ret = {}
    for group_df in groups:
        op = process_group(group_df, blosum62)
        ret.update(op)
    return ret


def partition_dataframe(df, chunk_size=100, list_size=5):
    all_chunked_dfs = []
    num_chunks = (len(df) + chunk_size - 1) // chunk_size

    for i in range(num_chunks):
        start_index = i * chunk_size
        end_index = min((i + 1) * chunk_size, len(df))
        all_chunked_dfs.append(df.iloc[start_index:end_index])

    final_partitions = []
    num_sublists = (len(all_chunked_dfs) + list_size - 1) // list_size

    for i in range(num_sublists):
        start_index_sublist = i * list_size
        end_index_sublist = min((i + 1) * list_size, len(all_chunked_dfs))
        final_partitions.append(all_chunked_dfs[start_index_sublist:end_index_sublist])

    return final_partitions


def process_data_cc(df, blosum62_mat):
    gps = partition_dataframe(df)
    ops = {}
    with ProcessPoolExecutor(10) as executor:
        results = executor.map(partial(process_groups, blosum62=blosum62_mat), gps)
        for result in results:
            ops.update(result)
            print("MADE: ", len(ops.keys()))
    return ops

# def process_data_cc(df, blosum62_mat):
#     gps = partition_dataframe(df)
#     ops = {}

#     total_groups = len(gps)
#     print(f"[process_data_cc] Total partitions to process: {total_groups}")

#     print("[process_data_cc] Starting multiprocessing")
#     with ProcessPoolExecutor(10) as executor:
#         # Wrap results in tqdm for progress bar
#         results = executor.map(partial(process_groups, blosum62=blosum62_mat), gps)

#         for idx, result in enumerate(tqdm(results, total=total_groups, desc="Processing groups", unit="group")):
#             ops.update(result)

#     print("[process_data_cc] All processing complete")
#     return ops



if __name__ == "__main__":
    df = pd.read_excel('/Users/haripat/Desktop/SF/protein/data/protein_constructs.xlsx')
    df['pdb_sequence_sanitized'] = df['pbd_id'].apply(sanitize_sequence_advanced)
    df['pdb_sequence_sanitized'] = df['pdb_sequence_sanitized'].str.replace("U", "C")
    df['uniprot_seq'] = df['uniprot_seq'].str.replace("U", "C")
    blosum62_matrix = substitution_matrices.load("BLOSUM62")

    # def mask_wrapper(row):
    #     return create_multi_class_mask(row['uniprot_seq'], row['pdb_sequence_sanitized'], blosum62_matrix)

    print("STARTING PROCESSING")
    df['label_mask'] = df.swifter.progress_bar(True).allow_dask_on_strings(True).apply(
        lambda row: create_multi_class_mask(row['uniprot_seq'], row['pdb_sequence_sanitized'], blosum62_matrix),
        axis=1
    )
    # df['label_mask'] = df.swifter.apply(lambda row: create_multi_class_mask(row['uniprot_seq'], row['pdb_sequence_sanitized'], blosum62_matrix), axis=1)
    # with ProcessPoolExecutor() as executor:
    #     df['label_mask'] = list(executor.map(mask_wrapper, [row for _, row in df.iterrows()]))
    # result = process_data_cc(df, blosum62_mat=blosum62_matrix)

    # df['label_mask'] = float('nan')

    # for row_index, value in result.items():
    #     if row_index in df.index:
    #         df.loc[row_index, 'label_mask'] = value

    # df.to_excel("/Users/haripat/Desktop/SF/protein/data/protein_constructs_w_label_masks.xlsx")