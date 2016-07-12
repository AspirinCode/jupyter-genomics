# third-party libraries
import pandas

# ccbb libraries
from ccbbucsd.utilities.bio_seq_utilities import trim_seq

__author__ = "Amanda Birmingham"
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"

_CONSTRUCT_ID = "CONSTRUCT_ID"
_GENE_A_SEQ = "GENE_A_SEQ"
_GENE_B_SEQ = "GENE_B_SEQ"
_CONSTRUCT_A_NAME = "CONSTRUCT_A_NAME"
_CONSTRUCT_B_NAME = "CONSTRUCT_B_NAME"


def get_construct_separator():
    return "__"


def extract_construct_and_grna_info(constructs_fp, column_indices):
    construct_table = _read_in_construct_table(constructs_fp, column_indices, rows_to_skip=1)
    seq_name_sets = _extract_grnas_from_construct_table(construct_table)
    grna_name_seq_pairs = _format_and_check_grnas_input(seq_name_sets)
    construct_names = construct_table[_CONSTRUCT_ID].unique().tolist()
    return construct_names, grna_name_seq_pairs


def _read_in_construct_table(constructs_fp, column_indices, rows_to_skip=1):
    result = pandas.read_table(constructs_fp, skiprows=rows_to_skip, header=None)
    result = _rename_columns(result, column_indices)
    return result


def _rename_columns(construct_table, column_indices):
    new_names = [_CONSTRUCT_ID, _GENE_A_SEQ, _GENE_B_SEQ]
    existing_names = list(construct_table.columns.values)

    existing_to_new_names = {}
    for curr_index in range(0, len(column_indices)):
        curr_col_index = column_indices[curr_index]
        curr_existing_name = existing_names[curr_col_index]
        existing_to_new_names[curr_existing_name] = new_names[curr_index]

    return construct_table.rename(columns=existing_to_new_names)


def _extract_grnas_from_construct_table(construct_table):
    grna_name_key = "grna_name"
    grna_seq_key = "grna_seq"

    # split the construct id in each row into two pieces--the two construct names--and put them into new columns
    construct_table[_CONSTRUCT_A_NAME], construct_table[_CONSTRUCT_B_NAME] = \
        zip(*construct_table[_CONSTRUCT_ID].str.split(get_construct_separator()).tolist())

    # get the gene/sequence pairs for each of the two genes and concatenate them
    gene_a_pairs = _extract_grna_name_and_seq_df(construct_table, _CONSTRUCT_A_NAME,
                                                 _GENE_A_SEQ, grna_name_key, grna_seq_key)
    gene_b_pairs = _extract_grna_name_and_seq_df(construct_table, _CONSTRUCT_B_NAME,
                                                 _GENE_B_SEQ, grna_name_key, grna_seq_key)
    gene_pairs = pandas.concat([gene_a_pairs, gene_b_pairs])
    gene_pairs[grna_seq_key] = gene_pairs[grna_seq_key].str.upper()  # upper-case all gRNA sequences

    # extract only the unique pairs
    grna_seq_name_groups = gene_pairs.groupby([grna_seq_key, grna_name_key]).groups
    result = [x for x in grna_seq_name_groups]
    return sorted(result)  # NB sort so that output order is predictable


def _extract_grna_name_and_seq_df(construct_table, name_key, seq_key, grna_name_key, grna_seq_key):
    name_and_seq_df = construct_table[[name_key, seq_key]]
    name_and_seq_df.rename(columns={name_key: grna_name_key, seq_key: grna_seq_key}, inplace=True)
    return name_and_seq_df


def trim_grnas(grnas_name_and_seq_list, retain_len):
    result = []
    for name_seq_tuple in grnas_name_and_seq_list:
        grna_name = name_seq_tuple[0]
        full_seq = name_seq_tuple[1]
        trimmed_seq = trim_seq(full_seq, retain_len, False)  # False = do not retain from 5p end but from 3p end
        result.append((grna_name, trimmed_seq))
    return result


def _read_in_grnas(grnas_fp):
    list_from_file = _read_grnas_input(grnas_fp)
    return _format_and_check_grnas_input(list_from_file)


def _read_grnas_input(grnas_fp, comment_prefix="#", delimiter="\t"):
    result = []
    with open(grnas_fp, 'r') as file_handle:
        for line in file_handle:
            if not line.startswith(comment_prefix):
                pieces = line.strip().split(delimiter)
                result.append(pieces)
    return result


def _format_and_check_grnas_input(grnas_seq_and_name_list):
    expected_num_pieces = 2
    seqs_by_names = {}
    names_by_seqs = {}
    result = []

    for curr_set in grnas_seq_and_name_list:
        if len(curr_set) != expected_num_pieces:
            raise ValueError(
                "input '{0}' has {1} pieces instead of the expected {2}".format(
                    curr_set, len(curr_set), expected_num_pieces
                ))
        curr_seq = curr_set[0]
        curr_name = curr_set[1]

        if curr_seq in names_by_seqs:
            raise ValueError(
                "sequence '{0}' associated with name '{1}' but was already associated with name '{2}'".format(
                    curr_seq, curr_name, names_by_seqs[curr_seq]
                ))

        if curr_name in seqs_by_names:
            raise ValueError(
                "name '{0}' associated with sequence '{1}' but was already associated with sequence '{2}'".format(
                    curr_name, curr_seq, seqs_by_names[curr_name]
                ))

        names_by_seqs[curr_seq] = curr_name
        seqs_by_names[curr_name] = curr_seq
        result.append((curr_name, curr_seq))
    # next pair in

    return result
