from itertools import combinations
import re
import glob
import sys
from preprocessing.antismash_fasta import read_fasta, write_fasta

def check_region(sequences: dict, region: str):
    """ Searches the coordinates of the region where the sequence is found.
        Input:
            sequences: dictionary from fasta file, where keys are sample ids and values are sequences as strings
            region: residues
        Output:
            integers (indeces of the start and end charactor of the input residues in the sequence) """
    pattern = re.compile(re.escape(region))

    region_coordinates = {}
    for id, sequence in sequences.items():
        matches = pattern.finditer(sequence)
        for match in matches:
            region_coordinates[id] = [match.start(), match.end()]

    return region_coordinates

def cut_alignment(file_name: str, fasta_dict: dict, positions):
    """ Takes the defined region from the alignment, writes it in a fasta file """
    for name, sequence in fasta_dict:
        fasta_dict[name] = sequence[positions]

    write_fasta(fasta_dict.keys(), fasta_dict.values(), file_name)

def crossval_preproc(fasta_file):
    """ Creates n fasta files with n-1 sample sizes from a fasta files with n samples.
        Input: fasta file
        Output: n fasta file, every file has one missing sample from the original fasta"""
    name_to_seq = read_fasta(fasta_file)

    name_to_seq = [list(element) for element in name_to_seq.items()]
    all_combinations = combinations(name_to_seq, len(name_to_seq)-1)
    all_ids = set(name for name, _ in name_to_seq)
    for combination in all_combinations:
        left_out = list(all_ids - set(id[0] for id in combination))[0].lstrip(">")
        with open(f'crossval_out_{left_out}.fasta', 'w') as crossval_file:
            for name, seq in combination:
                line = ">" + name + '\n' + seq + '\n'
                crossval_file.write(line)

def get_left_out_seq_results(folder, sequence_names):
    result_files = sorted(glob.glob(folder))
    csv_file = "left_out_seq_results.csv"

    sequences = {name.replace("\n","") for name in open(sequence_names, "r").readlines()}
    content = ""

    for sequence in sequences:
        result_file = [file for file in result_files if re.search(sequence, file)]
        if result_file:
            with open(result_file[0]) as left_out_seq_file:
                for line in left_out_seq_file:
                    f"^[^#Q].*{sequence}"
                    if sequence in line and not line.startswith("Query") and not line.startswith("#"):
                        content += line
                        break

    with open(csv_file, "w") as csv:
        csv.write(content)

def cut_sequence(file_name: str, motif: str, output_filename: str):
    """ Cut sequences from the search conserved residues
        Input: fasta file
        Output: dictionary of cut protein sequences with the name of the proteins as key"""
    sequences = read_fasta(file_name)
    proteins = dict(sequences.items())
    pattern = re.compile(motif)
    for protein in proteins:
        # maxsplit has to be set to one so it consideres the first occurence of the motif
        proteins[protein] = re.split(pattern, proteins[protein], maxsplit=1)[1]
    write_fasta(proteins.keys(), proteins.values(), output_filename)

crossval_preproc("HALOGENATION/DIMETAL-CARBOXYLATE/dimetal-carboxylate.fasta")

# FDHs
crossval_preproc("HALOGENATION/FLAVIN_DEPENDENT_HALOGENASES/general_profile/cut_all_general_profile.afa")
get_left_out_seq_results("HALOGENATION/FLAVIN_DEPENDENT_HALOGENASES/general_profile/cut/crossval_hmmsearch_res/crossval_out_*.out",
                         "HALOGENATION/FLAVIN_DEPENDENT_HALOGENASES/general_profile/cut/samples.txt")

# if __name__ == "__main__":
#     fasta_file_path = sys.argv[1]
#     # test_fasta = "/home/szenei/ENZYMES/FLAVIN_DEPENDENT_HALOGENASES/test/alkyl_halides_FDH.fasta"
#     # cut_sequence(test_fasta, "G.G..G", "/home/szenei/ENZYMES/FLAVIN_DEPENDENT_HALOGENASES/test/cut_alkyl_halides_FDH.fasta")

#     crossval_preproc(fasta_file_path)
#     # print(get_left_out_seq_results("/home/szenei/ENZYMES/HALOPEROXIDASES/VCPOs/crossval_hmmsearch_res/crossval_out_*.out",
#     #                                "/home/szenei/ENZYMES/HALOPEROXIDASES/VCPOs/samples.txt"))
