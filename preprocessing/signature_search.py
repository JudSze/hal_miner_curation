import sys
import re
from typing import List, Dict

from antismash.common.fasta import read_fasta
from operator import itemgetter


import logging
from typing import Optional

from antismash.common import subprocessing, utils
from antismash.config import build_config

build_config([])

ILLEGAL_CHARS = "!@#$%^&*(){}:\"<>?/.,';][`~1234567890*-+-=_\\|"
PROFILE = "/home/szenei/antismash/MAGic-MOLFUN/HALOGENASES/trp_6_7/trp_6_7_v2.hmm"
PROF_NAME = "trp_6_7"
# WxWxIx and Fx.Px.Sx.G motifs
GENERAL_FDH = {"second_motif": [206, 207, 208, 209, 210, 211],
               "last_motif": [280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304]}

PYRROLE_SIGNATURE = [110, 111, 318, 322, 362]
TYROSINE_Hpg_SIGNATURE = {"common": [58, 74, 89, 92, 99, 107, 149, 150, 152, 209, 215, 217, 219, 245, 267, 268, 282, 284, 289, 290, 293, 295, 305, 331, 357],
                          "glycopeptide": [34, 35, 36, 37],
                          "cyanobacterial": [24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
                          "Hpg": [66, 158, 196, 200, 246, 259]}

OTHER_PHENOLIC_SIGNATURE = [23, 27, 39, 40, 59, 74, 109, 113, 120, 124, 133, 165, 166, 168,
                            233, 284, 291, 303, 305, 306, 309, 311]

TRP_5_SIGNATURE = [33, 35, 41, 75, 77, 78, 80, 92, 101, 128, 165, 186, 187, 195, 302, 310, 342, 350, 400, 446, 450, 452, 482 ]
TRP_6_SIGNATURE = [19, 37, 45, 73, 75, 90, 129, 130, 142, 157, 181, 192, 194, 219, 221, 225, 227, 237, 287, 306, 337, 339, 350, 353, 356, 411, 462, 505]

# from antismash code
def verify_good_sequence(sequence: str) -> bool:
    """ Ensures a sequence is valid """
    for char in ILLEGAL_CHARS:
        if char in sequence:
            return False
    return True

# modified from antismash code
def search_signatures(hmm: str, hmm_name: str, sequence, protein_name, positions: list, max_evalue: float = 0.1) \
        -> tuple[Optional[str], Optional[str]]:
    """ Extract signatures from enzyme profiles """
    assert verify_good_sequence(sequence)

    args = ["-E", str(max_evalue)]
    results = subprocessing.hmmpfam.run_hmmpfam2(hmm, f">query\n{sequence}", extra_args=args)

    if not (results and results[0].hsps):
        logging.debug("no hits for query %s, is this a legitimate specific enzyme?", protein_name)
        return None, None

    found = False
    hit = None  # will be set in the loop or we abort anyway, just to make pylint happy
    for hit in results[0].hsps:
        if hit.hit_id == hmm_name:
            found = True
            break

    if not found:
        logging.debug(
            "no hits for the active site in %s, is this a legitimate A-domain?", protein_name)
        return None, None

    profile = hit.aln[1].seq
    query = hit.aln[0].seq
    offset = hit.hit_start

    sites = utils.extract_by_reference_positions(query, profile, [p - offset for p in positions])

    return sites, query

class FindConsensus:
    def __init__(self, sequences: str) -> None:
        self.sequences = read_fasta(sequences)
        self.translations = {protein: translation for protein, translation in self.sequences.items()}

    @staticmethod
    def find_conserved_residues(sequences):
        """ Find the positions of conserved residues in a given list of sequences"""
        conserved_residues = []
        residues_same_position = list(map(set, zip(*sequences)))
        for position in range(len(residues_same_position)):
            if len(list(residues_same_position[position])) == 1:
                conserved_residues.append(position)

        return conserved_residues

    @staticmethod
    def retrieve_signature(alignment: dict, conserved_residues: list):
        "Retrieve signature from the alignment"
        signature = dict()
        for id, seq in alignment.items():
            signature[id] = "".join(itemgetter(*conserved_residues)(seq))
        return signature

    @staticmethod
    def consensus_for_subgroups(sequences: dict, targets: list):
        common_positions = []

        translated_targets = {id: sequence for id, sequence in sequences.items() if id in targets}
        subset = {id: sequence for id, sequence in sequences.items() if id not in targets}

        targets_consensus = FindConsensus.find_conserved_residues(translated_targets.values())
        targets_signature = FindConsensus.retrieve_signature(translated_targets, targets_consensus)
        mapped_target_res_pos = dict(zip(targets_consensus, targets_signature[targets[0]]))

        complement_group_signature = FindConsensus.retrieve_signature(subset,
                                                                      targets_consensus)

        for protein in complement_group_signature.values():
            mapped_res_pos = dict(zip(targets_consensus, protein))
            for pos, residue in mapped_res_pos.items():
                if residue == mapped_target_res_pos[pos]:
                    common_positions.append(pos)

        return [pos for pos in targets_consensus if pos not in common_positions]

if __name__ == "__main__":
    if sys.argv[1] == "consensus":
        enzyme_group = FindConsensus(sys.argv[2])
        consensus_pos = FindConsensus.find_conserved_residues(enzyme_group.translations.values())
        signatures = FindConsensus.retrieve_signature(enzyme_group.translations, consensus_pos)
        print(consensus_pos)
        print(signatures)

    elif sys.argv[1] == "subgroup_consensus":
        enzyme_group = FindConsensus(sys.argv[2])
        target_consensus_pos = FindConsensus.consensus_for_subgroups(enzyme_group.sequences,
                                                                     sys.argv[3])
        signatures = FindConsensus.retrieve_signature(enzyme_group.translations,
                                                      target_consensus_pos)

        print(target_consensus_pos)
        print(signatures[target] for target in sys.argv[3])

    elif sys.argv[1] == "phmm_pos":
        # path to the pHMM profile in your local directory
        profile = sys.argv[2]
        # name of the pHMM (can be found in the pHMM file next to the 'NAME' label)
        # the translation cannot contain any separator (e.g. space, tab)
        profile_name = sys.argv[3]
        if sys.argv[5] == "manual":
            signature = [int(pos) for pos in sys.argv[6:]]
        elif sys.argv[5] == "range":
            signature = list(range(int(sys.argv[6]), int(sys.argv[7])))

        sequences = FindConsensus(sys.argv[4]) # sequences you would like to query

        for protein, translation in sequences.sequences.items():
            siganture_residues = search_signatures(profile, profile_name,
                                                   translation,
                                                   protein, signature)
            print(f'{protein}', siganture_residues[0])
            # if siganture_residues[0] == "RDAGDGSGGFDPFSGD":
            #     print(f"{protein}")

    elif sys.argv[1] == "motifs":
        # path to the pHMM profile in your local directory
        profile = sys.argv[2]
        # name of the pHMM (can be found in the pHMM file next to the 'NAME' label)
        # the translation cannot contain any separator (e.g. space, tab)
        profile_name = sys.argv[3]
        sequences = FindConsensus(sys.argv[4]) # sequences you would like to query

        signature = sys.argv[5]
        motif = sys.argv[6]

        for protein, translation in sequences.sequences:
            siganture_residues = search_signatures(profile, profile_name, translation,
                                                   protein, signature)
            print(f"{protein}", rf"{siganture_residues[0]}")
            for signature in siganture_residues:
                if not siganture_residues[0]:
                    print("None")
                elif not re.search(f"{motif}", siganture_residues[0]):
                    print("no motif found in" f'{protein}')
                else:
                    print("ok " f"{protein}")
