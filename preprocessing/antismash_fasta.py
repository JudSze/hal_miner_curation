# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions supporting the FASTA format
"""

from collections import OrderedDict
import logging
from typing import Dict, Iterator, List

def build_fasta(mapping: dict[str, str]) -> Iterator[str]:
    """ Generates a FASTA formatted string

        Arguments:
            mapping: a mapping of identifier to sequence

        Returns:
            an iterator with each identifier/sequence pair
    """
    for name, seq in mapping.items():
        yield f">{name}\n{seq}\n"


def write_fasta(names: List[str], seqs: List[str], filename: str) -> None:
    """ Writes name/sequence pairs to file in FASTA format

        Argumnets:
            names: a list of sequence identifiers
            seqs: a list of sequences as strings
            filename: the filename to write the FASTA formatted data to

        Returns:
            None
    """
    with open(filename, "w", encoding="utf-8") as out_file:
        for chunk in build_fasta(dict(zip(names, seqs))):
            out_file.write(chunk)


def read_fasta(filename: str) -> Dict[str, str]:
    """ Reads a fasta file into a dictionary

        Arguments:
            filename: the path to the FASTA file to read

        Returns:
            a dictionary mapping sequence ID to sequence

    """
    ids = []
    sequence_info = []
    with open(filename, "r", encoding="utf-8") as fasta:
        current_seq: List[str] = []
        for line in fasta:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':
                ids.append(line[1:].replace(" ", "_"))
                if current_seq:
                    sequence_info.append("".join(current_seq))
                    current_seq.clear()
            else:
                if not ids:
                    raise ValueError("Sequence before identifier in fasta file")
                if not line.replace("-", "z").isalpha():
                    raise ValueError("Sequence contains non-alphabetic characters")
                current_seq.append(line)
    if current_seq:
        sequence_info.append("".join(current_seq))
    if len(ids) != len(sequence_info):
        raise ValueError("Fasta files contains different counts of sequences and ids")
    if not ids:
        logging.debug("Fasta file %s contains no sequences", filename)
        raise ValueError("Fasta file contains no sequences")
    return OrderedDict(zip(ids, sequence_info))
