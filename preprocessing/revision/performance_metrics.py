# parse_loocv.py
import os, glob
import pandas as pd
from math import sqrt

def parse_hmmsearch_out(filepath):
    """Return set of sequence names that scored above --notextw default output."""
    hits = set()
    with open(filepath) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            # columns: target_name, accession, query_name, accession, ...
            # full-sequence E-value is col 4, bitscore is col 5
            try:
                seq_name = parts[8]
                bitscore = float(parts[1])
                hits.add((seq_name, bitscore))
            except (IndexError, ValueError):
                continue
    return hits

def loocv_metrics(loocv_dir, group_fasta, threshold):
    """
    loocv_dir   : directory containing one .out file per left-out sequence
    group_fasta : all positive sequences for this group (to know set size)
    threshold   : bitscore cutoff
    """
    # parse positive set names from fasta
    positives = set()
    with open(group_fasta) as fh:
        for line in fh:
            if line.startswith('>'):
                positives.add(line[1:].strip().split()[0])

    TP = FN = 0
    for out_file in glob.glob(os.path.join(loocv_dir, '*.out')):
        left_out = os.path.basename(out_file).replace('.out', '')
        hits = {seq.strip("|fluorinase") for seq, score in parse_hmmsearch_out(out_file)
                if score >= threshold}
        if left_out in hits:
            TP += 1
        else:
            print(left_out)
            FN += 1

    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    return {'n_positive': len(positives), 'TP': TP, 'FN': FN, 'Recall': recall}

loocv_metrics("HALOGENATION/NON_HEME_HALOGENASES/variant_A_unactivated_sp3/amino_acids/crossval_hmmsearch_res","HALOGENATION/NON_HEME_HALOGENASES/variant_A_unactivated_sp3/amino_acids/variant_A_small_amino_acids_copy.fasta", 0)