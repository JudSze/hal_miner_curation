import re
from Bio import SeqIO

apnu_homologs_path = "hal_miner_curation/copper-dependent_hits.fasta"

apnu_homologs = SeqIO.parse(open(apnu_homologs_path), "fasta")
apnu_2x_homologs = "/home/szenei/hal_miner_curation/uniref50_2x_homologs.fasta"
apnu_2x_homologs = open("uniref50_2x_homologs.fasta", "wb")

with open("uniref50_2x_homologs.fasta", "w") as output_handle:
    for homolog in apnu_homologs:
        name, seq = homolog.id, str(homolog.seq)
        r = re.findall(r"H..HC", seq)
        if len(r) == 2:
            print(name, r)
            SeqIO.write(homolog, output_handle, "fasta")
        else:
            print(name)

# apnu_homologs