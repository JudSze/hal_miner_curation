# Testing results
Files with **.out** format are the results of the hmmsearch against a pHMM.<br>
In each of these output file the out-commented lines in the very beginning of the file contain the<br>
filename of the pHMM that was searched against and the name of the fasta file in which the sequences were aligned against the pHMM.

Files retrieved from UniProt/SwissProt (these files were too large to be stored in the github repository):
* uniprotkb_ec_1_14_14_13_OR_ec_1_NOT_tax_2025_07_24.fasta
* uniprotkb_diiron_NOT_taxonomy_id_9606_A_2025_07_09.fasta
* uniprotkb_flavin_dependent_monooxygenas_2025_07_09.fasta

# Cross-validation workflow
For the crossalidation process, the following folders need to be in the folder of the specific group:
* crossval_fasta
* crossval_hmmsearch_res
* crossval_profile

The fasta files, which need to be in the **crossval_fasta** folder can be generated from an MSA file using the **crossval_preproc** function defined in the preprocessing/fasta_preproc.py file.

The preprocessing/crossval_hmmbuild_search.bash file runs the whole crossvalidation workflow. The first parameter for this command needs to be the folder of the specific group Example for the Trp-5 flavin-dependent halogenases:
1. Make sure you are in the hal_miner_curation folder (you can check it by running pwd in the command line)
2. bash preprocessing/crossval_hmmbuild_search.sh HALOGENATION/FLAVIN_DEPENDENT_HALOGENASES/MAFFT/trp_5 /home/szenei/hal_miner_curation/HALOGENATION/FLAVIN_DEPENDENT_HALOGENASES/MAFFT/trp_5/trp_5_with_ids.fasta
* the first parameter is the path to the folder
* the second parameter is the path to the test file