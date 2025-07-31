![image](https://github.com/JudSze/ENZYMES/assets/77547600/662840a0-6233-4ea8-a853-405847dde014)# Flavin-dependent Halogenases (FDH)

## The pHMM built for Trp-5 halogenases were built from 7 whole protein sequences.
Signature is used to differentiate Trp-5 halogenases from Trp-6/-7 halogenases.

Signature (Trp-5): 
```diff
+ VSILIREPGLPRGVPRAVLPGEA
+ LOW_QUALITY_CUTOFF: 350
+ HIGH_QUALITY_CUTOFF: 830
```
Positions in pHMM: 
```diff
+ 33, 35, 41, 75, 77, 78, 80, 92, 101, 128, 165, 186, 187, 195, 302, 310, 342, 350, 400, 446, 450, 452, 482 
```

## The pHMM built for Trp-6-7 halogenases were built from 11 whole protein sequences
A signature is used to differentiate Trp-6 halogenases from Trp-7 halogenases.
```diff
+ CUTOFF: 770
```

Signature (Trp-6):
```diff
+ TEGCAGFDAYHDRFGNADYGLSIIAKIL
```

Positions in pHMM: 
```
+ 19, 37, 45, 73, 75, 90, 129, 130, 142, 157, 181, 192, 194, 219, 221, 225, 227, 237, 287, 306, 337, 339, 350, 353, 356, 411, 462, 505 
```

## The pHMM for pyrrole-accepting halogenases were built from 15 sequences.
It includes:
* 5 sequences of mono/dihalogenating enzymes (Mpy16, Pyr29, Clz5, PltA, HrmQ) accepting substrate attached to ACP
* 5 sequences of tetrabrominating enzymes (Bmp2 and 4 homologs)
* 4 sequences of ambiguous, possibly monohalogenating enzymes (HalB from pentachloropseudilin biosynthesis, its 3 homologs and PrnC)

The sequences included in this pHMM can be found in the pyrrole/cut_pyrrole_accepting.fasta
```diff
- NOT WHOLE SEQUENCES (the part of the sequences before GxGxxG motif was cut, as well as the end after WI/WK/YR/FI/YT residues)
```

The whole sequences and further homologs (BASED ON BLAST SEARCH) of the enzymes can be found in the test folder.

In the preprocessing folder in the signature_search.py, you find a tutorial to check the signatures in the sequences.

Signature (pyrrole-accepting):
```diff
+ CUTOFF: 400
```
```diff
+ conventional mono/dihalogenases: (D)RSVW
+ tetrahalogenases: (R)RYFA
+ ambiguous mono/dihalogenases: (Y)RRNN
```
Positions in pHMM: 
```diff
+ 110, 111, 318, 322, 362 
+ position 111 (Arg) is specific to non-Trp accepting halogenases
+ position 318, 322, 362 are specific to halogenation activity in this group
+ position 110 was added as a "safety net" in case there would be deviation in the signature, but the bitscore would be high (well above 400)
```
## The pHMM for Tyr-like or Hpg-accepting halogenases were built from 12 sequences.
The sequences included in this pHMM can be found in the phenolic/substrate_based/tyrosine_ish/cut_tyrosine_ish_FDH_v2.fasta
```diff
- NOT WHOLE SEQUENCES (cut from the beginning until the end pf the GxGxxG motif and the end after KFI&SZF plus 5 residues)
```
Signature (Tyr-like/Hpg-accepting): 
```diff
+ LOW_QUALITY_CUTOFF: 300
+ HIGH_QUALITY_CUTOFF: 390
```
```diff
+ Hpg-accepting: SHCGMQ
+ Tyrosine and Hpg common: GFQRLGDAGLSGVPSYGADPSGLYW
```
Positions in pHMM: 
```diff
+ 66, 158, 196, 200, 246, 259 (for the Hpg-accepting ones) 
+ 58, 74, 89, 92, 99, 107, 149, 150, 152, 209, 215, 217, 219, 245, 267, 268, 282, 284, 289, 290, 293, 295, 305, 331, 357 (for Tyr and Hpg)
```
## The pHMM for other phenolic-accepting (orsellinic-acid-like) halogenases were built from 9 sequences.
The sequences included in this pHMM can be found in the phenolic/substrate_based/cycline_orsellinic_acid/cycline_orsellinic_FDH_v2.hmm
```diff
+ CUTOFF: 500
```
```diff
- WHOLE SEQUENCES
```
Signature (Tyr-like/Hpg-accepting): 
```diff
+ LGPRGGRDAGVDAGGYGFDPSG
```
Positions in pHMM: 
```diff
+ 23, 27, 39, 40, 59, 74, 109, 113, 120, 124, 133, 165, 166, 168, 233, 284, 291, 303, 305, 306, 309, 311
```
### [Conserved motifs](https://doi.org/10.1038/s41557-019-0349-z):
GxGxxG (present in flavin-dependent monooxygenases)

WxWxIP (absent in unusual FDHs)

Fx.Px.Sx.G (potentially the best for identification)

### Types of FDHs based on substrate-specificity
[Indole accepting halogenases](https://doi.org/10.1016/bs.enz.2020.05.009):
* PrnA [P95480](https://www.uniprot.org/uniprotkb/P95480/entry), [G0ZGJ1](https://www.uniprot.org/uniprotkb/G0ZGJ1/entry)
  * https://doi.org/10.1126%2Fscience.1116510
  * https://doi.org/10.1002/anie.201007896
  * https://doi.org/10.1016%2Fj.jmb.2009.06.008
  * https://doi.org/10.1002/anie.200802466
  * https://doi.org/10.1002/anie.201300762
* RebH [TRP7H_LENAE](https://www.uniprot.org/uniprotkb/Q8KHZ8/entry)
  * https://doi.org/10.1073/pnas.0500755102
  * https://doi.org/10.1002/anie.201300762
  * https://doi.org/10.1002/cbic.201300780
  * https://doi.org/10.1016/j.jmb.2010.01.020
  * http://dx.doi.org/10.1002/ange.201411901
  * https://doi.org/10.1039/C5SC04680G
* PyrH [A4D0H5](https://www.uniprot.org/uniprotkb/A4D0H5/entry), [W1J423](https://www.uniprot.org/uniprotkb/W1J423/entry), [K7QVV7](https://www.uniprot.org/uniprotkb/K7QVV7/entry)
  * https://doi.org/10.1016/j.jmb.2009.06.008
  * https://doi.org/10.1002/cbic.201600051
* Thal [A1E280](https://www.uniprot.org/uniprotkb/A1E280/entry)
  * https://doi.org/10.1002/pro.3739
  * https://doi.org/10.1074/jbc.ra118.005393
* SttH [E9P162](https://www.uniprot.org/uniprotkb/E9P162/entry)
* Th-Hal [A0A1L1QK36](https://www.uniprot.org/uniprotkb/A0A1L1QK36/entry)
  * https://doi.org/10.1039/C6OB01861K
* Tar14 [W5VG40_SACPI](https://www.uniprot.org/uniprotkb/W5VG40/entry)
  * https://doi.org/10.1002%2Fanie.201901571
  * https://doi.org/10.1021/acs.chemrev.6b00571
* BrvH
  * https://doi.org/10.1002/cbic.202000444
  * https://doi.org/10.1039/9781788017008-00001
* BorH [M9QSI0](https://www.uniprot.org/uniprotkb/M9QSI0/entry)
  * https://doi.org/10.1002/cbic.201900667
* MalA [L0E155](https://www.uniprot.org/uniprotkb/L0E155/entry)
  * http://dx.doi.org/10.1021/jacs.7b06773![image](https://github.com/JudSze/ENZYMES/assets/77547600/789852aa-7f3b-44fc-8e9c-ce8387e655a4)
* MibH [E2IHC5](https://www.uniprot.org/uniprotkb/E2IHC5/entry)
  * https://doi.org/10.1073/pnas.1008285107
  * https://doi.org/10.1021/acschembio.6b01031
* KtzR [A8CF74](https://www.uniprot.org/uniprotkb/A8CF74/entry)
  * https://doi.org/10.1021/ja806467a
  * https://doi.org/10.1007/s10529-011-0595-7
* KtzQ [A8CF75](https://www.uniprot.org/uniprotkb/A8CF75/entry)
  * https://doi.org/10.1021/ja806467a
  * https://doi.org/10.1007/s10529-011-0595-7
* aclH 
* ACRI_5336 (potentially)
* CmdE [Q0VZ69](https://www.uniprot.org/uniprotkb/Q0VZ69/entry)
  * 10.1016/j.chembiol.2006.06.002
* KrmI
  * https://doi.org/10.1021/acschembio.6b01115
* XszenFHal [W1J423](https://www.uniprot.org/uniprotkb/W1J423/entry)
* SpmH
  * https://doi.org/10.1039/C8OB02775G
* AdeV
  * https://doi.org/10.1002/anie.201914994
* SpH1
  * https://doi.org/10.3390/catal11040485
* SpH2
  * https://doi.org/10.3390/catal11040485
* AetA [A0A861B8S3](https://www.uniprot.org/uniprotkb/A0A861B8S3/entry)
  * https://doi.org/10.1021/jacs.1c12778
* AetF [A0A861B9Z9](https://www.uniprot.org/uniprotkb/A0A861B9Z9/entry)
  * https://doi.org/10.1021/jacs.1c12778
  * https://doi.org/10.1002/ange.202214610
* satH []()
  * https://doi.org/10.1002/cbic.201900723

Pyrrole accepting halogenases (https://doi.org/10.1016/bs.enz.2020.05.009, 10.1039/D0CS01551B):
* HalB [Q71ME2](https://www.uniprot.org/uniprotkb/Q71ME2/entry) (free pyrrole)
* PltA [B3G2A5](https://www.uniprot.org/uniprotkb/B3G2A5/entry) (tethered pyrrole)
* Bmp2 [U6BGC3](https://www.uniprot.org/uniprotkb/U6BGC3/entry) (tethered pyrrole)
* PrnC [P95482](https://www.uniprot.org/uniprotkb/P95482/entry) (free pyrrole)
* Clz5 [U6A3L4](https://www.uniprot.org/uniprotkb/U6A3L4/entry) 
* HrmQ [C1IHU5](https://www.uniprot.org/uniprotkb/C1IHU5/entry) (tethered pyrrole)
* Pyr29 [A3R4S0](https://www.uniprot.org/uniprotkb/A3R4S0/entry)
* Mpy16 [J7H1A1](https://www.uniprot.org/uniprotkb/J7H1A1/entry)

[Phenolic substrate accepting halogenases](https://doi.org/10.1016/bs.enz.2020.05.009):
* Bmp5
* Nat1 [F8QPG4](https://www.uniprot.org/uniprotkb/F8QPG4/entry)
* NapH2
* Ram20
* End30
* TiaM [E9LIP1](https://www.uniprot.org/uniprotkb/E9LIP1/entry)
* VhaA [G4V4R7](https://www.uniprot.org/uniprotkb/G4V4R7/entry)
* Gsfl [D7PI14](https://www.uniprot.org/uniprotkb/D7PI14/entry)
* GedL 
* Tcp21 [Q70AY7](https://www.uniprot.org/uniprotkb/Q70AY7/entry)
* CcxI [A0A4D6Q3Y9](https://www.uniprot.org/uniprotkb/A0A4D6Q3Y9/entry)
* AcoTAhal 
* McnD
* PtaM
* ChlA [Q54FI4](https://www.uniprot.org/uniprotkb/Q54FI4/entry)
* RadH
* Clo-hal
* ChIB4
* StaL
* CalO3 [Q8KND5](https://www.uniprot.org/uniprotkb/Q8KND5/entry)
* BhaA [O87676](https://www.uniprot.org/uniprotkb/O87676/entry)
* AviH 
* AscD
* Arm1
* CndH [B9ZUJ5](https://www.uniprot.org/uniprotkb/B9ZUJ5/entry)
* PltM

[Aliphatic substrate accepting halogenases](https://doi.org/10.1016/bs.enz.2020.05.009):
* AoiQ
* ClmS
* Hydrox

[Diverse substrate accepting halogenases](https://doi.org/10.1016/bs.enz.2020.05.009):
* VirX1 (from cyanophage)
