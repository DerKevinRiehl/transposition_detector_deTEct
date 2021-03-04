# transposition_detector_deTEct
Transposition event detection tool using NGS alignment data and SV calling outputs (VCF files) from PBSV or Sniffles.

Input: Structural variants (VCF file) of PBSV (on PBMM2 alignments) or Sniffles (on NGMLR alignments), transposon annotations (by resonaTE), reference genome (FASTA)

Output: Annotation and classification of detected transposition events (GFF3).

## Installation
Installation as [CondaPackage](https://anaconda.org/derkevinriehl/transposition_detector_detect):
```
 conda install -c derkevinriehl transposition_detector_detect 
```
*Note: Otherwise you can find all source codes in this Github repository.*

## Usage 
(using demo files of this repository, we use reference genome CB4856 and probe alignments SX3351)
```
transposition_deTEct -help

# demo for sniffles_ngmlr alignments
transposition_deTEct -seqHeadTXT demoFiles/sequence_heads.txt -transpGFF3 demoFiles/FinalAnnotations_Transposons.gff3 -assmFasta demoFiles/sequence_CB4856.fasta -svTool sniffles -svFile demoFiles/SX3351_addisababa.sniffles_ngmlr.vcf -outParsedFile demoFiles/sniffles_ngmlr/SX3351_addisababa.SV.vcf.gff3 -outResultFile demoFiles/sniffles_ngmlr/SX3351_addisababa.transpositionEvents.gff3

# demo for pbsv_pbsmm2 alignments
transposition_deTEct -seqHeadTXT demoFiles/sequence_heads.txt -transpGFF3 demoFiles/FinalAnnotations_Transposons.gff3 -assmFasta demoFiles/sequence_CB4856.fasta -svTool pbsv -svFile demoFiles/SX3351_addisababa.pbsv_pbmm2.vcf -outParsedFile demoFiles/pbsv_pbmm2/SX3351_addisababa.SV.vcf.gff3 -outResultFile demoFiles/pbsv_pbmm2/SX3351_addisababa.transpositionEvents.gff3
```

## Explanation of output files
SX3351_addisababa.SV.vcf.gff3.matches.gff3
```
seq1	PBSV	duplication	2909118	2910241	.	+	.	['23769', '23770', '23771'];Sseq1TYPE=DUP;END=2910240;Sseq1LEN=1122
seq1	PBSV	deletion	163800	164962	.	+	.	['1', '11827 '];Sseq1TYPE=DEL;END=164961;Sseq1LEN=-1161
seq1	PBSV	deletion	290360	290514	.	+	.	['11843 '];Sseq1TYPE=DEL;END=290513;Sseq1LEN=-153
seq1	PBSV	insertion	343890	344420	.	+	.	['538 '];merged;Sseq1TYPE=seq4NS;END=344424;Sseq1LEN=533
...
```

SX3351_addisababa.transpositionEvents.gff3
```
seq1	PBSV	deletion	290360	290514	.	+	.	Transposon=11843;Class=2/1/2(hAT,TIR,DNATransposon);Sseq1TYPE=DEL;END=290513;Sseq1LEN=-153
seq1	PBSV	insertion	610241	614786	.	+	.	Transposon=545;Class=2/1/3(CMC,TIR,DNATransposon);merged;merged;Sseq1TYPE=seq4NS;END=611763;Sseq1LEN=1521
seq1	PBSV	deletion	879772	884345	.	+	.	Transposon=556;Class=1/1/2(Gypsy,LTR,Retrotransposon);Sseq1TYPE=DEL;END=884344;Sseq1LEN=-4572
seq1	PBSV	insertion	1126531	1126860	.	+	.	Transposon=23592;Class=2/1/1(Tc1-Mariner,TIR,DNATransposon);Sseq1TYPE=seq4NS;END=1126859;Sseq1LEN=327
...
```

## Citations
Please cite our paper if you find transposition event detector "deTEct" useful:
(in progress)
