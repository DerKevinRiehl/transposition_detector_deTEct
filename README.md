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

## Citations
Please cite our paper if you find transposition event detector "deTEct" useful:
(in progress)
