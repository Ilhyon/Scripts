#!/bin/bash

for chr in `seq 1 22` X Y 
do 
	python ReplaceInformation.py --chromosome $chr
	### python g4_classifier_dev/fasta_enumerator.py --fasta-file ../data/chr$chr/HS_gene_unspliced_chr$chr_*.fas --nt-limit 2000000 --large --output ../data/chr$chr -e
	### python g4_classifier_dev/fasta_enumerator.py --fasta-file ../data/chr$chr/HS_transcript_unspliced_chr$chr_Sequence_*.txt --nt-limit 2000000 --large --output ../data/chr$chr -e
	python G4Annotation.py --chromosome $chr
done

## merged files ../results/perChromosome/HS_chr$chr_G4InTranscript.txt
## removed headers ---> ../results/all/HS_All_G4InTranscript.txt
python G4Calcul2.py --chromosome all -- choice Coding

## Mammouth
## run g4_classifier_dev/fasta_enum_qsub.sh 
### python g4_classifier_dev/fasta_enumerator.py --fasta-file ../data/chr$chr/HS_gene_unspliced_chr$chr_*.fas --nt-limit 2000000 --large --output ../data/chr$chr -e
### python g4_classifier_dev/fasta_enumerator.py --fasta-file ../data/chr$chr/HS_transcript_unspliced_chr$chr_Sequence_*.txt --nt-limit 2000000 --large --output ../data/chr$chr -e


## and G4 screener



