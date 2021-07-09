#!/bin/bash


#conda activate py27

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M
do
	echo $CHR
	python ~/local/G4Hunter/G4Hunter.py -i ~/QNAP2/iGenome/Homo_sapiens/NCBI/GRCh38/Sequence/Chromosomes/chr${CHR}.fa -o ~/work/G4Hunter/ -w 25 -s 1.2

done

#conda deactivate
