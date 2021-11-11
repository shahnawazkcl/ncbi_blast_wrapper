#!/bin/bash
mkdir data
for i in `ls -1 *.fasta`
do
	foldername=`basename ${i}`
	mkdir data/${foldername}
	awk -F '>' '/^>/ {F=sprintf("data/${foldername}%s.fasta", $2); print > F;next;} {print F; close(F)}' < ${i}
done