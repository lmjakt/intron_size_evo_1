#!/bin/bash
## get primary assemblies as these do not contain so much variance
## data that we run out of memory

dbs=`head -n 1 ../vertebrate_family_members_1_15.txt`

base=ftp://ftp.ensembl.org/pub/release-97/fasta/

for db in $dbs
do
    sp=`echo $db | cut -f 1,2 -d '_'`
    file=${base}${sp}/dna/*dna.primary_assembly.fa.gz
    echo "calling wget with $file"
    wget $file
done
