#!/bin/bash
## get primary assemblies as these do not contain so much variance
## data that we run out of memory

dbs=`head -n 1 ../family_members/vertebrate_family_members_1_130.txt`

base=ftp://ftp.ensembl.org/pub/release-98/fasta/

cd seq
echo changed dir to `pwd`

for db in $dbs
do
    sp=`echo $db | cut -f 1,2 -d '_'`
    file=${base}${sp}/dna/*dna.toplevel.fa.gz
    echo "calling wget with $file"
    wget -q $file
done
