#!/bin/bash

dbs=`cat ensembl_98_all_dbs | perl -e 'while(<STDIN>){ s/.+?\s+(\S+)$/$1/; print; }' | grep -E "^[^_]+_[^_]+_core"`

base=ftp://ftp.ensembl.org/pub/release-98/mysql/

for db in $dbs
do
    if [ -e $db ]
    then
	echo $db already exists. Nothing will be done.
	continue
    fi
    file=${base}${db}/*
    ##sp=`echo $db | cut -f 1,2 -d '_'`
    echo "wget $file"
    mkdir $db
    cd $db
    wget -q $file
    cd ..
    ./install_db.sh $db &
done
	 
