#!/bin/bash

dbs=`grep -v "#" db_list`

for db in $dbs
do
    echo $db
    cd $db
    gunzip *.gz
    cd ..
done
