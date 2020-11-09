#!/bin/bash

dbs=`grep -v "#" db_list`

for db in $dbs
do
    cd $db
    echo "creating $db"
    echo "drop database if exists $db" | mysql
    echo "create database $db" | mysql
    sqf=${db}.sql
    echo $sqf
    mysql $db < $sqf
    mysqlimport --fields-terminated-by='\t' --fields_escaped_by=\\ $db -L *.txt
    cd ..
done
