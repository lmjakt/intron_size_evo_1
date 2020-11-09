#!/bin/bash

db=$1

if [ -d $db ]
then
    cd $db
    echo "uncompressing files in $db"
    gunzip *.gz
    echo "creating and installing database $db"
    echo "drop database if exists $db" | mysql
    echo "create database $db" | mysql
    sqf=${db}.sql
    echo "sql file is $sqf"
    mysql $db < $sqf
    mysqlimport --fields-terminated-by='\t' --fields_escaped_by=\\ $db -L *.txt
else
    echo $db does not exist or is not a directory
fi
