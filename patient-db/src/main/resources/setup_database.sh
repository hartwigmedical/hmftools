#!/usr/bin/env bash

SCRIPT_EPOCH=$(date -r $1 '+%s') 

DB_EPOCH=$(mysql --defaults-file=~/mysql.login << HERE | sed -n '2p'
SELECT UNIX_TIMESTAMP(MAX(create_time)) db_creation FROM INFORMATION_SCHEMA.TABLES WHERE table_schema = "hmfpatients";
HERE
)

if [ $SCRIPT_EPOCH -gt $DB_EPOCH ];  
then
echo REBUILDING DATABASE
mysql --defaults-file=~/mysql.login < $1
fi

