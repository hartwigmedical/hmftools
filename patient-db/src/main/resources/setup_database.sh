#!/usr/bin/env bash

OS=$(uname)
if [ ${OS} = "Darwin" ];
then
    DARWIN_VERSION=$(uname -r | cut -d. -f1)
    if [ ${DARWIN_VERSION} -gt 15 ];
    then
        SCRIPT_EPOCH=$(date -r $1 '+%s')
    else
        SCRIPT_EPOCH=$(stat -t '%s' $1 | cut -d\" -f4)
    fi
else
    #Assume Linux with GNU date syntax
    SCRIPT_EPOCH=$(date -r $1 '+%s')
fi

DB_EPOCH=$(mysql --defaults-file=~/mysql.login << HERE | sed -n '2p'
SELECT UNIX_TIMESTAMP(MAX(create_time)) db_creation FROM INFORMATION_SCHEMA.TABLES WHERE table_schema = "hmfpatients";
HERE
)

if [ ${SCRIPT_EPOCH} -gt ${DB_EPOCH} ];
then
    echo REBUILDING DATABASE
    mysql --defaults-file=~/mysql.login < $1
fi
