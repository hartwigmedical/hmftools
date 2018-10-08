#!/usr/bin/env bash

db_generate_script=$1

OS=$(uname)
if [ ${OS} = "Darwin" ];
then
    script_epoch=$(stat -t '%s' ${db_generate_script} | cut -d\" -f4)
else
    # HEKE: Assume Linux with GNU date syntax
    script_epoch=$(date -r ${db_generate_script} '+%s')
fi

db_epoch=$(mysql --defaults-file=~/mysql.login << HERE | sed -n '2p'
SELECT UNIX_TIMESTAMP(MAX(create_time)) db_creation FROM INFORMATION_SCHEMA.TABLES WHERE table_schema = "hmfpatients_test";
HERE
)

echo "[INFO]: Script Epoch: ${script_epoch}"
echo "[INFO]: DB Epoch: ${db_epoch}"

if [ "$db_epoch" = "NULL" ] || [ ${script_epoch} -gt ${db_epoch} ];
then
    echo "[INFO] Rebuilding test database based on ${db_generate_script}"
    mysql --defaults-file=~/mysql.login < ${db_generate_script}
fi
