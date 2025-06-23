#!/bin/bash
set -euo pipefail

# Disable glob expansion. "*" in filenames are handled within hmftools java code
set -f

# Parse arguments
jvm_mem_opts=""
jvm_gen_opts=""
declare -a app_args
for arg in "$@"; do
  case ${arg} in
    '-Xm'*)
      jvm_mem_opts="${jvm_mem_opts} ${arg}"
      ;;
    '-D'* | '-XX'*)
      jvm_gen_opts="${jvm_gen_opts} ${arg}"
      ;;
    *)
      app_args+=("${arg}")
      ;;
  esac
done

# Form and execute command
if [[ ${app_args[0]:=} == com.hartwig.* ]]; then
  java ${jvm_mem_opts} ${jvm_gen_opts} -cp __JAR_PATH__ ${app_args[@]}
else
  java ${jvm_mem_opts} ${jvm_gen_opts} -jar __JAR_PATH__ ${app_args[@]}
fi

