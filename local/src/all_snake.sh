#!/bin/bash
shopt -s expand_aliases
alias lsnakemake='/home/egrassi/.local/bin/snakemake --log-handler-script /home/egrassi/sysadm/snakemake_slack.py'

dir=$1
target=$2
shift 2
for crc in "$@"; do
    echo $crc/$dir $target
    cd $crc/$dir
    #lsnakemake $target
done

