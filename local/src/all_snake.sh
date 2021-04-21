#!/bin/bash
shopt -s expand_aliases
alias lsnakemake='/home/egrassi/.local/bin/snakemake --log-handler-script /home/egrassi/sysadm/snakemake_slack.py'

N=3
dir=$1
target=$2
shift 2
for crc in "$@"; do
    ((i=i%N)); ((i++==0)) && wait;
    echo $crc/$dir $target;
    cd $crc/$dir;
    lsnakemake $target &
    cd ../..
done

