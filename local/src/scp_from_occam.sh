#!/bin/bash
mkdir -p qc
scp egrassi@occam.c3s.unito.it:/archive/home/egrassi/bit/prj/snakegatk/dataset/$1/mutect_paired/*pass.vcf.gz .
scp egrassi@occam.c3s.unito.it:/archive/home/egrassi/bit/prj/snakegatk/dataset/$1/sequenza/*/*segments.txt  .
scp egrassi@occam.c3s.unito.it:/archive/home/egrassi/bit/prj/snakegatk/dataset/$1/platypus/platypus_filtered.vcf.gz*  .
scp egrassi@occam.c3s.unito.it:/archive/home/egrassi/bit/prj/snakegatk/dataset/$1/depth/callable_covered.bed.gz  .
scp egrassi@occam.c3s.unito.it:/archive/home/egrassi/bit/prj/snakegatk/dataset/$1/align/*wgsmetrics  qc
scp egrassi@occam.c3s.unito.it:/archive/home/egrassi/bit/prj/snakegatk/dataset/$1/align/*flagstat  qc
scp egrassi@occam.c3s.unito.it:/archive/home/egrassi/bit/prj/snakegatk/dataset/$1/fastqc*/*html  qc
