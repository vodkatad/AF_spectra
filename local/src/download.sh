rsync -av -L mviviani@occam.c3s.unito.it:/archive/home/mviviani/snakegatk/dataset/$1/platypus/platypus_filtered.vcf.gz .
#rsync -av -L mviviani@occam.c3s.unito.it:/archive/home/mviviani/snakegatk/dataset/$1/depth/callable_covered.bed.gz .
rsync -av -L mviviani@occam.c3s.unito.it:/archive/home/mviviani/snakegatk/dataset/$1/mutect_paired/*pass.vcf.gz* .
#rsync -av -L mviviani@occam.c3s.unito.it:/archive/home/mviviani/snakegatk/dataset/$1/sequenza/*/*segments.txt .
