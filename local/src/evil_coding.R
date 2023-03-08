data <- read.table('/scratch/trcanmed/AF_spectra/datasetV2/CRC1502_clones_all/platypus_nobin_00/all.annovar.gz', header=FALSE, sep="\t", stringsAsFactors = F)
# WHY DO WE NEED UNIQUE?
data <- unique(data)
evil <- unique(data$V1[grepl('CRC1502-09C-2', data$V1)])
data$id <- paste0(data$V2,"_", data$V3)
evilmut <- setdiff(data[data$V1 %in% evil, 'id'], data[!data$V1 %in% evil, 'id'])

evilcodingmut <- data[data$id %in% evilmut & data$V7=="exonic" & (data$V10=="nonsynonymous SNV" | data$V10=="stopgain"),]


dft <- as.data.frame(table(evilcodingmut$id))

unique(data[data$id %in% dft[dft$Freq == 3,'Var1'],'V8'])

w <- c('ADAP2' , 'CCDC115' , 'CFAP54' , 'FOLH1B','MUC4' , 'SOX11')

data[data$V8 %in% w & (data$V10=="nonsynonymous SNV" | data$V10=="stopgain"),] # old data

data <- read.table('/scratch/trcanmed/AF_spectra/dataset/CRC1502_clones_all/platypus_nobin/all.annovar.gz', header=FALSE, sep="\t", stringsAsFactors = F)
# WHY DO WE NEED UNIQUE?
data <- unique(data)
evil <- unique(data$V1[grepl('CRC1502-09C-2', data$V1)])
data$id <- paste0(data$V2,"_", data$V3)
evilmut <- setdiff(data[data$V1 %in% evil, 'id'], data[!data$V1 %in% evil, 'id'])

evilcodingmut <- data[data$id %in% evilmut & data$V7=="exonic" & (data$V10=="nonsynonymous SNV" | data$V10=="stopgain"),]


dft <- as.data.frame(table(evilcodingmut$id))

data[data$id %in% dft[dft$Freq == 3,'Var1'],]


#egrassi@godot:/scratch/trcanmed/AF_spectra/datasetV2/CRC1502_clones_all/platypus_nobin_00$ zgrep SOX11 all.annovar.gz  | sort | uniq
#CRC1502-09C-2-1 2       5693358 5693358 G       A       exonic  SOX11   .       nonsynonymous SNV       SOX11:NM_003108:exon1:c.G637A:p.A213T   .       .
#CRC1502-09C-2-2 2       5693358 5693358 G       A       exonic  SOX11   .       nonsynonymous SNV       SOX11:NM_003108:exon1:c.G637A:p.A213T   .       .
#CRC1502-09C-2-3 2       5693358 5693358 G       A       exonic  SOX11   .       nonsynonymous SNV       SOX11:NM_003108:exon1:c.G637A:p.A213T   .       .

w <- c('ADAP2' , 'CCDC115' , 'CFAP54' , 'FOLH1B','MUC4' , 'SOX11')

unique(data[data$id %in% dft[dft$Freq == 3,'Var1'],'V8'])
data[data$id %in% dft[dft$Freq == 3,'Var1'],]

## 