library(annotatr)
#library(rtracklayer)
bed <- '/scratch/trcanmed/AF_spectra/dataset/all_vitro_all_merged_gained.bed'
#gained = read_regions(con = rtracklayer::import(bed, format="bed"), genome = 'hg38')
gained = read_regions(con = bed, genome = 'hg38', format='bed')
seqlevels(gained, pruning.mode="coarse") <- paste0('chr', seq(1,22))

annots = c('hg38_basicgenes', 'hg38_genes_intergenic')
           #'hg38_genes_intronexonboundaries')
annotations = build_annotations(genome = 'hg38', annotations = annots)

# Intersect the regions we read in with the annotations
gained_annotated = annotate_regions(
  regions = gained,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

ctheme <- theme_bw()+theme(text=element_text(size=10), axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), 
                axis.title.y=element_text(size=20), axis.text.y=element_text(size=15), 
                plot.title = element_text(face = "bold", size = 20, hjust = 0.5), legend.text = element_text(size = 15),
                legend.title = element_text(size = 15))





order <- c('intergenic','intron', 'exon','1 to 5kb', '5UTRs', '3UTRs')#,'intron/exon boundaries')
plot_annotations = plot_annotation(
  annotated_regions = gained_annotated,
  #annotation_order = annots,
  plot_title = '# of gained muts overall',
  x_label = 'knownGene Annotations',
  y_label = 'Count')
print(plot_annotations)#+current_theme)

####
dm_random_regions = randomize_regions(
  regions = gained,
  allow.overlaps = TRUE,
  per.chromosome = TRUE)

# Annotate the random regions using the same annotations as above
# These will be used in later functions
dm_random_annotated = annotate_regions(
  regions = dm_random_regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = TRUE)


rplot_annotations = plot_annotation(
  annotated_regions = dm_random_annotated,
  #annotation_order = order,
  plot_title = '# randomized muts',
  x_label = 'knownGene Annotations',
  y_label = 'Count')
print(rplot_annotations+current_theme)


dm_annsum = summarize_annotations(
annotated_regions = dm_random_annotated,
quiet = TRUE)

real_annsum <- summarize_annotations(annotated_regions=gained_annotated, quiet=T)

real <- as.data.frame(real_annsum)
random <- as.data.frame(dm_annsum)
real$SNV <- 'gained'
random$SNV <- 'randomized'

pd <- rbind(real, random)
pd$annot.type <- as.factor(pd$annot.type)
#> levels(pd$annot.type)
#[1] "hg38_genes_1to5kb"               "hg38_genes_3UTRs"                "hg38_genes_5UTRs"               
#[4] "hg38_genes_exons"                "hg38_genes_intergenic"           "hg38_genes_intronexonboundaries"
#[7] "hg38_genes_introns"              "hg38_genes_promoters"           
levels(pd$annot.type) <- c('1to5kb','3UTRs', '5UTRs','exons','intergenic','introns','promoters')
pd$annot.type <- factor(pd$annot.type, levels=c('introns','intergenic','1to5kb','promoters','3UTRs','5UTRs','exons'))
ggplot(data=pd, aes(x=annot.type,y=n,fill=SNV))+geom_bar(stat="identity", position='dodge')+ctheme+xlab('')+ylab('#')+scale_fill_manual(values=c('darkgreen', 'darkgoldenrod'))+theme()



read_annotations(con = '/scratch/trcanmed/AF_spectra/dataset/all_RepliSeq_median_early_merged.bed', genome = 'hg38', name = 'early', format = 'bed')
read_annotations(con = '/scratch/trcanmed/AF_spectra/dataset/all_RepliSeq_median_intermediate_merged.bed', genome = 'hg38', name = 'intermediate', format = 'bed')
read_annotations(con = '/scratch/trcanmed/AF_spectra/dataset/all_RepliSeq_median_late_merged.bed', genome = 'hg38', name = 'late', format = 'bed')

print(annotatr_cache$list_env())

library(reshape)
large <- cast(pd, value='n', formula = "annot.type ~ SNV")
rownames(large) <- large$annot.type
large$annot.type <- NULL
fisher.test(large)
# merge and plot together + chisq

library(maftools)
maff <- '/scratch/trcanmed/AF_spectra/dataset/pcgr/all_gained.pcgr_acmg.grch38.maf'

maff2 <- '/scratch/trcanmed/AF_spectra/dataset/pcgr/test.maf'
maf2 = read.maf(maf = maff2, useAll=TRUE)

titvd = titv(maf = maf2, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = titvd)

plotmafSummary(maf = maf2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

### unfruitful tried with maf from pcgr, ended up using vcf2maf from his annotated vcf
# had to add 'chr' though :/


maf = my.read.maf(maf = maff, useAll=TRUE, vc_nonSyn=c("SNP"))



my.read.maf = function(maf, clinicalData = NULL, removeDuplicatedVariants = TRUE, useAll = TRUE, gisticAllLesionsFile = NULL, gisticAmpGenesFile = NULL,
                    gisticDelGenesFile = NULL, gisticScoresFile = NULL, cnLevel = 'all', cnTable = NULL, isTCGA = FALSE, vc_nonSyn = NULL, verbose = TRUE){
  
  #1. Read MAF if its a file or convert to data.table if its data.frame
  start_time = proc.time()
  if (is.data.frame(x = maf)) {
    maf  = data.table::as.data.table(maf)
  } else{
    if (verbose) {
      cat('-Reading\n')
    }
    
    maf <-
      data.table::fread(
        file = maf,
        sep = "\t",
        stringsAsFactors = FALSE,
        verbose = FALSE,
        data.table = TRUE,
        showProgress = TRUE,
        header = TRUE,
        fill = TRUE,
        skip = "Hugo_Symbol",
        quote = ""
      )
    
    # if(as.logical(length(grep(pattern = 'gz$', x = maf, fixed = FALSE)))){
    #   #If system is Linux use fread, else use gz connection to read gz file.
    #   if(Sys.info()[['sysname']] == 'Windows'){
    #     maf.gz = gzfile(description = maf, open = 'r')
    #     suppressWarnings(maf <- data.table::as.data.table(read.csv(file = maf.gz, header = TRUE, sep = '\t', stringsAsFactors = FALSE, comment.char = "#")))
    #     close(maf.gz)
    #   } else{
    #     maf = suppressWarnings(data.table::fread(cmd = paste('zcat <', maf), sep = '\t', stringsAsFactors = FALSE, verbose = FALSE, data.table = TRUE, showProgress = TRUE, header = TRUE, fill = TRUE, skip = "Hugo_Symbol", quote = ""))
    #   }
    # } else{
    #   suppressWarnings(maf <- data.table::fread(input = maf, sep = "\t", stringsAsFactors = FALSE, verbose = FALSE, data.table = TRUE, showProgress = TRUE, header = TRUE, fill = TRUE, skip = "Hugo_Symbol", quote = ""))
    # }
  }
  
  #2. validate MAF file
  if(verbose){
    cat("-Validating\n")
  }
  #maf = validateMaf(maf = maf, isTCGA = isTCGA, rdup = removeDuplicatedVariants, chatty = verbose)
  
  #3. validation check for variants classified as Somatic in Mutation_Status field.
  if(!useAll){
    cat('--Using only `Somatic` variants from Mutation_Status. Set useAll = TRUE to include everything.')
    if(length(colnames(maf)[colnames(x = maf) %in% 'Mutation_Status']) > 0){
      maf = maf[Mutation_Status %in% "Somatic"]
      if(nrow(maf) == 0){
        stop('No more Somatic mutations left after filtering for Mutation_Status! Maybe set useAll to TRUE ?')
      }
    }else{
      cat('Mutation_Status not found. Assuming all variants are Somatic and validated\n')
    }
  }
  
  #4. Seperate synonymous variants from non-syn variants
  #Variant Classification with Low/Modifier variant consequences. http://asia.ensembl.org/Help/Glossary?id=535
  if(is.null(vc_nonSyn)){
    vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                     "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                     "In_Frame_Ins", "Missense_Mutation")
  }else{
    vc.nonSilent = vc_nonSyn
  }
  # silent = c("3'UTR", "5'UTR", "3'Flank", "Targeted_Region", "Silent", "Intron",
  #            "RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Start_Codon_Ins", "Start_Codon_SNP", "Stop_Codon_Del")
  #Variant Classification with High/Moderate variant consequences. http://asia.ensembl.org/Help/Glossary?id=535
  
  
  maf.silent = maf[!Variant_Classification %in% vc.nonSilent] #Silent variants
  if(nrow(maf.silent) > 0){
    maf.silent.vc = maf.silent[,.N, .(Tumor_Sample_Barcode, Variant_Classification)]
    maf.silent.vc.cast = data.table::dcast(data = maf.silent.vc, formula = Tumor_Sample_Barcode ~ Variant_Classification, fill = 0, value.var = 'N') #why dcast is not returning it as data.table ?
    summary.silent = data.table::data.table(ID = c('Samples',colnames(maf.silent.vc.cast)[2:ncol(maf.silent.vc.cast)]),
                                            N = c(nrow(maf.silent.vc.cast), colSums(maf.silent.vc.cast[,2:ncol(maf.silent.vc.cast), with = FALSE])))
    
    #maf = maf[Variant_Classification %in% vc.nonSilent] #Choose only non-silent variants from main table
    if(verbose){
      cat(paste0('-Silent variants: ', nrow(maf.silent)), '\n')
      #print(summary.silent)
    }
  }
  
  if(nrow(maf) == 0){
    stop("No non-synonymous mutations found\nCheck `vc_nonSyn`` argumet in `read.maf` for details")
  }
  
  #5. Process CN data if available.
  if(!is.null(gisticAllLesionsFile)){
    gisticIp = readGistic(gisticAllLesionsFile = gisticAllLesionsFile, gisticAmpGenesFile = gisticAmpGenesFile,
                          gisticDelGenesFile = gisticDelGenesFile, isTCGA = isTCGA, gisticScoresFile = gisticScoresFile, cnLevel = cnLevel, verbose = verbose)
    gisticIp = gisticIp@data
    
    suppressWarnings(gisticIp[, id := paste(Hugo_Symbol, Tumor_Sample_Barcode, sep=':')])
    gisticIp = gisticIp[!duplicated(id)]
    gisticIp[,id := NULL]
    
    maf = data.table::rbindlist(list(maf, gisticIp), fill = TRUE, use.names = TRUE)
    maf$Tumor_Sample_barcode = factor(x = maf$Tumor_Sample_barcode,
                                      levels = unique(c(levels(maf$Tumor_Sample_barcode), unique(as.character(gisticIp$Tumor_Sample_barcode)))))
    
    #oncomat = createOncoMatrix(maf, chatty = verbose)
  }else if(!is.null(cnTable)){
    if(verbose){
      cat('-Processing copy number data\n')
    }
    if(is.data.frame(cnTable)){
      cnDat = data.table::copy(cnTable)
      data.table::setDT(x = cnDat)
    }else{
      cnDat = data.table::fread(input = cnTable, sep = '\t', stringsAsFactors = FALSE, header = TRUE, colClasses = 'character')
    }
    colnames(cnDat) = c('Hugo_Symbol', 'Tumor_Sample_Barcode', 'Variant_Classification')
    if(isTCGA){
      cnDat[,Tumor_Sample_Barcode := substr(x = cnDat$Tumor_Sample_Barcode, start = 1, stop = 12)]
    }
    cnDat$Variant_Type = 'CNV'
    suppressWarnings(cnDat[, id := paste(Hugo_Symbol, Tumor_Sample_Barcode, sep=':')])
    cnDat = cnDat[!duplicated(id)]
    cnDat[,id := NULL]
    maf = data.table::rbindlist(l = list(maf, cnDat), fill = TRUE, use.names = TRUE)
    maf$Tumor_Sample_barcode = factor(x = maf$Tumor_Sample_barcode,
                                      levels = unique(c(levels(maf$Tumor_Sample_barcode), unique(as.character(cnDat$Tumor_Sample_barcode)))))
  }
  
  #Set factors
  maf$Variant_Classification = as.factor(as.character(maf$Variant_Classification))
  maf$Variant_Type = as.factor(as.character(maf$Variant_Type))
  
  if(verbose){
    cat('-Summarizing\n')
  }
  #mafSummary = summarizeMaf(maf = maf, anno = clinicalData, chatty = verbose)
  
  
  #7. Create MAF object
  m = MAF(data = maf)#, variants.per.sample = mafSummary$variants.per.sample, variant.type.summary = mafSummary$variant.type.summary,
          #variant.classification.summary = mafSummary$variant.classification.summary, gene.summary = mafSummary$gene.summary,
          #summary = mafSummary$summary, maf.silent = maf.silent, clinical.data = mafSummary$sample.anno)
  #m = mafSetKeys(maf = m)
  
  if(verbose){
    cat("-Finished in",data.table::timetaken(start_time),"\n")
  }
  
  return(m)
}
summarizeMaf = function(maf, anno = NULL, chatty = TRUE){
  
  if('NCBI_Build' %in% colnames(maf)){
    NCBI_Build = unique(maf[!Variant_Type %in% 'CNV', NCBI_Build])
    NCBI_Build = NCBI_Build[!is.na(NCBI_Build)]
    
    if(chatty){
      if(length(NCBI_Build) > 1){
        cat('--Mutiple reference builds found\n')
        NCBI_Build = do.call(paste, c(as.list(NCBI_Build), sep=";"))
        cat(NCBI_Build)
      }
    }
  }else{
    NCBI_Build = NA
  }
  
  if('Center' %in% colnames(maf)){
    Center = unique(maf[!Variant_Type %in% 'CNV', Center])
    #Center = Center[is.na(Center)]
    if(length(Center) > 1){
      Center = do.call(paste, c(as.list(Center), sep=";"))
      if(chatty){
        cat('--Mutiple centers found\n')
        cat(Center)
      }
    }
  }else{
    Center = NA
  }
  
  
  #nGenes
  nGenes = length(unique(maf[,Hugo_Symbol]))
  maf.tsbs = levels(maf[,Tumor_Sample_Barcode])
  nSamples = length(levels(maf$Tumor_Sample_Barcode))
  
  #Top 20 FLAGS - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267152/
  flags = flags(top = 20)
  
  #Variants per TSB
  tsb = maf[,.N, Tumor_Sample_Barcode]
  colnames(tsb)[2] = 'Variants'
  tsb = tsb[order(tsb$Variants, decreasing = TRUE),]
  
  #summarise and casting by 'Variant_Classification'
  vc = maf[,.N, .(Tumor_Sample_Barcode, Variant_Classification )]
  vc.cast = data.table::dcast(data = vc, formula = Tumor_Sample_Barcode ~ Variant_Classification, fill = 0, value.var = 'N', drop = FALSE)
  
  if(any(colnames(vc.cast) %in% c('Amp', 'Del'))){
    vc.cast.cnv = vc.cast[,c('Tumor_Sample_Barcode', colnames(vc.cast)[colnames(vc.cast) %in% c('Amp', 'Del')]), with =FALSE]
    vc.cast.cnv$CNV_total = rowSums(vc.cast.cnv[,2:ncol(vc.cast.cnv)], na.rm = TRUE)
    
    vc.cast = vc.cast[,!colnames(vc.cast)[colnames(vc.cast) %in% c('Amp', 'Del')], with =FALSE]
    vc.cast[,total:=rowSums(vc.cast[,2:ncol(vc.cast), with = FALSE])]
    
    vc.cast = merge(vc.cast, vc.cast.cnv, by = 'Tumor_Sample_Barcode', all = TRUE)[order(total, CNV_total, decreasing = TRUE)]
    
    vc.mean = as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, mean))))
    vc.median = as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, median))))
    
  }else{
    vc.cast = vc.cast[,total:=rowSums(vc.cast[,2:ncol(vc.cast), with = FALSE])][order(total, decreasing = TRUE)]
    
    vc.mean = round(as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, mean)))), 3)
    vc.median = round(as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, median)))), 3)
  }
  
  #summarise and casting by 'Variant_Type'
  vt = maf[,.N, .(Tumor_Sample_Barcode, Variant_Type )]
  vt.cast = data.table::dcast(data = vt, formula = Tumor_Sample_Barcode ~ Variant_Type, value.var = 'N', fill = 0, drop = FALSE)
  
  if(any(colnames(vt.cast) %in% c('CNV'))){
    vt.cast.cnv = vt.cast[,c('Tumor_Sample_Barcode', colnames(vt.cast)[colnames(vt.cast) %in% c('CNV')]), with =FALSE]
    
    vt.cast = vt.cast[,!colnames(vt.cast)[colnames(vt.cast) %in% c('CNV')], with =FALSE]
    vt.cast = vt.cast[,total:=rowSums(vt.cast[,2:ncol(vt.cast), with = FALSE])]
    
    vt.cast = merge(vt.cast, vt.cast.cnv, by = 'Tumor_Sample_Barcode', all = TRUE)[order(total, CNV, decreasing = TRUE)]
  }else{
    vt.cast = vt.cast[,total:=rowSums(vt.cast[,2:ncol(vt.cast), with = FALSE])][order(total, decreasing = TRUE)]
  }
  
  #summarise and casting by 'Hugo_Symbol'
  hs = maf[,.N, .(Hugo_Symbol, Variant_Classification)]
  hs.cast = data.table::dcast(data = hs, formula = Hugo_Symbol ~Variant_Classification, fill = 0, value.var = 'N')
  #----
  if(any(colnames(hs.cast) %in% c('Amp', 'Del'))){
    hs.cast.cnv = hs.cast[,c('Hugo_Symbol', colnames(hs.cast)[colnames(hs.cast) %in% c('Amp', 'Del')]), with = FALSE]
    hs.cast.cnv$CNV_total = rowSums(x = hs.cast.cnv[,2:ncol(hs.cast.cnv), with = FALSE], na.rm = TRUE)
    
    hs.cast = hs.cast[,!colnames(hs.cast)[colnames(hs.cast) %in% c('Amp', 'Del')], with = FALSE]
    hs.cast[,total:=rowSums(hs.cast[,2:ncol(hs.cast), with = FALSE], na.rm = TRUE)]
    
    hs.cast = merge(hs.cast, hs.cast.cnv, by = 'Hugo_Symbol', all = TRUE)[order(total, CNV_total, decreasing = TRUE)]
  }else{
    hs.cast[,total:=rowSums(hs.cast[,2:ncol(hs.cast), with = FALSE])]
    hs.cast = hs.cast[order(total, decreasing = TRUE)]
  }
  #----
  
  #Get in how many samples a gene ismutated
  numMutatedSamples = maf[!Variant_Type %in% 'CNV', .(MutatedSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
  numAlteredSamples = maf[, .(AlteredSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
  numAlteredSamples = merge(numMutatedSamples, numAlteredSamples, by = 'Hugo_Symbol', all = TRUE)
  #Merge and sort
  hs.cast = merge(hs.cast, numAlteredSamples, by = 'Hugo_Symbol', all = TRUE)[order(MutatedSamples, total, decreasing = TRUE)]
  #Replace NAs with 0
  hs.cast$AlteredSamples = ifelse(test = is.na(x = hs.cast$AlteredSamples), yes = 0, no = hs.cast$AlteredSamples)
  hs.cast$MutatedSamples = ifelse(test = is.na(x = hs.cast$MutatedSamples), yes = 0, no = hs.cast$MutatedSamples)
  #Make a summarized table
  summary = data.table::data.table(ID = c('NCBI_Build', 'Center','Samples', 'nGenes',colnames(vc.cast)[2:ncol(vc.cast)]),
                                   summary = c(NCBI_Build, Center, nSamples, nGenes, colSums(vc.cast[,2:ncol(vc.cast), with =FALSE])))
  summary[,Mean := vc.mean]
  summary[,Median := vc.median]
  
  # if(chatty){
  #   print(summary)
  #
  #   cat("Gene Summary:\n")
  #   print(hs.cast)
  # }
  
  #Check for flags.
  if(nrow(hs.cast) > 10){
    topten = as.character(hs.cast[1:10, Hugo_Symbol])
    topten = topten[topten %in% flags]
    if(chatty){
      if(length(topten) > 0){
        cat('--Possible FLAGS among top ten genes:\n')
        for(temp in topten){
          cat(paste0("  ", temp, "\n"))
        }
      }
    }
  }
  
  if(chatty){
    cat("-Processing clinical data\n")
  }
  
  if(is.null(anno)){
    if(chatty){
      cat("--Missing clinical data\n")
    }
    sample.anno = data.table::data.table(Tumor_Sample_Barcode = maf.tsbs)
  }else if(is.data.frame(x = anno)){
    sample.anno = data.table::copy(x = anno)
    data.table::setDT(sample.anno)
    if(!'Tumor_Sample_Barcode' %in% colnames(sample.anno)){
      message(paste0('Available fields in provided annotations..'))
      print(colnames(sample.anno))
      stop(paste0('Tumor_Sample_Barcode column not found in provided clinical data. Rename column containing sample names to Tumor_Sample_Barcode if necessary.'))
    }
  }else{
    if(file.exists(anno)){
      sample.anno = data.table::fread(anno, stringsAsFactors = FALSE)
      if(!'Tumor_Sample_Barcode' %in% colnames(sample.anno)){
        message(paste0('Available fields in ', basename(anno), '..'))
        print(colnames(sample.anno))
        stop(paste0('Tumor_Sample_Barcode column not found in provided clinical data. Rename column name containing sample names to Tumor_Sample_Barcode if necessary.'))
      }
    }
  }
  
  #clean up annotation data
  colnames(sample.anno) = gsub(pattern = ' ', replacement = '_', x = colnames(sample.anno), fixed = TRUE) #replace spaces in column names for annotation data
  sample.anno = as.data.frame(apply(sample.anno, 2, function(y) trimws(y))) #remove trailing whitespaces
  sample.anno[sample.anno == ""] = NA #Replace blanks with NA
  #sample.anno = as.data.frame(apply(sample.anno, 2, function(y) gsub(pattern = " ", replacement = "_", x = y))) #replace spaces with _
  data.table::setDT(x = sample.anno)
  if(ncol(sample.anno) == 1){
    colnames(sample.anno)[1] = c("Tumor_Sample_Barcode")
  }
  
  sample.anno = sample.anno[!duplicated(Tumor_Sample_Barcode)] #sample.anno[Tumor_Sample_Barcode %in% maf.tsbs]
  anno.tsbs = sample.anno[,Tumor_Sample_Barcode]
  
  if(!length(maf.tsbs[!maf.tsbs %in% anno.tsbs]) == 0){
    if(chatty){
      cat('--Annotation missing for below samples in MAF:\n')
      for(temp in maf.tsbs[!maf.tsbs %in% anno.tsbs]){
        cat(paste0("  ", temp, "\n"))
      }
    }
  }
  sample.anno = sample.anno[Tumor_Sample_Barcode %in% maf.tsbs]
  
  return(list(variants.per.sample = tsb, variant.type.summary = vt.cast, variant.classification.summary = vc.cast,
              gene.summary = hs.cast, summary = summary, sample.anno = sample.anno))
}

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

flags = function(top = NULL){
  top100flags = c("TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG", "MUC5B",
                  "DNAH17", "PLEC", "DST", "SYNE2", "NEB", "HSPG2", "LAMA5", "AHNAK",
                  "HMCN1", "USH2A", "DNAH11", "MACF1", "MUC17", "DNAH5", "GPR98",
                  "FAT1", "PKD1", "MDN1", "RNF213", "RYR1", "DNAH2", "DNAH3", "DNAH8",
                  "DNAH1", "DNAH9", "ABCA13", "APOB", "SRRM2", "CUBN", "SPTBN5",
                  "PKHD1", "LRP2", "FBN3", "CDH23", "DNAH10", "FAT4", "RYR3", "PKHD1L1",
                  "FAT2", "CSMD1", "PCNT", "COL6A3", "FRAS1", "FCGBP", "DNAH7",
                  "RP1L1", "PCLO", "ZFHX3", "COL7A1", "LRP1B", "FAT3", "EPPK1",
                  "VPS13C", "HRNR", "MKI67", "MYO15A", "STAB1", "ZAN", "UBR4",
                  "VPS13B", "LAMA1", "XIRP2", "BSN", "KMT2C", "ALMS1", "CELSR1",
                  "TG", "LAMA3", "DYNC2H1", "KMT2D", "BRCA2", "CMYA5", "SACS",
                  "STAB2", "AKAP13", "UTRN", "VWF", "VPS13D", "ANK3", "FREM2",
                  "PKD1L1", "LAMA2", "ABCA7", "LRP1", "ASPM", "MYOM2", "PDE4DIP",
                  "TACC2", "MUC2", "TEP1", "HELZ2", "HERC2", "ABCA4")
  
  if(is.null(top)){
    top100flags
  }else{
    top100flags[1:top]
  }
}
