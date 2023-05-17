library(AnnotationHub)
ah = AnnotationHub()
query(ah, c("EpigenomeRoadMap", "segmentations"))
wanted <- query(ah, c("EpigenomeRoadMap", "segmentations", 'E075'))
codes = c('ChromHMM' = 'AH46928')
# Fetch ah_codes from AnnotationHub and create annotations annotatr understands
build_ah_annots(genome = 'hg38', ah_codes = codes, annotation_class = 'E075_chromhmm')
