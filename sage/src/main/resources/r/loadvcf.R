library(VariantAnnotation)
library(tidyr)
library(dplyr)
library(GenomicRanges)

vcf = sample

rna_data_frame<- function(vcf) {
  vcf.info = info(vcf)
  vcf.alt = CharacterList(alt(vcf))
  vcf.alt[elementNROWS(vcf.alt) > 1L ] <- lapply(vcf.alt[elementNROWS(vcf.alt) > 1L ], paste0, collapse=",")
  
  rnaAD = unlist(lapply(vcf.info$RNA_AD, paste0, collapse=","))
  
  vcf.df = data.frame(
    chromosome = seqnames(vcf),
    pos = start(vcf),
    ref = as.character(ref(vcf)), 
    alt = as.character(vcf.alt),
    rnaDepth = vcf.info$RNA_DP,
    rnaAD = rnaAD,
    stringsAsFactors = F)
  
  vcf.df = vcf.df %>% 
    separate(rnaAD,c('rnaRefSupport', 'rnaAltSupport'), sep = ",") %>%
    mutate(rnaRefSupport = ifelse(rnaRefSupport == "NA", NA, as.integer(rnaRefSupport)), rnaAltSupport = ifelse(rnaAltSupport == "NA", NA, as.integer(rnaAltSupport)))

  vcf.df[is.na(vcf.df)] <- 0
  
  return (vcf.df)
}


dna_data_frame<- function(vcf) {
  vcf.rowRanges = rowRanges(vcf)
  vcf.info = info(vcf)
  vcf.alt = CharacterList(alt(vcf))
  vcf.qual = qual(vcf)
  vcf.alt[elementNROWS(vcf.alt) > 1L ] <- lapply(vcf.alt[elementNROWS(vcf.alt) > 1L ], paste0, collapse=",")
  vcf.geno = geno(vcf)
  vcf.ad = vcf.geno$AD
  vcf.fixed = fixed(vcf)
  
  normalAD = unlist(lapply(vcf.ad[ ,1], paste0, collapse=","))
  tumorAD = unlist(lapply(vcf.ad[ ,2], paste0, collapse=","))
  tumorIPC = vcf.geno$RC_IPC[ ,2]
  
  normalQual = vcf.geno$RC_QUAL[ ,1, ]
  tumorQual = vcf.geno$RC_QUAL[ ,2, ]
  
  normalRCC = vcf.geno$RC_CNT[ ,1, ]
  tumorRCC = vcf.geno$RC_CNT[ ,2, ]
  
  normalJIT = vcf.geno$RC_JIT[ ,1, ]
  tumorJIT = vcf.geno$RC_JIT[ ,2, ]
  
  vcf.df = data.frame(
    chromosome = seqnames(vcf),
    pos = start(vcf),
    ref = as.character(ref(vcf)), 
    alt = as.character(vcf.alt),
    filter = vcf.fixed$FILTER,
    
    #preStrelka = vcf.info$PRE_STRELKA,
    #postStrelka = vcf.info$POST_STRELKA,
    #highConfidence = vcf.info$HC,
    #mappability = vcf.info$MAPPABILITY,
    #ponCount = vcf.info$PON_COUNT,

        tier = vcf.info$TIER,
    normalAD = normalAD,
    tumorAD = tumorAD,
    tumorAF = vcf.info$AF,
    
    rawNormalAltBaseQual = vcf.geno$RABQ[, 1, 2],
    #rawTumorAltBaseQual = vcf.geno$RABQ[, 2, 2],
    
    rawNormalDP = vcf.geno$RDP[, 1],
    rawNormalAD0 = vcf.geno$RAD[, 1, 1],
    rawNormalAD1 = vcf.geno$RAD[, 1, 2],
    
    rawTumorDP = vcf.geno$RDP[, 2],
    rawTumorAD0 = vcf.geno$RAD[, 2, 1],
    rawTumorAD1 = vcf.geno$RAD[, 2, 2],

    tumorRCCFull = tumorRCC[, 1],
    tumorRCCPartial = tumorRCC[, 2],
    tumorRCCCore = tumorRCC[, 3],
    tumorRCCRealigned = tumorRCC[, 4],
    tumorRCCReference = tumorRCC[, 5],
    tumorRCCTotal = tumorRCC[, 6],
    tumorRCCShortened = tumorJIT[, 1],
    tumorRCCLengthened = tumorJIT[, 2],

    tumorQual = vcf.qual,
    tumorQualFull = tumorQual[, 1],
    tumorQualPartial = tumorQual[, 2],
    tumorQualCore = tumorQual[, 3],
    tumorQualRealigned = tumorQual[, 4],
    tumorQualReference = tumorQual[, 5],
    tumorQualTotal = tumorQual[, 6],
    tumorQualJitterPenalty = tumorJIT[, 3],

    normalRCCFull = normalRCC[, 1],
    normalRCCPartial = normalRCC[, 2],
    normalRCCCore = normalRCC[, 3],
    normalRCCRealigned = normalRCC[, 4],
    normalRCCReference = normalRCC[, 5],
    normalRCCTotal = normalRCC[, 6],
    normalRCCShortened = normalJIT[, 1],
    normalRCCLengthened = normalJIT[, 2],

    normalQualFull = normalQual[, 1],
    normalQualPartial = normalQual[, 2],
    normalQualCore = normalQual[, 3],
    normalQualRealigned = normalQual[, 4],
    normalQualReference = normalQual[, 5],
    normalQualTotal = normalQual[, 6],
    normalQualJitterPenalty = normalJIT[, 3],

    tumorDistanceFromRef = vcf.info$RC_DIS,
    tumorDiffFromRef = vcf.info$RC_DIF,
    microhomology = vcf.info$MH,
    readContextMicrohomology = vcf.info$RC_MH,
    readContext = vcf.info$RC,
    refContext = vcf.info$TNC,
    repeatSequence = vcf.info$REP_S,
    readContextRepeatSequence = vcf.info$RC_REPS,
    repeatCount = vcf.info$REP_C,
    readContextRepeatCount = vcf.info$RC_REPC,

    phase = vcf.info$LPS,
    tumorImproperPairCount = tumorIPC,
    stringsAsFactors = F)
  
  vcf.df[is.na(vcf.df)] <- 0
  
  vcf.df = vcf.df %>% 
    separate(normalAD,c('germlineRefSupport','germlineAltSupport'),sep=',') %>% 
    separate(tumorAD,c('tumorRefSupport','tumorAltSupport'),sep=',') %>%
    mutate(
      normalQual = normalQualFull + normalQualPartial + normalQualCore +  normalQualRealigned - normalQualJitterPenalty,
      normalAF = (normalQual) / normalQualTotal,
      tumorRefSupport=as.numeric(tumorRefSupport),
      tumorAltSupport=as.numeric(tumorAltSupport),
      germlineRefSupport=as.numeric(germlineRefSupport),
      germlineAltSupport=as.numeric(germlineAltSupport),
      readContextCentre = readContext)
  
  return (vcf.df)
}

sample = readVcf("/Users/jon/hmf/analysis/x/x.rna.final.vcf.gz")
dnaDF = dna_data_frame(sample)
rnaDf = rna_data_frame(sample)
combined = left_join(rnaDf, dnaDF, by = c("chromosome", "pos", "ref", "alt"))


