library(VariantAnnotation)
library(tidyr)
library(dplyr)
library(GenomicRanges)

snp_data_frame <- function(vcf) {
  vcf.info = info(vcf)
  vcf.alt = CharacterList(alt(vcf))
  vcf.alt[elementNROWS(vcf.alt) > 1L ] <- lapply(vcf.alt[elementNROWS(vcf.alt) > 1L ], paste0, collapse=",")

  vcf.df = data.frame(
    chromosome = seqnames(vcf),
    pos = start(vcf),
    ref = as.character(ref(vcf)), 
    alt = as.character(vcf.alt), stringsAsFactors = F)
  
  snpEffWorst = data.frame(as.matrix(vcf.info$SEW), stringsAsFactors = F)
  names(snpEffWorst) <- c("worstGene", "worstTranscript", "worstEffect", "worstCodingEffect", "genesAffected")
  snpEffWorst = snpEffWorst %>% mutate(genesAffected = as.integer(genesAffected))

  snpEffCanonical = data.frame(as.matrix(vcf.info$SEC), stringsAsFactors = F)
  names(snpEffCanonical) <- c("canonicalGene", "canonicalTranscript", "canonicalEffect", "canonicalCodingEffect", "canonicalHgvsCodingImpact", "canonicalHgvsProteinImpact")

  result = bind_cols(snpEffWorst, snpEffCanonical)
  result = bind_cols(vcf.df, result)
  return (result)
}

info_data_frame <- function(vcf) {
  vcf.info = info(vcf)
  vcf.fixed = fixed(vcf)
  vcf.alt = CharacterList(alt(vcf))
  vcf.alt[elementNROWS(vcf.alt) > 1L ] <- lapply(vcf.alt[elementNROWS(vcf.alt) > 1L ], paste0, collapse=",")
 
  vcf.df = data.frame(
    chromosome = seqnames(vcf),
    pos = start(vcf),
    ref = as.character(ref(vcf)), 
    alt = as.character(vcf.alt),
    filter = vcf.fixed$FILTER,
    qual = qual(vcf),
    tier = vcf.info$TIER,
    localPhaseSet = vcf.info$LPS,
    microhomology = vcf.info$MH,
    readContext = vcf.info$RC,
    readContextDistance = vcf.info$RC_DIS,
    readContextDifference = vcf.info$RC_DIF,
    readContextMicrohomology = vcf.info$RC_MH,
    readContextRepeatSequence = vcf.info$RC_REPS,
    readContextRepeatCount = vcf.info$RC_REPC,
    trinucleotideContext = vcf.info$TNC,
    referenceRepeatSequence = vcf.info$REP_S,
    referenceRepeatCount = vcf.info$REP_C,
    rightAlignedMicrohomology = vcf.info$RAM,
    phasedInframeIndel = vcf.info$PII,
    mixedSomaticGermline = vcf.info$MSG,
    
    #preStrelka = vcf.info$PRE_STRELKA,
    #postStrelka = vcf.info$POST_STRELKA,
    #mappability = vcf.info$MAPPABILITY,
    ponCount = vcf.info$PON_COUNT,
    stringsAsFactors = F) %>%
  mutate(
    ponCount = ifelse(is.na(ponCount), 0, ponCount),
    microhomology = ifelse(is.na(microhomology), "", microhomology),
    readContextMicrohomology = ifelse(is.na(readContextMicrohomology), "", readContextMicrohomology),
    readContextRepeatSequence = ifelse(is.na(readContextRepeatSequence), "", readContextRepeatSequence),
    readContextRepeatCount = ifelse(is.na(readContextRepeatCount), 0, readContextRepeatCount),
    referenceRepeatSequence = ifelse(is.na(referenceRepeatSequence), "", referenceRepeatSequence),
    referenceRepeatCount = ifelse(is.na(referenceRepeatCount), 0, referenceRepeatCount),
    localPhaseSet = ifelse(is.na(localPhaseSet), 0, localPhaseSet)
  )
}

geno_data_frame<- function(vcf, index, sample) {
  vcf.alt = CharacterList(alt(vcf))
  vcf.alt[elementNROWS(vcf.alt) > 1L ] <- lapply(vcf.alt[elementNROWS(vcf.alt) > 1L ], paste0, collapse=",")
  vcf.geno = geno(vcf)
  vcf.ad = vcf.geno$AD

  AD = unlist(lapply(vcf.ad[ ,index], paste0, collapse=","))
  IPC = vcf.geno$RC_IPC[ ,index]
  
  qual = vcf.geno$RC_QUAL[ ,index, ]
  RCC = vcf.geno$RC_CNT[ ,index, ]
  JIT = vcf.geno$RC_JIT[ ,index, ]

  vcf.df = data.frame(
    chromosome = seqnames(vcf),
    pos = start(vcf),
    ref = as.character(ref(vcf)), 
    alt = as.character(vcf.alt),

    AF = vcf.geno$AF[, index],
    AD = AD,
    Depth = vcf.geno$DP[, index],
    
    RawDP = vcf.geno$RDP[, index],
    RawAD0 = vcf.geno$RAD[, index, 1],
    RawAD1 = vcf.geno$RAD[, index, 2],
    RawBQ0 = vcf.geno$RABQ[, index, 1],
    RawBQ1 = vcf.geno$RABQ[, index, 2],
    
    RCCFull = RCC[, 1],
    RCCPartial = RCC[, 2],
    RCCCore = RCC[, 3],
    RCCRealigned = RCC[, 4],
    RCCReference = RCC[, 5],
    RCCTotal = RCC[, 6],
    RCCShortened = RCC[, 1],
    RCCLengthened = RCC[, 2],
    
    QualFull = qual[, 1],
    QualPartial = qual[, 2],
    QualCore = qual[, 3],
    QualRealigned = qual[, 4],
    QualReference = qual[, 5],
    QualTotal = qual[, 6],
    QualJitterPenalty = JIT[, 3],

    ImproperPairCount = IPC,
    
    stringsAsFactors = F)
  
  vcf.df[is.na(vcf.df)] <- 0
  
  vcf.df = vcf.df %>% 
    separate(AD,c('RefSupport','AltSupport'),sep=',') %>% 
    mutate(
      Qual = pmax(0, QualFull + QualPartial - QualJitterPenalty),
      RefSupport=as.numeric(RefSupport),
      AltSupport=as.numeric(AltSupport))
  
  names(vcf.df) <- paste0(sample, names(vcf.df))
  names(vcf.df)[1:4] <- c("chromosome", "pos", "ref", "alt")
  
  return (vcf.df)
}
