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
    localRealignSet = vcf.info$LRS,
    microhomology = vcf.info$MH,
    readContext = vcf.info$RC,
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


    QualFull = qual[, 1],
    QualPartial = qual[, 2],
    QualCore = qual[, 3],
    QualRealigned = qual[, 4],
    QualReference = qual[, 5],
    QualTotal = qual[, 6],

    JitterShortened = JIT[, 1],
    JitterLengthened = JIT[, 2],
    JitterQualPenalty = JIT[, 3],

    ImproperPairCount = IPC,

    stringsAsFactors = F)

  vcf.df[is.na(vcf.df)] <- 0

  vcf.df = vcf.df %>%
    separate(AD,c('RefSupport','AltSupport'),sep=',') %>%
    mutate(
      Qual = pmax(0, QualFull + QualPartial - JitterQualPenalty),
      RefSupport=as.numeric(RefSupport),
      AltSupport=as.numeric(AltSupport))

  names(vcf.df) <- paste0(sample, names(vcf.df))
  names(vcf.df)[1:4] <- c("chromosome", "pos", "ref", "alt")

  return (vcf.df)
}

##### GIAB
sample = readVcf("/Users/jon/hmf/tmp/GIABvsSELFv004.sage.pon.filtered.vcf.gz")
infoDf = info_data_frame(sample)
normalDf = geno_data_frame(sample, 1, "normal")
tumorDf = geno_data_frame(sample, 2, "tumor")
giab = left_join(infoDf, normalDf, by = c("chromosome", "pos", "ref", "alt")) %>%
  left_join(tumorDf, by = c("chromosome", "pos", "ref", "alt"))
save(giab, file = "/Users/jon/hmf/tmp/GIABvsSELFv004.sage.RData")
passing = giab %>% filter(filter == 'PASS')


##### COLO829
sample = readVcf("/Users/jon/hmf/tmp/COLO829v003.sage.pon.filtered.vcf.gz")
infoDf = info_data_frame(sample)
normalDf = geno_data_frame(sample, 1, "normal")
tumorDf = geno_data_frame(sample, 2, "tumor")
colo829 = left_join(infoDf, normalDf, by = c("chromosome", "pos", "ref", "alt")) %>%
  left_join(tumorDf, by = c("chromosome", "pos", "ref", "alt"))
save(colo829, file = "/Users/jon/hmf/tmp/COLO829v003.sage.RData")
passing = colo829 %>% filter(filter == 'PASS')

### COLO829 Multisample
sample = readVcf("/Users/jon/hmf/tmp/COLO829.multisample.sage.pon.filtered.vcf.gz")
info = info_data_frame(sample)
normal = geno_data_frame(sample, 1, "COLO829v003R")
tumor1 = geno_data_frame(sample, 2, "COLO829v001T")
tumor2 = geno_data_frame(sample, 3, "COLO829v002T")
tumor3 = geno_data_frame(sample, 4, "COLO829v003T")
tumor4 = geno_data_frame(sample, 5, "COLO829v004T")
multicolo829 = left_join(info, normal, by = c("chromosome", "pos", "ref", "alt")) %>%
  left_join(tumor1, by = c("chromosome", "pos", "ref", "alt")) %>%
  left_join(tumor2, by = c("chromosome", "pos", "ref", "alt")) %>%
  left_join(tumor3, by = c("chromosome", "pos", "ref", "alt")) %>%
  left_join(tumor4, by = c("chromosome", "pos", "ref", "alt"))
save(multicolo829, file = "~/hmf/tmp/multicolo829.RData")
