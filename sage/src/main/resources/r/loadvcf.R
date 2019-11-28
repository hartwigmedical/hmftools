library(VariantAnnotation)
library(tidyr)
library(dplyr)

vcf_data_frame<- function(vcf) {
    vcf.rowRanges = rowRanges(vcf)
    vcf.info = info(vcf)
    vcf.alt = CharacterList(alt(vcf))
    vcf.qual = qual(vcf)
    vcf.alt[elementNROWS(vcf.alt) > 1L ] <- lapply(vcf.alt[elementNROWS(vcf.alt) > 1L ], paste0, collapse=",")
    vcf.geno = geno(vcf)
    vcf.ad = vcf.geno$AD
    vcf.fixed = fixed(vcf)
    ann.df = as.data.frame(vcf.info$ANN) %>% group_by(group) %>% filter(row_number()== 1) %>% separate(value, sep = "\\|", into = c("alt", "impact","severity", "gene")) %>%
      select(group, impact, gene)
  
  
    normalAD = unlist(lapply(vcf.ad[ ,1], paste0, collapse=","))
    tumorAD = unlist(lapply(vcf.ad[ ,2], paste0, collapse=","))

    tumorIPC = vcf.geno$RC_IPC[ ,2]
  
    normalQual = vcf.geno$RC_QUAL[ ,1, ]
    tumorQual = vcf.geno$RC_QUAL[ ,2, ]
    
    normalRCC = vcf.geno$RC_CNT[ ,1, ]
    tumorRCC = vcf.geno$RC_CNT[ ,2, ]

    vcf.df = data.frame(
      chromosome = seqnames(vcf),
      pos = start(vcf),
      ref = as.character(ref(vcf)), 
      alt = as.character(vcf.alt),
      filter = vcf.fixed$FILTER,
    
      tier = vcf.info$TIER,
      normalAD = normalAD,
      normalDP = vcf.geno$DP[, 1],
      tumorAD = tumorAD,
      tumorDP = vcf.geno$DP[, 2],
      tumorAF = vcf.info$AF,
      #normalQual = vcf.qual
      normalQual = normalQual[, 1],
      normalBaseQual = normalQual[, 2],
      normalMapQual = normalQual[, 3],
      tumorQual = vcf.qual,
      tumorBaseQual = tumorQual[, 2],
      tumorMapQual = tumorQual[, 3],
      jitterPenalty = tumorQual[, 4],
      normalRCCFull = normalRCC[, 1],
      normalRCCPartial = normalRCC[, 2],
      normalRCCRealigned = normalRCC[, 3],
      normalRCCCoverage = normalRCC[, 6],
      tumorRCCFull = tumorRCC[, 1],
      tumorRCCPartial = tumorRCC[, 2],
      tumorRCCRealigned = tumorRCC[, 3],
      tumorRCCShortened = tumorRCC[, 4],
      tumorRCCLengthened = tumorRCC[, 5],
      tumorRCCCoverage = tumorRCC[, 6],
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
      #jitterPenalty = vcf.info$JITTER,
      stringsAsFactors = F)
    
    vcf.df[is.na(vcf.df)] <- 0

    vcf.df = vcf.df %>%
        separate(normalAD,c('germlineRefSupport','germlineAltSupport'),sep=',') %>%
        separate(tumorAD,c('tumorRefSupport','tumorAltSupport'),sep=',') %>%
        mutate(
        tumorRefSupport=as.numeric(tumorRefSupport),
        tumorAltSupport=as.numeric(tumorAltSupport),
        germlineRefSupport=as.numeric(germlineRefSupport),
        germlineAltSupport=as.numeric(germlineAltSupport),
        readContextCentre = ifelse(nchar(readContext) > 50, substring(readContext, 26, nchar(readContext) - 25), "")
        )

    result = vcf.df %>% mutate(group = row_number()) %>% left_join(ann.df, by = "group")
    
    return (result)
}

sagePanelResult = data.frame()
for (file in list.files(path = "/Users/jon/hmf/analysis/sagePanel/", pattern = "*vcf.gz$")) {
  cat (file, "\n")
  vcf = readVcf(paste0("~/hmf/analysis/sagePanel/", file))
  vcfDF = vcf_data_frame(vcf) %>% mutate(sample = substr(file,1, 13))    
  sagePanelResult = bind_rows(sagePanelResult, vcfDF)
}

save(sagePanelResult, file = "/Users/jon/hmf/analysis/sagePanel/sagePanelResult.RData")

