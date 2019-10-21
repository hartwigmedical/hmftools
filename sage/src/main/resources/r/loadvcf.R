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

    normalAD = unlist(lapply(vcf.ad[ ,1], paste0, collapse=","))
    tumorAD = unlist(lapply(vcf.ad[ ,2], paste0, collapse=","))

    normalQual = vcf.geno$QUAL[ ,1, ]
    tumorQual = vcf.geno$QUAL[ ,2, ]

    normalRCC = vcf.geno$RCC[ ,1, ]
    tumorRCC = vcf.geno$RCC[ ,2, ]

    normalRCQ = vcf.geno$RCQ[ ,1, ]
    tumorRCQ = vcf.geno$RCQ[ ,2, ]

    vcf.df = data.frame(
    chromosome = seqnames(vcf),
    pos = start(vcf),
    ref = as.character(ref(vcf)),
    alt = as.character(vcf.alt),
    preStrelka = vcf.info$PRE_STRELKA,
    postStrelka = vcf.info$POST_STRELKA,
    strelkaFilter = vcf.fixed$FILTER,
    mappability = vcf.info$MAPPABILITY,
    germlinePonCount = vcf.info$GERMLINE_PON_COUNT,
    somaticPonCount = vcf.info$SOMATIC_PON_COUNT,
    normalAD = normalAD,
    tumorAD = tumorAD,
    normalQual = normalQual[, 1],
    normalBaseQual = normalQual[, 2],
    normalMapQual = normalQual[, 3],
    tumorQual = tumorQual[, 1],
    tumorBaseQual = tumorQual[, 2],
    tumorMapQual = tumorQual[, 3],
    normalRCCFull = normalRCC[, 1],
    normalRCCPartial = normalRCC[, 2],
    normalRCCRealigned = normalRCC[, 3],
    tumorRCCFull = tumorRCC[, 1],
    tumorRCCPartial = tumorRCC[, 2],
    tumorRCCRealigned = tumorRCC[, 3],
    tumorRCCShortened = tumorRCC[, 4],
    tumorRCCLengthened = tumorRCC[, 5],
    tumorDistanceFromRef = vcf.info$RDIS,
    tumorDiffFromRef = vcf.info$RDIF,
    microhomology = vcf.info$MH,
    normalSubprimeDepth = vcf.geno$SDP[, 1],
    tumorSubprimeDepth = vcf.geno$SDP[, 2],
    readContext = vcf.info$RC,
    refContext = vcf.info$TNC,
    repeatSequence = vcf.info$REP_S,
    repeatCount = vcf.info$REP_C,
    highConfidence = vcf.info$HC,
    tumorRCImproperPair = tumorRCQ[, 1],
    tumorRCInconsistentChromosome = tumorRCQ[, 2],
    tumorRCExcesssiveSize = tumorRCQ[, 3],
    normalRCImproperPair = normalRCQ[, 1],
    normalRCInconsistentChromosome = normalRCQ[, 2],
    normalRCExcesssiveSize = normalRCQ[, 3],
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

    return (vcf.df)
}
