package com.hartwig.hmftools.gripss.filters;

public enum FilterType
{
    PASS,
    HARD_FILTERED,
    MIN_NORMAL_COVERAGE,
    MAX_NORMAL_RELATIVE_SUPPORT,
    MIN_TUMOR_AF,
    SHORT_DEL_INS_ARTIFACT,
    MIN_QUAL,
    IMPRECISE,
    MAX_POLY_G_LENGTH,
    MAX_POLY_A_HOM_LENGTH,
    MAX_HOM_LENGTH_SHORT_INV,
    SHORT_SR_SUPPORT,
    SHORT_SR_NORMAL,
    SHORT_STRAND_BIAS,
    SGL_STRAND_BIAS,
    SGL_INSERT_SEQ_MIN_LENGTH,
    MIN_LENGTH,
    DISCORDANT_PAIR_SUPPORT,
    DEDUP,
    PON;

    public static String vcfName(final FilterType filter)
    {
        switch(filter)
        {
            case PASS: return "PASS";
            case PON: return "PON";
            case DEDUP: return "dedup";
            case MIN_LENGTH: return "minLength";
            case MAX_HOM_LENGTH_SHORT_INV: return "maxHomLengthShortInv";
            case SHORT_SR_SUPPORT: return "shortSRTumorSupport";
            case SHORT_SR_NORMAL: return "shortSRNormalSupport";
            case DISCORDANT_PAIR_SUPPORT: return "discordantPairSupport";
            case MIN_TUMOR_AF: return "minTumorAF";
            case MAX_NORMAL_RELATIVE_SUPPORT: return "maxNormalRelativeSupport";
            case MIN_NORMAL_COVERAGE: return "minNormalCoverage";
            case SHORT_STRAND_BIAS: return "shortStrandBias";
            case SGL_STRAND_BIAS: return "sglStrandBias";
            case SGL_INSERT_SEQ_MIN_LENGTH: return "sglInsertSequenceMinLength";
            case MIN_QUAL: return "minQual";
            case MAX_POLY_G_LENGTH: return "maxPolyGLength";
            case MAX_POLY_A_HOM_LENGTH: return "maxPolyAHomLength";
            case IMPRECISE: return "imprecise";
            case SHORT_DEL_INS_ARTIFACT: return "smallDelInsertionArtifact";
            default: return "unknown";
        }
    }

    public static String vcfInfoString(final FilterType filter)
    {
        switch(filter)
        {
            case DEDUP:
                return "Event is duplicate of another";
            case MIN_LENGTH:
                return "Event is too short";
            case MAX_HOM_LENGTH_SHORT_INV:
                return "Short inversion with significant sequence homology";
            case SHORT_SR_SUPPORT:
                return "Short event not supported by any split reads either directly or via assembly";
            case SHORT_SR_NORMAL:
                return "Short event with split reads support in the normal sample";
            case DISCORDANT_PAIR_SUPPORT:
                return "Large event not supported by any read pairs either directly or via assembly";
            case PON:
                return "Found in panel of normals";
            case MIN_TUMOR_AF:
                return "Variant allele fraction too low";
            case MAX_NORMAL_RELATIVE_SUPPORT:
                return "Too many support reads from the normal sample relative to the tumor";
            case MIN_NORMAL_COVERAGE:
                return "Insufficient normal coverage to determine somatic status";
            case SHORT_STRAND_BIAS:
                return "Short event with excessive strand bias in split reads/soft clipped reads overlapping breakpoint";
            case SGL_STRAND_BIAS:
                return "Single breakend with excessive strand bias and not a mobile line insertion";
            case SGL_INSERT_SEQ_MIN_LENGTH:
                return "Single breakend with too short insert sequence";
            case MIN_QUAL:
                return "Insufficient quality";
            case MAX_POLY_G_LENGTH:
                return "Single breakend containing long polyC or polyG run. Likely to be an artifact";
            case MAX_POLY_A_HOM_LENGTH:
                return "Homology containing long polyA or polyT run";
            case IMPRECISE:
                return "Imprecise variant";
            case PASS:
                return "Variant passes all filters";
            case SHORT_DEL_INS_ARTIFACT:
                return "Filter any short DEL where the insert sequence length + 1 = deletion length. This is a known GRIDSS artifact";
            default: return "unknown";
        }
    }

}
