package com.hartwig.hmftools.gripss.filters;

public enum FilterType
{
    PASS("PASS", "Variant passes all filters"),
    HARD_FILTERED("hardFiltered", ""),
    MIN_NORMAL_COVERAGE("minNormalCoverage", "Insufficient normal coverage to determine somatic status"),
    MAX_NORMAL_RELATIVE_SUPPORT("maxNormalRelativeSupport", "Too many support reads from the normal sample relative to the tumor"),
    MIN_TUMOR_AF("minTumorAF", "Variant allele fraction too low"),
    SHORT_DEL_INS_ARTIFACT("smallDelInsertionArtifact", "Filter any short DEL where the insert sequence length + 1 = deletion length"),
    MIN_QUAL("minQual", "Insufficient quality"),
    IMPRECISE("imprecise", "Imprecise variant"),
    MAX_POLY_G_LENGTH("maxPolyGLength", "Single breakend containing long polyC or polyG run. Likely to be an artifact"),
    MAX_POLY_A_HOM_LENGTH("maxPolyAHomLength", "Homology containing long polyA or polyT run"),
    MAX_HOM_LENGTH_SHORT_INV("maxHomLengthShortInv", "Short inversion with significant sequence homology"),
    SHORT_SR_SUPPORT("shortSRTumorSupport", "Short event not supported by any split reads either directly or via assembly"),
    SHORT_SR_NORMAL("shortSRNormalSupport", "Short event with split reads support in the normal sample"),
    SHORT_STRAND_BIAS("shortStrandBias", "Short event with excessive strand bias in split reads/soft clipped reads overlapping breakpoint"),
    SGL_STRAND_BIAS("sglStrandBias", "Single breakend with excessive strand bias and not a mobile line insertion"),
    SGL_INSERT_SEQ_MIN_LENGTH("sglInsertSequenceMinLength", "Single breakend with too short insert sequence"),
    MIN_LENGTH("minLength", "Event is too short"),
    DISCORDANT_PAIR_SUPPORT("discordantPairSupport", "Large event not supported by any read pairs either directly or via assembly"),
    DEDUP("dedup", "Event is duplicate of another"),
    QUAL_PER_AD("qualPerAD", "Qual per AD below threshold"),
    MODIFIED_AF("modifiedAF", "Modified AF below threshol"),
    PON("PON", "Found in panel of normals");

    private final String mVcfTag;
    private final String mVcfDesc;

    FilterType(final String vcfTag, final String vcfDesc)
    {
        mVcfTag = vcfTag;
        mVcfDesc = vcfDesc;
    }

    public String vcfTag() { return mVcfTag; }
    public String vcfDesc() { return mVcfDesc; }
}
