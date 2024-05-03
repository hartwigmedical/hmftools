package com.hartwig.hmftools.esvee.common;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;

import java.util.Set;
import java.util.stream.Collectors;

public enum FilterType
{
    HARD_FILTERED("hardFiltered", ""),
    MIN_NORMAL_COVERAGE("minNormalCoverage", "Insufficient normal coverage to determine somatic status"),
    MAX_NORMAL_RELATIVE_SUPPORT("maxNormalRelativeSupport", "Too many support reads from the normal sample relative to the tumor"),
    MIN_TUMOR_AF("minTumorAF", "Variant allele fraction too low"),
    MIN_QUALITY("minQual", "Insufficient quality"),
    MAX_POLY_A_HOM_LENGTH("maxPolyAHomLength", "Homology containing long polyA or polyT run"),
    MAX_HOM_LENGTH_SHORT_INV("maxHomLengthShortInv", "Short inversion with significant sequence homology"),
    SHORT_SR_SUPPORT("shortSRTumorSupport", "Short event not supported by any split reads either directly or via assembly"),
    SHORT_SR_NORMAL("shortSRNormalSupport", "Short event with split reads support in the normal sample"),
    STRAND_BIAS("strandBias", "Outside valid strand bias range"),
    SHORT_STRAND_BIAS("shortStrandBias", "Short event with excessive strand bias in split reads/soft clipped reads overlapping breakpoint"),
    SGL_STRAND_BIAS("sglStrandBias", "Single breakend with excessive strand bias and not a mobile line insertion"),
    SGL_INSERT_SEQ_MIN_LENGTH("sglInsertSequenceMinLength", "Single breakend with too short insert sequence"),
    MIN_LENGTH("minLength", "Event is too short"),
    DUPLICATE("dedup", "Event is duplicate of another"),
    MODIFIED_AF("modifiedAF", "Modified AF below threshold"),
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

    public static String filtersAsStr(final Set<FilterType> filters)
    {
        if(filters.isEmpty())
            return PASS;

        return filters.stream().map(x -> x.vcfTag()).collect(Collectors.joining(";"));
    }
}
