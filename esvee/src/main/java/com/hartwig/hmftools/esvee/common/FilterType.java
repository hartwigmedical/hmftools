package com.hartwig.hmftools.esvee.common;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;

import java.util.Set;
import java.util.stream.Collectors;

public enum FilterType
{
    MIN_NORMAL_COVERAGE("minNormalCoverage", "Insufficient normal coverage to determine somatic status", true),
    MAX_NORMAL_RELATIVE_SUPPORT(
            "maxNormalRelativeSupport", "Too many support reads from the normal sample relative to the tumor", true),
    MIN_AF("minTumorAF", "Variant allele fraction too low", false),
    MIN_QUALITY("minQual", "Insufficient quality", false),
    MAX_POLY_A_HOM_LENGTH("maxPolyAHomLength", "Homology containing long polyA or polyT run", false),
    MAX_HOM_LENGTH_SHORT_INV("maxHomLengthShortInv", "Short inversion with significant sequence homology", false),
    SHORT_SR_SUPPORT("shortSRTumorSupport", "Short event not supported by any split reads either directly or via assembly", false),
    SHORT_SR_NORMAL("shortSRNormalSupport", "Short event with split reads support in the normal sample", false),
    STRAND_BIAS("strandBias", "Outside valid strand bias range", false),
    SHORT_STRAND_BIAS(
            "shortStrandBias",
            "Short event with excessive strand bias in split reads/soft clipped reads overlapping breakpoint", false),
    SGL_STRAND_BIAS("sglStrandBias", "Single breakend with excessive strand bias and not a mobile line insertion", false),
    SGL_INSERT_SEQ_MIN_LENGTH("sglInsertSequenceMinLength", "Single breakend with too short insert sequence", false),
    MIN_LENGTH("minLength", "Variant is too short", false),
    DUPLICATE("dedup", "Event is duplicate of another", false),
    MODIFIED_AF("modifiedAF", "Modified AF below threshold", false),
    PON("PON", "Found in panel of normals", true);

    private final String mVcfTag;
    private final String mVcfDesc;
    private final boolean mGermlineOnly;

    FilterType(final String vcfTag, final String vcfDesc, boolean germlineOnly)
    {
        mVcfTag = vcfTag;
        mVcfDesc = vcfDesc;
        mGermlineOnly = germlineOnly;
    }

    public String vcfTag() { return mVcfTag; }
    public String vcfDesc() { return mVcfDesc; }
    public boolean germlineOnly() { return mGermlineOnly; }

    public static String filtersAsStr(final Set<FilterType> filters)
    {
        if(filters.isEmpty())
            return PASS;

        return filters.stream().map(x -> x.vcfTag()).collect(Collectors.joining(";"));
    }
}
