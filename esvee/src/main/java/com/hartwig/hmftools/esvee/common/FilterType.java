package com.hartwig.hmftools.esvee.common;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;

import java.util.Set;
import java.util.stream.Collectors;

public enum FilterType
{
    MIN_AF("minAF", "Variant allele fraction too low", false),
    MIN_SUPPORT("minSupport", "Below minimum fragment support", false),
    MIN_QUALITY("minQual", "Insufficient quality", false),
    STRAND_BIAS("strandBias", "Outside valid strand bias range", false),
    SGL_STRAND_BIAS("sglStrandBias", "Single breakend with excessive strand bias and not a mobile line insertion", false),
    MIN_LENGTH("minLength", "Variant is too short", false),
    SHORT_FRAG_LENGTH("shortFrags", "Average variant fragment length is too short", false),
    DUPLICATE("dedup", "Event is duplicate of another", false),
    SGL("sgl", "SGLs filtered entirely", false),
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
