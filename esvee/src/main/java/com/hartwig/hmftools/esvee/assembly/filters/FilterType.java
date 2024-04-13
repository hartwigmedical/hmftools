package com.hartwig.hmftools.esvee.assembly.filters;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;

import java.util.Set;
import java.util.stream.Collectors;

public enum FilterType
{
    STRAND_BIAS("strandBias", "Outside valid strand bias range"),
    MIN_QUALITY("minQual", "Below min quality threshold"),
    DUPLICATE("dedup", "Duplicate of similar or matching breakend"),
    MIN_TUMOR_SUPPORT("minTumorSupport", "Below tumor min fragment support");

    private final String mFilter;
    private final String mVcfDescription;

    FilterType(final String filter, final String desc)
    {
        mFilter = filter;
        mVcfDescription = desc;
    }

    public String filterName() { return mFilter; }
    public String vcfDescription() { return mVcfDescription; }

    public String toString() { return mFilter; }

    public static String filtersAsStr(final Set<FilterType> filters)
    {
        if(filters.isEmpty())
            return PASS;

        return filters.stream().map(x -> x.filterName()).collect(Collectors.joining(";"));
    }
}

