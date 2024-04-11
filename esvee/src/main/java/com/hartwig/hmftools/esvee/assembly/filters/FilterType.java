package com.hartwig.hmftools.esvee.assembly.filters;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;

import java.util.Set;
import java.util.stream.Collectors;

public enum FilterType
{
    MULTIPLE_ASSEMBLIES("multipleAssemblies", "Variant is supported by multiple assemblies"),
    MIN_OVERHANG("minOverhang", "Insufficient read distance over junction"),
    MIN_QUALITY("minQuality", "Below min quality threshold"),
    MIN_SUPPORT("minSupport", "Below min fragment support");

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

