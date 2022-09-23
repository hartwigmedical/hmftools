package com.hartwig.hmftools.bammetrics;

public enum FilterType
{
    UNFILTERED("Unfiltered"),
    LOW_MAP_QUAL("LowMapQual"),
    DUPLICATE("Duplicate"),
    MATE_UNMAPPED("MateUnampped"),
    LOW_BASE_QUAL("LowBaseQual"),
    OVERLAPPED("Overlapped"),
    MAX_COVERAGE("MaxCoverage");

    private final String mDescription;

    FilterType(final String desc)
    {
        mDescription = desc;
    }

    public String description() { return mDescription; }
}
