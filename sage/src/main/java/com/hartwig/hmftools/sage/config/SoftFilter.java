package com.hartwig.hmftools.sage.config;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public enum SoftFilter
{
    MIN_TUMOR_QUAL("minTumorQual", true, false),
    MIN_TUMOR_VAF("minTumorVAF", true, false),
    MIN_GERMLINE_DEPTH("minGermlineDepth", false, true),
    MIN_GERMLINE_DEPTH_ALLOSOME("minGermlineDepthAllosome", false, true),
    MAX_GERMLINE_VAF("maxGermlineVAF", false, true),
    MAX_GERMLINE_REL_RAW_BASE_QUAL("maxGermlineRelRawBaseQual", false, true),
    MAX_GERMLINE_ALT_SUPPORT("maxGermlineAltSupport", false, true),
    STRAND_BIAS("strandBias", true, false);

    private static final Set<String> TUMOR_FILTERS = Sets.newHashSet();
    private static final Set<String> GERMLINE_FILTERS = Sets.newHashSet();

    static
    {
        for(SoftFilter softFilter : SoftFilter.values())
        {
            if(softFilter.mGermline)
            {
                GERMLINE_FILTERS.add(softFilter.toString());
            }
            if(softFilter.mTumor)
            {
                TUMOR_FILTERS.add(softFilter.toString());
            }
        }
    }

    private final String mFilter;
    private final boolean mTumor;
    private final boolean mGermline;

    SoftFilter(@NotNull final String filter, boolean tumor, boolean germline)
    {
        mFilter = filter;
        mTumor = tumor;
        mGermline = germline;
    }

    public static boolean isGermlineAndNotTumorFiltered(final Set<String> softFilters)
    {
        for(final String softFilter : softFilters)
        {
            if(TUMOR_FILTERS.contains(softFilter))
                return false;

            if(!GERMLINE_FILTERS.contains(softFilter))
                return false;
        }

        return !softFilters.isEmpty();
    }

    @Override
    public String toString()
    {
        return mFilter;
    }
}

