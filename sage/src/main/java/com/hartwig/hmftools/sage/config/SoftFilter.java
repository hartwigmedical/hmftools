package com.hartwig.hmftools.sage.config;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public enum SoftFilter
{
    MIN_TUMOR_QUAL("minTumorQual", "min_tumor_qual", true, false),
    MIN_TUMOR_VAF("minTumorVAF", "min_tumor_vaf", true, false),
    MIN_GERMLINE_DEPTH("minGermlineDepth", "min_germline_depth", false, true),
    MIN_GERMLINE_DEPTH_ALLOSOME("minGermlineDepthAllosome", "min_germline_depth_allosome",  false, true),
    MAX_GERMLINE_VAF("maxGermlineVAF", "max_germline_vaf", false, true),
    MAX_GERMLINE_REL_RAW_BASE_QUAL("maxGermlineRelRawBaseQual", "max_germline_rel_raw_base_qual", false, true),
    MAX_GERMLINE_ALT_SUPPORT("maxGermlineAltSupport", "max_germline_alt_support", false, true),
    STRAND_BIAS("strandBias", "", true, false);

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
    private final String mConfig;
    private final boolean mTumor;
    private final boolean mGermline;

    SoftFilter(final String filter, final String config, boolean tumor, boolean germline)
    {
        mFilter = filter;
        mConfig = config;
        mTumor = tumor;
        mGermline = germline;
    }

    public String config() { return mConfig; }

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

