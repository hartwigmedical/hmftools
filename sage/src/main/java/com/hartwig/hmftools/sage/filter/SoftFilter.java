package com.hartwig.hmftools.sage.filter;

import java.util.Set;

import com.google.common.collect.Sets;

public enum SoftFilter
{
    MIN_TUMOR_QUAL("minTumorQual", "min_tumor_qual", true, false, "Insufficient tumor quality"),
    MIN_TUMOR_VAF("minTumorVAF", "min_tumor_vaf", true, false, "Insufficient tumor VAF"),
    MIN_GERMLINE_DEPTH("minGermlineDepth", "min_germline_depth", false, true, "Insufficient germline depth"),
    MAX_GERMLINE_VAF("maxGermlineVAF", "max_germline_vaf", false, true, "Excess germline VAF"),
    MAX_GERMLINE_REL_RAW_BASE_QUAL(
            "maxGermlineRelRawBaseQual", "max_germline_rel_raw_base_qual", false, true,
            "Excess germline relative quality"),
    MAX_GERMLINE_ALT_SUPPORT(
            "maxGermlineAltSupport", "max_germline_alt_support", false, true, "Excess germline alt support"),
    MIN_AVG_BASE_QUALITY("minAvgBaseQual", "", true, false, "Variant average base quality below limit"),
    STRAND_BIAS("strandBias", "", true, false, "Variant exceeds strand bias limit"),
    FRAGMENT_COORDS("minFragmentCoords", "", true, false, "Insufficient fragment coordinate variation"),
    MAX_EDGE_DISTANCE("maxEdgeDistance", "", true, false, "Variant close to read edge"),
    JITTER("jitter", "", true, false, "Jitter filter");

    private static final Set<String> TUMOR_FILTERS = Sets.newHashSet();
    public static final Set<String> GERMLINE_FILTERS = Sets.newHashSet();

    static
    {
        for(SoftFilter softFilter : SoftFilter.values())
        {
            if(softFilter.mGermline)
            {
                GERMLINE_FILTERS.add(softFilter.filterName());
            }
            if(softFilter.mTumor)
            {
                TUMOR_FILTERS.add(softFilter.filterName());
            }
        }
    }

    private final String mFilter;
    private final String mConfig;
    private final boolean mTumor;
    private final boolean mGermline;
    private final String mVcfDescription;

    SoftFilter(final String filter, final String config, boolean tumor, boolean germline, final String desc)
    {
        mFilter = filter;
        mConfig = config;
        mTumor = tumor;
        mGermline = germline;
        mVcfDescription = desc;
    }

    public String configName() { return mConfig; }
    public String filterName() { return mFilter; }
    public String vcfDescription() { return mVcfDescription; }

    public String toString() { return mFilter; }

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

}

