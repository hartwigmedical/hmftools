package com.hartwig.hmftools.sage.filter;

import java.util.Set;

import com.google.common.collect.Sets;

public enum SoftFilter
{
    MIN_TUMOR_QUAL("minTumorQual", "min_tumor_qual", true, false, "Insufficient tumor quality"),
    MIN_TUMOR_VAF("minTumorVAF", "min_tumor_vaf", true, false, "Insufficient tumor VAF"),
    MIN_GERMLINE_DEPTH("minGermlineDepth", "min_germline_depth", false, true, "Insufficient germline depth"),
    MAX_GERMLINE_VAF("maxGermlineVAF", "max_germline_vaf", false, true, "Excess germline VAF"),
    MAX_GERMLINE_RELATIVE_QUAL(
            "maxGermlineRelQual", "max_germline_rel_qual", false, true, "Excess germline relative qual"),
    MAX_GERMLINE_ALT_SUPPORT(
            "maxGermlineAltSupport", "max_germline_alt_support", false, true, "Excess germline alt support"),
    MIN_AVG_BASE_QUALITY("minAvgBaseQual", "", true, false, "Variant average base quality below limit"),
    FRAGMENT_STRAND_BIAS("fragmentStrandBias", "", true, false, "Variant exceeds fragment strand bias limit"),
    READ_STRAND_BIAS("readStrandBias", "", true, false, "Variant exceeds read strand bias limit"),
    FRAGMENT_COORDS("minFragmentCoords", "", true, false, "Insufficient fragment coordinate variation"),
    MAX_EDGE_DISTANCE("maxEdgeDistance", "", true, false, "Variant close to read edge"),
    MAP_QUAL_REF_ALT_DIFFERENCE(
            "mapQualRefAltDiff", "", true, false, "Alt support map qual well below ref support"),
    JITTER("jitter", "", true, false, "Jitter filter"),
    DEDUP_MNV("dedupMnv", "", true, false, "Filter MNV duplicate"),
    DEDUP_MIXED_GERMLINE_SOMATIC(
            "dedupMixedGermlineSomatic", "", true, false, "Variant duplicate mixed somatic/germline"),
    DEDUP_SNV_MNV("dedupSnvMnv", "", true, false, "Variant duplicate MNV vs SNV"),
    DEDUP_INDEL("dedupIndel", "", true, false, "Variant duplicate SNV/MNV vs INDEL"),
    DEDUP_MATCH("dedupMatch", "", true, false, "Variant duplicate with different read contexts");

    private static final Set<SoftFilter> TUMOR_FILTERS = Sets.newHashSet();
    public static final Set<SoftFilter> GERMLINE_FILTERS = Sets.newHashSet();

    static
    {
        for(SoftFilter softFilter : SoftFilter.values())
        {
            if(softFilter.mGermline)
            {
                GERMLINE_FILTERS.add(softFilter);
            }
            if(softFilter.mTumor)
            {
                TUMOR_FILTERS.add(softFilter);
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

    public static boolean isGermlineAndNotTumorFiltered(final Set<SoftFilter> softFilters)
    {
        for(SoftFilter softFilter : softFilters)
        {
            if(TUMOR_FILTERS.contains(softFilter))
                return false;

            if(!GERMLINE_FILTERS.contains(softFilter))
                return false;
        }

        return !softFilters.isEmpty();
    }
}

