package com.hartwig.hmftools.sage.filter;

import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FILTERED_MAX_GERMLINE_ALT_SUPPORT;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_HARD_MIN_TUMOR_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_HARD_MIN_TUMOR_VAF;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_HIGH_CONFIDENCE_FILTER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_HOTSPOT_FILTER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_LOW_CONFIDENCE_FILTER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MIN_AVG_BASE_QUALITY;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MIN_AVG_BASE_QUALITY_HOTSPOT;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_PANEL_FILTER;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class FilterConfig
{
    public final boolean DisableHardFilter;
    public final boolean DisableSoftFilter;
    public final int HardMinTumorQual;
    public final double HardMinTumorVaf;
    public final int HardMinTumorRawAltSupport;
    public final int FilteredMaxGermlineAltSupport;
    public final int MinAvgBaseQual;
    public final int MinAvgBaseQualHotspot;
    public final int ReferenceSampleCount;

    public final SoftFilterConfig SoftHotspotFilter;
    public final SoftFilterConfig SoftPanelFilter;
    public final SoftFilterConfig SoftHighConfidenceFilter;
    public final SoftFilterConfig SoftLowConfidenceFilter;

    private static final String DISABLE_SOFT_FILTER = "disable_soft_filter";
    private static final String DISABLE_HARD_FILTER = "disable_hard_filter";

    private static final String HARD_MIN_TUMOR_QUAL = "hard_min_tumor_qual";
    private static final String HARD_MIN_TUMOR_VAF = "hard_min_tumor_vaf";
    private static final String HARD_MIN_TUMOR_RAW_ALT_SUPPORT = "hard_min_tumor_raw_alt_support";
    private static final String FILTERED_MAX_GERMLINE_ALT_SUPPORT = "filtered_max_germline_alt_support";
    private static final String MIN_AVG_BASE_QUAL = "min_avg_base_qual";
    private static final String MIN_AVG_BASE_QUAL_HOTSPOT = "min_avg_base_qual_hotspot";
    private static final String REF_SAMPLE_COUNT = "ref_sample_count";

    public FilterConfig(final ConfigBuilder configBuilder)
    {
        DisableSoftFilter = configBuilder.hasFlag(DISABLE_SOFT_FILTER);
        DisableHardFilter = configBuilder.hasFlag(DISABLE_HARD_FILTER);
        FilteredMaxGermlineAltSupport = configBuilder.getInteger(FILTERED_MAX_GERMLINE_ALT_SUPPORT);
        HardMinTumorQual = configBuilder.getInteger(HARD_MIN_TUMOR_QUAL);
        HardMinTumorVaf = configBuilder.getDecimal(HARD_MIN_TUMOR_VAF);
        HardMinTumorRawAltSupport = configBuilder.getInteger(HARD_MIN_TUMOR_RAW_ALT_SUPPORT);
        MinAvgBaseQual = configBuilder.getInteger(MIN_AVG_BASE_QUAL);
        MinAvgBaseQualHotspot = configBuilder.getInteger(MIN_AVG_BASE_QUAL_HOTSPOT);
        SoftHotspotFilter = new SoftFilterConfig(configBuilder, "hotspot", DEFAULT_HOTSPOT_FILTER);
        SoftPanelFilter = new SoftFilterConfig(configBuilder, "panel", DEFAULT_PANEL_FILTER);
        SoftHighConfidenceFilter = new SoftFilterConfig(configBuilder, "high_confidence", DEFAULT_HIGH_CONFIDENCE_FILTER);
        SoftLowConfidenceFilter = new SoftFilterConfig(configBuilder, "low_confidence", DEFAULT_LOW_CONFIDENCE_FILTER);
        ReferenceSampleCount = configBuilder.getInteger(REF_SAMPLE_COUNT);
    }

    public FilterConfig()
    {
        DisableHardFilter = true;
        DisableSoftFilter = false;
        HardMinTumorQual = DEFAULT_HARD_MIN_TUMOR_QUAL;
        HardMinTumorVaf = DEFAULT_HARD_MIN_TUMOR_VAF;
        HardMinTumorRawAltSupport = DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT;
        SoftHotspotFilter = DEFAULT_HOTSPOT_FILTER;
        SoftPanelFilter = DEFAULT_PANEL_FILTER;
        SoftHighConfidenceFilter = DEFAULT_HIGH_CONFIDENCE_FILTER;
        SoftLowConfidenceFilter = DEFAULT_LOW_CONFIDENCE_FILTER;
        FilteredMaxGermlineAltSupport = DEFAULT_FILTERED_MAX_GERMLINE_ALT_SUPPORT;
        MinAvgBaseQual = DEFAULT_MIN_AVG_BASE_QUALITY;
        MinAvgBaseQualHotspot = DEFAULT_MIN_AVG_BASE_QUALITY_HOTSPOT;
        ReferenceSampleCount = 1;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(DISABLE_SOFT_FILTER, "Disable soft filters");
        configBuilder.addFlag(DISABLE_HARD_FILTER, "Disable hard filters");

        configBuilder.addInteger(
                FILTERED_MAX_GERMLINE_ALT_SUPPORT, "Filtered max germline alt support", DEFAULT_FILTERED_MAX_GERMLINE_ALT_SUPPORT);

        configBuilder.addInteger(HARD_MIN_TUMOR_QUAL, "Hard minimum tumor quality", DEFAULT_HARD_MIN_TUMOR_QUAL);
        configBuilder.addDecimal(HARD_MIN_TUMOR_VAF, "Hard minimum tumor VAF", DEFAULT_HARD_MIN_TUMOR_VAF);

        configBuilder.addInteger(HARD_MIN_TUMOR_RAW_ALT_SUPPORT,
                "Hard minimum tumor raw alt support", DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT);

        configBuilder.addInteger(MIN_AVG_BASE_QUAL, "Min average base qua", DEFAULT_MIN_AVG_BASE_QUALITY);

        configBuilder.addInteger(
                MIN_AVG_BASE_QUAL_HOTSPOT,
                "Min average base qual for hotspots", DEFAULT_MIN_AVG_BASE_QUALITY_HOTSPOT);

        configBuilder.addInteger(
                REF_SAMPLE_COUNT,
                "Number of reference samples for applying tumor-reference filters, use 0 in 'germline' mode", 1);

        SoftFilterConfig.registerConfig(configBuilder, DEFAULT_HOTSPOT_FILTER);
        SoftFilterConfig.registerConfig(configBuilder, DEFAULT_PANEL_FILTER);
        SoftFilterConfig.registerConfig(configBuilder, DEFAULT_HIGH_CONFIDENCE_FILTER);
        SoftFilterConfig.registerConfig(configBuilder, DEFAULT_LOW_CONFIDENCE_FILTER);
    }
}
