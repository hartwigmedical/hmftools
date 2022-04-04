package com.hartwig.hmftools.sage.filter;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FILTERED_MAX_NORMAL_ALT_SUPPORT;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_HARD_MIN_TUMOR_BASE_QUALITY;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_HARD_MIN_TUMOR_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_HARD_MIN_TUMOR_VAF;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_HIGH_CONFIDENCE_FILTER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_HOTSPOT_FILTER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_LOW_CONFIDENCE_FILTER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_PANEL_FILTER;

import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.common.VariantTier;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class FilterConfig
{
    public final boolean HardFilter;
    public final boolean SoftFilter;
    public final boolean MnvFilter;
    public final int HardMinTumorQual;
    public final double HardMinTumorVaf;
    public final int HardMinTumorRawAltSupport;
    public final int HardMinTumorRawBaseQuality;
    public final int FilteredMaxNormalAltSupport;
    public final int ReferenceSampleCount;

    public final SoftFilterConfig SoftHotspotFilter;
    public final SoftFilterConfig SoftPanelFilter;
    public final SoftFilterConfig SoftHighConfidenceFilter;
    public final SoftFilterConfig SoftLowConfidenceFilter;

    private static final String SOFT_FILTER = "soft_filter_enabled";
    private static final String HARD_FILTER = "hard_filter_enabled";
    private static final String MNV_FILTER = "mnv_filter_enabled";

    private static final String HARD_MIN_TUMOR_QUAL = "hard_min_tumor_qual";
    private static final String HARD_MIN_TUMOR_VAF = "hard_min_tumor_vaf";
    private static final String HARD_MIN_TUMOR_RAW_ALT_SUPPORT = "hard_min_tumor_raw_alt_support";
    private static final String HARD_MIN_TUMOR_RAW_BASE_QUALITY = "hard_min_tumor_raw_base_quality";
    private static final String FILTERED_MAX_NORMAL_ALT_SUPPORT = "filtered_max_normal_alt_support";
    private static final String REF_SAMPLE_COUNT = "ref_sample_count";

    private static final boolean DEFAULT_SOFT_FILTER_ENABLED = true;
    private static final boolean DEFAULT_HARD_FILTER_ENABLED = false;
    private static final boolean DEFAULT_MNV_FILTER_ENABLED = true;

    public FilterConfig(final CommandLine cmd)
    {
        SoftFilter = getConfigValue(cmd, SOFT_FILTER, DEFAULT_SOFT_FILTER_ENABLED);
        HardFilter = getConfigValue(cmd, HARD_FILTER, DEFAULT_HARD_FILTER_ENABLED);
        MnvFilter = getConfigValue(cmd, MNV_FILTER, DEFAULT_MNV_FILTER_ENABLED);
        FilteredMaxNormalAltSupport = getConfigValue(cmd, FILTERED_MAX_NORMAL_ALT_SUPPORT, DEFAULT_FILTERED_MAX_NORMAL_ALT_SUPPORT);
        HardMinTumorQual = getConfigValue(cmd, HARD_MIN_TUMOR_QUAL, DEFAULT_HARD_MIN_TUMOR_QUAL);
        HardMinTumorVaf = getConfigValue(cmd, HARD_MIN_TUMOR_VAF, DEFAULT_HARD_MIN_TUMOR_VAF);
        HardMinTumorRawAltSupport = getConfigValue(cmd, HARD_MIN_TUMOR_RAW_ALT_SUPPORT, DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT);
        HardMinTumorRawBaseQuality = getConfigValue(cmd, HARD_MIN_TUMOR_RAW_BASE_QUALITY, DEFAULT_HARD_MIN_TUMOR_BASE_QUALITY);
        SoftHotspotFilter = new SoftFilterConfig(cmd, "hotspot", DEFAULT_HOTSPOT_FILTER);
        SoftPanelFilter = new SoftFilterConfig(cmd, "panel", DEFAULT_PANEL_FILTER);
        SoftHighConfidenceFilter = new SoftFilterConfig(cmd, "high_confidence", DEFAULT_HIGH_CONFIDENCE_FILTER);
        SoftLowConfidenceFilter = new SoftFilterConfig(cmd, "low_confidence", DEFAULT_LOW_CONFIDENCE_FILTER);
        ReferenceSampleCount = getConfigValue(cmd, REF_SAMPLE_COUNT, 1);
    }

    public FilterConfig()
    {
        HardFilter = false;
        SoftFilter = true;
        MnvFilter = true;
        HardMinTumorQual = DEFAULT_HARD_MIN_TUMOR_QUAL;
        HardMinTumorVaf = DEFAULT_HARD_MIN_TUMOR_VAF;
        HardMinTumorRawAltSupport = DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT;
        HardMinTumorRawBaseQuality = DEFAULT_HARD_MIN_TUMOR_BASE_QUALITY;
        SoftHotspotFilter = DEFAULT_HOTSPOT_FILTER;
        SoftPanelFilter = DEFAULT_PANEL_FILTER;
        SoftHighConfidenceFilter = DEFAULT_HIGH_CONFIDENCE_FILTER;
        SoftLowConfidenceFilter = DEFAULT_LOW_CONFIDENCE_FILTER;
        FilteredMaxNormalAltSupport = DEFAULT_FILTERED_MAX_NORMAL_ALT_SUPPORT;
        ReferenceSampleCount = 1;
    }

    public static Options createOptions()
    {
        final Options options = new Options();

        options.addOption(SOFT_FILTER, false, "Enable soft filters [" + DEFAULT_SOFT_FILTER_ENABLED + "]");
        options.addOption(HARD_FILTER, false, "All filters are hard [" + DEFAULT_HARD_FILTER_ENABLED + "]");
        options.addOption(MNV_FILTER, false, "Enable max_germline_alt_support mnv filter [" + DEFAULT_MNV_FILTER_ENABLED + "]");
        options.addOption(HARD_MIN_TUMOR_QUAL, true, "Hard minimum tumor quality [" + DEFAULT_HARD_MIN_TUMOR_QUAL + "]");
        options.addOption(HARD_MIN_TUMOR_VAF, true, "Hard minimum tumor VAF [" + DEFAULT_HARD_MIN_TUMOR_VAF + "]");
        options.addOption(HARD_MIN_TUMOR_RAW_ALT_SUPPORT, true,
                "Hard minimum tumor raw alt support [" + DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT + "]");
        options.addOption(HARD_MIN_TUMOR_RAW_BASE_QUALITY, true,
                "Hard minimum tumor raw base quality [" + DEFAULT_HARD_MIN_TUMOR_BASE_QUALITY + "]");
        options.addOption(REF_SAMPLE_COUNT, true, "Number of reference samples for applying turmor-reference filters (default=1)");

        SoftFilterConfig.createOptions("hotspot", DEFAULT_HOTSPOT_FILTER).getOptions().forEach(options::addOption);
        SoftFilterConfig.createOptions("panel", DEFAULT_PANEL_FILTER).getOptions().forEach(options::addOption);
        SoftFilterConfig.createOptions("high_confidence", DEFAULT_HIGH_CONFIDENCE_FILTER).getOptions().forEach(options::addOption);
        SoftFilterConfig.createOptions("low_confidence", DEFAULT_LOW_CONFIDENCE_FILTER).getOptions().forEach(options::addOption);

        return options;
    }
}
