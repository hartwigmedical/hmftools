package com.hartwig.hmftools.sage.config;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;

import java.util.function.Predicate;

import com.hartwig.hmftools.sage.candidate.AltContext;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.select.HotspotSelector;
import com.hartwig.hmftools.sage.variant.VariantTier;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class FilterConfig
{
    public final boolean HardFilter;
    public final boolean SoftFilter;
    public final boolean MnvFilter;
    public final int HardMinTumorQual;
    public final int HardMinTumorRawAltSupport;
    public final int HardMinTumorRawBaseQuality;
    public final int FilteredMaxNormalAltSupport;

    public final SoftFilterConfig SoftHotspotFilter;
    public final SoftFilterConfig SoftPanelFilter;
    public final SoftFilterConfig SoftHighConfidenceFilter;
    public final SoftFilterConfig SoftLowConfidenceFilter;

    public double hotspotMinTumorVafToSkipQualCheck()
    {
        return 0.08;
    }

    public int hotspotMinTumorAltSupportToSkipQualCheck() { return 8; }

    private static final String SOFT_FILTER = "soft_filter_enabled";
    private static final String HARD_FILTER = "hard_filter_enabled";
    private static final String MNV_FILTER = "mnv_filter_enabled";

    private static final String HARD_MIN_TUMOR_QUAL = "hard_min_tumor_qual";
    private static final String HARD_MIN_TUMOR_RAW_ALT_SUPPORT = "hard_min_tumor_raw_alt_support";
    private static final String HARD_MIN_TUMOR_RAW_BASE_QUALITY = "hard_min_tumor_raw_base_quality";
    private static final String FILTERED_MAX_NORMAL_ALT_SUPPORT = "filtered_max_normal_alt_support";

    private static final boolean DEFAULT_SOFT_FILTER_ENABLED = true;
    private static final boolean DEFAULT_HARD_FILTER_ENABLED = false;
    private static final boolean DEFAULT_MNV_FILTER_ENABLED = true;

    private static final int DEFAULT_HARD_MIN_TUMOR_BASE_QUALITY = 0;
    private static final int DEFAULT_HARD_MIN_TUMOR_QUAL = 30;
    private static final int DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT = 2;
    private static final int DEFAULT_FILTERED_MAX_NORMAL_ALT_SUPPORT = 3;

    public FilterConfig(final CommandLine cmd)
    {
        SoftFilter = getConfigValue(cmd, SOFT_FILTER, DEFAULT_SOFT_FILTER_ENABLED);
        HardFilter = getConfigValue(cmd, HARD_FILTER, DEFAULT_HARD_FILTER_ENABLED);
        MnvFilter = getConfigValue(cmd, MNV_FILTER, DEFAULT_MNV_FILTER_ENABLED);
        FilteredMaxNormalAltSupport = getConfigValue(cmd, FILTERED_MAX_NORMAL_ALT_SUPPORT, DEFAULT_FILTERED_MAX_NORMAL_ALT_SUPPORT);
        HardMinTumorQual = getConfigValue(cmd, HARD_MIN_TUMOR_QUAL, DEFAULT_HARD_MIN_TUMOR_QUAL);
        HardMinTumorRawAltSupport = getConfigValue(cmd, HARD_MIN_TUMOR_RAW_ALT_SUPPORT, DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT);
        HardMinTumorRawBaseQuality = getConfigValue(cmd, HARD_MIN_TUMOR_RAW_BASE_QUALITY, DEFAULT_HARD_MIN_TUMOR_BASE_QUALITY);
        SoftHotspotFilter = new SoftFilterConfig(cmd, "hotspot", DEFAULT_HOTSPOT_FILTER);
        SoftPanelFilter = new SoftFilterConfig(cmd, "panel", DEFAULT_PANEL_FILTER);
        SoftHighConfidenceFilter = new SoftFilterConfig(cmd, "high_confidence", DEFAULT_HIGH_CONFIDENCE_FILTER);
        SoftLowConfidenceFilter = new SoftFilterConfig(cmd, "low_confidence", DEFAULT_LOW_CONFIDENCE_FILTER);
    }

    public FilterConfig()
    {
        HardFilter = false;
        SoftFilter = true;
        MnvFilter = true;
        HardMinTumorQual = DEFAULT_HARD_MIN_TUMOR_QUAL;
        HardMinTumorRawAltSupport = DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT;
        HardMinTumorRawBaseQuality = DEFAULT_HARD_MIN_TUMOR_BASE_QUALITY;
        SoftHotspotFilter = DEFAULT_HOTSPOT_FILTER;
        SoftPanelFilter = DEFAULT_PANEL_FILTER;
        SoftHighConfidenceFilter = DEFAULT_HIGH_CONFIDENCE_FILTER;
        SoftLowConfidenceFilter = DEFAULT_LOW_CONFIDENCE_FILTER;
        FilteredMaxNormalAltSupport = DEFAULT_FILTERED_MAX_NORMAL_ALT_SUPPORT;
    }

    private static final SoftFilterConfig NO_FILTER = new SoftFilterConfig(
            0, 0, 0, 0,
            1d, 1d);

    private static final SoftFilterConfig DEFAULT_HOTSPOT_FILTER = new SoftFilterConfig(70, 0.005,
            0, 0, 0.1, 0.5);

    private static final SoftFilterConfig DEFAULT_PANEL_FILTER = new SoftFilterConfig(100, 0.015,
            0, 0, 0.04, 0.04);

    private static final SoftFilterConfig DEFAULT_HIGH_CONFIDENCE_FILTER = new SoftFilterConfig(160, 0.025,
            10, 6, 0.04, 0.04);

    private static final SoftFilterConfig DEFAULT_LOW_CONFIDENCE_FILTER = new SoftFilterConfig(240, 0.025,
            10, 6, 0.04, 0.04);

    public static Options createOptions()
    {
        final Options options = new Options();

        options.addOption(SOFT_FILTER, false, "Enable soft filters [" + DEFAULT_SOFT_FILTER_ENABLED + "]");
        options.addOption(HARD_FILTER, false, "All filters are hard [" + DEFAULT_HARD_FILTER_ENABLED + "]");
        options.addOption(MNV_FILTER, false, "Enable max_germline_alt_support mnv filter [" + DEFAULT_MNV_FILTER_ENABLED + "]");
        options.addOption(HARD_MIN_TUMOR_QUAL, true, "Hard minimum tumor quality [" + DEFAULT_HARD_MIN_TUMOR_QUAL + "]");
        options.addOption(HARD_MIN_TUMOR_RAW_ALT_SUPPORT, true,
                "Hard minimum tumor raw alt support [" + DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT + "]");
        options.addOption(HARD_MIN_TUMOR_RAW_BASE_QUALITY, true,
                "Hard minimum tumor raw base quality [" + DEFAULT_HARD_MIN_TUMOR_BASE_QUALITY + "]");

        SoftFilterConfig.createOptions("hotspot", DEFAULT_HOTSPOT_FILTER).getOptions().forEach(options::addOption);
        SoftFilterConfig.createOptions("panel", DEFAULT_PANEL_FILTER).getOptions().forEach(options::addOption);
        SoftFilterConfig.createOptions("high_confidence", DEFAULT_HIGH_CONFIDENCE_FILTER).getOptions().forEach(options::addOption);
        SoftFilterConfig.createOptions("low_confidence", DEFAULT_LOW_CONFIDENCE_FILTER).getOptions().forEach(options::addOption);

        return options;
    }


    public SoftFilterConfig softConfig(@NotNull final VariantTier tier)
    {
        switch(tier)
        {
            case HOTSPOT:
                return SoftHotspotFilter;
            case PANEL:
                return SoftPanelFilter;
            case HIGH_CONFIDENCE:
                return SoftHighConfidenceFilter;
            default:
                return SoftLowConfidenceFilter;
        }
    }

    public Predicate<AltContext> altContextFilter(final HotspotSelector hotspotSelector)
    {
        return altContext ->
        {
            if(hotspotSelector.isHotspot(altContext))
            {
                return true;
            }
            return altContext.rawAltBaseQuality() >= HardMinTumorRawBaseQuality
                    && altContext.rawAltSupport() >= HardMinTumorRawAltSupport;
        };
    }

    public Predicate<ReadContextCounter> readContextFilter()
    {
        return readContextCounter ->
        {
            if(readContextCounter.Tier.equals(VariantTier.HOTSPOT))
            {
                return true;
            }
            return readContextCounter.rawAltBaseQuality() >= HardMinTumorRawBaseQuality
                    && readContextCounter.rawAltSupport() >= HardMinTumorRawAltSupport
                    && readContextCounter.tumorQuality() >= HardMinTumorQual;
        };
    }
}
