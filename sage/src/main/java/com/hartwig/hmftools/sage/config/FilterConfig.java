package com.hartwig.hmftools.sage.config;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;

import java.util.function.Predicate;

import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.select.HotspotSelector;
import com.hartwig.hmftools.sage.variant.SageVariantTier;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface FilterConfig
{

    String SOFT_FILTER = "soft_filter_enabled";
    String HARD_FILTER = "hard_filter_enabled";
    String MNV_FILTER = "mnv_filter_enabled";

    String HARD_MIN_TUMOR_QUAL = "hard_min_tumor_qual";
    String HARD_MIN_TUMOR_RAW_ALT_SUPPORT = "hard_min_tumor_raw_alt_support";
    String HARD_MIN_TUMOR_RAW_BASE_QUALITY = "hard_min_tumor_raw_base_quality";
    String FILTERED_MAX_NORMAL_ALT_SUPPORT = "filtered_max_normal_alt_support";

    boolean DEFAULT_SOFT_FILTER_ENABLED = true;
    boolean DEFAULT_HARD_FILTER_ENABLED = false;
    boolean DEFAULT_MNV_FILTER_ENABLED = true;

    int DEFAULT_HARD_MIN_TUMOR_BASE_QUALITY = 0;
    int DEFAULT_HARD_MIN_TUMOR_QUAL = 30;
    int DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT = 2;
    int DEFAULT_FILTERED_MAX_NORMAL_ALT_SUPPORT = 3;

    SoftFilterConfig NO_FILTER = ImmutableSoftFilterConfig.builder()
            .minTumorQual(0)
            .minTumorVaf(0)
            .minGermlineReadContextCoverage(0)
            .minGermlineReadContextCoverageAllosome(0)
            .maxGermlineVaf(1d)
            .maxGermlineRelativeQual(1d)
            .build();

    SoftFilterConfig DEFAULT_HOTSPOT_FILTER = ImmutableSoftFilterConfig.builder()
            .from(NO_FILTER)
            .minTumorQual(70)
            .minTumorVaf(0.005)
            .maxGermlineVaf(0.1)
            .maxGermlineRelativeQual(0.5)
            .build();

    SoftFilterConfig DEFAULT_PANEL_FILTER = ImmutableSoftFilterConfig.builder()
            .from(NO_FILTER)
            .minTumorQual(100)
            .minTumorVaf(0.015)
            .maxGermlineVaf(0.04)
            .maxGermlineRelativeQual(0.04)
            .build();

    SoftFilterConfig DEFAULT_HIGH_CONFIDENCE_FILTER = ImmutableSoftFilterConfig.builder()
            .from(NO_FILTER)
            .minTumorQual(160)
            .minTumorVaf(0.025)
            .minGermlineReadContextCoverage(10)
            .minGermlineReadContextCoverageAllosome(6)
            .maxGermlineVaf(0.04)
            .maxGermlineRelativeQual(0.04)
            .build();

    SoftFilterConfig DEFAULT_LOW_CONFIDENCE_FILTER = ImmutableSoftFilterConfig.builder()
            .from(NO_FILTER)
            .minTumorQual(240)
            .minTumorVaf(0.025)
            .minGermlineReadContextCoverage(10)
            .minGermlineReadContextCoverageAllosome(6)
            .maxGermlineVaf(0.04)
            .maxGermlineRelativeQual(0.04)
            .build();

    boolean hardFilter();

    boolean softFilter();

    boolean mnvFilter();

    int hardMinTumorQual();

    int hardMinTumorRawAltSupport();

    int hardMinTumorRawBaseQuality();

    int filteredMaxNormalAltSupport();

    default double hotspotMinTumorVafToSkipQualCheck()
    {
        return 0.08;
    }

    default int hotspotMinTumorAltSupportToSkipQualCheck()
    {
        return 8;
    }

    @NotNull
    SoftFilterConfig softHotspotFilter();

    @NotNull
    SoftFilterConfig softPanelFilter();

    @NotNull
    SoftFilterConfig softHighConfidenceFilter();

    @NotNull
    SoftFilterConfig softLowConfidenceFilter();

    @NotNull
    static Options createOptions()
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

    @NotNull
    static FilterConfig createConfig(@NotNull final CommandLine cmd)
    {
        return ImmutableFilterConfig.builder()
                .softFilter(getConfigValue(cmd, SOFT_FILTER, DEFAULT_SOFT_FILTER_ENABLED))
                .hardFilter(getConfigValue(cmd, HARD_FILTER, DEFAULT_HARD_FILTER_ENABLED))
                .mnvFilter(getConfigValue(cmd, MNV_FILTER, DEFAULT_MNV_FILTER_ENABLED))
                .filteredMaxNormalAltSupport(getConfigValue(cmd, FILTERED_MAX_NORMAL_ALT_SUPPORT, DEFAULT_FILTERED_MAX_NORMAL_ALT_SUPPORT))
                .hardMinTumorQual(getConfigValue(cmd, HARD_MIN_TUMOR_QUAL, DEFAULT_HARD_MIN_TUMOR_QUAL))
                .hardMinTumorRawAltSupport(getConfigValue(cmd, HARD_MIN_TUMOR_RAW_ALT_SUPPORT, DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT))
                .hardMinTumorRawBaseQuality(getConfigValue(cmd, HARD_MIN_TUMOR_RAW_BASE_QUALITY, DEFAULT_HARD_MIN_TUMOR_BASE_QUALITY))
                .softHotspotFilter(SoftFilterConfig.createConfig(cmd, "hotspot", DEFAULT_HOTSPOT_FILTER))
                .softPanelFilter(SoftFilterConfig.createConfig(cmd, "panel", DEFAULT_PANEL_FILTER))
                .softHighConfidenceFilter(SoftFilterConfig.createConfig(cmd, "high_confidence", DEFAULT_HIGH_CONFIDENCE_FILTER))
                .softLowConfidenceFilter(SoftFilterConfig.createConfig(cmd, "low_confidence", DEFAULT_LOW_CONFIDENCE_FILTER))
                .build();
    }

    @NotNull
    default SoftFilterConfig softConfig(@NotNull final SageVariantTier tier)
    {
        switch(tier)
        {
            case HOTSPOT:
                return softHotspotFilter();
            case PANEL:
                return softPanelFilter();
            case HIGH_CONFIDENCE:
                return softHighConfidenceFilter();
            default:
                return softLowConfidenceFilter();
        }
    }

    @NotNull
    default Predicate<AltContext> altContextFilter(@NotNull final HotspotSelector hotspotSelector)
    {
        return altContext ->
        {
            if(hotspotSelector.isHotspot(altContext))
            {
                return true;
            }
            return altContext.rawAltBaseQuality() >= hardMinTumorRawBaseQuality()
                    && altContext.rawAltSupport() >= hardMinTumorRawAltSupport();
        };
    }

    @NotNull
    default Predicate<ReadContextCounter> readContextFilter()
    {
        return readContextCounter ->
        {
            if(readContextCounter.tier().equals(SageVariantTier.HOTSPOT))
            {
                return true;
            }
            return readContextCounter.rawAltBaseQuality() >= hardMinTumorRawBaseQuality()
                    && readContextCounter.rawAltSupport() >= hardMinTumorRawAltSupport()
                    && readContextCounter.tumorQuality() >= hardMinTumorQual();
        };
    }
}
