package com.hartwig.hmftools.sage.config;

import com.hartwig.hmftools.sage.SageConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface FilterConfig {

    String HARD_MIN_TUMOR_QUAL = "hard_min_tumor_qual";
    String HARD_MIN_TUMOR_ALT_SUPPORT = "hard_min_tumor_alt_support";

    int DEFAULT_HARD_MIN_TUMOR_QUAL = 1;
    int DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT = 3;

    SoftFilterConfig NO_FILTER = ImmutableSoftFilterConfig.builder()
            .minTumorQual(0)
            .minTumorVaf(0)
            .minGermlineDepth(0)
            .maxGermlineVaf(1d)
            .maxGermlineRelativeQual(1d)
            .maxGermlineRelativeReadContextCount(1d)
            .build();

    SoftFilterConfig DEFAULT_HOTSPOT_FILTER = ImmutableSoftFilterConfig.builder()
            .from(NO_FILTER)
            .minTumorQual(50)
            .minTumorVaf(0.05)
            .build();

    SoftFilterConfig DEFAULT_PANEL_FILTER = ImmutableSoftFilterConfig.builder()
            .from(NO_FILTER)
            .minTumorQual(100)
            .minTumorVaf(0.01)
            .maxGermlineVaf(0.05)
            .maxGermlineRelativeQual(0.05)
            .maxGermlineRelativeReadContextCount(0.05)
            .build();

    SoftFilterConfig DEFAULT_WIDE_FILTER = ImmutableSoftFilterConfig.builder()
            .from(NO_FILTER)
            .minTumorQual(150)
            .minTumorVaf(0.025)
            .minGermlineDepth(7)
            .maxGermlineVaf(0.05)
            .maxGermlineRelativeQual(0.05)
            .maxGermlineRelativeReadContextCount(0.05)
            .build();

    int hardMinTumorQual();

    int hardMinTumorAltSupport();

    @NotNull
    SoftFilterConfig softHotspotFilter();

    @NotNull
    SoftFilterConfig softPanelFilter();

    @NotNull
    SoftFilterConfig softWideFilter();

    @NotNull
    static Options createOptions() {
        final Options options = new Options();

        options.addOption(HARD_MIN_TUMOR_QUAL, true, "Hard minimum tumor quality [" + DEFAULT_HARD_MIN_TUMOR_QUAL + "]");
        options.addOption(HARD_MIN_TUMOR_ALT_SUPPORT, true, "Hard minimum tumor alt support [" + DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT + "]");

        SoftFilterConfig.createOptions("hotspot", DEFAULT_HOTSPOT_FILTER).getOptions().forEach(options::addOption);
        SoftFilterConfig.createOptions("panel", DEFAULT_PANEL_FILTER).getOptions().forEach(options::addOption);
        SoftFilterConfig.createOptions("wide", DEFAULT_WIDE_FILTER).getOptions().forEach(options::addOption);

        return options;
    }

    @NotNull
    static FilterConfig createConfig(@NotNull final CommandLine cmd) throws ParseException {
        return ImmutableFilterConfig.builder()
                .hardMinTumorQual(SageConfig.defaultIntValue(cmd, HARD_MIN_TUMOR_QUAL, DEFAULT_HARD_MIN_TUMOR_QUAL))
                .hardMinTumorAltSupport(SageConfig.defaultIntValue(cmd, HARD_MIN_TUMOR_ALT_SUPPORT, DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT))
                .softHotspotFilter(SoftFilterConfig.createConfig(cmd, "hotspot", DEFAULT_HOTSPOT_FILTER))
                .softPanelFilter(SoftFilterConfig.createConfig(cmd, "panel", DEFAULT_PANEL_FILTER))
                .softWideFilter(SoftFilterConfig.createConfig(cmd, "wide", DEFAULT_WIDE_FILTER))
                .build();

    }
}
