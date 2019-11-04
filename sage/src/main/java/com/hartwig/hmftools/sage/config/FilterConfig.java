package com.hartwig.hmftools.sage.config;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface FilterConfig {

    FilterTierConfig NO_FILTER = ImmutableFilterTierConfig.builder()
            .minTumorQual(0)
            .minTumorVaf(0)
            .minGermlineDepth(0)
            .maxGermlineVaf(1d)
            .maxGermlineRelativeQual(1d)
            .maxGermlineRelativeReadContextCount(1d)
            .build();

    FilterTierConfig DEFAULT_HARD_FILTER = ImmutableFilterTierConfig.builder()
            .from(NO_FILTER)
            .minTumorQual(1)
            .maxGermlineVaf(0.3)
            .build();

    FilterTierConfig DEFAULT_HOTSPOT_FILTER = ImmutableFilterTierConfig.builder()
            .from(NO_FILTER)
            .minTumorQual(50)
            .minTumorVaf(0.05)
            .build();

    FilterTierConfig DEFAULT_PANEL_FILTER = ImmutableFilterTierConfig.builder()
            .from(NO_FILTER)
            .minTumorQual(100)
            .minTumorVaf(0.01)
            .maxGermlineVaf(0.05)
            .maxGermlineRelativeQual(0.05)
            .maxGermlineRelativeReadContextCount(0.05)
            .build();

    FilterTierConfig DEFAULT_WIDE_FILTER = ImmutableFilterTierConfig.builder()
            .from(NO_FILTER)
            .minTumorQual(150)
            .minTumorVaf(0.025)
            .minGermlineDepth(7)
            .maxGermlineVaf(0.05)
            .maxGermlineRelativeQual(0.05)
            .maxGermlineRelativeReadContextCount(0.05)
            .build();

    @NotNull
    FilterTierConfig hardFilter();

    @NotNull
    FilterTierConfig hotspotFilter();

    @NotNull
    FilterTierConfig panelFilter();

    @NotNull
    FilterTierConfig wideFilter();

    @NotNull
    static Options createOptions() {
        final Options options = new Options();

        FilterTierConfig.createOptions("hard", DEFAULT_HARD_FILTER).getOptions().forEach(options::addOption);
        FilterTierConfig.createOptions("hotspot", DEFAULT_HOTSPOT_FILTER).getOptions().forEach(options::addOption);
        FilterTierConfig.createOptions("panel", DEFAULT_PANEL_FILTER).getOptions().forEach(options::addOption);
        FilterTierConfig.createOptions("wide", DEFAULT_WIDE_FILTER).getOptions().forEach(options::addOption);

        return options;
    }

    @NotNull
    static FilterConfig createConfig(@NotNull final CommandLine cmd) throws ParseException {
        return ImmutableFilterConfig.builder()
                .hardFilter(FilterTierConfig.createConfig(cmd, "hard", DEFAULT_HARD_FILTER))
                .hotspotFilter(FilterTierConfig.createConfig(cmd, "hotspot", DEFAULT_HOTSPOT_FILTER))
                .panelFilter(FilterTierConfig.createConfig(cmd, "panel", DEFAULT_PANEL_FILTER))
                .wideFilter(FilterTierConfig.createConfig(cmd, "wide", DEFAULT_WIDE_FILTER))
                .build();

    }
}
