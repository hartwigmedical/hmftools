package com.hartwig.hmftools.sage.config;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface SoftFilterConfig {

    String MIN_TUMOR_QUAL = "min_tumor_qual";
    String MIN_TUMOR_VAF = "min_tumor_vaf";
    String MIN_GERMLINE_DEPTH = "min_germline_depth";
    String MAX_GERMLINE_VAF = "max_germline_vaf";
    String MAX_GERMLINE_REL_QUAL = "max_germline_rel_qual";
    String MAX_GERMLINE_REL_RCC = "max_germline_rel_rcc";

    int minTumorQual();

    double minTumorVaf();

    int minGermlineDepth();

    double maxGermlineVaf();

    double maxGermlineRelativeQual();

    double maxGermlineRelativeReadContextCount();

    @NotNull
    static Options createOptions(@NotNull final String prefix, @NotNull final SoftFilterConfig defaultValue) {
        final Options options = new Options();
        options.addOption(prefix + "_" + MIN_TUMOR_QUAL,true,"Minimum " + prefix + " tumor quality [" + defaultValue.minTumorQual() + "]");
        options.addOption(prefix + "_" + MIN_TUMOR_VAF, true, "Minimum " + prefix + " tumor VAF [" + defaultValue.minTumorVaf() + "]");
        options.addOption(prefix + "_" + MIN_GERMLINE_DEPTH, true, "Minimum " + prefix + " germline depth [" + defaultValue.minGermlineDepth() + "]");

        options.addOption(prefix + "_" + MAX_GERMLINE_VAF,true,"Maximum " + prefix + " germline VAF [" + defaultValue.maxGermlineVaf() + "]");
        options.addOption(prefix + "_" + MAX_GERMLINE_REL_QUAL, true, "Maximum " + prefix + " germline relative quality [" + defaultValue.maxGermlineRelativeQual() + "]");
        options.addOption(prefix + "_" + MAX_GERMLINE_REL_RCC, true, "Maximum " + prefix + " germline relative read context count [" + defaultValue.maxGermlineRelativeReadContextCount() + "]");

        return options;
    }

    @NotNull
    static SoftFilterConfig createConfig(@NotNull final CommandLine cmd, @NotNull final String prefix, @NotNull final SoftFilterConfig defaultValue) throws ParseException {
        final int minTumorQual = SageConfig.defaultIntValue(cmd, prefix + "_" + MIN_TUMOR_QUAL, defaultValue.minTumorQual());
        final double minTumorVaf = SageConfig.defaultValue(cmd, prefix + "_" + MIN_TUMOR_VAF, defaultValue.minTumorVaf());
        final int minGermlineDepth = SageConfig.defaultIntValue(cmd, prefix + "_" + MIN_GERMLINE_DEPTH, defaultValue.minGermlineDepth());

        final double maxGermlineVaf = SageConfig.defaultValue(cmd, prefix + "_" + MAX_GERMLINE_VAF, defaultValue.maxGermlineVaf());
        final double maxGermlineRelativeQual = SageConfig.defaultValue(cmd, prefix + "_" + MAX_GERMLINE_REL_QUAL, defaultValue.maxGermlineRelativeQual());
        final double maxGermlineRelativeReadContextCount = SageConfig.defaultValue(cmd, prefix + "_" + MAX_GERMLINE_REL_RCC, defaultValue.maxGermlineRelativeReadContextCount());

        return ImmutableSoftFilterConfig.builder()
                .minTumorQual(minTumorQual)
                .minTumorVaf(minTumorVaf)
                .minGermlineDepth(minGermlineDepth)
                .maxGermlineVaf(maxGermlineVaf)
                .maxGermlineRelativeQual(maxGermlineRelativeQual)
                .maxGermlineRelativeReadContextCount(maxGermlineRelativeReadContextCount)
                .build();
    }

}
