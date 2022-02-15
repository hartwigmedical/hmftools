package com.hartwig.hmftools.sage.config;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.sage.config.SoftFilter.MAX_GERMLINE_REL_RAW_BASE_QUAL;
import static com.hartwig.hmftools.sage.config.SoftFilter.MAX_GERMLINE_VAF;
import static com.hartwig.hmftools.sage.config.SoftFilter.MIN_GERMLINE_DEPTH;
import static com.hartwig.hmftools.sage.config.SoftFilter.MIN_GERMLINE_DEPTH_ALLOSOME;
import static com.hartwig.hmftools.sage.config.SoftFilter.MIN_TUMOR_QUAL;
import static com.hartwig.hmftools.sage.config.SoftFilter.MIN_TUMOR_VAF;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class SoftFilterConfig
{
    public final int MinTumorQual;
    public final double MinTumorVaf;
    public final int MinGermlineReadContextCoverage;
    public final int MinGermlineReadContextCoverageAllosome;
    public final double MaxGermlineVaf;
    public final double MaxGermlineRelativeQual;

    public SoftFilterConfig(final CommandLine cmd, final String prefix, final SoftFilterConfig defaultValue)
    {
        MinTumorQual = getConfigValue(cmd, prefix + "_" + MIN_TUMOR_QUAL, defaultValue.MinTumorQual);
        MinTumorVaf = getConfigValue(cmd, prefix + "_" + MIN_TUMOR_VAF, defaultValue.MinTumorVaf);

        MinGermlineReadContextCoverage = getConfigValue(
                cmd, prefix + "_" + MIN_GERMLINE_DEPTH, defaultValue.MinGermlineReadContextCoverage);

        MinGermlineReadContextCoverageAllosome =
                getConfigValue(cmd, prefix + "_" + MIN_GERMLINE_DEPTH_ALLOSOME, defaultValue.MinGermlineReadContextCoverageAllosome);

        MaxGermlineVaf = getConfigValue(cmd, prefix + "_" + MAX_GERMLINE_VAF, defaultValue.MaxGermlineVaf);

        MaxGermlineRelativeQual =
                getConfigValue(cmd, prefix + "_" + MAX_GERMLINE_REL_RAW_BASE_QUAL, defaultValue.MaxGermlineRelativeQual);
    }

    public SoftFilterConfig(
            final int minTumorQual, final double minTumorVaf, final int minGermlineReadContextCoverage,
            final int minGermlineReadContextCoverageAllosome, final double maxGermlineVaf, final double maxGermlineRelativeQual)
    {
        MinTumorQual = minTumorQual;
        MinTumorVaf = minTumorVaf;
        MinGermlineReadContextCoverage = minGermlineReadContextCoverage;
        MinGermlineReadContextCoverageAllosome = minGermlineReadContextCoverageAllosome;
        MaxGermlineVaf = maxGermlineVaf;
        MaxGermlineRelativeQual = maxGermlineRelativeQual;
    }

    @NotNull
    public static Options createOptions(@NotNull final String prefix, @NotNull final SoftFilterConfig defaultValue)
    {
        final Options options = new Options();

        options.addOption(
                prefix + "_" + MIN_TUMOR_QUAL.config(), true, "Minimum " + prefix + " tumor quality [" + defaultValue.MinTumorQual + "]");
        options.addOption(prefix + "_" + MIN_TUMOR_VAF.config(), true, "Minimum " + prefix + " tumor VAF [" + defaultValue.MinTumorVaf + "]");
        options.addOption(
                prefix + "_" + MIN_GERMLINE_DEPTH.config(), true,
                "Minimum " + prefix + " germline depth [" + defaultValue.MinGermlineReadContextCoverage + "]");
        options.addOption(
                prefix + "_" + MIN_GERMLINE_DEPTH_ALLOSOME.config(), true,
                "Minimum " + prefix + " germline depth [" + defaultValue.MinGermlineReadContextCoverageAllosome + "]");

        options.addOption(
                prefix + "_" + MAX_GERMLINE_VAF.config(), true, "Maximum " + prefix + " germline VAF [" + defaultValue.MaxGermlineVaf + "]");
        options.addOption(
                prefix + "_" + MAX_GERMLINE_REL_RAW_BASE_QUAL.config(), true,
                "Maximum " + prefix + " germline relative quality [" + defaultValue.MaxGermlineRelativeQual + "]");

        return options;
    }
}
