package com.hartwig.hmftools.sage.filter;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_REL_RAW_BASE_QUAL;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_VAF;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MIN_GERMLINE_DEPTH;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MIN_GERMLINE_DEPTH_ALLOSOME;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MIN_TUMOR_QUAL;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MIN_TUMOR_VAF;

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
        MinTumorQual = getConfigValue(cmd, prefix + "_" + MIN_TUMOR_QUAL.configName(), defaultValue.MinTumorQual);
        MinTumorVaf = getConfigValue(cmd, prefix + "_" + MIN_TUMOR_VAF.configName(), defaultValue.MinTumorVaf);

        MinGermlineReadContextCoverage = getConfigValue(
                cmd, prefix + "_" + MIN_GERMLINE_DEPTH.configName(), defaultValue.MinGermlineReadContextCoverage);

        MinGermlineReadContextCoverageAllosome =
                getConfigValue(cmd, prefix + "_" + MIN_GERMLINE_DEPTH_ALLOSOME.configName(), defaultValue.MinGermlineReadContextCoverageAllosome);

        MaxGermlineVaf = getConfigValue(cmd, prefix + "_" + MAX_GERMLINE_VAF.configName(), defaultValue.MaxGermlineVaf);

        MaxGermlineRelativeQual =
                getConfigValue(cmd, prefix + "_" + MAX_GERMLINE_REL_RAW_BASE_QUAL.configName(), defaultValue.MaxGermlineRelativeQual);
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
                prefix + "_" + MIN_TUMOR_QUAL.configName(), true, "Minimum " + prefix + " tumor quality [" + defaultValue.MinTumorQual + "]");
        options.addOption(prefix + "_" + MIN_TUMOR_VAF.configName(), true, "Minimum " + prefix + " tumor VAF [" + defaultValue.MinTumorVaf + "]");
        options.addOption(
                prefix + "_" + MIN_GERMLINE_DEPTH.configName(), true,
                "Minimum " + prefix + " germline depth [" + defaultValue.MinGermlineReadContextCoverage + "]");
        options.addOption(
                prefix + "_" + MIN_GERMLINE_DEPTH_ALLOSOME.configName(), true,
                "Minimum " + prefix + " germline depth [" + defaultValue.MinGermlineReadContextCoverageAllosome + "]");

        options.addOption(
                prefix + "_" + MAX_GERMLINE_VAF.configName(), true, "Maximum " + prefix + " germline VAF [" + defaultValue.MaxGermlineVaf + "]");
        options.addOption(
                prefix + "_" + MAX_GERMLINE_REL_RAW_BASE_QUAL.configName(), true,
                "Maximum " + prefix + " germline relative quality [" + defaultValue.MaxGermlineRelativeQual + "]");

        return options;
    }
}
