package com.hartwig.hmftools.sage.filter;

import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_REL_RAW_BASE_QUAL;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_VAF;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MIN_GERMLINE_DEPTH;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MIN_TUMOR_QUAL;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MIN_TUMOR_VAF;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class SoftFilterConfig
{
    public final String Name;
    public final int MinTumorQual;
    public final double MinTumorVaf;
    public final int MinGermlineCoverage;
    public final int MinGermlineCoverageLongInsert;
    public final int MinGermlineCoverageAllosome;
    public final int MinGermlineCoverageAllosomeLongInsert;
    public final double MaxGermlineVaf;
    public final double MaxGermlineRelativeQual;

    public SoftFilterConfig(final ConfigBuilder configBuilder, final String prefix, final SoftFilterConfig defaultValue)
    {
        Name = defaultValue.Name;
        MinTumorQual = configBuilder.getInteger(prefix + "_" + MIN_TUMOR_QUAL.configName());
        MinTumorVaf = configBuilder.getDecimal(prefix + "_" + MIN_TUMOR_VAF.configName());

        MinGermlineCoverage = configBuilder.getInteger(prefix + "_" + MIN_GERMLINE_DEPTH.configName());
        MinGermlineCoverageLongInsert = configBuilder.getInteger(prefix + "_" + MIN_GERMLINE_DEPTH.configName());
        MinGermlineCoverageAllosome = defaultValue.MinGermlineCoverageAllosome;
        MinGermlineCoverageAllosomeLongInsert = defaultValue.MinGermlineCoverageAllosomeLongInsert;
        MaxGermlineVaf = configBuilder.getDecimal(prefix + "_" + MAX_GERMLINE_VAF.configName());
        MaxGermlineRelativeQual = configBuilder.getDecimal(prefix + "_" + MAX_GERMLINE_REL_RAW_BASE_QUAL.configName());
    }

    public SoftFilterConfig(
            final String name, final int minTumorQual, final double minTumorVaf, final int minGermlineCoverage,
            final int minGermlineCoverageLongInsert, final int minGermlineCoverageAllosome, final int minGermlineCoverageAllosomeLongInsert,
            final double maxGermlineVaf, final double maxGermlineRelativeQual)
    {
        Name = name;
        MinTumorQual = minTumorQual;
        MinTumorVaf = minTumorVaf;
        MinGermlineCoverage = minGermlineCoverage;
        MinGermlineCoverageLongInsert = minGermlineCoverageLongInsert;
        MinGermlineCoverageAllosome = minGermlineCoverageAllosome;
        MinGermlineCoverageAllosomeLongInsert = minGermlineCoverageAllosomeLongInsert;
        MaxGermlineVaf = maxGermlineVaf;
        MaxGermlineRelativeQual = maxGermlineRelativeQual;
    }

    public static void registerConfig(final ConfigBuilder configBuilder, final SoftFilterConfig defaultConfig)
    {
        String prefix = defaultConfig.Name;

        configBuilder.addInteger(
                prefix + "_" + MIN_TUMOR_QUAL.configName(),
                "Minimum " + prefix + " tumor quality", defaultConfig.MinTumorQual);

        configBuilder.addDecimal(
                prefix + "_" + MIN_TUMOR_VAF.configName(), "Minimum " + prefix + " tumor VAF",defaultConfig.MinTumorVaf);

        configBuilder.addInteger(
                prefix + "_" + MIN_GERMLINE_DEPTH.configName(),
                "Minimum " + prefix + " germline depth", defaultConfig.MinGermlineCoverage);

        configBuilder.addDecimal(
                prefix + "_" + MAX_GERMLINE_VAF.configName(), "Maximum " + prefix + " germline VAF", defaultConfig.MaxGermlineVaf);

        configBuilder.addDecimal(
                prefix + "_" + MAX_GERMLINE_REL_RAW_BASE_QUAL.configName(),
                "Maximum " + prefix + " germline relative qualit", defaultConfig.MaxGermlineRelativeQual);
    }
}
