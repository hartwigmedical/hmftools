package com.hartwig.hmftools.qsee.coverage;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.QseeConfig;

public class ExonExclusionsConfig
{
    public final QseeConfig Qsee;

    public final double MinLowCoverageSampleFraction;

    private static final String MIN_LOW_COVERAGE_SAMPLE_FRACTION = "min_low_coverage_sample_fraction";
    private static final String MIN_LOW_COVERAGE_SAMPLE_FRACTION_DESC = "Min fraction of tumor samples with exon status WARN for an exon to be excluded";
    private static final double DEFAULT_MIN_LOW_COVERAGE_SAMPLE_FRACTION = 0.95;

    public ExonExclusionsConfig(final ConfigBuilder configBuilder)
    {
        Qsee = new QseeConfig(configBuilder);

        MinLowCoverageSampleFraction = configBuilder.getDecimal(MIN_LOW_COVERAGE_SAMPLE_FRACTION);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        QseeConfig.registerConfig(configBuilder);

        configBuilder.addDecimal(MIN_LOW_COVERAGE_SAMPLE_FRACTION, MIN_LOW_COVERAGE_SAMPLE_FRACTION_DESC, DEFAULT_MIN_LOW_COVERAGE_SAMPLE_FRACTION);
    }
}
