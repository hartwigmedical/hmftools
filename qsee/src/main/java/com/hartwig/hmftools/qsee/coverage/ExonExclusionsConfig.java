package com.hartwig.hmftools.qsee.coverage;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.QseeConfig;

public class ExonExclusionsConfig
{
    public final QseeConfig Qsee;

    public final double MinLowCoverageLodVaf;
    public final double MinLowCoverageCohortFraction;

    private static final String MIN_LOW_COVERAGE_LOD_VAF = "min_low_coverage_lod_vaf";
    private static final String MIN_LOW_COVERAGE_LOD_VAF_DESC = "Min LoD VAF for an exon to count as low coverage in a sample";
    private static final double DEFAULT_MIN_LOW_COVERAGE_LOD_VAF = 0.07;

    private static final String MIN_LOW_COVERAGE_COHORT_FRACTION = "min_low_coverage_cohort_fraction";
    private static final String MIN_LOW_COVERAGE_COHORT_FRACTION_DESC = "Min fraction of samples with low coverage for an exon to be excluded";
    private static final double DEFAULT_MIN_LOW_COVERAGE_COHORT_FRACTION = 0.3;

    public ExonExclusionsConfig(final ConfigBuilder configBuilder)
    {
        Qsee = new QseeConfig(configBuilder);

        MinLowCoverageLodVaf = configBuilder.getDecimal(MIN_LOW_COVERAGE_LOD_VAF);
        MinLowCoverageCohortFraction = configBuilder.getDecimal(MIN_LOW_COVERAGE_COHORT_FRACTION);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        QseeConfig.registerConfig(configBuilder);

        configBuilder.addDecimal(MIN_LOW_COVERAGE_LOD_VAF, MIN_LOW_COVERAGE_LOD_VAF_DESC, DEFAULT_MIN_LOW_COVERAGE_LOD_VAF);
        configBuilder.addDecimal(MIN_LOW_COVERAGE_COHORT_FRACTION, MIN_LOW_COVERAGE_COHORT_FRACTION_DESC, DEFAULT_MIN_LOW_COVERAGE_COHORT_FRACTION);
    }
}
