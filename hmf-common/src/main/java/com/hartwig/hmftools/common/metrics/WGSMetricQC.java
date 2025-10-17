package com.hartwig.hmftools.common.metrics;

public final class WGSMetricQC
{
    public static final double MIN_REF_10X_COVERAGE = 0.9;
    public static final double MIN_REF_20X_COVERAGE = 0.7;
    public static final double MIN_TUMOR_30X_COVERAGE = 0.8;
    public static final double MIN_TUMOR_60X_COVERAGE = 0.65;

    public static final int TUMOR_COVERAGE_LEVEL_30x = 30;
    public static final int TUMOR_COVERAGE_LEVEL_60x = 60;
    public static final int REF_COVERAGE_LEVEL_10x = 10;
    public static final int REF_COVERAGE_LEVEL_20x = 20;

    public static boolean hasSufficientCoverage(final BamMetricSummary tumorMetrics, final BamMetricSummary refMetrics)
    {
        boolean wgsQCRef10 = refMetrics.coveragePercent(REF_COVERAGE_LEVEL_10x) >= MIN_REF_10X_COVERAGE;
        boolean wgsQCRef20 = refMetrics.coveragePercent(REF_COVERAGE_LEVEL_20x) >= MIN_REF_20X_COVERAGE;
        boolean wgsQCTumor30 = tumorMetrics.coveragePercent(TUMOR_COVERAGE_LEVEL_30x) >= MIN_TUMOR_30X_COVERAGE;
        boolean wgsQCTumor60 = tumorMetrics.coveragePercent(TUMOR_COVERAGE_LEVEL_60x) >= MIN_TUMOR_60X_COVERAGE;

        return wgsQCRef10 && wgsQCRef20 && wgsQCTumor30 && wgsQCTumor60;
    }

    public static final double MIN_MAPPED_PROPORTION = 0.95;

    public static boolean hasSufficientMappedProportion(final BamFlagStats flagstat)
    {
        return hasSufficientMappedProportion(flagstat.mappedProportion());
    }

    public static boolean hasSufficientMappedProportion(double mappedProportion)
    {
        return mappedProportion >= MIN_MAPPED_PROPORTION;
    }

}
