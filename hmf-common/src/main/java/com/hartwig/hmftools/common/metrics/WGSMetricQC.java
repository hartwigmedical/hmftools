package com.hartwig.hmftools.common.metrics;

import com.hartwig.hmftools.common.healthchecker.HealthCheckEvaluation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class WGSMetricQC {

    private static final Logger LOGGER = LogManager.getLogger(WGSMetricQC.class);

    private WGSMetricQC() {
    }

    public static ImmutableWGSMetricWithQC checkQCMetric(@NotNull WGSMetrics metrics) {

        // We only write metrics for somatic runs.
        assert metrics.tumor30xCoveragePercentage() != null;
        assert metrics.tumor60xCoveragePercentage() != null;
        assert metrics.tumorMeanCoverage() != null;

        if (metrics.tumor30xCoveragePercentage() == null || metrics.tumor60xCoveragePercentage() == null) {
            LOGGER.error("No somatic run");
        }

        boolean WGSqcRef10 = HealthCheckEvaluation.succeedCoverage(Double.toString(metrics.ref10xCoveragePercentage()),
                "Ref 10x",
                HealthCheckEvaluation.MIN_REF_10X_COVERAGE);
        boolean WGSqcRef20 = HealthCheckEvaluation.succeedCoverage(Double.toString(metrics.ref20xCoveragePercentage()),
                "Ref 20x",
                HealthCheckEvaluation.MIN_REF_20X_COVERAGE);
        boolean WGSqcTumor30 = HealthCheckEvaluation.succeedCoverage(Double.toString(metrics.tumor30xCoveragePercentage()),
                "Tumor 30x",
                HealthCheckEvaluation.MIN_TUMOR_30X_COVERAGE);
        boolean WGSqcTumor60 = HealthCheckEvaluation.succeedCoverage(Double.toString(metrics.tumor60xCoveragePercentage()),
                "Tumor 60x",
                HealthCheckEvaluation.MIN_TUMOR_60X_COVERAGE);

        boolean WGSqc = true;
        if (!WGSqcRef10 || !WGSqcRef20 || !WGSqcTumor30 || !WGSqcTumor60) {
            WGSqc = false;
        }

        if (WGSqc) {
            LOGGER.info("PASS - The run has enough coverage ");
        } else {
            LOGGER.info("FAIL - The run has not enough coverage ");
        }


        return ImmutableWGSMetricWithQC.builder().wgsMetrics(metrics).qcMetric(WGSqc).build();
    }
}
