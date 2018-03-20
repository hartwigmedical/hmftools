package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.result.MultiValueResult;
import com.hartwig.hmftools.healthchecker.result.PatientResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CoverageChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(CoverageChecker.class);

    public CoverageChecker() {
    }

    @NotNull
    public BaseResult run(@NotNull final RunContext runContext) throws IOException {
        WGSMetrics metrics = extractMetrics(runContext);
        final List<HealthCheck> refChecks = Lists.newArrayList(new HealthCheck(runContext.refSample(),
                        CoverageCheck.COVERAGE_10X.toString(),
                        String.valueOf(metrics.ref10xCoveragePercentage())),
                new HealthCheck(runContext.refSample(),
                        CoverageCheck.COVERAGE_20X.toString(),
                        String.valueOf(metrics.ref20xCoveragePercentage())));

        if (runContext.isSomaticRun()) {
            assert metrics.tumor30xCoveragePercentage() != null;
            assert metrics.tumor60xCoveragePercentage() != null;

            final List<HealthCheck> tumorChecks = Lists.newArrayList(new HealthCheck(runContext.tumorSample(),
                            CoverageCheck.COVERAGE_30X.toString(),
                            String.valueOf(metrics.tumor30xCoveragePercentage())),
                    new HealthCheck(runContext.tumorSample(),
                            CoverageCheck.COVERAGE_60X.toString(),
                            String.valueOf(metrics.tumor60xCoveragePercentage())));

            return toPatientResult(refChecks, tumorChecks);
        } else {
            return toMultiValueResult(refChecks);
        }
    }

    @NotNull
    private static WGSMetrics extractMetrics(@NotNull RunContext runContext) throws IOException {
        String refFile = WGSMetricsFile.generateFilename(runContext.runDirectory(), runContext.refSample());
        if (runContext.isSomaticRun()) {
            String tumorFile = WGSMetricsFile.generateFilename(runContext.runDirectory(), runContext.tumorSample());
            return WGSMetricsFile.read(refFile, tumorFile);
        } else {
            return WGSMetricsFile.read(refFile);
        }
    }

    @NotNull
    private BaseResult toPatientResult(@NotNull final List<HealthCheck> refChecks, @NotNull final List<HealthCheck> tumorChecks) {
        HealthCheck.log(LOGGER, refChecks);
        HealthCheck.log(LOGGER, tumorChecks);

        return new PatientResult(CheckType.COVERAGE, refChecks, tumorChecks);
    }

    @NotNull
    private static BaseResult toMultiValueResult(@NotNull final List<HealthCheck> checks) {
        HealthCheck.log(LOGGER, checks);

        return new MultiValueResult(CheckType.COVERAGE, checks);
    }
}
