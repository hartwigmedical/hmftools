package com.hartwig.hmftools.virusinterpreter.coverages;

import java.io.IOException;

import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.common.purple.PurityContext;

import org.jetbrains.annotations.NotNull;

public final class CoveragesAnalyzer {

    private CoveragesAnalyzer() {
    }

    @NotNull
    public static CoveragesAnalysis run(@NotNull PurityContext purityContext,
            @NotNull String tumorSampleWGSMetricsFile) throws IOException {
        return ImmutableCoveragesAnalysis.builder()
                .expectedClonalCoverage(calculateExpectedClonalCoverage(purityContext, tumorSampleWGSMetricsFile))
                .build();
    }

    static Double calculateExpectedClonalCoverage(@NotNull PurityContext purityContext,
            @NotNull String tumorSampleWGSMetricsFile) throws IOException {
        double ploidy = purityContext.bestFit().ploidy();
        double purity = purityContext.bestFit().purity();

        WGSMetrics metrics = WGSMetricsFile.read(tumorSampleWGSMetricsFile);
        double tumorMeanCoverage = metrics.meanCoverage();
        return (tumorMeanCoverage * purity) / ploidy;
    }
}