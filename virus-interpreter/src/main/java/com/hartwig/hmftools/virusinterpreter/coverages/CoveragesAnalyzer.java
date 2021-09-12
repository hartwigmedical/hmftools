package com.hartwig.hmftools.virusinterpreter.coverages;

import java.io.IOException;

import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;

import org.jetbrains.annotations.NotNull;

public final class CoveragesAnalyzer {

    private CoveragesAnalyzer() {
    }

    @NotNull
    public static CoveragesAnalysis run(@NotNull String purplePurityTsv, @NotNull String purpleQcFile,
            @NotNull String tumorSampleWGSMetricsFile) throws IOException {
        return ImmutableCoveragesAnalysis.builder()
                .expectedClonalCoverage(calculateExpectedClonalCoverage(purplePurityTsv, purpleQcFile, tumorSampleWGSMetricsFile))
                .build();
    }

    static double calculateExpectedClonalCoverage(@NotNull String purplePurityTsv, @NotNull String purpleQcFile,
            @NotNull String tumorSampleWGSMetricsFile) throws IOException {
        PurityContext purityContext = PurityContextFile.readWithQC(purpleQcFile, purplePurityTsv);
        double ploidy = purityContext.bestFit().ploidy();
        double purity = purityContext.bestFit().purity();

        WGSMetrics metrics = WGSMetricsFile.read(tumorSampleWGSMetricsFile);
        double tumorMeanCoverage = metrics.meanCoverage();
        return tumorMeanCoverage * purity / ploidy;
    }
}