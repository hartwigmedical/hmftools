package com.hartwig.hmftools.virusinterpreter.coverages;

import java.io.IOException;

import com.hartwig.hmftools.common.metrics.BamMetricsSummary;
import com.hartwig.hmftools.common.purple.PurityContext;

public final class CoveragesAnalyzer
{
    public static CoveragesAnalysis run(final PurityContext purityContext, final BamMetricsSummary tumorMetrics) throws IOException
    {
        return new CoveragesAnalysis(calculateExpectedClonalCoverage(purityContext, tumorMetrics));
    }

    static Double calculateExpectedClonalCoverage(final PurityContext purityContext, final BamMetricsSummary tumorMetrics) throws IOException
    {
        double ploidy = purityContext.bestFit().ploidy();
        double purity = purityContext.bestFit().purity();

        double tumorMeanCoverage = tumorMetrics.meanCoverage();
        return (tumorMeanCoverage * purity) / ploidy;
    }
}