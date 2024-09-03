package com.hartwig.hmftools.virusinterpreter.coverages;

import java.io.IOException;

import com.hartwig.hmftools.common.metrics.OldWGSMetrics;
import com.hartwig.hmftools.common.metrics.OldWGSMetricsFile;
import com.hartwig.hmftools.common.purple.PurityContext;

public final class CoveragesAnalyzer
{
    public static CoveragesAnalysis run(final PurityContext purityContext, final String tumorSampleWGSMetricsFile) throws IOException
    {
        return new CoveragesAnalysis(calculateExpectedClonalCoverage(purityContext, tumorSampleWGSMetricsFile));
    }

    static Double calculateExpectedClonalCoverage(final PurityContext purityContext, final String tumorSampleWGSMetricsFile) throws IOException
    {
        double ploidy = purityContext.bestFit().ploidy();
        double purity = purityContext.bestFit().purity();

        OldWGSMetrics metrics = OldWGSMetricsFile.read(tumorSampleWGSMetricsFile);
        double tumorMeanCoverage = metrics.meanCoverage();
        return (tumorMeanCoverage * purity) / ploidy;
    }
}