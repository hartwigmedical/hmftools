package com.hartwig.hmftools.viridian.reporting;

import java.util.OptionalDouble;

import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.purple.PurityContext;

public final class ClonalCoverage
{
    // Expected virus genome read depth from one virus genome copy carried in each tumor cell.
    // Divide a virus genome's mean depth by this number to get its copies per tumor cell.
    // Empty when Purple's fit cannot support the calculation.
    public static OptionalDouble expectedViralDepthPerCopy(PurityContext purityContext, BamMetricSummary tumorBamMetrics)
    {
        if(purityFitFailed(purityContext))
        {
            return OptionalDouble.empty();
        }
        double tumorPurity = purityContext.bestFit().purity();
        double tumorPloidy = purityContext.bestFit().ploidy();
        double hostMeanDepth = tumorBamMetrics.meanCoverage();
        return OptionalDouble.of(hostMeanDepth * tumorPurity / (tumorPurity * tumorPloidy + 2 * (1 - tumorPurity)));
    }

    private static boolean purityFitFailed(PurityContext purity)
    {
        return purity.qc().status().stream().anyMatch(status -> status.name().startsWith("FAIL"));
    }
}
