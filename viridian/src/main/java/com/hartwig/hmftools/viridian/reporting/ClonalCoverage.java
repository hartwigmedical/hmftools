package com.hartwig.hmftools.viridian.reporting;

import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.purple.PurityContext;

import org.jetbrains.annotations.Nullable;

public final class ClonalCoverage
{
    // Expected virus genome read depth from one virus genome copy carried in each tumor cell.
    // Divide a virus genome's mean depth by this number to get its copies per tumor cell.
    // Null when Purple's fit cannot support the calculation.
    @Nullable
    public static Double expectedViralDepthPerCopy(PurityContext purityContext, BamMetricSummary tumorBamMetrics)
    {
        if(purityFitFailed(purityContext))
        {
            return null;
        }
        double tumorPurity = purityContext.bestFit().purity();
        double tumorPloidy = purityContext.bestFit().ploidy();
        double hostMeanDepth = tumorBamMetrics.meanCoverage();
        return hostMeanDepth * tumorPurity / (tumorPurity * tumorPloidy + 2 * (1 - tumorPurity));
    }

    private static boolean purityFitFailed(PurityContext purity)
    {
        return purity.qc().status().stream().anyMatch(status -> status.name().startsWith("FAIL"));
    }
}
