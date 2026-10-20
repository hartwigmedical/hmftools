package com.hartwig.hmftools.viridian.reporting;

import static com.hartwig.hmftools.common.purple.PurpleQCStatus.FAIL_CONTAMINATION;
import static com.hartwig.hmftools.common.purple.PurpleQCStatus.FAIL_NO_TUMOR;

import java.util.Set;

import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;

import org.jetbrains.annotations.Nullable;

public record ClonalityInputs(
        // Expected virus genome read depth from one virus genome copy carried in each tumor cell.
        // Divide a virus genome's mean depth by this number to get its copies per tumor cell.
        @Nullable Double expectedViralDepthPerCopy,
        // FAIL_NO_TUMOR or FAIL_CONTAMINATION, passed through to the output for visibility.
        @Nullable PurpleQCStatus purpleFailure
)
{
    public ClonalityInputs
    {
        if((expectedViralDepthPerCopy == null) == (purpleFailure == null))
        {
            throw new IllegalArgumentException("Must have either expectedViralDepthPerCopy or purpleFailure");
        }
        if(purpleFailure != null && purpleFailure != FAIL_NO_TUMOR && purpleFailure != FAIL_CONTAMINATION)
        {
            throw new IllegalArgumentException("Not a clonality-blocking Purple failure: " + purpleFailure);
        }
    }

    public static ClonalityInputs from(PurityContext purityContext, BamMetricSummary tumorBamMetrics)
    {
        Set<PurpleQCStatus> statuses = purityContext.qc().status();
        if(statuses.contains(FAIL_NO_TUMOR))
        {
            return new ClonalityInputs(null, FAIL_NO_TUMOR);
        }
        else if(statuses.contains(FAIL_CONTAMINATION))
        {
            return new ClonalityInputs(null, FAIL_CONTAMINATION);
        }
        else
        {
            double tumorPurity = purityContext.bestFit().purity();
            double tumorPloidy = purityContext.bestFit().ploidy();
            double hostMeanDepth = tumorBamMetrics.meanCoverage();
            double expectedViralDepthPerCopy =
                    hostMeanDepth * tumorPurity / (tumorPurity * tumorPloidy + 2 * (1 - tumorPurity));
            return new ClonalityInputs(expectedViralDepthPerCopy, null);
        }
    }
}
