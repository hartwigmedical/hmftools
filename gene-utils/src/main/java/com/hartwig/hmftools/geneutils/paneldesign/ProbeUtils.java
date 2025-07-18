package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.min;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_COVERAGE_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_OVERLAP_MAX;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Miscellaneous probe maths and utilities.
public class ProbeUtils
{
    public static CandidateProbe probeCenteredAt(final String chromosome, int centrePosition, final CandidateProbeContext context)
    {
        return probeStartingAt(chromosome, centrePosition - PROBE_LENGTH / 2, context);
    }

    public static CandidateProbe probeStartingAt(final String chromosome, int startPosition, final CandidateProbeContext context)
    {
        return context.createProbe(new ChrBaseRegion(chromosome, startPosition, startPosition + PROBE_LENGTH - 1));
    }

    // Calculates the minimum probe starting position such that the specified position is contained within the probe.
    public static int minProbeStartContaining(int targetPosition)
    {
        // start + PROBE_LENGTH - 1 >= targetPosition
        return targetPosition - PROBE_LENGTH + 1;
    }

    // Calculates the maximum probe ending position such that the specified position is contained within the probe.
    public static int maxProbeEndContaining(int targetPosition)
    {
        // end - PROBE_LENGTH + 1 <= targetPosition
        return targetPosition + PROBE_LENGTH - 1;
    }

    // TODO: unit test

    // minProbeStartContaining() and maxProbeEndContaining() in one call.
    public static BaseRegion probeBoundsContaining(int targetPosition)
    {
        return new BaseRegion(minProbeStartContaining(targetPosition), maxProbeEndContaining(targetPosition));
    }

    // Calculates the minimum probe starting position such that the probe sufficiently overlaps the target region.
    public static int minProbeStartCovering(final BaseRegion targetRegion)
    {
        // Must handle the case where the region is smaller than the required coverage.
        int end = targetRegion.start() + min(targetRegion.baseLength(), PROBE_COVERAGE_MIN) - 1;
        return end - PROBE_LENGTH + 1;
    }

    // Calculates the maximum probe starting position such that the probe sufficiently overlaps the target region.
    public static int maxProbeStartCovering(final BaseRegion targetRegion)
    {
        // Must handle the case where the region is smaller than the required coverage.
        return targetRegion.end() - min(targetRegion.baseLength(), PROBE_COVERAGE_MIN) + 1;
    }

    // Calculates the next position a probe may start at, respecting the probe overlap constraint.
    public static int nextProbeStartPosition(int prevProbeEnd)
    {
        // prevProbeEnd - nextStart + 1 < PROBE_OVERLAP_MAX
        return prevProbeEnd - PROBE_OVERLAP_MAX + 1;
    }
}
