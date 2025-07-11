package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_OVERLAP_MAX;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.UNCOVERED_BASES_MAX;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class ProbeUtils
{
    public static CandidateProbe probeCenteredAt(final String chromosome, int centrePosition, final ProbeFactory factory)
    {
        return probeStartingAt(chromosome, centrePosition - PROBE_LENGTH / 2, factory);
    }

    public static CandidateProbe probeStartingAt(final String chromosome, int startPosition, final ProbeFactory factory)
    {
        return factory.create(new ChrBaseRegion(chromosome, startPosition, startPosition + PROBE_LENGTH - 1));
    }

    public static boolean canCoverRegionWithOneProbe(final BaseRegion target)
    {
        return target.baseLength() - UNCOVERED_BASES_MAX <= PROBE_LENGTH;
    }

    public static int minProbeStartCovering(final BaseRegion target)
    {
        // Need to take care of the case where target is smaller than a probe: don't allow the probe to end before the target start.
        int minProbeEnd = max(target.start(), target.end() - UNCOVERED_BASES_MAX);
        return minProbeEnd - PROBE_LENGTH + 1;
    }

    public static int maxProbeEndCovering(final BaseRegion target)
    {
        // Need to take care of the case where target is smaller than a probe: don't allow the probe to start after the target end.
        int maxProbeStart = min(target.end(), target.start() + UNCOVERED_BASES_MAX);
        return maxProbeStart + PROBE_LENGTH - 1;
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

    // Calculates the next position a probe may start at, respecting the probe overlap constraint.
    public static int nextProbeStartPosition(int prevProbeEnd)
    {
        // prevProbeEnd - nextStart + 1 < PROBE_OVERLAP_MAX
        return prevProbeEnd - PROBE_OVERLAP_MAX + 1;
    }
}
