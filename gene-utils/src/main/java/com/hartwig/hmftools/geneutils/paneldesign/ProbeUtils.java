package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_OVERLAP_MAX;

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
