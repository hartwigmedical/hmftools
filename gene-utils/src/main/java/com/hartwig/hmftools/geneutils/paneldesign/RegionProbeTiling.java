package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;

import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class RegionProbeTiling
{
    public static ProbeGenerationResult fillRegionWithProbes(final ChrBaseRegion region, final ProbeSourceInfo source,
            final ProbeEvaluator probeEvaluator)
    {
        // TODO
    }

    // Yields consecutive nonoverlapping probes starting at a position.
    public static Stream<BaseRegion> tileBaseRegionsFrom(int startPosition)
    {
        return IntStream.iterate(startPosition, start -> start + PROBE_LENGTH)
                .mapToObj(start -> new BaseRegion(start, start + PROBE_LENGTH - 1));
    }
}
