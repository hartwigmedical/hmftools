package com.hartwig.hmftools.geneutils.paneldesign;

import static java.util.Collections.emptyList;

import java.util.List;
import java.util.stream.Stream;

public record ProbeGenerationResult(
        List<TargetRegion> targetRegions,
        List<Probe> probes,
        List<RejectedRegion> rejectedRegions
)
{
    public ProbeGenerationResult()
    {
        this(emptyList(), emptyList(), emptyList());
    }

    public ProbeGenerationResult add(ProbeGenerationResult other)
    {
        // Naive merge, does not try to deduplicate data.
        return new ProbeGenerationResult(
                Stream.concat(targetRegions.stream(), other.targetRegions.stream()).toList(),
                Stream.concat(probes.stream(), other.probes.stream()).toList(),
                Stream.concat(rejectedRegions.stream(), other.rejectedRegions.stream()).toList()
        );
    }
}
